"""
Opt-in safety limits of the CDS API route (``--cdsapi``) for controlled runs
such as the D-Bench benchmark. Nothing changes unless ``--cds-max-seconds`` is
given: the normal route keeps ``cdsapi.Client(timeout=600, retry_max=500)``.

Guarded mode (``--cds-max-seconds`` with ``--cds-ledger``):

- one global monotonic deadline per case, from the start of the acquisition
  through the queue, the transfers and the assembly of the file; no request
  is submitted once it has passed, and the wait for a queued request stops
  at it;
- few HTTP retries (``--cds-retries``) and a bounded network timeout
  (``--cds-timeout``); a request is never resubmitted;
- one active CDS submission at a time per ledger *and* per CDS account on
  this host (two exclusive locks taken without waiting: one next to the
  ledger, one per credential in ``LECTK_CDS_LOCK_DIR`` or ``~/.cache/lectk``,
  named after a truncated hash of the credential); runs on other hosts are
  not covered;
- a private ledger (JSON lines, mode 0600, outside Git) with an entry
  written *before* each submission and the job id as soon as the API gives
  it (``null`` if it never does);
- stop policy: after any failure, deadline or interruption with a request
  submitted, the ledger blocks every later submission until an operator has
  checked the account and appended a ``resolved`` entry
  (``python -m src.utils.cds_guard resolve LEDGER --confirm-no-pending --note ...``;
  refused while an acquisition holds the ledger lock); ``resolve`` is a human
  attestation, it does not query the CDS;
- an explicit credential is required (``CDSAPI_URL`` and ``CDSAPI_KEY``, or a
  private ``CDSAPI_RC``), so that ``~/.cdsapirc`` is never used silently; the
  key is never read, logged or written by this module;
- only the current CDS API (ecmwf-datastores-client behind cdsapi >= 0.7) is
  supported: it gives the job id at submission and lets the wait be bounded.

What the limits do NOT guarantee: stopping the client (deadline, error,
Ctrl+C, killed process) does not cancel the remote request, which may stay
queued or running and use the account's quota until it ends; its state is
recorded as uncertain. With ``--cds-delete-on-stop`` the guard also sends,
after a deadline or an interruption, a best-effort ``DELETE`` of *that* job
only: the remote object returned by this process's own submission, never an
id taken from the ledger or from a listing, and only while its status is
``accepted`` or ``running`` (a completed job and its results are never
deleted). It reads the status again and records whether the API confirmed it
(``dismissed``/``deleted``/not found).
An unconfirmed or failed delete changes nothing: the ledger blocks new
submissions in every case, and a killed process cannot send it at all. Before any new submission the account must be checked
for pending requests. A single blocking download can overrun the deadline by
its own duration (bounded by the network timeout only while stalled); the
overrun is detected and the case is reported as over the deadline, never as
a success.
"""

import argparse
import contextlib
import datetime
import fcntl
import hashlib
import json
import os
import sys
import time

DATASET = "reanalysis-era5-pressure-levels"

# Ledger states. A submission ends in one of TERMINAL; BLOCKING states (and a
# submission left without a terminal state) stop every later submission until
# a "resolved" entry is appended.
TERMINAL = {"downloaded", "failed", "abandoned", "interrupted"}
BLOCKING = {"failed", "abandoned", "interrupted"}

DEFAULT_RETRIES = 3
DEFAULT_TIMEOUT = 120
DEFAULT_SLEEP_MAX = 60


class CdsGuardError(RuntimeError):
    """Guarded CDS acquisition refused or stopped."""


class CdsDeadlineExceeded(CdsGuardError):
    """The global deadline of the case passed; the case is censored, not a success."""


class CdsSubmissionBlocked(CdsGuardError):
    """A new submission is refused (active submission, or unresolved stop in the ledger)."""


def _utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _category(exc):
    """Coarse class of a failure, for the ledger (no credential is ever included)."""
    status = getattr(getattr(exc, "response", None), "status_code", None)
    if status in (401, 403):
        return "credential"
    if status == 429:
        return "rate_limit"
    text = f"{type(exc).__name__} {exc}".lower()
    if any(word in text for word in ("licen", "quota", "limit", "policy", "forbidden")):
        return "policy_or_quota"
    if type(exc).__name__ == "ProcessingFailedError":
        return "remote_failed"
    if isinstance(exc, (TimeoutError,)) or "timeout" in text or "timed out" in text:
        return "timeout"
    return "error"


def _error(exc):
    return {"type": type(exc).__name__, "message": str(exc)[:500], "category": _category(exc)}


def read_ledger(path):
    """Entries of a ledger (empty if it does not exist)."""
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def blocking_reasons(entries):
    """
    Why a new submission must be refused, from the entries after the last
    ``resolved`` one: a blocking state, or a submission with no terminal state
    (its remote state is unknown, e.g. a killed process).
    """
    last_resolved = max((i for i, e in enumerate(entries) if e.get("state") == "resolved"), default=-1)
    reasons, last = [], {}
    for entry in entries[last_resolved + 1:]:
        if entry.get("state") in BLOCKING:
            reasons.append(f"{entry['state']}: case {entry.get('case_id')} day {entry.get('day')} "
                           f"job {entry.get('job_id')}")
        key = (entry.get("case_id"), entry.get("day"), entry.get("attempt"))
        if entry.get("state") in TERMINAL | {"submitting", "submitted"}:
            last[key] = entry
    for (case, day, _), entry in last.items():
        if entry["state"] in ("submitting", "submitted"):
            reasons.append(f"no final state for case {case} day {day} job {entry.get('job_id')} "
                           f"(remote state unknown)")
    return reasons


def check_credentials(environ=None):
    """
    Require an explicit credential for guarded mode: CDSAPI_URL and CDSAPI_KEY,
    or CDSAPI_RC naming a file readable by its owner only. Values are not read.
    """
    environ = os.environ if environ is None else environ
    if environ.get("CDSAPI_URL") and environ.get("CDSAPI_KEY"):
        return "environment"
    rc = environ.get("CDSAPI_RC")
    if rc:
        if os.path.realpath(rc) == os.path.realpath(os.path.expanduser("~/.cdsapirc")):
            raise CdsGuardError("CDSAPI_RC names the global ~/.cdsapirc; guarded mode needs an isolated credential")
        if not os.path.isfile(rc):
            raise CdsGuardError("CDSAPI_RC does not name an existing file")
        if os.stat(rc).st_mode & 0o077:
            raise CdsGuardError("the file named by CDSAPI_RC is readable by other users; use mode 600")
        return "CDSAPI_RC"
    raise CdsGuardError("guarded CDS mode needs an explicit credential (CDSAPI_URL and CDSAPI_KEY, or "
                        "CDSAPI_RC); ~/.cdsapirc is not used in this mode")


def _credential_identity(environ=None):
    """Truncated SHA-256 of the explicit credential (url and key), to name the account lock."""
    environ = os.environ if environ is None else environ
    url, key = environ.get("CDSAPI_URL"), environ.get("CDSAPI_KEY")
    rc = environ.get("CDSAPI_RC")
    if (not url or not key) and rc and os.path.isfile(rc):
        with open(rc, encoding="utf-8") as handle:
            for line in handle:
                name, _, value = line.partition(":")
                if name.strip() == "url" and not url:
                    url = value.strip()
                elif name.strip() == "key" and not key:
                    key = value.strip()
    if not url or not key:
        raise CdsGuardError("guarded CDS mode needs an explicit credential to lock its account")
    return hashlib.sha256(f"{url}\n{key}".encode("utf-8")).hexdigest()[:16]


def account_lock_path(environ=None):
    environ = os.environ if environ is None else environ
    directory = environ.get("LECTK_CDS_LOCK_DIR") or os.path.expanduser("~/.cache/lectk")
    os.makedirs(directory, mode=0o700, exist_ok=True)
    return os.path.join(directory, f"cds-account-{_credential_identity(environ)}.lock")


def _lock(path, message):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        os.close(fd)
        raise CdsSubmissionBlocked(message)
    return fd


def _unlock(fd):
    with contextlib.suppress(OSError):
        fcntl.flock(fd, fcntl.LOCK_UN)
    os.close(fd)


def check_client(client):
    """Guarded mode supports only the current CDS API client (ecmwf-datastores-client)."""
    module = type(client).__module__ or ""
    if not module.startswith("ecmwf.datastores"):
        raise CdsGuardError(
            f"guarded CDS mode needs the current CDS API client (ecmwf-datastores-client behind cdsapi >= 0.7), "
            f"got {type(client).__module__}.{type(client).__name__}; legacy 'UID:KEY' credentials are not supported"
        )


class CdsGuard:
    """
    Limits and ledger of one guarded CDS acquisition (one case).

    Args:
        max_seconds (float): Global deadline of the case, in seconds.
        ledger (str): Private ledger file (created with mode 0600).
        case_id (str): Identifier of the case in the ledger.
        retries (int): HTTP tries of the CDS client (``maximum_tries``).
        timeout (float): Network timeout of the CDS client, in seconds.
        sleep_max (float): Longest pause between two status checks.
        delete_on_stop (bool): After a deadline or an interruption, send a
            best-effort delete of the job just submitted (never of another).
        clock, sleep: Monotonic clock and sleep function (replaced in tests).
    """

    def __init__(self, max_seconds, ledger, case_id, retries=DEFAULT_RETRIES, timeout=DEFAULT_TIMEOUT,
                 sleep_max=DEFAULT_SLEEP_MAX, delete_on_stop=False, clock=time.monotonic, sleep=time.sleep):
        if max_seconds is None or max_seconds <= 0:
            raise CdsGuardError("--cds-max-seconds must be positive")
        if not ledger:
            raise CdsGuardError("guarded CDS mode needs --cds-ledger")
        if retries < 1:
            raise CdsGuardError("--cds-retries must be at least 1")
        self.max_seconds = float(max_seconds)
        self.ledger = str(ledger)
        self.case_id = str(case_id)
        self.retries = int(retries)
        self.timeout = float(timeout)
        self.sleep_max = float(sleep_max)
        self.delete_on_stop = bool(delete_on_stop)
        self._clock = clock
        self._sleep = sleep
        self._start = clock()
        self._locks = []

    @classmethod
    def from_args(cls, args):
        """The guard of a command line, or None when guarded mode is not requested."""
        max_seconds = getattr(args, "cds_max_seconds", None)
        if max_seconds is None:
            return None
        case_id = getattr(args, "cds_case_id", None) or os.path.basename(str(args.infile))
        return cls(max_seconds, getattr(args, "cds_ledger", None), case_id,
                   retries=getattr(args, "cds_retries", None) or DEFAULT_RETRIES,
                   timeout=getattr(args, "cds_timeout", None) or DEFAULT_TIMEOUT,
                   delete_on_stop=bool(getattr(args, "cds_delete_on_stop", False)))

    # ------------------------------------------------------------- client
    def client_options(self):
        """Keyword arguments of cdsapi.Client in guarded mode."""
        return {"timeout": self.timeout, "retry_max": self.retries, "sleep_max": self.sleep_max,
                "wait_until_complete": False}

    # -------------------------------------------------------------- time
    def elapsed(self):
        return self._clock() - self._start

    def remaining(self):
        return self.max_seconds - self.elapsed()

    def check_time(self, moment):
        """Stop the case when the deadline has passed (nothing is submitted after it)."""
        if self.remaining() <= 0:
            raise CdsDeadlineExceeded(
                f"CDS deadline of {self.max_seconds:g} s for case {self.case_id} passed {moment} "
                f"(after {self.elapsed():.1f} s); the case is censored"
            )

    # ------------------------------------------------------------ ledger
    def _append(self, state, **fields):
        entry = {"utc": _utc(), "case_id": self.case_id, "state": state,
                 "elapsed_s": round(self.elapsed(), 3), "pid": os.getpid()}
        entry.update(fields)
        fd = os.open(self.ledger, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        try:
            if os.fstat(fd).st_mode & 0o077:
                raise CdsGuardError(f"the ledger {self.ledger} is readable by other users; use mode 600")
            os.write(fd, (json.dumps(entry, sort_keys=True) + "\n").encode("utf-8"))
            os.fsync(fd)
        finally:
            os.close(fd)
        return entry

    def check_not_blocked(self):
        reasons = blocking_reasons(read_ledger(self.ledger))
        if reasons:
            raise CdsSubmissionBlocked(
                "new CDS submissions are blocked by the ledger until the account has been checked for pending "
                "requests and a 'resolved' entry appended: " + "; ".join(reasons)
            )

    # -------------------------------------------------------------- lock
    def acquire(self, environ=None):
        """Take the ledger lock and the account lock: one active CDS submission at a time."""
        self.check_not_blocked()
        account = account_lock_path(environ)
        self._locks.append(_lock(f"{self.ledger}.lock", "another guarded CDS acquisition holds the ledger lock "
                                                        "(one active CDS submission at a time)"))
        try:
            self._locks.append(_lock(account, "another guarded CDS acquisition is active on this CDS account "
                                              "(one active CDS submission at a time)"))
        except BaseException:
            self.release()
            raise

    def release(self):
        while self._locks:
            _unlock(self._locks.pop())

    # ------------------------------------------------------------ delete
    def _delete_own_job(self, remote, job_id, base):
        """Best-effort delete of the job just submitted, then its status again; recorded, never assumed."""
        if not (self.delete_on_stop and job_id) or getattr(remote, "request_id", None) != job_id:
            return  # only the job returned by this process's own submission
        outcome = {"job_id": job_id, "reply_status": None, "status_after": None}
        try:
            outcome["status_before"] = remote.status
        except Exception as exc:
            outcome["status_before"] = None
            outcome["status_error"] = _error(exc)
        if outcome["status_before"] not in ("accepted", "running"):
            # completed, failed or unknown: never delete results or act blindly
            self._append("delete_skipped", remote_state=outcome["status_before"] or "unknown", **outcome, **base,
                         note="delete only while accepted or running; the ledger still blocks new submissions")
            return
        try:
            reply = remote.delete()
            outcome["reply_status"] = (reply or {}).get("status") if isinstance(reply, dict) else None
        except Exception as exc:
            outcome["delete_error"] = _error(exc)
        try:
            outcome["status_after"] = remote.status
        except Exception as exc:
            code = getattr(getattr(exc, "response", None), "status_code", None)
            outcome["status_after"] = "not_found" if code == 404 else None
            outcome["status_error"] = _error(exc)
        outcome["confirmed_by_api"] = outcome["status_after"] in ("dismissed", "deleted", "not_found")
        self._append("delete_requested", remote_state=outcome["status_after"] or "unknown", **outcome, **base,
                     note="best-effort delete of this job only; the ledger still blocks new submissions")

    # ----------------------------------------------------------- request
    def retrieve(self, client, request, target, day):
        """
        Submit one request, wait for it within the deadline and download it.

        Returns:
            str: ``target``.
        """
        self.check_time(f"before submitting day {day}")
        self.check_not_blocked()
        entries = read_ledger(self.ledger)
        attempt = 1 + sum(e.get("state") == "submitting" and e.get("case_id") == self.case_id
                          and e.get("day") == day for e in entries)
        base = {"day": day, "attempt": attempt}
        self._append("submitting", job_id=None,
                     request={k: request[k] for k in ("year", "month", "day", "time", "area") if k in request},
                     **base)
        job_id = None
        try:
            try:
                remote = client.retrieve(DATASET, request)
            except Exception as exc:
                # the request may or may not exist remotely
                self._append("failed", job_id=None, error=_error(exc), remote_state="unknown", **base)
                raise
            job_id = getattr(remote, "request_id", None)
            self._append("submitted", job_id=job_id, **base)

            delay = 1.0
            while True:
                if self.remaining() <= 0:
                    self._append("abandoned", job_id=job_id, remote_state="unknown", **base,
                                 note="deadline passed while waiting; stopping the client does not cancel "
                                      "the remote request, which may still be queued or running")
                    self._delete_own_job(remote, job_id, base)
                    self.check_time(f"while waiting for job {job_id}")
                try:
                    ready = remote.results_ready
                except Exception as exc:
                    self._append("failed", job_id=job_id, error=_error(exc), **base)
                    raise
                if ready:
                    break
                self._sleep(max(0.0, min(delay, self.remaining())))
                delay = min(delay * 1.5, self.sleep_max)

            try:
                path = remote.get_results().download(target)
            except Exception as exc:
                self._append("failed", job_id=job_id, error=_error(exc), remote_state="successful", **base)
                raise
        except BaseException as exc:
            if not isinstance(exc, Exception):  # Ctrl+C, SystemExit: remote state unknown
                self._append("interrupted", job_id=job_id, remote_state="unknown",
                             error={"type": type(exc).__name__, "message": str(exc)[:500]}, **base)
                if job_id is not None:
                    self._delete_own_job(remote, job_id, base)
            raise
        size = os.path.getsize(path) if os.path.exists(path) else None
        self._append("downloaded", job_id=job_id, bytes=size, **base)
        self.check_time(f"after downloading day {day}")
        return target


# ------------------------------------------------------------------ command
def main(argv=None):
    """``status`` and ``resolve`` of a ledger (no network access)."""
    parser = argparse.ArgumentParser(prog="python -m src.utils.cds_guard",
                                     description="Ledger of the guarded CDS mode (no network access).")
    sub = parser.add_subparsers(dest="command", required=True)
    status = sub.add_parser("status", help="say whether new submissions are blocked")
    status.add_argument("ledger")
    resolve = sub.add_parser("resolve", help="unblock after checking the account for pending requests")
    resolve.add_argument("ledger")
    resolve.add_argument("--note", required=True, help="what was checked (no credentials)")
    resolve.add_argument("--confirm-no-pending", action="store_true",
                         help="attest that the account was checked and holds no pending request of the campaign")
    args = parser.parse_args(argv)

    entries = read_ledger(args.ledger)
    reasons = blocking_reasons(entries)
    if args.command == "status":
        print(f"{len(entries)} entries; " + ("BLOCKED: " + "; ".join(reasons) if reasons else "not blocked"))
        return 1 if reasons else 0
    if not args.confirm_no_pending or not args.note.strip():
        parser.error("resolve needs --confirm-no-pending and a non-empty --note")
    if not reasons:
        print("not blocked; nothing to resolve")
        return 0
    lock = _lock(f"{args.ledger}.lock", "an acquisition holds the ledger lock; resolve refused")
    try:
        guard = CdsGuard(1, args.ledger, "-")
        guard._append("resolved", note=args.note, cleared=reasons)
    finally:
        _unlock(lock)
    print(f"resolved {len(reasons)} blocking entr{'y' if len(reasons) == 1 else 'ies'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

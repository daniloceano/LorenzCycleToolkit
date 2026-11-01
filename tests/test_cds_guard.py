"""
Guarded CDS mode (src/utils/cds_guard.py), offline.

The CDS client is replaced by fakes that follow the interface of the current
CDS API client (ecmwf-datastores-client behind cdsapi >= 0.7): ``retrieve``
without a target returns a remote request with ``request_id``,
``results_ready`` and ``get_results().download(target)``. No connection is
opened and no credential is loaded: the credential variables hold dummy
values. They check the global deadline (before submitting, while queued,
after the download), the bounded client options, the ledger entry written
before each submission and the job id (known or not), the stop policy after
failures, deadlines and interruptions, the one-submission-at-a-time lock,
the credential and client checks, and that the normal --cdsapi route is
unchanged.
"""

import argparse
import json
import logging
import os
import stat
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from lorenzcycletoolkit import create_arg_parser
from src.utils import cds_guard as cg
from src.utils import preprocessing
from src.utils.cds_guard import CdsDeadlineExceeded, CdsGuard, CdsGuardError, CdsSubmissionBlocked

LOGGER = logging.getLogger("test_cds_guard")
DUMMY_KEY = "dummy-not-a-real-key-0000"
REQUEST = {"year": "1980", "month": "03", "day": "01", "time": ["00:00"], "area": [-10, -70, -60, -10],
           "variable": ["temperature"], "pressure_level": ["850"], "product_type": "reanalysis"}


class Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds


class HTTPError(Exception):
    def __init__(self, status):
        super().__init__(f"{status} Client Error")
        self.response = argparse.Namespace(status_code=status)


class ProcessingFailedError(Exception):
    pass


class FakeRemote:
    def __init__(self, statuses, request_id="job-1", clock=None, download_seconds=0.0, fail=None,
                 after_delete="dismissed", delete_error=None, status_before="accepted"):
        self.statuses = list(statuses)
        self.after_delete = after_delete
        self.status_before = status_before
        self.delete_error = delete_error
        self.deletes = 0
        if request_id is not None:
            self.request_id = request_id
        self.clock = clock
        self.download_seconds = download_seconds
        self.fail = fail

    @property
    def results_ready(self):
        status = self.statuses.pop(0) if len(self.statuses) > 1 else self.statuses[0]
        if isinstance(status, BaseException):
            raise status
        if status in ("failed", "rejected"):
            raise ProcessingFailedError("request failed")
        return status == "successful"

    @property
    def status(self):
        if not self.deletes:
            if isinstance(self.status_before, BaseException):
                raise self.status_before
            return self.status_before
        if isinstance(self.after_delete, BaseException):
            raise self.after_delete
        return self.after_delete

    def delete(self):
        self.delete_calls = getattr(self, "delete_calls", 0) + 1
        if self.delete_error:
            raise self.delete_error
        self.deletes += 1
        return {"jobID": getattr(self, "request_id", None), "status": "dismissed"}

    def get_results(self):
        return self

    def download(self, target):
        if self.fail:
            raise self.fail
        if self.clock:
            self.clock.now += self.download_seconds
        xr.Dataset({"t": ("valid_time", np.arange(2.0))},
                   coords={"valid_time": pd.date_range("1980-03-01", periods=2, freq="3h")}).to_netcdf(target)
        return target


class FakeClient:
    """Stands for ecmwf.datastores.legacy_client.LegacyClient."""

    def __init__(self, remote=None, ledger=None, submit_error=None, **options):
        self.options = options
        self.remote = remote
        self.ledger = ledger
        self.submit_error = submit_error
        self.calls = []

    def retrieve(self, name, request, target=None):
        # the ledger entry must exist before the request leaves the process
        if self.ledger:
            assert cg.read_ledger(self.ledger)[-1]["state"] == "submitting"
        self.calls.append((name, request, target))
        if self.submit_error:
            raise self.submit_error
        return self.remote


FakeClient.__module__ = "ecmwf.datastores.legacy_client"


@pytest.fixture
def clock():
    return Clock()


@pytest.fixture
def ledger(tmp_path):
    return str(tmp_path / "private" / "ledger.jsonl")


@pytest.fixture(autouse=True)
def private_dir(tmp_path, monkeypatch):
    (tmp_path / "private").mkdir(mode=0o700)
    # dummy credential and private lock directory: never a real key, never ~/.cache
    monkeypatch.setenv("CDSAPI_URL", "https://example.invalid/api")
    monkeypatch.setenv("CDSAPI_KEY", DUMMY_KEY)
    monkeypatch.delenv("CDSAPI_RC", raising=False)
    monkeypatch.setenv("LECTK_CDS_LOCK_DIR", str(tmp_path / "locks"))


def _guard(ledger, clock, max_seconds=100.0, case="19800296-cds-r1", delete_on_stop=False):
    return CdsGuard(max_seconds, ledger, case, clock=clock, sleep=clock.sleep, delete_on_stop=delete_on_stop)


def _states(ledger):
    return [e["state"] for e in cg.read_ledger(ledger)]


# ------------------------------------------------------------ retrieve
def test_success_writes_ledger_before_submission_and_job_id(tmp_path, ledger, clock):
    remote = FakeRemote(["accepted", "running", "successful"], clock=clock)
    client = FakeClient(remote, ledger=ledger)
    guard = _guard(ledger, clock)
    guard.retrieve(client, REQUEST, str(tmp_path / "day.nc"), day="1980-03-01")
    entries = cg.read_ledger(ledger)
    assert [e["state"] for e in entries] == ["submitting", "submitted", "downloaded"]
    assert entries[0]["job_id"] is None and entries[1]["job_id"] == "job-1"
    assert entries[2]["bytes"] > 0 and entries[0]["attempt"] == 1
    assert client.calls == [(cg.DATASET, REQUEST, None)]  # submitted without waiting (wait_until_complete=False)
    assert stat.S_IMODE(os.stat(ledger).st_mode) == 0o600
    assert DUMMY_KEY not in open(ledger).read()
    assert not cg.blocking_reasons(entries)


def test_unknown_job_id_is_null(tmp_path, ledger, clock):
    guard = _guard(ledger, clock)
    guard.retrieve(FakeClient(FakeRemote(["successful"], request_id=None)), REQUEST, str(tmp_path / "d.nc"), day="d")
    assert [e["job_id"] for e in cg.read_ledger(ledger)] == [None, None, None]


def test_deadline_while_queued_abandons_and_blocks_later_submissions(tmp_path, ledger, clock):
    guard = _guard(ledger, clock, max_seconds=30)
    remote = FakeRemote(["accepted"])
    with pytest.raises(CdsDeadlineExceeded, match="censored"):
        guard.retrieve(FakeClient(remote), REQUEST, str(tmp_path / "d.nc"), day="1980-03-01")
    entries = cg.read_ledger(ledger)
    assert [e["state"] for e in entries] == ["submitting", "submitted", "abandoned"]
    assert entries[-1]["remote_state"] == "unknown" and "does not cancel" in entries[-1]["note"]
    assert 30 <= entries[-1]["elapsed_s"] <= 31  # the wait stops at the deadline, not much later
    assert not (tmp_path / "d.nc").exists()

    # stop policy: a new case sharing the ledger submits nothing
    other = FakeClient(FakeRemote(["successful"]))
    with pytest.raises(CdsSubmissionBlocked, match="abandoned"):
        _guard(ledger, clock, case="19800318-cds-r1").retrieve(other, REQUEST, str(tmp_path / "e.nc"), day="x")
    assert other.calls == []


def test_resolved_entry_unblocks(tmp_path, ledger, clock):
    with pytest.raises(CdsDeadlineExceeded):
        _guard(ledger, clock, max_seconds=5).retrieve(FakeClient(FakeRemote(["accepted"])), REQUEST,
                                                      str(tmp_path / "d.nc"), day="d")
    assert cg.main(["status", ledger]) == 1
    assert cg.main(["resolve", ledger, "--confirm-no-pending", "--note", "account checked: no pending request"]) == 0
    assert cg.main(["status", ledger]) == 0
    _guard(ledger, clock).retrieve(FakeClient(FakeRemote(["successful"])), REQUEST, str(tmp_path / "e.nc"), day="d")
    assert _states(ledger)[-1] == "downloaded"


def test_deadline_before_submission_submits_nothing(tmp_path, ledger, clock):
    guard = _guard(ledger, clock, max_seconds=10)
    clock.now += 11
    client = FakeClient(FakeRemote(["successful"]))
    with pytest.raises(CdsDeadlineExceeded):
        guard.retrieve(client, REQUEST, str(tmp_path / "d.nc"), day="d")
    assert client.calls == [] and _states(ledger) == []


def test_download_over_deadline_is_censored_not_blocking(tmp_path, ledger, clock):
    remote = FakeRemote(["successful"], clock=clock, download_seconds=50)
    with pytest.raises(CdsDeadlineExceeded, match="after downloading"):
        _guard(ledger, clock, max_seconds=20).retrieve(FakeClient(remote), REQUEST, str(tmp_path / "d.nc"), day="d")
    assert _states(ledger) == ["submitting", "submitted", "downloaded"]
    assert not cg.blocking_reasons(cg.read_ledger(ledger))  # the remote request is complete


@pytest.mark.parametrize("error, category", [
    (HTTPError(401), "credential"),
    (HTTPError(403), "credential"),
    (HTTPError(429), "rate_limit"),
    (RuntimeError("required licences not accepted"), "policy_or_quota"),
    (RuntimeError("Your request is over the cost limit"), "policy_or_quota"),
    (TimeoutError("read timed out"), "timeout"),
])
def test_submission_errors_stop_and_block(tmp_path, ledger, clock, error, category):
    client = FakeClient(submit_error=error)
    with pytest.raises(type(error)):
        _guard(ledger, clock).retrieve(client, REQUEST, str(tmp_path / "d.nc"), day="d")
    last = cg.read_ledger(ledger)[-1]
    assert (last["state"], last["error"]["category"], last["remote_state"]) == ("failed", category, "unknown")
    with pytest.raises(CdsSubmissionBlocked):
        _guard(ledger, clock).retrieve(FakeClient(FakeRemote(["successful"])), REQUEST, str(tmp_path / "e"), day="e")
    assert len(client.calls) == 1  # never resubmitted


def test_remote_failure_is_recorded_with_job_id(tmp_path, ledger, clock):
    with pytest.raises(ProcessingFailedError):
        _guard(ledger, clock).retrieve(FakeClient(FakeRemote(["running", "failed"])), REQUEST,
                                       str(tmp_path / "d.nc"), day="d")
    last = cg.read_ledger(ledger)[-1]
    assert (last["state"], last["job_id"], last["error"]["category"]) == ("failed", "job-1", "remote_failed")


def test_download_failure_is_recorded(tmp_path, ledger, clock):
    remote = FakeRemote(["successful"], fail=OSError("connection reset"))
    with pytest.raises(OSError):
        _guard(ledger, clock).retrieve(FakeClient(remote), REQUEST, str(tmp_path / "d.nc"), day="d")
    last = cg.read_ledger(ledger)[-1]
    assert (last["state"], last["remote_state"]) == ("failed", "successful")


def test_interruption_while_queued_is_recorded_as_unknown(tmp_path, ledger, clock):
    remote = FakeRemote(["accepted", KeyboardInterrupt()])
    with pytest.raises(KeyboardInterrupt):
        _guard(ledger, clock).retrieve(FakeClient(remote), REQUEST, str(tmp_path / "d.nc"), day="d")
    last = cg.read_ledger(ledger)[-1]
    assert (last["state"], last["job_id"], last["remote_state"]) == ("interrupted", "job-1", "unknown")
    assert cg.blocking_reasons(cg.read_ledger(ledger))


def test_submission_without_final_state_blocks(ledger, clock):
    guard = _guard(ledger, clock)
    guard._append("submitting", day="d", attempt=1, job_id=None)
    guard._append("submitted", day="d", attempt=1, job_id="job-9")  # e.g. the process was killed here
    reasons = cg.blocking_reasons(cg.read_ledger(ledger))
    assert reasons and "job-9" in reasons[0] and "unknown" in reasons[0]


def test_attempts_are_counted_per_case_and_day(tmp_path, ledger, clock):
    _guard(ledger, clock).retrieve(FakeClient(FakeRemote(["successful"])), REQUEST, str(tmp_path / "a"), day="d")
    _guard(ledger, clock).retrieve(FakeClient(FakeRemote(["successful"])), REQUEST, str(tmp_path / "b"), day="d")
    assert [e["attempt"] for e in cg.read_ledger(ledger) if e["state"] == "submitting"] == [1, 2]


# ------------------------------------------------------ delete on stop
def _stop_with_delete(ledger, clock, remote, tmp_path):
    return _guard(ledger, clock, max_seconds=5, delete_on_stop=True).retrieve(
        FakeClient(remote), REQUEST, str(tmp_path / "d"), day="d")


def test_deadline_without_delete_option_sends_no_delete(tmp_path, ledger, clock):
    remote = FakeRemote(["accepted"])
    with pytest.raises(CdsDeadlineExceeded):
        _guard(ledger, clock, max_seconds=5).retrieve(FakeClient(remote), REQUEST, str(tmp_path / "d"), day="d")
    assert not hasattr(remote, "delete_calls") and "delete_requested" not in _states(ledger)


def test_deadline_with_delete_option_deletes_own_job_and_records_confirmation(tmp_path, ledger, clock):
    remote = FakeRemote(["accepted"], request_id="job-7")
    with pytest.raises(CdsDeadlineExceeded):
        _stop_with_delete(ledger, clock, remote, tmp_path)
    assert remote.deletes == 1
    entry = cg.read_ledger(ledger)[-1]
    assert (entry["state"], entry["job_id"], entry["status_after"], entry["confirmed_by_api"]) == \
        ("delete_requested", "job-7", "dismissed", True)
    assert cg.blocking_reasons(cg.read_ledger(ledger))  # still blocked: confirmation is not a guarantee


@pytest.mark.parametrize("remote_kwargs, status_after, confirmed", [
    (dict(after_delete="running"), "running", False),                          # not stopped
    (dict(delete_error=HTTPError(500)), "accepted", False),                    # delete refused
    (dict(after_delete=HTTPError(404)), "not_found", True),                    # job gone
    (dict(after_delete=OSError("down")), None, False),                         # state unknown after delete
])
def test_delete_outcomes_are_recorded_honestly(tmp_path, ledger, clock, remote_kwargs, status_after, confirmed):
    remote = FakeRemote(["accepted"], **remote_kwargs)
    with pytest.raises(CdsDeadlineExceeded):
        _stop_with_delete(ledger, clock, remote, tmp_path)
    entry = cg.read_ledger(ledger)[-1]
    assert (entry["status_after"], entry["confirmed_by_api"]) == (status_after, confirmed)
    assert cg.blocking_reasons(cg.read_ledger(ledger))


def test_interruption_with_delete_option_deletes_own_job(tmp_path, ledger, clock):
    remote = FakeRemote(["accepted", KeyboardInterrupt()])
    with pytest.raises(KeyboardInterrupt):
        _guard(ledger, clock, delete_on_stop=True).retrieve(FakeClient(remote), REQUEST, str(tmp_path / "d"), day="d")
    assert remote.delete_calls == 1 and _states(ledger)[-2:] == ["interrupted", "delete_requested"]


def test_no_delete_without_job_id_or_after_remote_failure(tmp_path, ledger, clock):
    remote = FakeRemote(["accepted"], request_id=None)
    with pytest.raises(CdsDeadlineExceeded):
        _stop_with_delete(ledger, clock, remote, tmp_path)
    failed = FakeRemote(["failed"])
    with pytest.raises(ProcessingFailedError):
        _guard(str(tmp_path / "private" / "l2"), clock, delete_on_stop=True).retrieve(
            FakeClient(failed), REQUEST, str(tmp_path / "e"), day="e")
    assert not hasattr(remote, "delete_calls") and not hasattr(failed, "delete_calls")


def test_delete_skipped_when_job_already_completed(tmp_path, ledger, clock):
    remote = FakeRemote(["accepted"], status_before="successful")
    with pytest.raises(CdsDeadlineExceeded):
        _stop_with_delete(ledger, clock, remote, tmp_path)
    entry = cg.read_ledger(ledger)[-1]
    assert not hasattr(remote, "delete_calls")  # completed results are never deleted
    assert (entry["state"], entry["status_before"]) == ("delete_skipped", "successful")
    assert cg.blocking_reasons(cg.read_ledger(ledger))


def test_delete_skipped_when_status_unknown(tmp_path, ledger, clock):
    remote = FakeRemote(["accepted"], status_before=OSError("down"))
    with pytest.raises(CdsDeadlineExceeded):
        _stop_with_delete(ledger, clock, remote, tmp_path)
    assert not hasattr(remote, "delete_calls") and _states(ledger)[-1] == "delete_skipped"


def test_delete_only_for_the_job_of_this_submission(ledger, clock):
    guard = _guard(ledger, clock, delete_on_stop=True)
    other = FakeRemote(["accepted"], request_id="someone-else")
    guard._delete_own_job(other, "job-1", {"day": "d", "attempt": 1})
    assert not hasattr(other, "delete_calls") and _states(ledger) == []


# ------------------------------------------------------------ resolve
def test_resolve_needs_attestation_and_note(ledger, clock):
    _guard(ledger, clock)._append("abandoned", day="d", attempt=1, job_id="j")
    for argv in (["resolve", ledger, "--note", "x"], ["resolve", ledger, "--confirm-no-pending", "--note", " "]):
        with pytest.raises(SystemExit):
            cg.main(argv)
    assert cg.main(["status", ledger]) == 1


def test_resolve_refused_during_an_active_acquisition(ledger, clock):
    _guard(ledger, clock)._append("abandoned", day="d", attempt=1, job_id="j")
    active = _guard(ledger, clock, case="active")
    active._locks.append(cg._lock(f"{ledger}.lock", "held"))  # an acquisition in progress
    try:
        with pytest.raises(CdsSubmissionBlocked, match="resolve refused"):
            cg.main(["resolve", ledger, "--confirm-no-pending", "--note", "checked"])
    finally:
        active.release()
    assert cg.main(["status", ledger]) == 1


def test_resolve_without_block_writes_nothing(ledger):
    assert cg.main(["resolve", ledger, "--confirm-no-pending", "--note", "checked"]) == 0
    assert not os.path.exists(ledger)


# ------------------------------------------------- lock, ledger, checks
def test_one_active_submission_per_account_across_ledgers(tmp_path, clock):
    first = _guard(str(tmp_path / "private" / "a.jsonl"), clock)
    first.acquire()
    try:
        with pytest.raises(CdsSubmissionBlocked, match="this CDS account"):
            _guard(str(tmp_path / "private" / "b.jsonl"), clock, case="other").acquire()
    finally:
        first.release()


def test_other_account_is_not_blocked(tmp_path, clock, monkeypatch):
    first = _guard(str(tmp_path / "private" / "a.jsonl"), clock)
    first.acquire()
    try:
        monkeypatch.setenv("CDSAPI_KEY", "another-dummy-key")
        second = _guard(str(tmp_path / "private" / "b.jsonl"), clock, case="other")
        second.acquire()
        second.release()
    finally:
        first.release()


def test_account_lock_name_does_not_reveal_the_key(tmp_path):
    path = cg.account_lock_path({"CDSAPI_URL": "https://example.invalid/api", "CDSAPI_KEY": DUMMY_KEY,
                                 "LECTK_CDS_LOCK_DIR": str(tmp_path / "locks")})
    assert DUMMY_KEY not in path and os.path.basename(path).startswith("cds-account-")
    assert oct(os.stat(tmp_path / "locks").st_mode & 0o777) == "0o700"


def test_global_cdsapirc_is_refused_as_cdsapi_rc(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    (home / ".cdsapirc").write_text("url: u\nkey: k\n")
    os.chmod(home / ".cdsapirc", 0o600)
    monkeypatch.setenv("HOME", str(home))
    with pytest.raises(CdsGuardError, match="global ~/.cdsapirc"):
        cg.check_credentials({"CDSAPI_RC": str(home / ".cdsapirc")})


def test_one_active_submission_per_ledger(ledger, clock):
    first = _guard(ledger, clock)
    first.acquire()
    try:
        with pytest.raises(CdsSubmissionBlocked, match="one active"):
            _guard(ledger, clock, case="other").acquire()
    finally:
        first.release()
    second = _guard(ledger, clock, case="other")
    second.acquire()
    second.release()


def test_ledger_readable_by_others_is_refused(ledger, clock):
    with open(ledger, "w"):
        pass
    os.chmod(ledger, 0o644)
    with pytest.raises(CdsGuardError, match="readable by other users"):
        _guard(ledger, clock)._append("submitting")


def test_credentials_must_be_explicit(tmp_path):
    with pytest.raises(CdsGuardError, match="explicit credential"):
        cg.check_credentials({"HOME": str(tmp_path)})
    assert cg.check_credentials({"CDSAPI_URL": "https://example.invalid/api", "CDSAPI_KEY": DUMMY_KEY}) == "environment"
    rc = tmp_path / "rc"
    rc.write_text("url: https://example.invalid/api\nkey: dummy\n")
    os.chmod(rc, 0o644)
    with pytest.raises(CdsGuardError, match="mode 600"):
        cg.check_credentials({"CDSAPI_RC": str(rc)})
    os.chmod(rc, 0o600)
    assert cg.check_credentials({"CDSAPI_RC": str(rc)}) == "CDSAPI_RC"


def test_legacy_client_is_refused():
    class LegacyResultClient:
        pass
    with pytest.raises(CdsGuardError, match="current CDS API client"):
        cg.check_client(LegacyResultClient())
    cg.check_client(FakeClient())


def test_invalid_limits_are_refused(ledger):
    for kwargs in (dict(max_seconds=0), dict(max_seconds=-1), dict(max_seconds=None)):
        with pytest.raises(CdsGuardError):
            CdsGuard(ledger=ledger, case_id="c", **kwargs)
    with pytest.raises(CdsGuardError):
        CdsGuard(10, "", "c")
    with pytest.raises(CdsGuardError):
        CdsGuard(10, ledger, "c", retries=0)


def test_client_options_bound_retries_and_do_not_wait(ledger, clock):
    guard = CdsGuard(60, ledger, "c", retries=2, timeout=30, clock=clock, sleep=clock.sleep)
    assert guard.client_options() == {"timeout": 30.0, "retry_max": 2, "sleep_max": cg.DEFAULT_SLEEP_MAX,
                                      "wait_until_complete": False}


# ------------------------------------------------ get_cdsapi_data / CLI
@pytest.fixture
def track():
    index = pd.to_datetime(["1980-03-01 00:00", "1980-03-01 03:00"])
    return pd.DataFrame({"Lat": [-30.0, -31.0], "Lon": [-50.0, -49.0]}, index=index)


def _cli_args(tmp_path, ledger, **extra):
    argv = [str(tmp_path / "case.nc"), "-t", "--cdsapi", "--cds-max-seconds", "600", "--cds-ledger", ledger]
    for key, value in extra.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    return create_arg_parser().parse_args(argv)


def test_get_cdsapi_data_guarded(tmp_path, ledger, track, monkeypatch):
    from src.utils.tools import get_cdsapi_data

    monkeypatch.setenv("CDSAPI_URL", "https://example.invalid/api")
    monkeypatch.setenv("CDSAPI_KEY", DUMMY_KEY)
    made = []

    def client_factory(**options):
        made.append(options)
        return FakeClient(FakeRemote(["accepted", "successful"]), ledger=ledger, **options)

    args = _cli_args(tmp_path, ledger, cds_retries=2, cds_case_id="19800296-cds-r1")
    with patch("cdsapi.Client", side_effect=client_factory):
        assert get_cdsapi_data(args, track, LOGGER) == args.infile
    assert made == [{"timeout": 120.0, "retry_max": 2, "sleep_max": cg.DEFAULT_SLEEP_MAX,
                     "wait_until_complete": False}]
    assert os.path.exists(args.infile)
    entries = cg.read_ledger(ledger)
    assert [e["state"] for e in entries] == ["submitting", "submitted", "downloaded"]
    assert {e["case_id"] for e in entries} == {"19800296-cds-r1"}
    assert DUMMY_KEY not in open(ledger).read()
    assert not os.path.exists(f"{ledger}.lock") or not cg.blocking_reasons(entries)


def test_get_cdsapi_data_guarded_refuses_implicit_credential(tmp_path, ledger, track, monkeypatch):
    from src.utils.tools import get_cdsapi_data

    for name in ("CDSAPI_URL", "CDSAPI_KEY", "CDSAPI_RC"):
        monkeypatch.delenv(name, raising=False)
    with patch("cdsapi.Client") as client, pytest.raises(CdsGuardError, match="explicit credential"):
        get_cdsapi_data(_cli_args(tmp_path, ledger), track, LOGGER)
    client.assert_not_called()
    assert not os.path.exists(tmp_path / "case.nc")


def test_normal_cdsapi_route_is_unchanged(tmp_path, track):
    from src.utils.tools import get_cdsapi_data

    args = create_arg_parser().parse_args([str(tmp_path / "case.nc"), "-t", "--cdsapi"])
    with patch("cdsapi.Client") as client, patch("xarray.open_dataset"), patch("xarray.concat") as concat:
        client.return_value.retrieve.side_effect = lambda name, request, target: open(target, "w").close()
        concat.return_value.to_netcdf.side_effect = lambda path: open(path, "w").close()
        get_cdsapi_data(args, track, LOGGER)
    client.assert_called_once_with(timeout=600, retry_max=500)
    name, request, target = client.return_value.retrieve.call_args.args
    assert name == cg.DATASET and target.endswith(".nc")


def test_guarded_mode_never_reuses_an_existing_file(tmp_path, ledger):
    existing = tmp_path / "case.nc"
    existing.write_text("old")
    args = _cli_args(tmp_path, ledger)
    with patch.object(preprocessing, "get_cdsapi_data") as download, pytest.raises(CdsGuardError, match="not reused"):
        preprocessing.get_data(args, LOGGER)
    download.assert_not_called()


@pytest.mark.parametrize("argv, message", [
    (["--cds-max-seconds", "60", "--cds-ledger", "l"], "need --cdsapi"),
    (["--cdsapi", "--cds-ledger", "l"], "need --cds-max-seconds"),
    (["--cdsapi", "--cds-max-seconds", "60"], "needs --cds-ledger"),
    (["--cdsapi", "--cds-max-seconds", "0", "--cds-ledger", "l"], "must be positive"),
    (["--cdsapi", "--cds-max-seconds", "60", "--cds-ledger", "l", "--cds-retries", "0"], "at least 1"),
    (["--ncar", "--cds-max-seconds", "60", "--cds-ledger", "l"], "need --cdsapi"),
    (["--cdsapi", "--cds-delete-on-stop"], "need --cds-max-seconds"),
])
def test_cli_rejects_inconsistent_guard_options(monkeypatch, capsys, argv, message):
    import lorenzcycletoolkit

    monkeypatch.setattr("sys.argv", ["lorenzcycletoolkit.py", "case.nc", "-t", *argv])
    with pytest.raises(SystemExit) as stop:
        lorenzcycletoolkit.main()
    assert stop.value.code == 2
    assert message in capsys.readouterr().err


def test_ledger_is_valid_json_lines(tmp_path, ledger, clock):
    _guard(ledger, clock).retrieve(FakeClient(FakeRemote(["successful"])), REQUEST, str(tmp_path / "d"), day="d")
    for line in open(ledger):
        entry = json.loads(line)
        assert {"utc", "case_id", "state", "elapsed_s", "pid"} <= set(entry)

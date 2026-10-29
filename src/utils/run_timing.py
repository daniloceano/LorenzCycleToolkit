"""
Wall-clock timing of one LorenzCycleToolkit run, reported in ``run_status.json``
(key ``timing``) and in the closing lines of the command line.

Every duration is measured with ``time.perf_counter`` (monotonic), never with
the system clock. The phases do not overlap:

    acquisition   the call that obtains the input file: get_cdsapi_data
                  (a --cdsapi download) or get_ncar_data (a --ncar download,
                  or the provenance check that lets an existing file be
                  reused). ``null`` when no such call runs: a local input file,
                  or an existing --cdsapi file (reused without any check).
    preparation   prepare_data without the acquisition: namelist, opening,
                  validation, processing and slicing of the data.
    lec           the Lorenz Energy Cycle framework (with --plots, the figures
                  drawn by the framework are included; see ``plots_requested``).
    total         from the start of the command line to the end of the run.

A phase that was not reached stays ``null``; a phase interrupted by an error
keeps the time measured until the error, with ``completed: false``. Nothing
is estimated.

``input_data`` says where the input file came from: ``downloaded`` (written
by this run), ``reused`` (an existing file accepted without downloading) or
``local`` (a file given without --cdsapi/--ncar). With --ncar the transfer
and assembly times of the provenance record are copied only when this run
wrote that record; they are sub-intervals of the acquisition, not extra time.
"""

import contextlib
import datetime
import json
import os
import time

CLOCK = "time.perf_counter (monotonic)"
PHASES = ("acquisition", "preparation", "lec")


def utc_now():
    """Current UTC time, truncated to the second like the NCAR provenance record."""
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0)


class RunTiming:
    """
    Durations of the phases of a run.

    Args:
        start (float): ``time.perf_counter()`` value of the start of the run
            (defaults to now).
    """

    def __init__(self, start=None):
        self._start = time.perf_counter() if start is None else start
        self._raw = {}
        self._completed = {}
        self.total_seconds = None
        self.provider = None
        self.input_data = None
        self.ncar = None
        self.plots_requested = None

    @contextlib.contextmanager
    def phase(self, name):
        """Measure a phase; the time until an exception is kept, marked not completed."""
        started = time.perf_counter()
        self._completed[name] = False
        try:
            yield self
            self._completed[name] = True
        finally:
            self._raw[name] = time.perf_counter() - started

    def set_source(self, provider, input_data):
        """Record the data source ('cds', 'ncar', 'local') and the origin of the file."""
        self.provider = provider
        self.input_data = input_data

    def record_ncar(self, path, retrieval_started_utc):
        """
        Record the --ncar origin of ``path`` from its provenance record: the file
        was downloaded by this run only if the record was started at or after
        ``retrieval_started_utc`` (another run may have written it meanwhile).
        Its transfer/assembly times are kept only in that case.
        """
        record = None
        with contextlib.suppress(OSError, ValueError):
            with open(f"{path}.provenance.json", encoding="utf-8") as handle:
                record = json.load(handle)
        timing = (record or {}).get("timing") or {}
        try:
            started = datetime.datetime.strptime(timing["started_utc"], "%Y-%m-%dT%H:%M:%SZ")
            started = started.replace(tzinfo=datetime.timezone.utc)
        except (KeyError, TypeError, ValueError):
            started = None
        if started is None:
            self.set_source("ncar", None)  # origin not verifiable: not guessed
            return
        if started < retrieval_started_utc:
            self.set_source("ncar", "reused")
            return
        self.set_source("ncar", "downloaded")
        totals = record.get("totals") or {}
        self.ncar = {
            "transfer_seconds": timing.get("transfer_seconds"),
            "assembly_seconds": timing.get("assembly_seconds"),
            "workers": timing.get("workers"),
            "downloaded_responses": totals.get("downloaded_responses"),
            "reused_responses": totals.get("reused_responses"),
            "downloaded_bytes": totals.get("downloaded_bytes"),
            "source": os.path.basename(f"{path}.provenance.json"),
        }

    def stop(self):
        """Close the total (only once)."""
        if self.total_seconds is None:
            self.total_seconds = time.perf_counter() - self._start

    def seconds(self, name):
        """Duration of a phase in seconds, or None if it was not measured."""
        if name not in self._raw:
            return None
        value = self._raw[name]
        if name == "preparation":
            # prepare_data contains the acquisition; report it only once
            value -= self._raw.get("acquisition", 0.0)
        return max(0.0, value)

    def to_dict(self):
        def rounded(value):
            return None if value is None else round(value, 3)

        return {
            "clock": CLOCK,
            "provider": self.provider,
            "input_data": self.input_data,
            "plots_requested": self.plots_requested,
            "phases": {
                name: {"seconds": rounded(self.seconds(name)),
                       "completed": self._completed.get(name) if name in self._raw else None}
                for name in PHASES
            },
            "total_seconds": rounded(self.total_seconds),
            "ncar_provenance": self.ncar,
        }

    def summary(self):
        """One line for the end of the command line."""
        def show(value):
            return "n/a" if value is None else f"{value:.2f} s"

        source = self.provider or "unknown"
        origin = self.input_data or "unknown"
        parts = [f"⏱️ Timing ({source}, {origin} data): "
                 + ", ".join(f"{name} {show(self.seconds(name))}" for name in PHASES)
                 + f", total {show(self.total_seconds)}"]
        if self.ncar:
            parts.append(f"(NCAR transfer {show(self.ncar.get('transfer_seconds'))}, "
                         f"assembly {show(self.ncar.get('assembly_seconds'))})")
        return " ".join(parts)

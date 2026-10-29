"""
Timing of a run (run_status.json ``timing`` and the closing CLI line).

They check that the phases are measured with a monotonic clock and do not
overlap; that the origin of the input file is reported for a local file,
a --cdsapi download, a reused --cdsapi file, a --ncar download, a reused
--ncar file and a --ncar file written by another run in the same second
(through the real get_ncar_data); that the NCAR transfer/assembly times
come only from this run's own retrieval; and that a failure keeps the time measured until the error and
leaves the unreached phases ``null``. Nothing is downloaded: the CDS and
NCAR calls are replaced by mocks, and the CLI test uses the offline launcher.
"""

import argparse
import contextlib
import datetime
import json
import logging
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
import xarray as xr

from src.utils import ncar_gdex, preprocessing
from src.utils import run_timing as rt
from src.utils.run_status import RunStatus
from src.utils.run_timing import RunTiming
from tests import output_checks as oc
from tests.test_cli_exit_status import BLOCKED, FIXED_RESULTS, TRACK_RESULTS, cli

LOGGER = logging.getLogger("test_run_timing")
REPO_ROOT = Path(__file__).resolve().parents[1]


class FakeClock:
    """Replaces time.perf_counter in run_timing with a controlled clock."""

    def __init__(self):
        self.now = 100.0

    def __call__(self):
        return self.now


@pytest.fixture
def clock(monkeypatch):
    fake = FakeClock()
    monkeypatch.setattr(rt.time, "perf_counter", fake)
    return fake


def _netcdf(path):
    xr.Dataset({"t": ("time", np.arange(2.0))}).to_netcdf(path)
    return path


def _args(tmp_path, infile, **flags):
    track = tmp_path / "track"
    track.write_text("time;Lat;Lon\n1980-01-01-0000;-30.0;-50.0\n1980-01-01-0300;-31.0;-49.0\n")
    options = dict(infile=str(infile), cdsapi=False, ncar=False, trackfile=str(track), time_resolution=3)
    return argparse.Namespace(**dict(options, **flags))


# ------------------------------------------------------------- RunTiming
def test_phases_do_not_overlap_and_unmeasured_are_null(clock):
    timing = RunTiming()
    with timing.phase("preparation"):
        clock.now += 1.0
        with timing.phase("acquisition"):
            clock.now += 5.0
        clock.now += 2.0
    clock.now += 0.5
    timing.stop()
    clock.now += 99.0
    timing.stop()  # the total is closed only once
    data = timing.to_dict()
    assert data["clock"] == rt.CLOCK
    assert data["phases"]["acquisition"] == {"seconds": 5.0, "completed": True}
    assert data["phases"]["preparation"] == {"seconds": 3.0, "completed": True}
    assert data["phases"]["lec"] == {"seconds": None, "completed": None}
    assert data["total_seconds"] == 8.5
    assert data["ncar_provenance"] is None
    assert "lec n/a" in timing.summary() and "total 8.50 s" in timing.summary()


def test_failed_phase_keeps_time_until_error(clock):
    timing = RunTiming()
    with pytest.raises(RuntimeError):
        with timing.phase("lec"):
            clock.now += 4.0
            raise RuntimeError("boom")
    assert timing.to_dict()["phases"]["lec"] == {"seconds": 4.0, "completed": False}


def test_timing_uses_monotonic_clock_not_system_time():
    # time.time() moved backwards must not affect the durations
    timing = RunTiming()
    with patch("time.time", return_value=0.0), timing.phase("lec"):
        pass
    assert timing.seconds("lec") >= 0.0


# --------------------------------------------------------- data sources
def test_local_file(tmp_path):
    args = _args(tmp_path, _netcdf(tmp_path / "local.nc"))
    timing = RunTiming()
    with patch.object(preprocessing, "get_cdsapi_data") as cds, patch.object(preprocessing, "get_ncar_data") as ncar:
        preprocessing.get_data(args, LOGGER, timing=timing)
    cds.assert_not_called()
    ncar.assert_not_called()
    data = timing.to_dict()
    assert (data["provider"], data["input_data"]) == ("local", "local")
    assert data["phases"]["acquisition"]["seconds"] is None


def test_get_data_without_timing_is_unchanged(tmp_path):
    args = _args(tmp_path, _netcdf(tmp_path / "local.nc"))
    assert "t" in preprocessing.get_data(args, LOGGER)


def test_cds_new_download(tmp_path, clock):
    infile = tmp_path / "cds.nc"
    args = _args(tmp_path, infile, cdsapi=True)

    def download(args, track, logger):
        clock.now += 30.0
        _netcdf(infile)
        return str(infile)

    timing = RunTiming()
    with patch.object(preprocessing, "get_cdsapi_data", side_effect=download) as cds:
        preprocessing.get_data(args, LOGGER, timing=timing)
    cds.assert_called_once()
    data = timing.to_dict()
    assert (data["provider"], data["input_data"]) == ("cds", "downloaded")
    assert data["phases"]["acquisition"] == {"seconds": 30.0, "completed": True}


def test_cds_reused_file(tmp_path):
    args = _args(tmp_path, _netcdf(tmp_path / "cds.nc"), cdsapi=True)
    timing = RunTiming()
    with patch.object(preprocessing, "get_cdsapi_data") as cds:
        preprocessing.get_data(args, LOGGER, timing=timing)
    cds.assert_not_called()
    data = timing.to_dict()
    assert (data["provider"], data["input_data"]) == ("cds", "reused")
    assert data["phases"]["acquisition"]["seconds"] is None  # no acquisition call ran


def test_cds_failure(tmp_path, clock):
    args = _args(tmp_path, tmp_path / "cds.nc", cdsapi=True)

    def fail(*_):
        clock.now += 7.0
        raise RuntimeError("CDS request rejected")

    timing = RunTiming()
    with patch.object(preprocessing, "get_cdsapi_data", side_effect=fail), pytest.raises(RuntimeError):
        preprocessing.get_data(args, LOGGER, timing=timing)
    data = timing.to_dict()
    assert data["input_data"] is None  # attempted, never completed: not "downloaded"
    assert data["phases"]["acquisition"] == {"seconds": 7.0, "completed": False}


NAMELIST = str(REPO_ROOT / "inputs" / "namelist_ERA5")


def _record(started_utc, transfer=12.5, assembly=1.25):
    return {
        "timing": {"started_utc": f"{started_utc:%Y-%m-%dT%H:%M:%SZ}", "transfer_seconds": transfer,
                   "assembly_seconds": assembly, "workers": 4},
        "totals": {"downloaded_responses": 3, "reused_responses": 1, "downloaded_bytes": 1000},
    }


def _write_ncar(path, record):
    _netcdf(path)
    with open(f"{path}.provenance.json", "w", encoding="utf-8") as handle:
        json.dump(record, handle)


def _get_ncar(args, timing):
    # the real get_ncar_data; only its retrieval/lock/verification are replaced
    return preprocessing.get_data(args, LOGGER, varlist=NAMELIST, timing=timing)


def test_ncar_new_download_copies_provenance_times(tmp_path, clock):
    infile = tmp_path / "ncar.nc"
    args = _args(tmp_path, infile, ncar=True)

    def retrieve(plan, path, names, **_):
        clock.now += 40.0
        record = _record(datetime.datetime.now(datetime.timezone.utc))
        _write_ncar(path, record)
        return record

    timing = RunTiming()
    with patch.object(ncar_gdex, "_retrieve", side_effect=retrieve) as retrieved:
        _get_ncar(args, timing)
    retrieved.assert_called_once()
    data = timing.to_dict()
    assert (data["provider"], data["input_data"]) == ("ncar", "downloaded")
    assert data["phases"]["acquisition"] == {"seconds": 40.0, "completed": True}
    ncar = data["ncar_provenance"]
    assert ncar["transfer_seconds"] == 12.5 and ncar["assembly_seconds"] == 1.25
    assert (ncar["downloaded_responses"], ncar["reused_responses"]) == (3, 1)
    assert "NCAR transfer 12.50 s" in timing.summary()


def test_ncar_reused_file(tmp_path, clock):
    infile = tmp_path / "ncar.nc"
    _write_ncar(infile, _record(datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(days=3)))
    args = _args(tmp_path, infile, ncar=True)

    def verify(*_, **__):
        clock.now += 2.0  # SHA-256 check of the existing file
        return {}

    timing = RunTiming()
    with patch.object(ncar_gdex, "check_existing", side_effect=verify), \
            patch.object(ncar_gdex, "_retrieve", side_effect=AssertionError("must not retrieve")):
        _get_ncar(args, timing)
    data = timing.to_dict()
    assert (data["provider"], data["input_data"]) == ("ncar", "reused")
    assert data["phases"]["acquisition"] == {"seconds": 2.0, "completed": True}
    assert data["ncar_provenance"] is None  # earlier transfer times are not this run's


def test_ncar_file_written_by_another_run_in_the_same_second_is_reused(tmp_path):
    """
    Race of get_ncar_data: the file is absent when this run starts; while it
    waits for the lock another run writes it, with a provenance record
    started in the same second. This run only verifies the file: it must be
    reported as reused, without the other run's transfer/assembly times.
    """
    infile = tmp_path / "ncar.nc"
    args = _args(tmp_path, infile, ncar=True)
    same_second = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0)

    @contextlib.contextmanager
    def lock_released_by_other_run(path):
        _write_ncar(path, _record(same_second))
        yield

    timing = RunTiming()
    with patch.object(ncar_gdex, "_exclusive", side_effect=lock_released_by_other_run), \
            patch.object(ncar_gdex, "check_existing", return_value={}) as verified, \
            patch.object(ncar_gdex, "_retrieve", side_effect=AssertionError("must not retrieve")):
        _get_ncar(args, timing)
    verified.assert_called_once()
    data = timing.to_dict()
    assert (data["provider"], data["input_data"]) == ("ncar", "reused")
    assert data["ncar_provenance"] is None
    assert "NCAR transfer" not in timing.summary()


def test_ncar_failure(tmp_path):
    args = _args(tmp_path, tmp_path / "ncar.nc", ncar=True)
    timing = RunTiming()
    with patch.object(ncar_gdex, "_retrieve", side_effect=OSError("network down")), pytest.raises(OSError):
        _get_ncar(args, timing)
    data = timing.to_dict()
    assert data["provider"] == "ncar" and data["input_data"] is None
    assert data["phases"]["acquisition"]["completed"] is False
    assert data["ncar_provenance"] is None


def test_ncar_reuse_refused_is_not_reported_as_reused(tmp_path):
    infile = _netcdf(tmp_path / "ncar.nc")  # no provenance record: get_ncar_data refuses it
    args = _args(tmp_path, infile, ncar=True)
    timing = RunTiming()
    with pytest.raises(ncar_gdex.NcarReuseError):
        _get_ncar(args, timing)
    assert timing.to_dict()["input_data"] is None


# ------------------------------------------------------- run_status.json
def test_failed_run_writes_timing(tmp_path, clock):
    args = argparse.Namespace(infile="x.nc", fixed=True, plots=False)
    timing = RunTiming()
    status = RunStatus(args, tmp_path, argv=["lorenzcycletoolkit.py"], timing=timing)
    with pytest.raises(RuntimeError), status.guard(), timing.phase("preparation"):
        with timing.phase("acquisition"):
            clock.now += 3.0
            raise RuntimeError("download failed")
    data = json.loads((tmp_path / "run_status.json").read_text())["timing"]
    assert data["phases"]["acquisition"] == {"seconds": 3.0, "completed": False}
    assert data["phases"]["preparation"] == {"seconds": 0.0, "completed": False}
    assert data["phases"]["lec"]["seconds"] is None
    assert data["total_seconds"] == 3.0
    assert data["plots_requested"] is False


def test_cli_local_run_reports_timing(lectk_workdir):
    lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    run = cli("samples/testdata_ERA5.nc", "-f", "-r")
    assert run.returncode == 0, run.stderr[-2000:]
    status = oc.check_run_status(FIXED_RESULTS, "success", 0, oc.products("fixed", figures="not_requested"))
    timing = status["timing"]
    assert (timing["provider"], timing["input_data"]) == ("local", "local")
    phases = timing["phases"]
    assert phases["acquisition"] == {"seconds": None, "completed": None}
    for name in ("preparation", "lec"):
        assert phases[name]["completed"] is True and phases[name]["seconds"] > 0
    assert timing["total_seconds"] >= phases["preparation"]["seconds"] + phases["lec"]["seconds"]
    log = (FIXED_RESULTS / "log.testdata_ERA5").read_text(encoding="utf-8")
    assert "⏱️ Timing (local, local data): acquisition n/a, preparation" in log
    assert BLOCKED not in run.stderr


def _check_track_timing(status, complete=True):
    timing = status["timing"]
    assert (timing["provider"], timing["input_data"]) == ("local", "local")
    assert status["command"]["method"] == "track"
    return timing["phases"]


def test_cli_track_run_reports_timing(lectk_workdir):
    """Moving framework (--track), the mode of the D-Bench cases; offline, no plots."""
    lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    run = cli("samples/testdata_ERA5.nc", "-t", "-r")
    assert run.returncode == 0, run.stderr[-2000:]
    status = oc.check_run_status(TRACK_RESULTS, "success", 0, oc.products("track", figures="not_requested"))
    phases = _check_track_timing(status)
    assert phases["acquisition"] == {"seconds": None, "completed": None}
    for name in ("preparation", "lec"):
        assert phases[name]["completed"] is True and phases[name]["seconds"] > 0
    assert status["timing"]["total_seconds"] >= phases["preparation"]["seconds"] + phases["lec"]["seconds"]
    log = (TRACK_RESULTS / "log.testdata_ERA5").read_text(encoding="utf-8")
    assert "⏱️ Timing (local, local data): acquisition n/a, preparation" in log
    assert BLOCKED not in run.stderr


def test_cli_track_failure_keeps_timing(lectk_workdir):
    """A --track run whose track lies outside the data fails in the LEC phase; its time is kept."""
    lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    lines = Path("inputs/track").read_text().splitlines()
    shifted = [lines[0]] + [f"{t};-80.0;100.0" for t in (line.split(";")[0] for line in lines[1:])]
    Path("inputs/track_outside").write_text("\n".join(shifted) + "\n")
    run = cli("samples/testdata_ERA5.nc", "-t", "-r", "--trackfile", "inputs/track_outside")
    assert run.returncode == 1
    status = oc.load_run_status(TRACK_RESULTS)
    assert status["status"] == "failed"
    phases = _check_track_timing(status)
    assert phases["preparation"]["completed"] is True
    assert phases["lec"]["completed"] is False and phases["lec"]["seconds"] is not None
    assert status["timing"]["total_seconds"] is not None
    assert BLOCKED not in run.stderr

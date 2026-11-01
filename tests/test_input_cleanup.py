"""
Input NetCDF files are kept after a run; ``--delete-downloaded-input``
deletes only a file this run downloaded, after a successful run, unchanged.

Unit tests of src/utils/input_cleanup.py, in-process runs of the command line
with the downloads replaced (CDS and NCAR; nothing is downloaded), and
command-line runs in a subprocess with the network blocked, a private HOME
and no CDS credential variables, so that a missing input never reaches the
CDS (regression of the deleted-input incident of the moving framework).
"""

import json
import os
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest
import xarray as xr

import lorenzcycletoolkit
from src.utils import preprocessing
from src.utils.input_cleanup import cleanup_input, downloaded_input, fingerprint
from src.utils.run_status import RunStatus
from src.utils.run_timing import RunTiming

REPO_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = REPO_ROOT / "tests" / "run_cli_offline.py"
TRACK_RESULTS = Path("LEC_Results/case_track")
CDS_NAMES = {"T": "t", "Z": "z", "W": "w", "U": "u", "V": "v", "time": "valid_time", "level": "pressure_level"}


def _cds_format(path):
    """The ERA5 sample with the variable and coordinate names of a CDS download."""
    with xr.open_dataset(REPO_ROOT / "samples" / "testdata_ERA5.nc") as ds:
        ds.rename(CDS_NAMES).to_netcdf(path)
    return str(path)


# ------------------------------------------------------------------- units
class Status:
    def __init__(self, status="success"):
        self.status = status
        self.input_cleanup = None
        self.writes = 0

    def write(self):
        self.writes += 1


def _args(path, **flags):
    import argparse

    return argparse.Namespace(**dict(dict(infile=str(path), cdsapi=True, ncar=False, delete_downloaded_input=False),
                                     **flags))


def _downloaded_timing():
    timing = RunTiming()
    timing.set_source("cds", "downloaded")
    return timing


def test_kept_without_flag(tmp_path):
    path = tmp_path / "d.nc"
    path.write_text("x")
    status = Status()
    record = cleanup_input(_args(path), status, fingerprint(path))
    assert path.exists() and record["deleted"] is False and "never deleted" in record["reason"]
    assert status.writes == 1 and status.input_cleanup == record


def test_deleted_with_flag_when_downloaded_successful_and_unchanged(tmp_path):
    path = tmp_path / "d.nc"
    path.write_text("x")
    record = cleanup_input(_args(path, delete_downloaded_input=True), Status(), fingerprint(path))
    assert not path.exists() and record["deleted"] is True


@pytest.mark.parametrize("state", ["partial", "failed", "running"])
def test_kept_when_run_not_successful(tmp_path, state):
    path = tmp_path / "d.nc"
    path.write_text("x")
    record = cleanup_input(_args(path, delete_downloaded_input=True), Status(state), fingerprint(path))
    assert path.exists() and state in record["reason"]


def test_kept_when_not_downloaded_by_this_run(tmp_path):
    path = tmp_path / "d.nc"
    path.write_text("x")
    record = cleanup_input(_args(path, delete_downloaded_input=True), Status(), None)
    assert path.exists() and "not downloaded by this run" in record["reason"]


def test_kept_when_changed_after_download(tmp_path):
    path = tmp_path / "d.nc"
    path.write_text("x")
    before = fingerprint(path)
    path.write_text("changed by someone else")
    record = cleanup_input(_args(path, delete_downloaded_input=True), Status(), before)
    assert path.exists() and "changed" in record["reason"]


def test_deletion_error_keeps_file_and_does_not_fail(tmp_path):
    path = tmp_path / "d.nc"
    path.write_text("x")
    with patch("src.utils.input_cleanup.os.remove", side_effect=PermissionError("denied")):
        record = cleanup_input(_args(path, delete_downloaded_input=True), Status(), fingerprint(path))
    assert path.exists() and record["deleted"] is False and "PermissionError" in record["reason"]


def test_downloaded_input_rules(tmp_path):
    path = tmp_path / "d.nc"
    path.write_text("x")
    assert downloaded_input(_args(path), False, _downloaded_timing()) == fingerprint(path)
    assert downloaded_input(_args(path), True, _downloaded_timing()) is None  # existed before the run
    assert downloaded_input(_args(path, cdsapi=False), False, _downloaded_timing()) is None  # local file
    reused = RunTiming()
    reused.set_source("ncar", "reused")
    assert downloaded_input(_args(path, cdsapi=False, ncar=True), False, reused) is None
    unknown = RunTiming()
    unknown.set_source("ncar", None)
    assert downloaded_input(_args(path, cdsapi=False, ncar=True), False, unknown) is None


def test_run_status_reports_input_cleanup(tmp_path):
    import argparse

    status = RunStatus(argparse.Namespace(infile="x.nc", track=True), tmp_path, argv=["lectk"])
    assert status.to_dict()["input_cleanup"] is None
    status.input_cleanup = {"requested": False, "path": "x.nc", "deleted": False, "reason": "kept"}
    assert status.to_dict()["input_cleanup"]["reason"] == "kept"


# ----------------------------------------------------- in-process command
@pytest.fixture
def cds_workdir(lectk_workdir):
    workdir = lectk_workdir(namelist="namelist_ERA5-cdsapi", track="track_testdata_ERA5")
    (workdir / "inputs" / "namelist_ERA5-cdsapi").write_bytes((workdir / "inputs" / "namelist").read_bytes())
    return workdir


def _main(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["lorenzcycletoolkit.py", "case.nc", "-t", "-r", "--trackfile", "inputs/track",
                                      "-o", "case", *argv])
    return lorenzcycletoolkit.main()


def _status():
    return json.loads((TRACK_RESULTS / "run_status.json").read_text())


def _fake_cds_download(args, track, logger):
    return _cds_format(args.infile)


def test_cds_download_is_kept_by_default(cds_workdir, monkeypatch):
    with patch.object(preprocessing, "get_cdsapi_data", side_effect=_fake_cds_download) as download:
        assert _main(monkeypatch, "--cdsapi") == 0
    download.assert_called_once()
    assert (cds_workdir / "case.nc").exists()  # formerly deleted by the moving framework
    status = _status()
    assert status["timing"]["input_data"] == "downloaded"
    assert status["input_cleanup"]["deleted"] is False


def test_cds_download_deleted_with_flag_after_success(cds_workdir, monkeypatch):
    with patch.object(preprocessing, "get_cdsapi_data", side_effect=_fake_cds_download):
        assert _main(monkeypatch, "--cdsapi", "--delete-downloaded-input") == 0
    assert not (cds_workdir / "case.nc").exists()
    assert _status()["input_cleanup"]["deleted"] is True
    assert (TRACK_RESULTS / "case_track_results.csv").exists()


def test_reused_cds_file_is_kept_even_with_flag(cds_workdir, monkeypatch):
    _cds_format(cds_workdir / "case.nc")
    with patch.object(preprocessing, "get_cdsapi_data") as download:
        assert _main(monkeypatch, "--cdsapi", "--delete-downloaded-input") == 0
    download.assert_not_called()
    assert (cds_workdir / "case.nc").exists()
    assert "not downloaded by this run" in _status()["input_cleanup"]["reason"]


def test_failed_run_keeps_downloaded_file_even_with_flag(cds_workdir, monkeypatch):
    with patch.object(preprocessing, "get_cdsapi_data", side_effect=_fake_cds_download), \
            patch.object(lorenzcycletoolkit, "lec_moving", side_effect=RuntimeError("computation failed")), \
            pytest.raises(RuntimeError):
        _main(monkeypatch, "--cdsapi", "--delete-downloaded-input")
    assert (cds_workdir / "case.nc").exists()
    status = _status()
    assert status["status"] == "failed" and status["input_cleanup"] is None


def test_ncar_download_deleted_with_flag_keeps_provenance(cds_workdir, monkeypatch):
    def fake_ncar(args, track, logger, namelist, on_retrieved):
        _cds_format(args.infile)
        Path(f"{args.infile}.provenance.json").write_text("{}")
        on_retrieved({"timing": {}, "totals": {}})
        return args.infile

    with patch.object(preprocessing, "get_ncar_data", side_effect=fake_ncar):
        assert _main(monkeypatch, "--ncar", "--delete-downloaded-input") == 0
    assert not (cds_workdir / "case.nc").exists()
    assert (cds_workdir / "case.nc.provenance.json").exists()


def test_flag_needs_a_download_source(cds_workdir, monkeypatch, capsys):
    with pytest.raises(SystemExit) as stop:
        _main(monkeypatch, "--delete-downloaded-input")
    assert stop.value.code == 2 and "local input files are never deleted" in capsys.readouterr().err


# ------------------------------------------- subprocess, network blocked
def _offline_cli(workdir, *argv):
    env = {k: v for k, v in os.environ.items() if not k.startswith("CDSAPI_")}
    env.update(HOME=str(workdir / "home"), MPLBACKEND="Agg")  # no ~/.cdsapirc reachable
    (workdir / "home").mkdir(exist_ok=True)
    flags = ["-I"] if sys.flags.isolated else []
    return subprocess.run([sys.executable, *flags, "-B", str(LAUNCHER), str(REPO_ROOT / "lorenzcycletoolkit.py"),
                           "missing.nc", "-t", "-r", "--trackfile", "inputs/track", "-o", "case", *argv],
                          cwd=workdir, env=env, capture_output=True, text=True, timeout=600)


def test_missing_input_with_cdsapi_never_reaches_the_cds(cds_workdir):
    run = _offline_cli(cds_workdir, "--cdsapi")
    assert run.returncode == 1
    assert "Request ID is" not in run.stderr + run.stdout  # nothing was submitted
    status = _status()
    assert status["status"] == "failed" and status["timing"]["phases"]["acquisition"]["completed"] is False
    assert not (cds_workdir / "missing.nc").exists()


def test_guarded_mode_without_explicit_credential_makes_no_network_attempt(cds_workdir):
    ledger = cds_workdir / "private"
    ledger.mkdir(mode=0o700)
    run = _offline_cli(cds_workdir, "--cdsapi", "--cds-max-seconds", "60", "--cds-ledger", str(ledger / "l.jsonl"))
    assert run.returncode == 1
    assert "network access blocked" not in run.stderr  # refused before any connection
    assert _status()["error"]["type"] == "CdsGuardError"

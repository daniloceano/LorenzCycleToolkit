"""
-o/--outname and the protection of earlier results (Wave 3: E0-12, DEC-07, DEC-08).

--outname NAME replaces the input-file stem everywhere, as the documentation
says: LEC_Results/NAME_<method>/ with NAME_<method>_results.csv,
NAME_<method>_trackfile and log.NAME. A run never overwrites a results
directory that already holds files unless --overwrite is given; with it, files
are rewritten in place and nothing is deleted. Each test runs the real command
line in a subprocess with the network blocked.
"""

import hashlib
import json
from pathlib import Path

import pytest

from src.utils.output_names import validate_outname
from tests.test_cli_exit_status import cli


def _tree(root):
    """{relative path: (sha256, mtime_ns)} of every file under root."""
    root = Path(root)
    return {p.relative_to(root).as_posix(): (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
            for p in sorted(root.rglob("*")) if p.is_file()}


def test_outname_fixed(lectk_workdir):
    lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    run = cli("samples/testdata_ERA5.nc", "-f", "-r", "-o", "case_A")
    assert run.returncode == 0, run.stderr[-2000:]
    results = Path("LEC_Results/case_A_fixed")
    for name in ("case_A_fixed_results.csv", "log.case_A", "run_status.json",
                 "lec_metadata.json", "domain_selection.csv"):
        assert (results / name).is_file(), name
    assert not Path("LEC_Results/testdata_ERA5_fixed").exists()
    assert sorted(p.name for p in Path("LEC_Results").iterdir()) == ["case_A_fixed"]
    status = json.loads((results / "run_status.json").read_text())
    assert status["log_file"] == "log.case_A"
    meta = json.loads((results / "lec_metadata.json").read_text())
    assert meta["main_results"]["file"] == "case_A_fixed_results.csv"
    assert meta["command"]["outname"] == "case_A"


def test_outname_track(lectk_workdir):
    lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    run = cli("samples/testdata_ERA5.nc", "-t", "-r", "--outname", "case_B")
    assert run.returncode == 0, run.stderr[-2000:]
    results = Path("LEC_Results/case_B_track")
    for name in ("case_B_track_results.csv", "case_B_track_trackfile", "log.case_B"):
        assert (results / name).is_file(), name
    assert sorted(p.name for p in Path("LEC_Results").iterdir()) == ["case_B_track"]


def test_default_names_are_unchanged(lectk_workdir):
    lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    run = cli("samples/testdata_ERA5.nc", "-f", "-r")
    assert run.returncode == 0, run.stderr[-2000:]
    results = Path("LEC_Results/testdata_ERA5_fixed")
    assert (results / "testdata_ERA5_fixed_results.csv").is_file()
    assert (results / "log.testdata_ERA5").is_file()


@pytest.mark.parametrize("name", ["", "   ", ".", "..", "a/b", "../up", "a\\b"])
def test_invalid_outname_exits_2_without_writing(lectk_workdir, name):
    lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    run = cli("samples/testdata_ERA5.nc", "-f", "-r", "-o", name)
    assert run.returncode == 2, (name, run.returncode, run.stderr[-500:])
    assert "outname" in run.stderr
    assert not Path("LEC_Results").exists()


@pytest.mark.parametrize("name", ["case", "case.v2", "Catarina_2004-03", "ação"])
def test_valid_outnames(name):
    assert validate_outname(name) == name


def test_existing_results_are_protected(lectk_workdir):
    """C1: a second run in the same directory stops before writing anything."""
    lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    first = cli("samples/testdata_ERA5.nc", "-f", "-r")
    assert first.returncode == 0, first.stderr[-2000:]
    results = Path("LEC_Results/testdata_ERA5_fixed")
    (results / "my_notes.txt").write_text("user file\n")
    before = _tree(results)
    second = cli("samples/testdata_ERA5.nc", "-f", "-r")
    assert second.returncode == 2
    assert "--overwrite" in second.stderr and "already contains" in second.stderr
    assert _tree(results) == before


def test_overwrite_rewrites_in_place_and_deletes_nothing(lectk_workdir):
    """C1 with --overwrite: files rewritten, older products kept and reported, user file intact."""
    lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    first = cli("samples/testdata_ERA5.nc", "-f", "-r")
    assert first.returncode == 0, first.stderr[-2000:]
    results = Path("LEC_Results/testdata_ERA5_fixed")
    (results / "my_notes.txt").write_text("user file\n")
    stale = results / "Figures" / "timeseries" / "old_figure.png"
    stale.parent.mkdir(parents=True)
    stale.write_bytes(b"old")
    before = _tree(results)
    second = cli("samples/testdata_ERA5.nc", "-f", "-r", "--overwrite")
    assert second.returncode == 0, second.stderr[-2000:]
    after = _tree(results)
    assert set(before) <= set(after), f"files deleted: {sorted(set(before) - set(after))}"
    assert after["my_notes.txt"] == before["my_notes.txt"]
    assert after[stale.relative_to(results).as_posix()] == before[stale.relative_to(results).as_posix()]
    status = json.loads((results / "run_status.json").read_text())
    previous = status["previous_results"]
    assert previous["detected"] is True and previous["overwrite"] is True
    assert previous["n_files_before"] == len(before)
    assert sorted(previous["stale_files"]) == ["Figures/timeseries/old_figure.png", "my_notes.txt"]
    log = (results / "log.testdata_ERA5").read_text()
    assert "--overwrite: reusing" in log


def test_empty_results_directory_is_not_an_earlier_run(lectk_workdir):
    lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    Path("LEC_Results/testdata_ERA5_fixed/Figures").mkdir(parents=True)
    run = cli("samples/testdata_ERA5.nc", "-f", "-r")
    assert run.returncode == 0, run.stderr[-2000:]
    status = json.loads(Path("LEC_Results/testdata_ERA5_fixed/run_status.json").read_text())
    assert status["previous_results"] is None

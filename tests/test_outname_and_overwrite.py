"""
-o/--outname (Wave 3: E0-12, DEC-07).

--outname NAME replaces the input-file stem everywhere, as the documentation
says: LEC_Results/NAME_<method>/ with NAME_<method>_results.csv,
NAME_<method>_trackfile and log.NAME. Each test runs the real command line in
a subprocess with the network blocked.
"""

import json
from pathlib import Path

import pytest

from src.utils.output_names import validate_outname
from tests.test_cli_exit_status import cli


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

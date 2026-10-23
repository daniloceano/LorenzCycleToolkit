"""
Negative controls for tests/output_checks.py.

The reference outputs are copied to a temporary directory and deliberately
damaged (a file removed, a time step dropped, a term renamed, a value made
non-finite or shifted beyond the registered tolerance). Each damaged output
must be rejected, and the intact reference must be accepted, so the
integration tests cannot pass on missing or corrupted results.
"""

import csv
import shutil

import pytest

from tests import output_checks as oc

FIXED = oc.REFERENCE_DIR / "ERA5_fixed"
TRACK = oc.REFERENCE_DIR / "ERA5_track"
FIXED_STEM = "testdata_ERA5_fixed"
TRACK_STEM = "testdata_ERA5_track"
TIMES = [f"2005-08-09 0{h}:00:00" for h in range(5)]


@pytest.fixture
def fixed_copy(tmp_path):
    target = tmp_path / "ERA5_fixed"
    shutil.copytree(FIXED, target)
    return target


@pytest.fixture
def track_copy(tmp_path):
    target = tmp_path / "ERA5_track"
    shutil.copytree(TRACK, target)
    return target


def _edit_cell(path, row, column, value, delimiter=","):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter=delimiter))
    rows[row][column] = value
    with open(path, "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter=delimiter, lineterminator="\n").writerows(rows)


def _drop_last_row(path):
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    path.write_text("".join(lines[:-1]), encoding="utf-8")


def test_intact_reference_is_accepted(fixed_copy, track_copy):
    oc.check_main_results(fixed_copy, FIXED_STEM, TIMES)
    oc.check_vertical_levels(fixed_copy, "level", "time", TIMES, "fixed")
    oc.assert_matches_reference(fixed_copy, FIXED, FIXED_STEM, "level")
    oc.check_main_results(track_copy, TRACK_STEM, TIMES)
    oc.check_trackfile(track_copy, TRACK_STEM, TIMES)
    oc.check_vertical_levels(track_copy, "level", "time", TIMES, "track")
    oc.assert_matches_reference(track_copy, TRACK, TRACK_STEM, "level")


def test_missing_main_results_is_detected(fixed_copy):
    (fixed_copy / f"{FIXED_STEM}_results.csv").unlink()
    with pytest.raises(AssertionError, match="main results file is missing"):
        oc.check_main_results(fixed_copy, FIXED_STEM, TIMES)


def test_missing_vertical_file_is_detected(fixed_copy):
    (fixed_copy / "results_vertical_levels" / "Ck_3_level.csv").unlink()
    with pytest.raises(AssertionError, match="per-level file is missing: Ck_3_level.csv"):
        oc.check_vertical_levels(fixed_copy, "level", "time", TIMES, "fixed")


def test_missing_trackfile_is_detected(track_copy):
    (track_copy / f"{TRACK_STEM}_trackfile").unlink()
    with pytest.raises(AssertionError, match="trackfile is missing"):
        oc.check_trackfile(track_copy, TRACK_STEM, TIMES)


def test_missing_time_step_is_detected(fixed_copy):
    _drop_last_row(fixed_copy / f"{FIXED_STEM}_results.csv")
    with pytest.raises(AssertionError, match="4 time steps, expected 5"):
        oc.check_main_results(fixed_copy, FIXED_STEM, TIMES)


def test_missing_time_step_in_vertical_file_is_detected(track_copy):
    _drop_last_row(track_copy / "results_vertical_levels" / "Kz_level.csv")
    with pytest.raises(AssertionError, match="Kz_level.csv: 4 rows"):
        oc.check_vertical_levels(track_copy, "level", "time", TIMES, "track")


def test_missing_term_is_detected(fixed_copy):
    _edit_cell(fixed_copy / f"{FIXED_STEM}_results.csv", 0, 3, "Kz_renamed")
    with pytest.raises(AssertionError, match=r"missing \['Kz'\]"):
        oc.check_main_results(fixed_copy, FIXED_STEM, TIMES)


@pytest.mark.parametrize("bad_value", ["nan", "inf", "-inf", ""])
def test_non_finite_value_is_detected(fixed_copy, bad_value):
    _edit_cell(fixed_copy / f"{FIXED_STEM}_results.csv", 3, 1, bad_value)
    with pytest.raises(AssertionError, match="not (finite|numeric)"):
        oc.check_main_results(fixed_copy, FIXED_STEM, TIMES)


def test_non_finite_vertical_value_is_detected(track_copy):
    _edit_cell(track_copy / "results_vertical_levels" / "Ae_level.csv", 2, 2, "nan")
    with pytest.raises(AssertionError, match="not finite"):
        oc.check_vertical_levels(track_copy, "level", "time", TIMES, "track")


def test_value_beyond_tolerance_is_detected(fixed_copy):
    # Kz at the first time step multiplied by (1 + 1e-6): about 0.058 J m-2,
    # far above the registered tolerance (1e-6 J m-2 + 1e-10 relative).
    path = fixed_copy / f"{FIXED_STEM}_results.csv"
    with open(path, newline="", encoding="utf-8") as handle:
        original = float(list(csv.reader(handle))[1][3])
    _edit_cell(path, 1, 3, repr(original * (1 + 1e-6)))
    with pytest.raises(AssertionError, match="values differ from the reference"):
        oc.assert_matches_reference(fixed_copy, FIXED, FIXED_STEM, "level")


def test_vertical_value_beyond_tolerance_is_detected(track_copy):
    # The per-level Cz value here is about -9.1e-5 W m-2 Pa-1; a relative change
    # of 1e-6 (about 9e-11) is well above the registered 1e-12 absolute tolerance.
    path = track_copy / "results_vertical_levels" / "Cz_level.csv"
    with open(path, newline="", encoding="utf-8") as handle:
        original = float(list(csv.reader(handle))[2][3])
    _edit_cell(path, 2, 3, repr(original * (1 + 1e-6)))
    with pytest.raises(AssertionError, match="Cz_level.csv: 1 values differ"):
        oc.assert_matches_reference(track_copy, TRACK, TRACK_STEM, "level")


def test_trackfile_coordinate_change_is_detected(track_copy):
    path = track_copy / f"{TRACK_STEM}_trackfile"
    _edit_cell(path, 1, 1, "-22.75", delimiter=";")
    with pytest.raises(AssertionError, match="values differ from the reference"):
        oc.assert_matches_reference(track_copy, TRACK, TRACK_STEM, "level")


def test_value_within_registered_tolerance_is_accepted(fixed_copy):
    # A change of 1e-7 J m-2 in Kz is inside the registered tolerance.
    path = fixed_copy / f"{FIXED_STEM}_results.csv"
    with open(path, newline="", encoding="utf-8") as handle:
        original = float(list(csv.reader(handle))[1][3])
    _edit_cell(path, 1, 3, repr(original + 1e-7))
    oc.assert_matches_reference(fixed_copy, FIXED, FIXED_STEM, "level")


def _write_log(directory, lines):
    path = directory / "log.testdata_ERA5"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


ARGS_LINE = "INFO - Command line arguments: Namespace(infile='{}', residuals=True)"


def test_isolated_log_is_accepted(tmp_path):
    _write_log(tmp_path, ["INFO - Starting LEC analysis", ARGS_LINE.format("samples/testdata_ERA5.nc")])
    oc.check_log(tmp_path, "testdata_ERA5")


def test_log_of_another_run_is_detected(tmp_path):
    _write_log(tmp_path, [
        "INFO - Starting LEC analysis", ARGS_LINE.format("samples/testdata_ERA5.nc"),
        "INFO - Starting LEC analysis", ARGS_LINE.format("samples/testdata_NCEP-R2.nc"),
    ])
    with pytest.raises(AssertionError, match="2 runs logged"):
        oc.check_log(tmp_path, "testdata_ERA5")


def test_missing_log_is_detected(tmp_path):
    with pytest.raises(AssertionError, match="log file is missing"):
        oc.check_log(tmp_path, "testdata_ERA5")

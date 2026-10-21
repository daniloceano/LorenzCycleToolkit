"""
Common vertical control volume of the moving framework (Wave 4: DEC-10 B, E0-21a).

Expectations from the pre-registration of the wave, written before the
implementation: a level incomplete inside the box of any time step is
excluded from every time step, so the run equals, byte for byte, the run of
the same input without that level (main results with tendencies and
residuals, per-level files, trackfile, domain selection); the per-level
files list the levels used and each value sits under its own level; an
internal gap or fewer than two levels stop the run explicitly. The input
NetCDF of the repository is never modified: each case writes its own copy.
"""

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from lorenzcycletoolkit import (create_arg_parser, initialize_logging,
                                prepare_data, run_lec_analysis,
                                setup_results_directory)
from src.utils.term_output import align_to_header
from src.utils.vertical_harmonization import VerticalHarmonizationError

LEVELS = [70000.0, 85000.0, 92500.0, 100000.0]
TIMES = ["2005-08-09 00:00:00", "2005-08-09 01:00:00", "2005-08-09 02:00:00",
         "2005-08-09 03:00:00", "2005-08-09 04:00:00"]


def _run(monkeypatch, infile, *flags):
    argv = ["lorenzcycletoolkit.py", infile, "-r", *flags]
    monkeypatch.setattr("sys.argv", argv)
    args = create_arg_parser().parse_args()
    method = "fixed" if args.fixed else ("track" if args.track else "choose")
    results, figures, levels = setup_results_directory(args, method)
    logger = initialize_logging(results, args)
    logger.info(f"Command line arguments: {args}")
    data = prepare_data(args, "inputs/namelist", logger)
    run_lec_analysis(data, args, results, figures, levels, logger)
    return Path(results)


def _sample():
    with xr.open_dataset("samples/testdata_ERA5.nc") as source:
        data = source.load()
    for name in data.variables:
        data[name].encoding = {}
    return data


def _write(data, path):
    data.to_netcdf(path)
    return path.name


def _omega_missing(workdir, name, level_hpa, times, point=None):
    """Copy of the sample with omega missing at one level (whole field or one point) at some times."""
    data = _sample()
    where = {"level": level_hpa}
    if point is not None:
        lat, lon = point
        where.update(latitude=lat, longitude=lon % 360)
    for it in times:
        data["W"].loc[{**where, "time": data.time.values[it]}] = np.nan
    return _write(data, workdir / name)


def _without_level(workdir, name, level_hpa):
    data = _sample()
    return _write(data.drop_sel(level=level_hpa), workdir / name)


def _outputs(results):
    """Scientific outputs of a run (CSV files, trackfile), keyed by their name without the stem."""
    stem = results.name
    return {
        path.relative_to(results).as_posix().replace(stem, "STEM"): path.read_bytes()
        for path in sorted(results.rglob("*"))
        if path.is_file() and path.suffix != ".json" and not path.name.startswith("log.")
    }


def _rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle))


def _harmonization(results):
    return json.loads((results / "lec_metadata.json").read_text(encoding="utf-8"))["vertical"]


@pytest.mark.parametrize(
    "times", [[2], [1, 3], [0], [4], [0, 1, 2, 3, 4]],
    ids=["middle", "several", "start", "end", "all"],
)
def test_level_missing_at_some_times_equals_input_without_it(monkeypatch, lectk_workdir, times):
    """C-A, C-B: omega missing at 700 hPa at some times gives the run of the input without 700 hPa."""
    workdir = lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    case = _run(monkeypatch, _omega_missing(workdir, "case.nc", 700.0, times), "-t")
    reference = _run(monkeypatch, _without_level(workdir, "no700.nc", 700.0), "-t")

    new, ref = _outputs(case), _outputs(reference)
    assert sorted(new) == sorted(ref)
    assert len(new) == 26  # main results, trackfile, domain selection, 23 per-level files
    assert [name for name in new if new[name] != ref[name]] == []

    vertical = _harmonization(case)
    harmonization = vertical["harmonization"]
    assert harmonization["available_levels_Pa"] == LEVELS
    assert harmonization["used_levels_Pa"] == LEVELS[1:]
    assert harmonization["excluded_levels_Pa"] == [70000.0]
    assert harmonization["excluded_because_incomplete_at"] == [
        {"level_Pa": 70000.0, "times": [TIMES[i] for i in times]}
    ]
    assert harmonization["input_layer_Pa"] == {"top": 70000.0, "bottom": 100000.0}
    assert harmonization["effective_layer_Pa"] == {"top": 85000.0, "bottom": 100000.0}
    assert vertical["levels_used_vary_with_time"] is False
    assert vertical["per_level_file_header_levels_Pa"] == LEVELS[1:]

    log = next(case.glob("log.*")).read_text(encoding="utf-8")
    warnings = [line for line in log.splitlines() if " - WARNING - " in line and "Vertical harmonization" in line]
    assert len(warnings) == 1
    assert "70000 Pa (incomplete at " + ", ".join(TIMES[i] for i in times) + ")" in warnings[0]


def test_bottom_level_missing_thins_the_layer(monkeypatch, lectk_workdir):
    """C-A: omega missing at 1000 hPa at one time gives the run of the input without 1000 hPa."""
    workdir = lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    case = _run(monkeypatch, _omega_missing(workdir, "case.nc", 1000.0, [2]), "-t")
    reference = _run(monkeypatch, _without_level(workdir, "no1000.nc", 1000.0), "-t")
    new, ref = _outputs(case), _outputs(reference)
    assert sorted(new) == sorted(ref)
    assert [name for name in new if new[name] != ref[name]] == []
    harmonization = _harmonization(case)["harmonization"]
    assert harmonization["used_levels_Pa"] == LEVELS[:3]
    assert harmonization["effective_layer_Pa"] == {"top": 70000.0, "bottom": 92500.0}


def test_per_level_values_sit_under_their_own_level(monkeypatch, lectk_workdir):
    """C-C (E0-21a): the value of 850 hPa is written under 85000.0, never under another level."""
    workdir = lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    clean = _run(monkeypatch, "samples/testdata_ERA5.nc", "-t")
    case = _run(monkeypatch, _omega_missing(workdir, "case.nc", 700.0, [2]), "-t")

    for path in sorted((case / "results_vertical_levels").glob("*.csv")):
        rows = _rows(path)
        assert [float(v) for v in rows[0][1:]] == LEVELS[1:], path.name
        assert all(len(row) == 1 + len(LEVELS[1:]) for row in rows[1:]), path.name

    # Kz and Ke per level depend on the fields of that level only.
    for term in ("Kz", "Ke"):
        reference = _rows(clean / f"results_vertical_levels/{term}_level.csv")
        new = _rows(case / f"results_vertical_levels/{term}_level.csv")
        ref_header = [float(v) for v in reference[0][1:]]
        new_header = [float(v) for v in new[0][1:]]
        for ref_row, new_row in zip(reference[1:], new[1:]):
            assert ref_row[0] == new_row[0]
            for level, text in zip(new_header, new_row[1:]):
                value = float(text)
                own = float(ref_row[1 + ref_header.index(level)])
                assert abs(value - own) <= 1e-14 * abs(own), (term, new_row[0], level)
                for other in ref_header:
                    if other != level:
                        neighbour = float(ref_row[1 + ref_header.index(other)])
                        assert abs(value - neighbour) > 1e-3 * abs(neighbour), (term, new_row[0], level, other)


def test_internal_gap_stops_the_run(monkeypatch, lectk_workdir):
    """C-D: a level missing between the top and the bottom of the levels kept stops the run (DEC-10-L)."""
    workdir = lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    infile = _omega_missing(workdir, "case.nc", 850.0, [2])
    with pytest.raises(VerticalHarmonizationError, match="internal gap") as excinfo:
        _run(monkeypatch, infile, "-t")
    assert "85000 Pa" in str(excinfo.value) and "DEC-10-L" in str(excinfo.value)
    results = workdir / "LEC_Results" / "case_track"
    assert not list(results.glob("*_results.csv"))
    assert not list((results / "results_vertical_levels").glob("*.csv"))
    status = json.loads((results / "run_status.json").read_text(encoding="utf-8"))
    assert status["status"] != "success"


def test_fewer_than_two_levels_stop_the_run(monkeypatch, lectk_workdir):
    """C-D: when only one level is complete at every time the run stops."""
    workdir = lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    data = _sample()
    for level, it in ((700.0, 1), (850.0, 2), (925.0, 3)):
        data["W"].loc[{"level": level, "time": data.time.values[it]}] = np.nan
    infile = _write(data, workdir / "case.nc")
    with pytest.raises(VerticalHarmonizationError, match="at least 2") as excinfo:
        _run(monkeypatch, infile, "-t")
    assert "100000 Pa" in str(excinfo.value)
    assert not list((workdir / "LEC_Results" / "case_track").glob("*_results.csv"))


def _variable_track(workdir, large_at=3):
    """Boxes of 10 degrees centred at (-25, -45), 16 degrees at one time step."""
    lines = ["time;Lat;Lon;width;length"]
    for i in range(5):
        size = 16 if i == large_at else 10
        lines.append(f"2005-08-09-0{i}00;-25.0;-45.0;{size};{size}")
    (workdir / "inputs" / "track").write_text("\n".join(lines) + "\n")


POINT = (-32.0, -52.0)  # inside the 16-degree box only


def test_level_incomplete_only_inside_a_larger_box(monkeypatch, lectk_workdir):
    """C-A, C-B: the boxes of each time step decide; a gap outside the box of a time step does not count."""
    workdir = lectk_workdir(namelist="namelist_ERA5")
    _variable_track(workdir)
    case = _run(monkeypatch, _omega_missing(workdir, "case.nc", 700.0, range(5), point=POINT), "-t")
    reference = _run(monkeypatch, _without_level(workdir, "no700.nc", 700.0), "-t")
    new, ref = _outputs(case), _outputs(reference)
    assert [name for name in new if new[name] != ref[name]] == []
    harmonization = _harmonization(case)["harmonization"]
    assert harmonization["excluded_because_incomplete_at"] == [{"level_Pa": 70000.0, "times": [TIMES[3]]}]


def test_gap_outside_every_box_changes_nothing(monkeypatch, lectk_workdir):
    """C-E: missing values that no box contains at their time leave the run unchanged."""
    workdir = lectk_workdir(namelist="namelist_ERA5")
    _variable_track(workdir)
    case = _run(monkeypatch, _omega_missing(workdir, "case.nc", 700.0, [0, 1, 2, 4], point=POINT), "-t")
    clean = _run(monkeypatch, _write(_sample(), workdir / "clean.nc"), "-t")
    new, ref = _outputs(case), _outputs(clean)
    assert [name for name in new if new[name] != ref[name]] == []
    assert _harmonization(case)["harmonization"]["excluded_levels_Pa"] == []


def test_choose_draws_every_box_before_computing(monkeypatch, lectk_workdir):
    """--choose: all boxes are chosen before the first term is computed, and the levels are harmonised."""
    import src.frameworks.lec_moving_framework as moving
    import src.utils.select_area as select_area

    workdir = lectk_workdir(namelist="namelist_ERA5")
    infile = _omega_missing(workdir, "case.nc", 700.0, [2])
    events = []
    box = {"min_lon": -52.5, "max_lon": -37.5, "min_lat": -30.0, "max_lat": -15.0}
    monkeypatch.setattr(select_area, "initial_domain",
                        lambda zeta, lat, lon: {"min_lon": -65.0, "max_lon": -25.0, "min_lat": -45.0, "max_lat": -15.0})

    def draw(*args):
        events.append("draw")
        return dict(box)

    real_box = moving.BoxData

    def box_data(*args, **kwargs):
        events.append("box")
        return real_box(*args, **kwargs)

    monkeypatch.setattr(moving, "draw_box_map", draw)
    monkeypatch.setattr(moving, "BoxData", box_data)
    results = _run(monkeypatch, infile, "--choose")
    assert events == ["draw"] * 5 + ["box"] * 5
    harmonization = _harmonization(results)["harmonization"]
    assert harmonization["excluded_because_incomplete_at"] == [{"level_Pa": 70000.0, "times": [TIMES[2]]}]


def test_align_to_header_puts_each_value_under_its_level(tmp_path):
    """E0-21a at unit level: a profile without one level leaves that cell empty."""
    path = tmp_path / "Kz_level.csv"
    pd.DataFrame(columns=["time"] + LEVELS).to_csv(path, index=None)
    frame = pd.DataFrame([[2.0, 3.0, 4.0]], index=["t2"], columns=[85000.0, 92500.0, 100000.0])
    frame.to_csv(path, mode="a", header=None)
    align_to_header(frame, path, "time").to_csv(path, mode="a", header=None)
    rows = _rows(path)
    assert rows[1] == ["t2", "2.0", "3.0", "4.0"]          # what was written before: shifted
    assert rows[2] == ["t2", "", "2.0", "3.0", "4.0"]      # aligned: 850 hPa under 85000.0


def test_align_to_header_refuses_a_level_not_in_the_header(tmp_path):
    path = tmp_path / "Kz_level.csv"
    pd.DataFrame(columns=["time"] + LEVELS[1:]).to_csv(path, index=None)
    frame = pd.DataFrame([[1.0, 2.0]], index=["t0"], columns=[70000.0, 85000.0])
    with pytest.raises(ValueError, match="not in the file header"):
        align_to_header(frame, path, "time")


def test_align_to_header_leaves_files_without_a_header(tmp_path):
    """Files not created by a framework (absent, or rows appended only) are written as before."""
    frame = pd.DataFrame([[1.0, 2.0]], index=["t0"], columns=[70000.0, 85000.0])
    path = tmp_path / "Kz_level.csv"
    assert align_to_header(frame, path, "time") is frame
    frame.to_csv(path, mode="a", header=None)
    assert align_to_header(frame, path, "time") is frame

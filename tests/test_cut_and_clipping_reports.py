"""
Transparency of the cut of the levels above 10 hPa (DEC-09) and of a box
requested beyond the data (DEC-11).

Expectations from the pre-registration, written before the implementation:
nothing in the computation changes; the cut keeps 1000 Pa (10 hPa) and
discards lower pressures, as before; a box requested beyond the data is cut at
the data edge, as before, and the run goes on. The terminal and the log warn,
and lec_metadata.json and domain_selection.csv record, what was discarded or
cut, with the limits of the fields actually used.
"""

import csv
import json

import numpy as np
import pytest
import xarray as xr

from tests.test_vertical_harmonization import (LABELS, LEVELS, MODES, TIMES, _mode_workdir, _outputs, _run,
                                               _sample, _write)


def _with_upper_levels(levels_hpa):
    """Sample with extra (complete) levels above 700 hPa, copies of the 700 hPa fields."""
    data = _sample()
    top = data.sel(level=[700.0])
    extra = xr.concat([top.assign_coords(level=[p]) for p in levels_hpa], dim="level")
    merged = xr.concat([extra, data], dim="level", combine_attrs="override").sortby("level")
    merged["level"].attrs = data["level"].attrs
    return merged


def _log(results):
    return next(results.glob("log.*")).read_text(encoding="utf-8")


def _vertical(results):
    return json.loads((results / "lec_metadata.json").read_text(encoding="utf-8"))["vertical"]


def _horizontal(results):
    return json.loads((results / "lec_metadata.json").read_text(encoding="utf-8"))["horizontal"]


def _domain_rows(results):
    with open(results / "domain_selection.csv", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


# ------------------------------------------------------------------ DEC-09
def test_levels_above_10_hpa_are_reported_and_equal_the_input_without_them(monkeypatch, lectk_workdir):
    """D9-1: 5 and 7 hPa are discarded as before; the warning gives the levels, the rule and both grids."""
    workdir = _mode_workdir(lectk_workdir, monkeypatch, "fixed")
    case = _run(monkeypatch, _write(_with_upper_levels([5.0, 7.0]), workdir / "case.nc"), "-f")
    reference = _run(monkeypatch, _write(_sample(), workdir / "plain.nc"), "-f")
    new, ref = _outputs(case), _outputs(reference)
    assert sorted(new) == sorted(ref)
    assert [name for name in new if new[name] != ref[name]] == []

    warnings = [line for line in _log(case).splitlines() if "Vertical levels above 10 hPa discarded" in line]
    assert len(warnings) == 1 and " - WARNING - " in warnings[0]
    assert "2 of 6 input levels have a pressure lower than 1000 Pa (10 hPa)" in warnings[0]
    assert "1000 Pa itself is kept (DEC-09)" in warnings[0]
    assert "Discarded: 500 Pa, 700 Pa." in warnings[0]
    assert "Input grid: 6 levels, 500 Pa, 700 Pa, 70000 Pa, 85000 Pa, 92500 Pa, 100000 Pa" in warnings[0]
    assert "grid after preprocessing: 4 levels, 70000 Pa, 85000 Pa, 92500 Pa, 100000 Pa." in warnings[0]
    cut = _vertical(case)["top_cut"]
    assert cut["input_grid_before_cut_Pa"] == [500.0, 700.0] + LEVELS
    assert cut["discarded_levels_Pa"] == [500.0, 700.0]
    assert _vertical(case)["input_levels_Pa"] == LEVELS


def test_no_level_discarded_gives_no_warning(monkeypatch, lectk_workdir):
    """D9-2: the sample (700 to 1000 hPa) loses nothing and says nothing about the cut."""
    _mode_workdir(lectk_workdir, monkeypatch, "fixed")
    results = _run(monkeypatch, "samples/testdata_ERA5.nc", "-f")
    assert "above 10 hPa discarded" not in _log(results)
    cut = _vertical(results)["top_cut"]
    assert cut["input_grid_before_cut_Pa"] == LEVELS and cut["discarded_levels_Pa"] == []


def test_10_hpa_is_kept_and_a_lower_pressure_is_discarded(monkeypatch, lectk_workdir):
    """D9-3: 1000 Pa is kept (the cut is inclusive); 9.99 hPa = 999 Pa alone is discarded."""
    workdir = _mode_workdir(lectk_workdir, monkeypatch, "fixed")
    kept = _run(monkeypatch, _write(_with_upper_levels([10.0, 20.0]), workdir / "ten.nc"), "-f")
    assert "above 10 hPa discarded" not in _log(kept)
    assert _vertical(kept)["top_cut"]["discarded_levels_Pa"] == []
    assert _vertical(kept)["input_levels_Pa"][:2] == [1000.0, 2000.0]

    below = _run(monkeypatch, _write(_with_upper_levels([9.99, 10.0]), workdir / "below.nc"), "-f")
    cut = _vertical(below)["top_cut"]
    assert cut["discarded_levels_Pa"] == pytest.approx([999.0])
    assert _vertical(below)["input_levels_Pa"][0] == 1000.0
    assert "Discarded: 999 Pa." in _log(below)


def test_cut_is_the_same_with_levels_in_pa(monkeypatch, lectk_workdir):
    """D9-4: the same grid in Pa gives the same cut, the same warning and the same results."""
    workdir = _mode_workdir(lectk_workdir, monkeypatch, "fixed")
    hpa = _with_upper_levels([5.0])
    pa = hpa.assign_coords(level=hpa.level.values * 100.0)
    pa["level"].attrs = {**hpa["level"].attrs, "units": "Pa"}
    in_hpa = _run(monkeypatch, _write(hpa, workdir / "hpa.nc"), "-f")
    in_pa = _run(monkeypatch, _write(pa, workdir / "pa.nc"), "-f")
    assert _vertical(in_pa)["top_cut"] == _vertical(in_hpa)["top_cut"]
    assert _vertical(in_pa)["top_cut"]["discarded_levels_Pa"] == [500.0]
    new, ref = _outputs(in_pa), _outputs(in_hpa)
    assert [name for name in new if new[name] != ref[name]] == []
    line = [x for x in _log(in_pa).splitlines() if "above 10 hPa discarded" in x]
    assert len(line) == 1 and "Discarded: 500 Pa." in line[0]


# ------------------------------------------------------------------ DEC-11
def _box_limits(workdir, north):
    (workdir / "inputs" / "box_limits").write_text(f"min_lon;-60\nmax_lon;-30\nmin_lat;-42.5\nmax_lat;{north}\n")


def _track_at(workdir, lat):
    lines = ["time;Lat;Lon"] + [f"2005-08-09-0{i}00;{lat:.2f};-45.00" for i in range(5)]
    (workdir / "inputs" / "track").write_text("\n".join(lines) + "\n")


BEYOND = {"min_lon": -52.5, "max_lon": -37.5, "min_lat": -27.0, "max_lat": -12.0}


def _clipping_case(lectk_workdir, monkeypatch, mode, beyond):
    """Run with the northern edge requested beyond the data (-15 degrees) or inside it."""
    boxes = {t: BEYOND for t in TIMES} if mode == "choose" and beyond else None
    workdir = _mode_workdir(lectk_workdir, monkeypatch, mode, boxes=boxes)
    if mode == "fixed":
        _box_limits(workdir, -10.0 if beyond else -17.5)
    elif mode == "track":
        _track_at(workdir, -20.0 if beyond else -25.0)
    return workdir, _run(monkeypatch, "samples/testdata_ERA5.nc", MODES[mode])


@pytest.mark.parametrize("mode", MODES)
def test_box_beyond_the_data_is_cut_and_reported(monkeypatch, lectk_workdir, mode):
    """D11-2: the run goes on; each affected box is reported with the limits requested and used."""
    _, results = _clipping_case(lectk_workdir, monkeypatch, mode, beyond=True)
    rows = _domain_rows(results)
    assert rows and all(r["requested_outside_input"] == "north" for r in rows)
    assert all(float(r["selected_north"]) == -15.0 for r in rows)
    warnings = [line for line in _log(results).splitlines() if "Domain clipping" in line]
    assert len(warnings) == len(rows)
    for line, row in zip(warnings, rows):
        assert " - WARNING - " in line and f"Domain clipping ({LABELS[mode]}" in line
        assert "on the north edge(s) and is cut at the data edge (DEC-11)" in line
        assert f"lat {float(row['requested_south']):g} to {float(row['requested_north']):g}; used:" in line
        assert (f"used: lon {float(row['selected_west']):g} to {float(row['selected_east']):g}, lat "
                f"{float(row['selected_south']):g} to {float(row['selected_north']):g}") in line
        assert f"({row['n_lon']} x {row['n_lat']} grid points)" in line
    cut = _horizontal(results)["boxes_cut_at_input_edge"]
    assert cut == [{"time": r["time"], "edges": ["north"]} for r in rows]


@pytest.mark.parametrize("mode", MODES)
def test_box_inside_the_data_is_not_reported(monkeypatch, lectk_workdir, mode):
    """D11-1: a box inside the data gives no clipping warning and no record."""
    _, results = _clipping_case(lectk_workdir, monkeypatch, mode, beyond=False)
    assert "Domain clipping" not in _log(results)
    assert all(r["requested_outside_input"] == "none" for r in _domain_rows(results))
    assert _horizontal(results)["boxes_cut_at_input_edge"] == []


def test_clipped_box_computes_with_the_box_actually_used(monkeypatch, lectk_workdir):
    """D11-3: requesting north = -10 gives the results of requesting the data edge, -15, directly."""
    workdir = _mode_workdir(lectk_workdir, monkeypatch, "fixed")
    _box_limits(workdir, -10.0)
    beyond = _run(monkeypatch, _write(_sample(), workdir / "beyond.nc"), "-f")
    _box_limits(workdir, -15.0)
    edge = _run(monkeypatch, _write(_sample(), workdir / "edge.nc"), "-f")
    new, ref = _outputs(beyond), _outputs(edge)
    differ = [name for name in new if new[name] != ref[name]]
    assert differ == ["domain_selection.csv"]  # only the requested limits differ
    assert np.isclose(float(_domain_rows(beyond)[0]["selected_north"]), -15.0)

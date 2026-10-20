"""
lec_metadata.json and domain_selection.csv (Wave 3: E0-17, E0-18, E0-19, E0-22, G3).

The expectations below are written from the pre-registration of the wave, not
from the implementation: the levels and grid of the ERA5 sample, the units of
each column and the coordinates predicted for the three reference cases by the
complementary Phase 2D report (section C). Each case runs the real pipeline
without figures in a private working directory.
"""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from lorenzcycletoolkit import (create_arg_parser, initialize_logging,
                                prepare_data, run_lec_analysis,
                                setup_results_directory)
from src.utils.output_spec import VERTICAL_PROFILE_SPEC
from tests import output_checks as oc

LEVELS = [70000.0, 85000.0, 92500.0, 100000.0]
INPUT_EXTENT = {"west": -65.0, "east": -25.0, "south": -45.0, "north": -15.0}
ENERGY = {"Az", "Ae", "Kz", "Ke"}
TIMES = ["2005-08-09 00:00:00", "2005-08-09 01:00:00", "2005-08-09 02:00:00",
         "2005-08-09 03:00:00", "2005-08-09 04:00:00"]
MOVING_TRACK = [(-23.00, -45.00), (-23.25, -44.75), (-23.50, -44.50), (-23.75, -44.25), (-24.00, -44.00)]


def _run(monkeypatch, infile, *flags):
    argv = ["lorenzcycletoolkit.py", infile, "-r", *flags]
    monkeypatch.setattr("sys.argv", argv)
    args = create_arg_parser().parse_args()
    method = "fixed" if args.fixed else "track"
    results, figures, levels = setup_results_directory(args, method)
    logger = initialize_logging(results, args)
    logger.info("Starting LEC analysis")
    logger.info(f"Command line arguments: {args}")
    data = prepare_data(args, "inputs/namelist", logger)
    status = run_lec_analysis(data, args, results, figures, levels, logger)
    return Path(results), status


def _write_track(workdir, centres, name="track"):
    lines = ["time;Lat;Lon"] + [f"2005-08-09-0{i}00;{lat:.2f};{lon:.2f}" for i, (lat, lon) in enumerate(centres)]
    (workdir / "inputs" / name).write_text("\n".join(lines) + "\n")


def _domain_rows(results):
    with open(results / "domain_selection.csv", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _metadata(results):
    return json.loads((results / "lec_metadata.json").read_text(encoding="utf-8"))


def _check_box(row, west, east, south, north, n_lon, n_lat, at_edge=(), outside="none"):
    for edge, value in zip(("west", "east", "south", "north"), (west, east, south, north)):
        assert float(row[f"selected_{edge}"]) == value, (edge, row)
    assert (int(row["n_lon"]), int(row["n_lat"])) == (n_lon, n_lat)
    for edge in ("west", "east", "south", "north"):
        assert float(row[f"input_{edge}"]) == INPUT_EXTENT[edge]
        assert row[f"{edge}_at_input_edge"] == ("true" if edge in at_edge else "false"), (edge, row)
    assert row["requested_outside_input"] == outside
    assert row["fields_share_selection"] == "true"


def _check_metadata_common(results, method, layout_cz1):
    meta = _metadata(results)
    assert meta["schema_version"] == 1
    assert meta["command"]["method"] == method
    digest = hashlib.sha256(Path("samples/testdata_ERA5.nc").read_bytes()).hexdigest()
    assert meta["input"]["sha256"] == digest
    assert meta["time"] == {"n": 5, "first": TIMES[0], "last": TIMES[-1], "distinct_steps_s": [3600]}
    vertical = meta["vertical"]
    assert vertical["levels_used_Pa"] == LEVELS
    assert vertical["dropped_levels_Pa"] == []
    assert (vertical["p_top_Pa"], vertical["p_bottom_Pa"]) == (70000.0, 100000.0)
    assert vertical["levels_used_vary_with_time"] is False
    assert "retained pressure layer" in vertical["column_statement"]
    columns = {c["name"]: c["units"] for c in meta["main_results"]["columns"]}
    header = next(csv.reader(open(results / meta["main_results"]["file"], encoding="utf-8")))
    assert header[0] == "" and list(columns) == header[1:]
    for name, unit in columns.items():
        expected = "J/m**2" if name in ENERGY else "kg/m**2/s" if name == "M" else "W/m**2"
        assert unit == expected, (name, unit)
    files = {e["term"]: e for e in meta["vertical_level_files"]}
    assert len(files) == 23
    for term, entry in files.items():
        spec = VERTICAL_PROFILE_SPEC[term]
        assert entry["units"] == spec["units"] and entry["role"] == spec["role"], term
        assert entry["column_factor"] == spec.get("column_factor"), term
        expected_layout = layout_cz1 if term in ("Cz_1", "Ce_1") else "time_by_level"
        assert entry["layout"] == expected_layout, (term, entry["layout"])
    return meta


def test_fixed_metadata_and_domain(monkeypatch, lectk_workdir):
    lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    results, status = _run(monkeypatch, "samples/testdata_ERA5.nc", "-f")
    assert status.status == "success"
    _check_metadata_common(results, "fixed", "time_invariant_profile")
    rows = _domain_rows(results)
    assert [r["time"] for r in rows] == ["all"]
    _check_box(rows[0], -60.0, -30.0, -42.5, -17.5, 121, 101)
    assert all(float(rows[0][f"offset_{e}"]) == 0.0 for e in ("west", "east", "south", "north"))
    oc.check_run_status(results, "success", 0, oc.products("fixed", figures="not_requested"))


def test_track_metadata_and_domain(monkeypatch, lectk_workdir):
    lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    results, status = _run(monkeypatch, "samples/testdata_ERA5.nc", "-t")
    assert status.status == "success"
    meta = _check_metadata_common(results, "track", "time_by_level")
    assert [r["time"] for r in meta["vertical"]["by_time"]] == TIMES
    rows = _domain_rows(results)
    assert [r["time"] for r in rows] == TIMES
    for row in rows:
        # the northern edge of the box is the northern edge of the data
        _check_box(row, -52.5, -37.5, -30.0, -15.0, 61, 61, at_edge=("north",))
    oc.check_run_status(results, "success", 0, oc.products("track", figures="not_requested"))


def test_moving_track_domain(monkeypatch, lectk_workdir):
    workdir = lectk_workdir(namelist="namelist_ERA5")
    _write_track(workdir, MOVING_TRACK)
    results, _ = _run(monkeypatch, "samples/testdata_ERA5.nc", "-t")
    rows = _domain_rows(results)
    assert len(rows) == 5
    for row, (lat, lon) in zip(rows, MOVING_TRACK):
        _check_box(row, lon - 7.5, lon + 7.5, lat - 7.5, lat + 7.5, 61, 61)
    _check_box(rows[0], -52.5, -37.5, -30.5, -15.5, 61, 61)
    _check_box(rows[-1], -51.5, -36.5, -31.5, -16.5, 61, 61)


def test_box_beyond_the_data_is_recorded(monkeypatch, lectk_workdir):
    """D2 / E0-23: a box reaching beyond the data is cut at the data edge; the cut is recorded."""
    workdir = lectk_workdir(namelist="namelist_ERA5")
    _write_track(workdir, [(-20.0, -45.0)] * 5)
    results, _ = _run(monkeypatch, "samples/testdata_ERA5.nc", "-t")
    for row in _domain_rows(results):
        assert float(row["requested_north"]) == -12.5
        _check_box(row, -52.5, -37.5, -27.5, -15.0, 61, 51, at_edge=("north",), outside="north")
        assert float(row["offset_north"]) == -2.5


def test_level_dropped_at_one_time_is_recorded(monkeypatch, lectk_workdir):
    """M2: omega missing at 700 hPa at one time only drops 700 hPa at that time only."""
    workdir = lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    with xr.open_dataset("samples/testdata_ERA5.nc") as source:
        data = source.load()
    data["W"].loc[{"time": data.time.values[2], "level": 700.0}] = np.nan
    for name in data.variables:
        data[name].encoding = {}
    data.to_netcdf(workdir / "case.nc")
    results, _ = _run(monkeypatch, "case.nc", "-t")
    vertical = _metadata(results)["vertical"]
    assert vertical["levels_used_vary_with_time"] is True
    assert vertical["levels_used_Pa"] is None and vertical["p_top_Pa"] is None
    assert vertical["levels_used_at_every_time_Pa"] == LEVELS[1:]
    by_time = {r["time"]: r for r in vertical["by_time"]}
    for time in TIMES:
        dropped = [70000.0] if time == TIMES[2] else []
        assert by_time[time]["dropped_levels_Pa"] == dropped, time
        assert by_time[time]["levels_used_Pa"] == [p for p in LEVELS if p not in dropped], time
    # the per-level files keep the header of the input levels (E0-21 is examined in Wave 4)
    assert vertical["per_level_file_header_levels_Pa"] == LEVELS


@pytest.mark.parametrize("edge", ["west", "east", "south", "north"])
def test_record_box_offsets_and_edges(edge):
    """D2 at unit level: nearest-neighbour selection away from the request is reported per edge."""
    from types import SimpleNamespace

    from src.utils.run_metadata import RunMetadata

    lon = np.arange(-65.0, -24.9, 0.25)
    lat = np.arange(-45.0, -14.9, 0.25)
    field = xr.DataArray(np.zeros((lat.size, lon.size)), coords={"latitude": lat, "longitude": lon},
                         dims=("latitude", "longitude"))
    requested = {"west": -60.0, "east": -30.0, "south": -42.5, "north": -17.5}
    sel = dict(requested)
    if edge in ("west", "south"):
        requested[edge] = INPUT_EXTENT[edge] - 1.0
        sel[edge] = INPUT_EXTENT[edge]
    else:
        requested[edge] = INPUT_EXTENT[edge] + 1.0
        sel[edge] = INPUT_EXTENT[edge]
    box = field.sel(latitude=slice(sel["south"], sel["north"]), longitude=slice(sel["west"], sel["east"]))
    box_obj = SimpleNamespace(LonIndexer="longitude", LatIndexer="latitude", tair=box, u=box, v=box,
                              omega=box, geopt=box, valid_levels=LEVELS, dropped_levels=[])
    meta = RunMetadata(SimpleNamespace(), "fixed", None, dict(INPUT_EXTENT))
    meta.record_box(box_obj, requested, "all", LEVELS)
    row = meta.domain_rows[0]
    assert row[f"selected_{edge}"] == INPUT_EXTENT[edge]
    assert row[f"offset_{edge}"] == pytest.approx(1.0 if edge in ("west", "south") else -1.0)
    assert row[f"{edge}_at_input_edge"] is True
    assert row["requested_outside_input"] == edge
    others = [e for e in ("west", "east", "south", "north") if e != edge]
    assert all(row[f"offset_{e}"] == 0.0 and row[f"{e}_at_input_edge"] is False for e in others)

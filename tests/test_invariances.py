"""
Properties of the equations that the outputs must keep (Wave 6, extended validation).

These are consequences of the formulation, not physical validation:

* homogeneity in the winds: with u, v and omega multiplied by ``a`` (T and
  Phi unchanged) each term is multiplied by ``a**k``, with ``k`` the number
  of wind factors in its equation (``docs/source/math.rst``); for ``a = -1``
  and ``a = 2`` the products are exact in floating point;
* the fixed and the moving frameworks give the same values for the same box;
* each column value is the column factor times the trapezoidal integral of
  its per-level profile (``src/utils/output_spec.py``).

Gz, Ge and the residuals mix powers of the winds (the diabatic heating
contains dT/dt and the advection) and are not homogeneous.
"""

import csv
import json
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from metpy.constants import g

from lorenzcycletoolkit import (create_arg_parser, initialize_logging,
                                prepare_data, run_lec_analysis,
                                setup_results_directory)

POWER = {
    "Az": 0, "Ae": 0, "∂Az/∂t (finite diff.)": 0, "∂Ae/∂t (finite diff.)": 0,
    "Cz": 1, "Ce": 1, "Ca": 1, "C_overturning": 1, "BAz": 1, "BAe": 1, "BΦZ": 1, "BΦE": 1, "M": 1,
    "Kz": 2, "Ke": 2, "∂Kz/∂t (finite diff.)": 2, "∂Ke/∂t (finite diff.)": 2,
    "Ck": 3, "BKz": 3, "BKe": 3,
}
FACTORS = {"1": 1.0, "1/(2g)": 1.0 / (2 * float(g.m)), "1/g": 1.0 / float(g.m)}


def _run(monkeypatch, infile, *flags):
    argv = ["lorenzcycletoolkit.py", infile, "-r", *flags]
    monkeypatch.setattr("sys.argv", argv)
    args = create_arg_parser().parse_args()
    results, figures, levels = setup_results_directory(args, "fixed" if args.fixed else "track")
    logger = initialize_logging(results, args)
    logger.info("Starting LEC analysis")
    logger.info(f"Command line arguments: {args}")
    data = prepare_data(args, "inputs/namelist", logger)
    run_lec_analysis(data, args, results, figures, levels, logger)
    return Path(results)


def _with_winds_times(workdir, factor):
    with xr.open_dataset("samples/testdata_ERA5.nc") as source:
        data = source.load()
    for name in data.variables:
        data[name].encoding = {}
    for name in ("U", "V", "W"):
        data[name] = data[name] * factor
    name = f"winds_{factor:+g}.nc".replace("+", "p").replace("-", "m")
    data.to_netcdf(workdir / name)
    return name


def _table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    return rows[0], {row[0]: dict(zip(rows[0][1:], map(float, row[1:]))) for row in rows[1:]}


def _main(results):
    return _table(next(results.glob("*_results.csv")))


@pytest.mark.parametrize("factor", [-1.0, 2.0])
def test_each_term_scales_with_its_power_of_the_winds(monkeypatch, lectk_workdir, factor):
    workdir = lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    _, base = _main(_run(monkeypatch, _with_winds_times(workdir, 1.0), "-f"))
    _, scaled = _main(_run(monkeypatch, _with_winds_times(workdir, factor), "-f"))
    for time, row in base.items():
        for column, power in POWER.items():
            assert scaled[time][column] == factor ** power * row[column], (time, column)


def test_fixed_and_moving_frameworks_agree_for_the_same_box(monkeypatch, lectk_workdir):
    workdir = lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1")
    # box_limits_Reg1: -60 to -30 E, -42.5 to -17.5 N
    lines = ["time;Lat;Lon;width;length"] + [f"2005-08-09-0{i}00;-30.0;-45.0;30;25" for i in range(5)]
    (workdir / "inputs" / "track").write_text("\n".join(lines) + "\n")
    fixed = _run(monkeypatch, "samples/testdata_ERA5.nc", "-f")
    moving = _run(monkeypatch, "samples/testdata_ERA5.nc", "-t")
    header, fixed_rows = _main(fixed)
    _, moving_rows = _main(moving)
    assert moving_rows == fixed_rows
    for path in sorted((fixed / "results_vertical_levels").glob("*.csv")):
        rows = list(csv.reader(path.open(newline="", encoding="utf-8")))
        if len(rows) > 1 and rows[1][0] == "level":
            continue  # time-invariant profiles of the fixed framework (Cz_1, Ce_1)
        assert _table(path) == _table(moving / "results_vertical_levels" / path.name), path.name


@pytest.mark.parametrize("flag", ["-f", "-t"])
def test_column_is_the_factor_times_the_profile_integral(monkeypatch, lectk_workdir, flag):
    lectk_workdir(namelist="namelist_ERA5", box_limits="box_limits_Reg1", track="track_testdata_ERA5")
    results = _run(monkeypatch, "samples/testdata_ERA5.nc", flag)
    _, main = _main(results)
    metadata = json.loads((results / "lec_metadata.json").read_text(encoding="utf-8"))
    checked = 0
    for entry in metadata["vertical_level_files"]:
        if entry["role"] not in ("integrand", "scaled_integrand") or entry["layout"] != "time_by_level":
            continue
        header, profiles = _table(results / entry["file"])
        levels = np.array([float(p) for p in header[1:]])
        column = entry["column_term"]
        scale = max(abs(row[column]) for row in main.values())
        for time, profile in profiles.items():
            integral = FACTORS[entry["column_factor"]] * np.trapezoid([profile[h] for h in header[1:]], levels)
            assert abs(integral - main[time][column]) <= 1e-12 * scale, (entry["term"], time)
        checked += 1
    assert checked == 11  # 8 integrands, Kz and Ke (1/(2g)), Ck (1/g)

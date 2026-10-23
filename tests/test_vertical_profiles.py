#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests of the per-level output convention recorded in ``src/utils/output_spec.py``.

The units of every profile are checked with Pint at the moment the profile is
written (synthetic box), the profile-to-column factors and the combination
of the split terms are checked on the reference outputs, and the column
relation of the mass-continuity residual is checked against a closed form.
These tests read the convention; none of them changes how a term is computed.
"""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pint
import pytest
import xarray as xr
from metpy.constants import g
from metpy.units import units

from src.analysis.conversion_terms import ConversionTerms
from src.analysis.energy_contents import EnergyContents
from src.analysis.generation_and_dissipation_terms import GenerationDissipationTerms
from src.analysis.mass_continuity import MassContinuity
from src.utils import term_output
from src.utils.calc_budget_and_residual import calc_budget_diff
from src.utils.output_spec import (PROVISIONAL_PROFILE_SPEC, SPLIT_TERM_RELATIONS,
                                   TENDENCY_UNITS, VERTICAL_PROFILE_SPEC,
                                   column_factor)
from tests import synthetic as syn
from tests.output_checks import VERTICAL_TERMS

REFERENCE_DIR = Path(__file__).resolve().parent / "data" / "reference"
CASES = {"ERA5_fixed": "testdata_ERA5_fixed", "ERA5_track": "testdata_ERA5_track"}
G = float(g.magnitude)


def _read_profile(case, term):
    """Return (levels in Pa, {row label: values}) of a per-level reference file."""
    path = REFERENCE_DIR / case / "results_vertical_levels" / f"{term}_level.csv"
    rows = list(csv.reader(path.open(encoding="utf-8")))
    levels = [float(x) for x in rows[0][1:]]
    return levels, {r[0]: [float(x) for x in r[1:]] for r in rows[1:]}


def _read_main(case):
    path = REFERENCE_DIR / case / f"{CASES[case]}_results.csv"
    rows = list(csv.reader(path.open(encoding="utf-8")))
    header = rows[0]
    return {r[0]: {h: float(x) for h, x in zip(header[1:], r[1:])} for r in rows[1:]}


# ---------------------------------------------------------------- the spec

def test_spec_covers_every_vertical_file_and_nothing_else():
    assert sorted(VERTICAL_PROFILE_SPEC) == sorted(VERTICAL_TERMS)
    for case in CASES:
        written = sorted(p.name[: -len("_level.csv")]
                         for p in (REFERENCE_DIR / case / "results_vertical_levels").glob("*_level.csv"))
        assert written == sorted(VERTICAL_PROFILE_SPEC), case


def test_dz_de_are_kept_apart_and_marked_not_verified():
    assert set(PROVISIONAL_PROFILE_SPEC) == {"Dz", "De"}
    assert not set(PROVISIONAL_PROFILE_SPEC) & set(VERTICAL_PROFILE_SPEC)
    for entry in PROVISIONAL_PROFILE_SPEC.values():
        assert entry["status"] == "NOT VERIFIED"
        # No unit, expression or factor is presented for them.
        assert not {"units", "expression", "column_factor"} & set(entry)


def test_spec_units_are_valid_pint_units():
    ureg = pint.UnitRegistry()
    for term, spec in VERTICAL_PROFILE_SPEC.items():
        ureg.Quantity(1, spec["units"])
        if "column_units" in spec:
            ureg.Quantity(1, spec["column_units"])
        assert spec["verification"].startswith("verified"), term


# ------------------------------------------------- T1: units when written

@pytest.fixture
def written_profiles(tmp_path, monkeypatch):
    """Run every term of the synthetic box and record what each profile write receives."""
    box = syn.make_box(tmp_path)
    captured = {}
    original = term_output.TermOutputMixin._save_vertical_levels

    def spy(self, function, variable_name):
        captured.setdefault(variable_name, []).append(function.data)
        return original(self, function, variable_name)

    monkeypatch.setattr(term_output.TermOutputMixin, "_save_vertical_levels", spy)
    log = syn.SilentLogger()
    ec = EnergyContents(box, "fixed", log)
    ec.calc_az(), ec.calc_ae(), ec.calc_kz(), ec.calc_ke()
    ct = ConversionTerms(box, "fixed", log)
    ct.calc_cz(), ct.calc_ca(), ct.calc_ce(), ct.calc_ck(), ct.calc_c_overturning()
    gd = GenerationDissipationTerms(box, "fixed", log)
    gd.calc_gz(), gd.calc_ge()
    MassContinuity(box, "fixed", log).calc_mass_residual()
    return captured


def test_every_written_profile_has_the_units_of_the_spec(written_profiles):
    assert sorted(written_profiles) == sorted(VERTICAL_PROFILE_SPEC)
    for term, arrays in written_profiles.items():
        expected = VERTICAL_PROFILE_SPEC[term]["units"]
        for data in arrays:
            assert isinstance(data, units.Quantity), f"{term} is written without units"
            scale = (1 * data.units).to(expected).magnitude
            assert scale == pytest.approx(1.0, rel=1e-12, abs=0), (
                f"{term}: written in {data.units}, the spec says {expected}"
            )


# ----------------------------------------- T2: profile-to-column factors

@pytest.mark.parametrize("case", sorted(CASES))
def test_column_is_factor_times_pressure_integral_of_profile(case):
    main = _read_main(case)
    checked = 0
    for term, spec in VERTICAL_PROFILE_SPEC.items():
        factor = column_factor(term)
        if factor is None:
            continue
        levels, rows = _read_profile(case, term)
        assert levels == sorted(levels)
        for time, values in rows.items():
            column = main[time][spec["column_term"]]
            got = factor * np.trapezoid(values, levels)
            assert got == pytest.approx(column, rel=1e-12, abs=0), (term, time)
            checked += 1
    assert checked == 11 * 5


def test_factors_of_kz_ke_and_ck_are_not_one():
    assert column_factor("Kz") == column_factor("Ke") == 1.0 / (2.0 * G)
    assert column_factor("Ck") == 1.0 / G
    assert column_factor("M") is None
    assert column_factor("Ck_1") is None


# ------------------------------------------- T3: column relation of M

def _mass_box(tmp_path, omega_of_p):
    """Synthetic box with no horizontal divergence and a horizontally uniform omega(p)."""
    from src.utils.box_data import BoxData

    ds = syn.make_dataset()
    ds["U"] = xr.full_like(ds["U"], 10.0)
    ds["V"] = xr.full_like(ds["V"], 0.0)
    p = np.asarray(ds["level"].values, dtype=float)
    ds["W"] = xr.zeros_like(ds["W"]) + xr.DataArray(omega_of_p(p), dims="level")
    vlevels = tmp_path / "vlevels"
    vlevels.mkdir(parents=True)
    box = BoxData(
        data=ds, variable_list_df=syn.variable_list(),
        western_limit=syn.LON_MIN, eastern_limit=syn.LON_MAX,
        southern_limit=syn.LAT_MIN, northern_limit=syn.LAT_MAX,
        args=syn.default_args(), results_subdirectory=str(tmp_path),
        results_subdirectory_vertical_levels=str(vlevels),
    )
    return box, np.asarray(box.PressureData.metpy.dequantify().values, dtype=float)


def _profile_and_column(box):
    mc = MassContinuity(box, "fixed", syn.SilentLogger())
    column = np.asarray(mc.calc_mass_residual().metpy.dequantify().values, dtype=float)
    profile = mc.residual_profile.metpy.dequantify().transpose("time", box.VerticalCoordIndexer).values
    return np.asarray(profile, dtype=float), column


def test_mass_column_uses_the_boundary_difference_of_omega(tmp_path):
    """
    With no horizontal divergence and omega = c2 p^2, the column is exactly
    (omega_b - omega_t) / g. The profile holds d(omega)/dp by finite
    differences: exact at interior levels for a quadratic, b (p0 + p1) at the
    first-order edges, so integral(M dp) / g differs from the column by an
    amount fixed in closed form.
    """
    c2 = 4.0e-12  # Pa^-1 s^-1
    box, p = _mass_box(tmp_path, lambda p: c2 * p**2)
    profile, column = _profile_and_column(box)

    expected_column = c2 * (p[-1] ** 2 - p[0] ** 2) / G
    assert np.allclose(column, expected_column, rtol=1e-12, atol=0)

    dodp = 2.0 * c2 * p
    dodp[0] = c2 * (p[0] + p[1])
    dodp[-1] = c2 * (p[-2] + p[-1])
    for row in profile:
        assert np.allclose(row, dodp, rtol=1e-9, atol=0)
    expected_trapz = np.trapezoid(dodp, p) / G
    got_trapz = np.trapezoid(profile, p, axis=1) / G
    assert np.allclose(got_trapz, expected_trapz, rtol=1e-9, atol=0)
    # The two differ by a non-negligible amount for this field.
    assert abs(expected_trapz - expected_column) > 1e-3 * abs(expected_column)


def test_mass_column_and_profile_integral_agree_for_linear_omega(tmp_path):
    box, p = _mass_box(tmp_path, lambda p: -0.05 + 2.0e-7 * p)
    profile, column = _profile_and_column(box)
    assert np.allclose(np.trapezoid(profile, p, axis=1) / G, column, rtol=1e-12, atol=0)


def test_spec_of_m_is_not_a_constant_factor():
    assert VERTICAL_PROFILE_SPEC["M"]["column_factor"] == "special"
    assert "omega_bottom - omega_top" in VERTICAL_PROFILE_SPEC["M"]["column_relation"]


# ------------------------------------------------- T4: split terms

@pytest.mark.parametrize("case", sorted(CASES))
def test_split_terms_combine_as_the_spec_says(case):
    for term, rel in SPLIT_TERM_RELATIONS.items():
        _, total = _read_profile(case, term)
        parts = [_read_profile(case, name)[1] for name in rel["terms"]]
        for time, values in total.items():
            # The pressure factor of Cz/Ce has no time axis in the fixed
            # framework (written as rows "level" and the term, E0-17).
            rows = [p[time] if time in p else p[name] for p, name in zip(parts, rel["terms"])]
            rows = np.array(rows)
            combined = rows.sum(axis=0) if rel["operation"] == "sum" else rows.prod(axis=0)
            got = rel["sign"] * combined
            scale = np.max(np.abs(values))
            assert np.allclose(got, values, rtol=0, atol=1e-12 * scale), (term, time)
    for term, rel in SPLIT_TERM_RELATIONS.items():
        for name in rel["terms"]:
            assert VERTICAL_PROFILE_SPEC[name]["part_of"] == term


# ------------------------------------------ T7: column figure labels

TERM_DETAILS_SHA256_81FD7D2 = "d698d6dd2f0792ef44f05b4e9c90fba823fc63eedceba3cc0429828351137153"


def test_term_details_of_the_column_figures_are_unchanged():
    """TERM_DETAILS labels the column (main results) figures, where its units are right."""
    import src.plots.utils as plot_utils

    text = json.dumps(plot_utils.TERM_DETAILS, sort_keys=True, ensure_ascii=False)
    assert hashlib.sha256(text.encode("utf-8")).hexdigest() == TERM_DETAILS_SHA256_81FD7D2


# --------------------------------------------------- T9: tendencies

def test_tendencies_are_in_watts_per_square_metre():
    ureg = pint.UnitRegistry()
    for name, unit in TENDENCY_UNITS.items():
        assert (ureg.Quantity(1, "J/m**2") / ureg.Quantity(1, "s")).to(unit).magnitude == pytest.approx(1.0)
    assert sorted(TENDENCY_UNITS) == sorted(f"∂{t}/∂t (finite diff.)" for t in ("Az", "Ae", "Kz", "Ke"))


def test_tendency_of_a_linear_series_is_its_slope_for_irregular_times():
    dates = pd.to_datetime(["2020-01-01 00:00", "2020-01-01 01:00", "2020-01-01 04:00",
                            "2020-01-01 06:00", "2020-01-01 12:00"])
    seconds = (dates - dates[0]).total_seconds().to_numpy()
    slopes = {"Az": 3.0, "Ae": -1.5, "Kz": 0.25, "Ke": 2.0}  # W m-2
    df = pd.DataFrame({t: 1.0e5 + s * seconds for t, s in slopes.items()}, index=dates)
    out = calc_budget_diff(df, dates, syn.SilentLogger())
    for t, s in slopes.items():
        assert np.allclose(out[f"∂{t}/∂t (finite diff.)"], s, rtol=1e-12, atol=0)

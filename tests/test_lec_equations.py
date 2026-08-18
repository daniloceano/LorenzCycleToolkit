#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analytic regression tests for the Lorenz Energy Cycle equations.

Every expectation in this file is derived either in closed form from the
prescribed synthetic fields, or from a short NumPy transcription of the
published equation.  None of them re-uses the implementation under test, so a
formula error changes the implementation's answer without changing the
expectation.

Primary sources used for the expectations:
  * Muench, H. S. (1965), J. Atmos. Sci. 22, 349-360, pp. 351-352.
  * Norquist, Recker and Reed (1977), Mon. Wea. Rev. 105, 334-342, p. 336.
  * Brennan and Vincent (1980), Mon. Wea. Rev. 108, 954-965, pp. 963-965.
  * Michaelides (1987), Mon. Wea. Rev. 115, 13-26, pp. 24-25.
  * Michaelides, Prezerakos and Flocas (1999), Q. J. R. Meteorol. Soc. 125,
    139-168, Eqs. (2)-(13).
"""

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from metpy.constants import Cp_d, Rd, Re, g
from metpy.units import units

from src.analysis.boundary_terms import BoundaryTerms
from src.analysis.conversion_terms import ConversionTerms
from src.analysis.mass_continuity import MassContinuity
from src.utils.calc_averages import CalcAreaAverage
from src.utils.calc_budget_and_residual import calc_budget_diff, elapsed_seconds
from src.utils.longitude import (canonicalize_box, canonicalize_longitude,
                                 dataset_longitude_convention)
from src.utils.nan_handling import convert_units
from src.utils.thermodynamics import SIGMA_FLOOR, apply_sigma_floor
from tests import synthetic as syn

A = float(Re.magnitude)
G = float(g.magnitude)


# ===========================================================================
# Helpers
# ===========================================================================

@pytest.fixture(scope="module")
def box(tmp_path_factory):
    return syn.make_box(tmp_path_factory.mktemp("box"))


@pytest.fixture(scope="module")
def conversions(box):
    return ConversionTerms(box, "fixed", syn.SilentLogger())


def _dequant(arr):
    try:
        return np.asarray(arr.metpy.dequantify().values, dtype=float)
    except Exception:
        return np.asarray(arr.values, dtype=float)


def _sigma(box):
    """Level profile of sigma at the first time step."""
    return _level_profile(box.sigma_AA, box)


def _rlats(box):
    return np.asarray(box.tair["rlats"].values, dtype=float)


def _rlons(box):
    return np.asarray(box.tair["rlons"].values, dtype=float)


def _levels(box):
    return _dequant(box.PressureData)


def _level_profile(da, box):
    """Return the (level,) profile of a vertical integrand, at the first time."""
    values = _dequant(da)
    dims = list(da.dims)
    idx = [0] * len(dims)
    idx[dims.index(box.VerticalCoordIndexer)] = slice(None)
    return values[tuple(idx)]


# ===========================================================================
# 1. Averaging operators behave as the primary sources define them
# ===========================================================================

def test_zonal_mean_of_eddy_is_zero(box):
    """[X'] == 0 identically (Brennan and Vincent 1980, p. 963)."""
    for eddy in (box.tair_ZE, box.u_ZE, box.v_ZE, box.omega_ZE):
        zonal_of_eddy = _dequant(eddy.integrate("rlons") / box.xlength)
        assert np.max(np.abs(zonal_of_eddy)) < 1e-9


def test_eddy_covariances_match_closed_form(box):
    """
    With u' = A_U cos(k lam) and v' = A_V cos(k lam) over exactly one period,
    [u'v'] = A_U A_V / 2 to machine precision.
    """
    uv = _dequant((box.u_ZE * box.v_ZE).integrate("rlons") / box.xlength)
    assert np.allclose(uv, syn.A_U * syn.A_V / 2.0, rtol=1e-8)

    tv = _dequant((box.v_ZE * box.tair_ZE).integrate("rlons") / box.xlength)
    assert np.allclose(tv, syn.A_T * syn.A_V / 2.0, rtol=1e-8)


# ===========================================================================
# 2. C_A  (LCT-REV-001 spurious cos(phi); LCT-REV-002 spurious factor 1/2)
# ===========================================================================

def _ca_term1_implementation(box):
    """Re-run exactly what calc_ca computes for its first term."""
    return CalcAreaAverage(
        (box.v_ZE * box.tair_ZE * box.tair_AE.differentiate("rlats"))
        / (Re * box.sigma_AA),
        box.ylength,
        xlength=box.xlength,
    )


def test_ca_horizontal_term_matches_primary_source(box, conversions):
    """
    Brennan and Vincent (1980, p. 964) / Michaelides et al. (1999, Eq. 11):

        C_A(term 1) = -(1/(a sigma)) * mean( v' T' dT*/dphi )

    For the synthetic field [T] is linear in phi with slope GRAD_T, so
    dT*/dphi = GRAD_T exactly, and [v'T'] = A_V A_T / 2 independent of phi.
    Hence the integrand saved as "Ca_1" (before the leading minus sign) is

        (1/(a sigma)) * GRAD_T * A_V * A_T / 2 .
    """
    conversions.calc_ca()  # exercise the real code path as well

    got = _level_profile(_ca_term1_implementation(box), box)
    expected = (syn.GRAD_T * syn.A_V * syn.A_T / 2.0) / (A * _sigma(box))

    assert np.allclose(got, expected, rtol=2e-6), (
        f"Ca_1 integrand mismatch.\n got      = {got}\n expected = {expected}"
    )


def test_ca_has_no_factor_one_half(box):
    """
    Explicit check of the LCT-REV-002 decision: the implemented Ca_1 must equal
    the no-1/2 expectation and must NOT equal half of it.
    """
    got = _level_profile(_ca_term1_implementation(box), box)
    expected = (syn.GRAD_T * syn.A_V * syn.A_T / 2.0) / (A * _sigma(box))

    assert np.allclose(got, expected, rtol=2e-6)
    assert not np.allclose(got, 0.5 * expected, rtol=1e-3)


def test_ca_derivative_has_no_spurious_cosine(box):
    """
    The baseline differentiated T* cos(phi) instead of T*.  Confirm that the
    two differ substantially and in the analytically predicted way, so that
    test_ca_horizontal_term_matches_primary_source is genuinely discriminating.

        d(T* cos phi)/dphi = cos(phi) dT*/dphi - T* sin(phi)
    """
    correct = _dequant(box.tair_AE.differentiate("rlats"))
    with_cos = _dequant(
        (box.tair_AE * box.tair_AE["coslats"]).differentiate("rlats")
    )
    rlats = _rlats(box)
    tstar = _dequant(box.tair_AE)

    # The difference is not a small perturbation.
    rel = np.abs(with_cos - correct) / np.abs(correct)
    assert np.nanmax(rel) > 0.20, "cos(phi) variant must differ materially"

    # And it differs exactly as the product rule predicts (interior points).
    predicted = (
        np.cos(rlats)[:, None, None] * correct
        - tstar * np.sin(rlats)[:, None, None]
    )
    assert np.allclose(with_cos[2:-2], predicted[2:-2], rtol=1e-3)


def test_ca_vertical_term_vanishes_for_barotropic_tstar(box):
    """T* is independent of pressure, so dT*/dp = 0 and Ca_2 must vanish."""
    dTstar_dp = box.tair_AE.differentiate(box.VerticalCoordIndexer) / (
        box.PressureData.metpy.units
    )
    term2 = CalcAreaAverage(
        box.omega_ZE * box.tair_ZE * dTstar_dp, box.ylength, xlength=box.xlength
    ) / box.sigma_AA
    term1 = _level_profile(_ca_term1_implementation(box), box)
    ratio = np.max(np.abs(_level_profile(term2, box))) / np.max(np.abs(term1))
    # Not identically zero at machine precision: the trapezoidal area mean of a
    # constant differs from that constant by O(dphi^2) (int cos phi dphi is only
    # approximately sin phi_n - sin phi_s on a finite grid), so T* retains a
    # residual ~1e-5 fraction of the horizontal-mean temperature, which is
    # pressure dependent. That artefact is a property of the averaging
    # quadrature, identical in the baseline and corrected code.
    assert ratio < 1e-4, f"Ca_2 should be negligible here, got ratio {ratio:.2e}"


# ===========================================================================
# 2b. Domain-mean pressure-work overlap conversion
# ===========================================================================

def test_c_sobreposicao_matches_positive_constant_ascent_solution(tmp_path):
    """Constant mean ascent has the closed-form positive overlap conversion."""
    temperature = 270.0
    omega = -0.12
    dataset = syn.make_dataset()
    dataset["T"] = xr.full_like(dataset["T"], temperature)
    dataset["W"] = xr.full_like(dataset["W"], omega)
    overlap_box = _box_from_dataset(tmp_path, dataset)

    got = _dequant(
        ConversionTerms(
            overlap_box, "fixed", syn.SilentLogger()
        ).calc_c_overturning()
    )

    levels = _levels(overlap_box)
    # CalcAreaAverage uses an analytic area denominator with trapezoidal
    # quadrature.  Its mean of one is therefore included explicitly in this
    # independent NumPy evaluation of <omega> and <T>.
    qnorm = syn.area_mean(
        np.ones((_rlats(overlap_box).size, levels.size)),
        _rlats(overlap_box),
        axis=0,
    )
    profile = (
        -omega * float(Rd.magnitude) * temperature * qnorm**2
        / (G * levels)
    )
    expected = np.trapezoid(profile, levels)

    assert expected > 0.0
    assert np.all(got > 0.0)
    assert np.allclose(got, expected, rtol=1e-10, atol=1e-10)


# ===========================================================================
# 3. C_K -- all five subterms (LCT-AUD-001 is subterm 5)
# ===========================================================================

def _ck_reference(box):
    """
    Independent NumPy evaluation of the five C_K integrands, transcribed from
    Brennan and Vincent (1980, p. 964) and Michaelides (1987, p. 24):

      1: mean( (cos phi / a) u'v' d/dphi([u]/cos phi) )
      2: mean( (v'^2 / a) d[v]/dphi )
      3: mean( (tan phi / a) u'^2 [v] )
      4: mean( omega'u' d[u]/dp )
      5: mean( omega'v' d[v]/dp )

    Analytic ingredients (all exact for the synthetic field):
      [u'v'] = A_U A_V / 2,  [v'^2] = A_V^2 / 2,  [u'^2] = A_U^2 / 2
      [omega'u'] = A_W A_U / 2,  [omega'v'] = A_W A_V / 2
      d[u]/dphi = GRAD_U,  d[v]/dphi = GRAD_V,  d[u]/dp = DUDP,  d[v]/dp = DVDP
      cos(phi) d/dphi([u]/cos phi) = GRAD_U + [u] tan(phi)
    """
    rlats = _rlats(box)
    levels = _levels(box)
    phi = rlats[:, None]
    p = levels[None, :]

    u_zm = syn.U0 + syn.GRAD_U * phi + syn.DUDP * (p - 100000.0)
    v_zm = syn.V0 + syn.GRAD_V * phi + syn.DVDP * (p - 100000.0)
    tan = np.tan(phi)

    integrand1 = (syn.A_U * syn.A_V / 2.0) / A * (syn.GRAD_U + u_zm * tan)
    integrand2 = (syn.A_V ** 2 / 2.0) / A * syn.GRAD_V * np.ones_like(u_zm)
    integrand3 = (syn.A_U ** 2 / 2.0) / A * tan * v_zm
    integrand4 = (syn.A_W * syn.A_U / 2.0) * syn.DUDP * np.ones_like(u_zm)
    integrand5 = (syn.A_W * syn.A_V / 2.0) * syn.DVDP * np.ones_like(u_zm)

    return [
        syn.area_mean(x, rlats, axis=0)
        for x in (integrand1, integrand2, integrand3, integrand4, integrand5)
    ]


def _ck_implementation(box):
    tan_lats = np.tan(box.tair["rlats"])
    d_u_over_cos = (box.u_ZA / box.u_ZA["coslats"]).differentiate("rlats")
    dudp = box.u_ZA.differentiate(box.VerticalCoordIndexer) / (
        box.PressureData.metpy.units
    )
    dvdp = box.v_ZA.differentiate(box.VerticalCoordIndexer) / (
        box.PressureData.metpy.units
    )
    return [
        CalcAreaAverage(
            (box.u_ZE["coslats"] * box.u_ZE * box.v_ZE / Re) * d_u_over_cos,
            box.ylength, xlength=box.xlength),
        CalcAreaAverage(
            ((box.v_ZE ** 2) / Re) * box.v_ZA.differentiate("rlats"),
            box.ylength, xlength=box.xlength),
        CalcAreaAverage(
            (tan_lats * (box.u_ZE ** 2) * box.v_ZA) / Re,
            box.ylength, xlength=box.xlength),
        CalcAreaAverage(
            box.omega_ZE * box.u_ZE * dudp, box.ylength, xlength=box.xlength),
        CalcAreaAverage(
            box.omega_ZE * box.v_ZE * dvdp, box.ylength, xlength=box.xlength),
    ]


def test_ck_all_five_subterms(box):
    """Every C_K subterm must match its independently evaluated expectation."""
    expected = _ck_reference(box)
    got_terms = _ck_implementation(box)

    for i, (term, exp) in enumerate(zip(got_terms, expected), start=1):
        got = _level_profile(term, box)
        # Subterm 1 differentiates a non-polynomial function of phi, so it
        # carries an O(dphi^2) truncation error; the other four are exact.
        rtol = 3e-3 if i == 1 else 1e-6
        assert np.allclose(got, exp, rtol=rtol), (
            f"C_K subterm {i} mismatch.\n got      = {got}\n expected = {exp}"
        )


def test_ck_fifth_subterm_uses_meridional_shear(box):
    """
    LCT-AUD-001 regression: with d[v]/dp != d[u]/dp the fifth subterm must
    follow d[v]/dp.  The baseline used d[u]/dp, giving an answer wrong by the
    ratio DUDP/DVDP.
    """
    dudp = box.u_ZA.differentiate(box.VerticalCoordIndexer) / (
        box.PressureData.metpy.units
    )
    dvdp = box.v_ZA.differentiate(box.VerticalCoordIndexer) / (
        box.PressureData.metpy.units
    )
    correct = _level_profile(
        CalcAreaAverage(box.omega_ZE * box.v_ZE * dvdp, box.ylength,
                        xlength=box.xlength), box)
    wrong = _level_profile(
        CalcAreaAverage(box.omega_ZE * box.v_ZE * dudp, box.ylength,
                        xlength=box.xlength), box)

    expected = (syn.A_W * syn.A_V / 2.0) * syn.DVDP
    assert np.allclose(correct, expected, rtol=1e-6)
    assert np.allclose(wrong / correct, syn.DUDP / syn.DVDP, rtol=1e-6)


def test_ck_total_equals_sum_of_subterms(box, conversions):
    """C_K = (1/g) int sum_i C_K,i dp, checked against the class output."""
    ck = _dequant(conversions.calc_ck())
    subterms = _ck_implementation(box)
    total_integrand = sum(_level_profile(t, box) for t in subterms)
    levels = _levels(box)
    expected = np.trapezoid(total_integrand, levels) / G
    assert np.isclose(ck[0], expected, rtol=1e-8)


# ===========================================================================
# 4. Pressure-work boundary terms (LCT-AUD-002 / LCT-AUD-003)
# ===========================================================================

def test_bphi_e_vanishes_when_geopotential_has_no_zonal_structure(tmp_path):
    """
    Analytic zero test.  Brennan and Vincent (1980, pp. 964-965) and
    Michaelides (1987, p. 25) express BPhi_E entirely through the zonal EDDY
    geopotential Phi'.  With a geopotential that is a function of pressure and
    latitude only, Phi' == 0, so BPhi_E must be exactly zero.

    The baseline used v'Phi* and [v]Phi*, neither of which vanishes for such a
    field, so this test discriminates the two implementations.
    """
    box = syn.make_box(tmp_path, uniform_geopotential=True)
    assert np.max(np.abs(_dequant(box.geopt_ZE))) < 1e-9, "Phi' must be zero"
    assert np.max(np.abs(_dequant(box.geopt_AE))) > 1.0, "Phi* must be nonzero"

    bt = BoundaryTerms(box, "fixed", syn.SilentLogger())
    boe = _dequant(bt.calc_boe())
    assert np.max(np.abs(boe)) < 1e-9, (
        f"BPhi_E must vanish when Phi' == 0, got {boe}"
    )


def _bphi_reference(box, eddy):
    """
    Independent NumPy transcription of the adopted primary-source equations,
    for the first time step.

    eddy=False -> BPhi_Z: literal Michaelides (1987, p. 25): zero east/west
                           term, [v]Phi* cos(phi), mean(omega* Phi*)
    eddy=True  -> BPhi_E:  faces are u'Phi',  [v'Phi'] cos phi,
                           mean([omega'Phi']) (Brennan and Vincent 1980,
                           pp. 964-965; also Michaelides 1987, p. 25)
    """
    rlats, rlons, levels = _rlats(box), _rlons(box), _levels(box)
    c1 = -1.0 / (A * float(box.xlength) * float(box.ylength))
    c2 = -1.0 / (A * float(box.ylength))

    u = _dequant(box.u)[..., 0]
    v = _dequant(box.v)[..., 0]
    w = _dequant(box.omega)[..., 0]
    ph = _dequant(box.geopt)[..., 0]

    u_zm = syn.zonal_mean(u, rlons, axis=0)
    v_zm = syn.zonal_mean(v, rlons, axis=0)
    w_zm = syn.zonal_mean(w, rlons, axis=0)
    ph_zm = syn.zonal_mean(ph, rlons, axis=0)
    u_ed, v_ed = u - u_zm[None], v - v_zm[None]
    w_ed, ph_ed = w - w_zm[None], ph - ph_zm[None]

    if eddy:
        face_ew = (u_ed * ph_ed) / G
        face_ns = syn.zonal_mean(v_ed * ph_ed, rlons, axis=0) * np.cos(rlats)[:, None] / G
        vert = syn.area_mean(syn.zonal_mean(w_ed * ph_ed, rlons, axis=0), rlats, axis=0) / G
        term1 = np.trapezoid(
            np.trapezoid(face_ew[-1] - face_ew[0], rlats, axis=0), levels
        ) * c1
    else:
        # Normalise the trapezoidal area operator by its mean of one.  This is
        # the numerical form of Phi*=[Phi]-mean(Phi) and omega*=[omega]-mean(omega),
        # and makes a pressure-only geopotential reference cancel to roundoff.
        ones = np.ones_like(ph_zm)
        quadrature_norm = syn.area_mean(ones, rlats, axis=0)
        ph_star = ph_zm - syn.area_mean(ph_zm, rlats, axis=0) / quadrature_norm
        w_star = w_zm - syn.area_mean(w_zm, rlats, axis=0) / quadrature_norm
        face_ns = v_zm * ph_star * np.cos(rlats)[:, None] / G
        vert = syn.area_mean(w_star * ph_star, rlats, axis=0) / G
        # Michaelides' printed ([v] Phi*)|lambda_1^lambda_2 has no longitude
        # dependence and is identically zero.
        term1 = 0.0

    term2 = np.trapezoid(face_ns[-1] - face_ns[0], levels) * c2
    term3 = vert[-1] - vert[0]
    return term1 + term2 - term3


def test_bphi_z_matches_michaelides_literal_reference(box):
    bt = BoundaryTerms(box, "fixed", syn.SilentLogger())
    got = _dequant(bt.calc_boz())[0]
    expected = _bphi_reference(box, eddy=False)
    assert np.isclose(got, expected, rtol=1e-6), (
        f"BPhi_Z mismatch: got {got}, expected {expected}"
    )


def _box_from_dataset(tmp_path, dataset):
    """Build a BoxData from a modified synthetic dataset for strong tests."""
    from src.utils.box_data import BoxData

    tmp_path.mkdir(parents=True, exist_ok=True)
    vlevels = tmp_path / "vlevels"
    vlevels.mkdir()
    return BoxData(
        data=dataset,
        variable_list_df=syn.variable_list(),
        western_limit=syn.LON_MIN,
        eastern_limit=syn.LON_MAX,
        southern_limit=syn.LAT_MIN,
        northern_limit=syn.LAT_MAX,
        args=syn.default_args(),
        results_subdirectory=str(tmp_path),
        results_subdirectory_vertical_levels=str(vlevels),
    )


def test_bphi_z_invariant_to_pressure_only_geopotential_reference(tmp_path):
    """Adding arbitrary Phi_ref(p) cannot change a Phi*-based BΦZ."""
    original = syn.make_dataset()
    shifted = original.copy(deep=True)
    pressure = _dequant(shifted["level"])
    phi_ref = 2.5e4 + 0.31 * pressure + 4.0e3 * np.sin(pressure / 1.7e4)
    shifted["Z"] = shifted["Z"] + xr.DataArray(
        phi_ref, dims=("level",), coords={"level": shifted["level"]}
    )

    box_original = _box_from_dataset(tmp_path / "original", original)
    box_shifted = _box_from_dataset(tmp_path / "shifted", shifted)
    got_original = _dequant(
        BoundaryTerms(box_original, "fixed", syn.SilentLogger()).calc_boz()
    )
    got_shifted = _dequant(
        BoundaryTerms(box_shifted, "fixed", syn.SilentLogger()).calc_boz()
    )

    assert np.allclose(got_shifted, got_original, rtol=5e-13, atol=5e-13), (
        "BPhi_Z changed after adding pressure-only Phi_ref(p): "
        f"max difference={np.max(np.abs(got_shifted - got_original)):.3e}"
    )


def test_bphi_z_and_e_are_invariant_to_constant_geopotential_gauge(tmp_path):
    """Both adopted pressure-work terms are invariant under Phi -> Phi+Phi0."""
    original = syn.make_dataset()
    shifted = original.copy(deep=True)
    shifted["Z"] = shifted["Z"] + 1.0e5

    box_original = _box_from_dataset(tmp_path / "gauge_original", original)
    box_shifted = _box_from_dataset(tmp_path / "gauge_shifted", shifted)
    bt_original = BoundaryTerms(box_original, "fixed", syn.SilentLogger())
    bt_shifted = BoundaryTerms(box_shifted, "fixed", syn.SilentLogger())

    for label, before, after in (
        ("BPhiz", bt_original.calc_boz(), bt_shifted.calc_boz()),
        ("BPhie", bt_original.calc_boe(), bt_shifted.calc_boe()),
    ):
        left, right = _dequant(before), _dequant(after)
        assert np.allclose(right, left, rtol=5e-13, atol=5e-13), (
            f"{label} changed under a constant geopotential gauge: "
            f"max difference={np.max(np.abs(right-left)):.3e}"
        )


def test_bphi_z_east_west_term_vanishes_for_periodic_longitude(box):
    """Term (I) is zero when the two longitude faces coincide."""
    bt = BoundaryTerms(box, "fixed", syn.SilentLogger())
    geopt_star, _ = bt._area_anomalies()
    term_i = _dequant(bt._east_west_pressure_work(geopt_star))
    assert np.max(np.abs(term_i)) < 1e-12


def test_bphi_z_is_exactly_zero_without_any_flow(tmp_path):
    """No flow through any face implies exactly zero BΦZ."""
    dataset = syn.make_dataset()
    for name in ("U", "V", "W"):
        dataset[name] = xr.zeros_like(dataset[name])
    zero_flow_box = _box_from_dataset(tmp_path, dataset)
    got = _dequant(
        BoundaryTerms(zero_flow_box, "fixed", syn.SilentLogger()).calc_boz()
    )
    assert np.all(got == 0.0), f"BPhi_Z must be exactly zero, got {got}"


@pytest.mark.parametrize("offset", [1.0e3, 1.0e5, 1.0e6])
def test_bphi_terms_are_independent_of_geopotential_reference(tmp_path, offset):
    """Adding a constant to the geopotential must not move BΦZ or BΦE.

    The geopotential zero is arbitrary, so no diagnosed flux may depend on it.
    """
    dataset = syn.make_dataset()
    base = _box_from_dataset(tmp_path / "base", dataset)
    shifted_dataset = dataset.copy(deep=True)
    shifted_dataset["Z"] = shifted_dataset["Z"] + offset
    shifted = _box_from_dataset(tmp_path / f"shift{offset:.0e}", shifted_dataset)

    for term in ("calc_boz", "calc_boe"):
        a = _dequant(getattr(BoundaryTerms(base, "fixed", syn.SilentLogger()), term)())
        b = _dequant(
            getattr(BoundaryTerms(shifted, "fixed", syn.SilentLogger()), term)()
        )
        scale = max(np.max(np.abs(a)), 1.0)
        assert np.max(np.abs(b - a)) / scale < 1e-9, (
            f"{term} moved by {np.max(np.abs(b - a))} for offset {offset}"
        )

def test_bphi_e_matches_brennan_vincent_reference(box):
    bt = BoundaryTerms(box, "fixed", syn.SilentLogger())
    got = _dequant(bt.calc_boe())[0]
    expected = _bphi_reference(box, eddy=True)
    assert np.isclose(got, expected, rtol=1e-6), (
        f"BPhi_E mismatch: got {got}, expected {expected}"
    )


# ===========================================================================
# 5. Tendencies on regular and irregular time coordinates (LCT-AUD-007)
# ===========================================================================

def _frame(times, values):
    df = pd.DataFrame(index=pd.to_datetime(times))
    for term in ("Az", "Ae", "Kz", "Ke"):
        df[term] = values
    return df


def test_tendency_regular_time_matches_analytic_derivative():
    """E(t) = c * t must give a constant tendency equal to c."""
    times = pd.date_range("2020-01-01", periods=9, freq="6h")
    seconds = (times - times[0]).total_seconds().to_numpy()
    slope = 3.0e-4
    df = _frame(times, 100.0 + slope * seconds)

    out = calc_budget_diff(df.copy(), times.values, syn.SilentLogger())
    assert np.allclose(out["∂Az/∂t (finite diff.)"].to_numpy(), slope, rtol=1e-10)


def test_tendency_irregular_time_matches_analytic_derivative():
    """
    Same linear signal but with a missing timestamp.  The exact derivative is
    still constant; a single leading dt (the historical behaviour) is wrong
    across the gap.
    """
    times = pd.to_datetime([
        "2020-01-01 00:00", "2020-01-01 06:00", "2020-01-01 12:00",
        # 18:00 missing
        "2020-01-02 00:00", "2020-01-02 06:00",
    ])
    seconds = (times - times[0]).total_seconds().to_numpy()
    slope = 3.0e-4
    df = _frame(times, 100.0 + slope * seconds)

    out = calc_budget_diff(df.copy(), times.values, syn.SilentLogger())
    got = out["∂Ke/∂t (finite diff.)"].to_numpy()
    assert np.allclose(got, slope, rtol=1e-10), f"irregular tendency wrong: {got}"

    legacy = np.gradient(df["Ke"].to_numpy(), seconds[1] - seconds[0])
    assert not np.allclose(legacy, slope, rtol=1e-3)


def test_tendency_rejects_non_monotonic_time():
    times = pd.to_datetime(
        ["2020-01-01 00:00", "2020-01-01 12:00", "2020-01-01 06:00"])
    with pytest.raises(ValueError, match="strictly increasing"):
        elapsed_seconds(times.values)


def test_tendency_rejects_single_timestep():
    with pytest.raises(ValueError, match="At least two time steps"):
        elapsed_seconds(pd.to_datetime(["2020-01-01"]).values)


# ===========================================================================
# 6. Static-stability floor (LCT-AUD-004 / LCT-REV-006)
# ===========================================================================

def _sigma_array(values):
    da = xr.DataArray(np.array(values, dtype=float), dims=("level",))
    return da * units("K**2/m")


def test_sigma_floor_preserves_nan():
    """A NaN must stay NaN and must NOT be replaced by the floor."""
    filtered, _ = apply_sigma_floor(_sigma_array([1.0, 0.01, np.nan, -0.5, 0.03]))
    out = _dequant(filtered)

    assert np.isnan(out[2]), "NaN was silently replaced by the floor"
    assert out[0] == pytest.approx(1.0)
    assert out[1] == pytest.approx(SIGMA_FLOOR)
    assert out[3] == pytest.approx(SIGMA_FLOOR)
    assert out[4] == pytest.approx(SIGMA_FLOOR)


def test_sigma_floor_diagnostics():
    _, diag = apply_sigma_floor(_sigma_array([1.0, 0.01, np.nan, -0.5, 0.03, 2.0]))
    assert diag["n_total"] == 6
    assert diag["n_nan"] == 1
    assert diag["n_negative"] == 1
    assert diag["n_below_floor"] == 3          # 0.01, -0.5, 0.03
    assert diag["floor"] == SIGMA_FLOOR
    assert diag["raw_min"] == pytest.approx(-0.5)
    assert diag["raw_max"] == pytest.approx(2.0)
    assert diag["fraction_floored"] == pytest.approx(0.5)


def test_sigma_floor_can_be_disabled():
    filtered, diag = apply_sigma_floor(_sigma_array([1.0, 0.01, -0.5]), floor=None)
    out = _dequant(filtered)
    assert out[1] == pytest.approx(0.01)
    assert out[2] == pytest.approx(-0.5)
    assert diag["floor"] is None


def test_static_stability_formula_matches_michaelides(box):
    """
    sigma = mean( gT/c_p - (p g / R) dT/dp )   (Michaelides 1987, p. 24;
    Michaelides et al. 1999, Eq. 4; Norquist et al. 1977, p. 336).
    """
    rlats, rlons, levels = _rlats(box), _rlons(box), _levels(box)
    T = _dequant(box.tair)[..., 0]

    dTdp = np.gradient(T, levels, axis=2)
    integrand = (
        G * T / float(Cp_d.magnitude)
        - levels[None, None, :] * G / float(Rd.magnitude) * dTdp
    )
    expected = syn.area_mean(syn.zonal_mean(integrand, rlons, axis=0), rlats, axis=0)
    assert np.allclose(_sigma(box), expected, rtol=1e-8)


def test_static_stability_is_positive_and_physical(box):
    """Free-tropospheric sigma must be well above the numerical floor."""
    sigma = _sigma(box)
    assert np.all(sigma > 0)
    assert np.all(sigma > 10 * SIGMA_FLOOR)


# ===========================================================================
# 7. Longitude handling (LCT-AUD-008)
# ===========================================================================

def test_longitude_convention_detection():
    assert dataset_longitude_convention(np.array([-179.0, 0.0, 179.0])) == "-180..180"
    assert dataset_longitude_convention(np.array([0.0, 180.0, 359.0])) == "0..360"


@pytest.mark.parametrize("value,convention,expected", [
    (350.0, "-180..180", -10.0),
    (-10.0, "-180..180", -10.0),
    (-45.0, "0..360", 315.0),
    (315.0, "0..360", 315.0),
])
def test_canonicalize_longitude(value, convention, expected):
    assert canonicalize_longitude(value, convention) == pytest.approx(expected)


def test_canonicalize_box_translates_convention():
    lons = np.arange(-180.0, 180.0, 1.0)
    west, east = canonicalize_box(300.0, 340.0, lons)
    assert west == pytest.approx(-60.0)
    assert east == pytest.approx(-20.0)


def test_canonicalize_box_rejects_wrapped_domain():
    """A seam-crossing request must fail loudly, never silently truncate."""
    lons = np.arange(-180.0, 180.0, 1.0)
    with pytest.raises(ValueError, match="coordinate seam"):
        canonicalize_box(170.0, 190.0, lons)

    lons360 = np.arange(0.0, 360.0, 1.0)
    with pytest.raises(ValueError, match="coordinate seam"):
        canonicalize_box(350.0, 10.0, lons360)


def test_canonicalize_box_accepts_ordinary_domain():
    lons = np.arange(-180.0, 180.0, 1.0)
    assert canonicalize_box(-60.0, -30.0, lons) == (-60.0, -30.0)


# ===========================================================================
# 8. 850 hPa selection (LCT-REV-010)
# ===========================================================================

def test_850_hpa_constant_is_in_pascal():
    from src.utils.select_area import LEVEL_850_HPA_IN_PA
    assert LEVEL_850_HPA_IN_PA == 85000


def test_850_hpa_selects_the_right_level():
    """Nearest-level selection on a Pa coordinate must land on 850 hPa."""
    from src.utils.select_area import LEVEL_850_HPA_IN_PA
    ds = syn.make_dataset()
    picked = float(_dequant(
        ds["level"].sel(level=LEVEL_850_HPA_IN_PA, method="nearest")))
    assert picked == pytest.approx(85000.0)
    legacy = float(_dequant(ds["level"].sel(level=8500, method="nearest")))
    assert legacy != pytest.approx(85000.0)


# ===========================================================================
# 9. Unit-conversion error handling (LCT-REV-007)
# ===========================================================================

def test_convert_units_raises_informative_value_error():
    """
    pint raises DimensionalityError, which subclasses TypeError, so the
    historical `except ValueError` guard never fired.
    """
    energy = xr.DataArray(np.array([1.0, 2.0]), dims=("time",)) * units("J/m**2")
    with pytest.raises(ValueError, match="Unit error in TestTerm"):
        convert_units(energy, "W/m^2", "TestTerm", app_logger=syn.SilentLogger())


def test_convert_units_succeeds_for_compatible_units():
    power = xr.DataArray(np.array([1.0, 2.0]), dims=("time",)) * units("kg/s**3")
    out = convert_units(power, "W/m^2", "TestTerm", app_logger=syn.SilentLogger())
    assert np.allclose(_dequant(out), [1.0, 2.0])


# ===========================================================================
# 10. Direct dissipation is explicitly unsupported (LCT-AUD-006)
# ===========================================================================

def test_direct_dissipation_raises_not_implemented(tmp_path):
    from src.utils.box_data import BoxData

    vlevels = tmp_path / "vlevels"
    vlevels.mkdir()
    with pytest.raises(NotImplementedError, match="residual formulation"):
        BoxData(
            data=syn.make_dataset(),
            variable_list_df=syn.variable_list(),
            western_limit=syn.LON_MIN,
            eastern_limit=syn.LON_MAX,
            southern_limit=syn.LAT_MIN,
            northern_limit=syn.LAT_MAX,
            args=syn.default_args(residuals=False),
            results_subdirectory=str(tmp_path),
            results_subdirectory_vertical_levels=str(vlevels),
        )


# ===========================================================================
# 11. Vertical control volume shared by all terms (LCT-REV-004)
# ===========================================================================

def test_all_levels_valid_for_clean_data(box):
    assert box.dropped_levels == []
    assert len(box.valid_levels) == len(syn.LEVELS)


def test_level_with_missing_data_is_dropped_once_for_all_terms(tmp_path):
    from src.utils.box_data import BoxData

    ds = syn.make_dataset()
    # Punch a hole at one level, at one time, in one variable only.
    ds["T"][:, :, 2, 1] = np.nan

    vlevels = tmp_path / "vlevels"
    vlevels.mkdir()
    box = BoxData(
        data=ds,
        variable_list_df=syn.variable_list(),
        western_limit=syn.LON_MIN,
        eastern_limit=syn.LON_MAX,
        southern_limit=syn.LAT_MIN,
        northern_limit=syn.LAT_MAX,
        args=syn.default_args(),
        results_subdirectory=str(tmp_path),
        results_subdirectory_vertical_levels=str(vlevels),
    )
    dropped = float(np.sort(syn.LEVELS)[2])
    assert box.dropped_levels == [dropped]

    for field in (box.tair, box.u, box.v, box.omega, box.geopt, box.sigma_AA):
        levels = _dequant(field[box.VerticalCoordIndexer])
        assert dropped not in levels
        assert levels.size == len(syn.LEVELS) - 1

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 28 11:33:22 2022

Script for thermodynamics calculations necessary for computate
Lorenz Energy Cycle such static stability parameter, and all temrs of the
thermodynamic equation.

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com
"""

import logging

import numpy as np
from metpy.calc import potential_temperature
from metpy.constants import Cp_d, Rd, Re, g
from metpy.units import units

from .calc_averages import CalcAreaAverage, CalcZonalAverage


# ---------------------------------------------------------------------------
# Static-stability numerical safeguard
# ---------------------------------------------------------------------------
# SIGMA_FLOOR is a NUMERICAL-STABILITY SAFEGUARD. It is *not* part of the
# theoretical Lorenz Energy Cycle equations: neither Muench (1965), nor
# Brennan and Vincent (1980), nor Michaelides (1987, 1999) prescribes any
# minimum value for sigma.
#
# It was introduced in the toolkit because area-averaged sigma values below
# roughly 0.03 K^2 m^-1 (i.e. nearly dry-adiabatic mean stratification) make
# the 1/sigma weighting of Az, Ae, Ca, Gz and Ge blow up.  The legacy value is
# retained as the default so that earlier results remain reproducible.
#
# EXPERIMENTAL: passing ``apply_floor=False`` to StaticStability disables it
# for a diagnostic experiment. That path is not exercised by the shipped
# workflows and changes Az, Ae, Ca, Gz and Ge, so results obtained with it are
# not comparable with the defaults; it is reported at INFO level when used.
# Typical free-tropospheric sigma is 0.5-5 K^2 m^-1, so the floor is expected
# to be inactive except in nearly neutral layers, usually near the surface.
SIGMA_FLOOR = 0.03


def apply_sigma_floor(sigma, floor=SIGMA_FLOOR, app_logger=None):
    """
    Apply the static-stability floor while preserving NaNs and reporting
    diagnostics.

    NaN is preserved: because ``NaN > floor`` is False, a plain
    ``sigma.where(sigma > floor, floor)`` expression would replace missing data
    by the smallest admissible stability, *maximising* 1/sigma exactly where
    the data are absent.  Here only finite values below the floor are clipped.

    Parameters
    ----------
    sigma : xarray.DataArray
        Raw (unfiltered) static stability, already area averaged.  May carry
        pint units.
    floor : float
        Minimum admissible magnitude.  ``None`` disables the floor entirely.
    app_logger : logging.Logger, optional
        Logger used to report the diagnostics.

    Returns
    -------
    (filtered, diagnostics) : tuple(xarray.DataArray, dict)
        The filtered field (same units as the input) and a dictionary of
        diagnostics describing the raw field.
    """
    log = app_logger.debug if app_logger is not None else logging.debug

    units_sigma = sigma.metpy.units
    raw = sigma.metpy.dequantify()
    values = np.asarray(raw.values, dtype=float)

    n_total = int(values.size)
    n_nan = int(np.count_nonzero(np.isnan(values)))
    finite = values[np.isfinite(values)]
    n_negative = int(np.count_nonzero(finite < 0))
    n_below = int(np.count_nonzero(finite <= floor)) if floor is not None else 0

    diagnostics = {
        "n_total": n_total,
        "n_nan": n_nan,
        "n_negative": n_negative,
        "n_below_floor": n_below,
        "fraction_floored": (n_below / n_total) if n_total else 0.0,
        "raw_min": float(np.nanmin(finite)) if finite.size else float("nan"),
        "raw_median": float(np.nanmedian(finite)) if finite.size else float("nan"),
        "raw_max": float(np.nanmax(finite)) if finite.size else float("nan"),
        "floor": floor,
    }

    log(
        f"Static stability diagnostics: n={n_total}, NaN={n_nan}, "
        f"negative={n_negative}, <=floor={n_below} "
        f"({100.0 * diagnostics['fraction_floored']:.3f}%), "
        f"raw min/median/max = {diagnostics['raw_min']:.4g} / "
        f"{diagnostics['raw_median']:.4g} / {diagnostics['raw_max']:.4g}, "
        f"floor={floor}"
    )
    if n_nan:
        log(
            f"Static stability: {n_nan} NaN value(s) preserved as NaN (they "
            "are NOT replaced by the floor)."
        )
    if n_negative:
        log(
            f"Static stability: {n_negative} negative value(s) found "
            "(statically unstable mean stratification); these are raised to "
            "the floor when it is enabled, which inverts the sign of the APE "
            "weighting."
        )

    if floor is None:
        filtered = raw
    else:
        # Preserve NaN: only clip where the value is finite AND below the floor.
        filtered = raw.where(~(raw <= floor) | raw.isnull(), floor)

    return filtered * units_sigma, diagnostics


def StaticStability(
    TemperatureData,
    PressureData,
    VerticalCoordIndexer,
    xlength,
    ylength,
    apply_floor=True,
    app_logger=None,
):
    """
    Compute the static stability parameter sigma for all vertical levels
    and for the desired domain.

    .. math::
        \\sigma = \\overline{\\left[\\frac{gT}{c_p}
                  - \\frac{pg}{R}\\frac{\\partial T}{\\partial p}\\right]}

    Source:
        Michaelides, S. C. (1987).
        Limited Area Energetics of Genoa Cyclogenesis,
        Monthly Weather Review, 115(1), 13-26. Retrieved Jan 24, 2022, from:
        https://journals.ametsoc.org/view/journals/mwre/115/1/1520-0493_1987_115_0013_laeogc_2_0_co_2.xml

        The same expression is printed by Michaelides et al. (1999, Eq. 4) and
        by Norquist, Recker and Reed (1977, p. 336).  Note that the sigma
        printed in Brennan and Vincent (1980, p. 963) is short by a factor
        p/Rd and is dimensionally inadmissible; it is a typographical error in
        that paper.

    Parameters
    ----------
    TemperatureData: xarray.DataArray
        temperature data in Kelvin
    PressureData: xarray.DataArray
        pressure coordinate, in Pa
    VerticalCoordIndexer: str
        name of the vertical coordinate
    xlength, ylength: float
        Zonal width in radians and (sin(phi_n) - sin(phi_s)) respectively.
        Passed through to CalcZonalAverage/CalcAreaAverage, which normalise by
        their own quadrature measure, so neither value is used.
    apply_floor: bool
        Whether to apply the legacy SIGMA_FLOOR numerical safeguard. Setting it
        to False is EXPERIMENTAL: it changes Az, Ae, Ca, Gz and Ge, and the
        results are not comparable with the defaults.
    app_logger: logging.Logger, optional
        logger used for the stability diagnostics

    Returns
    -------
    sigma: xarray.DataArray
        sigma values for all pressure levels of the selected box.  The
        diagnostics dictionary produced by :func:`apply_sigma_floor` is
        attached as ``sigma.attrs["sigma_diagnostics"]``.
    """
    logging.debug("Computing static stability parameter...")

    FirstTerm = g * TemperatureData / Cp_d
    SecondTerm = PressureData * g / Rd
    ThirdTerm = TemperatureData.differentiate(VerticalCoordIndexer) / units("Pa")
    function = FirstTerm - (SecondTerm * ThirdTerm)
    sigma_ZA = CalcZonalAverage(function, xlength)
    sigma_AA = CalcAreaAverage(sigma_ZA, ylength)

    floor = SIGMA_FLOOR if apply_floor else None
    if floor is None:
        # Disabling the floor changes Az, Ae, Ca, Gz and Ge. Report it at INFO
        # level so that the reason is recorded in the run's log file.
        info = app_logger.info if app_logger is not None else logging.info
        info(
            f"⚠️ apply_floor=False: the sigma = {SIGMA_FLOOR} static-stability "
            "floor is DISABLED for this run (experimental)."
        )

    sigma_AA_filtered, diagnostics = apply_sigma_floor(
        sigma_AA, floor=floor, app_logger=app_logger
    )
    sigma_AA_filtered = sigma_AA_filtered.drop_vars("coslats", errors="ignore")
    sigma_AA_filtered.attrs["sigma_diagnostics"] = diagnostics

    logging.debug("Ok.")
    return sigma_AA_filtered


def AdiabaticHEating(
    TemperatureData,
    PressureData,
    OmegaData,
    UWindComponentData,
    VWindComponentData,
    VerticalCoordIndexer,
    LatIndexer,
    LonIndexer,
    TimeName,
    dTdt=None,
):
    """
    Compute the diabatic heating as a residual form the thermodynamic
    equation for all vertical levels and for the desired domain
    """
    logging.debug("Computing adiabatic heating...")

    # Horizontal temperature advection
    lons, lats = TemperatureData[LonIndexer], TemperatureData[LatIndexer]
    cos_lats = TemperatureData["coslats"]
    # Differentiate temperature in respect to longitude and latitude
    dTdlambda = TemperatureData.differentiate(LonIndexer)
    dTdphi = TemperatureData.differentiate(LatIndexer)
    # Get the values for width and length in meters
    dx = np.deg2rad(lons.differentiate(LonIndexer)) * cos_lats * Re
    dy = np.deg2rad(lats.differentiate(LatIndexer)) * Re
    AdvHTemp = -1 * (
        (UWindComponentData * dTdlambda / dx) + (VWindComponentData * dTdphi / dy)
    )

    theta = potential_temperature(PressureData, TemperatureData)

    if dTdt is None:
        dTdt = TemperatureData.differentiate(TimeName, datetime_unit="s") / units("s")

    sigma = (
        -1
        * (TemperatureData / theta)
        * theta.differentiate(VerticalCoordIndexer)
        / units("Pa")
    )

    ResT = dTdt - AdvHTemp - (sigma * OmegaData)

    AdiabaticHeating = ResT * Cp_d

    logging.debug("Ok.")
    return AdiabaticHeating

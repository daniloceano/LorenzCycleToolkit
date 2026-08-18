#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for computing zonal and area averages in a sphere.

Source:
    Brennan, F. E., & Vincent, D. G. (1980).
    Zonal and Eddy Components of the Synoptic-Scale Energy Budget
    during Intensification of Hurricane Carmen (1974),
    Monthly Weather Review, 108(7), 954-965. Retrieved Jan 25, 2022, from:
    https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_0954_zaecot_2_0_co_2.xml

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com

"""

import xarray as xr


def CalcZonalAverage(VariableData, xlength):
    """
    Computates variable zonal average of some variable, for all z
    levels and time steps.

    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated. Requires dimension rlons
        (longitude in radians)
    xlength: float
        Length (in radians), of the data as eastern limit minus western limit.
        Kept for backwards compatibility; the average is normalised by the
        quadrature's own measure, so this value is not used.

    Returns
    -------
    ZA: xarray.Dataset
        Arrays of zonal avreages for all longitudes from the passed Dataset

    Notes
    -----
    Numerator and denominator use the same trapezoidal rule, so the zonal
    average of a constant is that constant on any grid and at any coordinate
    precision.  Normalising by the analytic ``xlength`` instead would leave a
    residual that, applied to the geopotential, injects a spurious constant
    into the departure fields.
    """
    measure = xr.ones_like(VariableData["rlons"]).integrate("rlons")
    return VariableData.integrate("rlons") / measure


def CalcAreaAverage(VariableData, ylength, xlength=False):
    """
    Computates the Area Average of a function.

    The default is to computate the zonal average and then a meridional average.
    If the xlength is provided, it will firstly compute the zonal average and
    then, the area average.

    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated
    ylength: float
        Length (in radians), of the data as northern limit minus southern limit.
        Kept for backwards compatibility; the average is normalised by the
        quadrature's own measure, so this value is not used.
    xlength: float (optional)
        Length (in radians), of the data as eastern limit minus western limit.
        If passed, it will first compute zonal averages

    Returns
    -------
    AA: xarray.Dataset
        Arrays of area avreages for all latitudes and longitudes from
        the passed Dataset

    Notes
    -----
    The cosine weight is integrated by the same trapezoidal rule used for the
    field itself, rather than replaced by the analytic
    ``sin(phi_n) - sin(phi_s)``.  The two differ by O(dphi^2) -- about
    1.6e-4 on a 2.5 degree grid and 1.6e-6 on a 0.25 degree grid -- which is
    negligible for most terms but not for departures from the area mean: with
    a geopotential of order 1e5 m2/s2 the mismatch injects a constant of a few
    m2/s2 into ``Phi*``, and makes the diagnosed fluxes depend on the arbitrary
    reference level of the geopotential.  Matching the two quadratures makes
    the area average of a constant exact by construction, so the departures
    satisfy their defining identity and every derived flux is independent of
    that reference level.
    """
    # Compute zonal average if requested
    if xlength:
        ZA = CalcZonalAverage(VariableData, xlength)
    else:
        ZA = VariableData
    weight = ZA["coslats"]
    measure = weight.integrate("rlats")
    return ((ZA * weight).integrate("rlats") / measure).drop_vars(
        "coslats", errors="ignore"
    )

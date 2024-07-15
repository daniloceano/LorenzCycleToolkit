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

import numpy as np


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

    Returns
    -------
    ZA: xarray.Dataset
        Arrays of zonal avreages for all longitudes from the passed Dataset
    """
    return VariableData.integrate("rlons") / xlength


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
    xlength: float (optional)
        Length (in radians), of the data as eastern limit minus western limit.
        If passed, it will first compute zonal averages

    Returns
    -------
    AA: xarray.Dataset
        Arrays of area avreages for all latitudes and longitudes from
        the passed Dataset
    """
    # Compute zonal average if requested
    if xlength:
        ZA = CalcZonalAverage(VariableData, xlength)
    else:
        ZA = VariableData
    ylength = np.sin(VariableData["rlats"][-1]) - np.sin(VariableData["rlats"][0])
    return ((ZA * ZA["coslats"]).integrate("rlats") / ylength).drop_vars(
        "coslats", errors="ignore"
    )

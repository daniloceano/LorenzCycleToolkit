#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared NaN and unit-conversion policy for the Lorenz Energy Cycle terms.

Two invariants of the toolkit live here:

1. **The vertical control volume is fixed once**, in
   :class:`~src.utils.box_data.BoxData`, and every term inherits it.  Terms
   never drop pressure levels of their own: :func:`handle_nans` interpolates
   interior gaps along the vertical, never extrapolates, and preserves any NaN
   that survives, so a term whose integration depth would change becomes NaN
   instead of silently changing depth.  All terms of a budget are therefore
   integrated over the same levels and the budget identities hold.

2. **Unit errors are reported, not swallowed.** :func:`convert_units` catches
   pint's ``DimensionalityError`` (a ``TypeError`` subclass) as well as
   ``ValueError``, and re-raises a ``ValueError`` naming the variable, its
   current units and the target.

"""

import logging

import numpy as np

try:  # pragma: no cover - import shape depends on the pint version
    from pint.errors import DimensionalityError, UndefinedUnitError

    _UNIT_ERRORS = (DimensionalityError, UndefinedUnitError, ValueError, TypeError)
except Exception:  # pragma: no cover
    _UNIT_ERRORS = (ValueError, TypeError)


def handle_nans(function, vertical_coord, variable_name="", app_logger=None):
    """
    Interpolate interior NaNs along the vertical coordinate, with logging.

    This function never drops pressure levels:
    the vertical control volume is fixed once for the whole budget (see
    :meth:`src.utils.box_data.BoxData._restrict_to_valid_levels`).  It also
    never extrapolates, so values below terrain are not invented; any NaN that
    survives interpolation is preserved and reported, making the integral NaN
    rather than quietly changing the integration depth.

    Parameters
    ----------
    function : xarray.DataArray
        The integrand.
    vertical_coord : str
        Name of the vertical coordinate.
    variable_name : str
        Label used in log messages.
    app_logger : logging.Logger, optional

    Returns
    -------
    xarray.DataArray
    """
    log = app_logger if app_logger is not None else logging.getLogger(
        "lorenzcycletoolkit"
    )

    if vertical_coord not in function.dims:
        return function

    n_nan = int(np.count_nonzero(np.isnan(np.asarray(function.values, dtype=float))))
    if n_nan == 0:
        return function

    units = None
    try:
        units = function.metpy.units
    except Exception:  # pragma: no cover - not all arrays are quantified
        units = None

    # ``interpolate_na`` without ``fill_value`` does not extrapolate, so points
    # outside the range of valid data (e.g. below terrain) stay NaN.
    interpolated = function.interpolate_na(dim=vertical_coord, use_coordinate=True)
    if units is not None:
        interpolated = interpolated * units

    n_remaining = int(
        np.count_nonzero(np.isnan(np.asarray(interpolated.values, dtype=float)))
    )
    n_filled = n_nan - n_remaining

    log.warning(
        "⚠️ %s: %d NaN value(s) found along %s; %d filled by interior "
        "interpolation, %d preserved as NaN (no extrapolation, no level "
        "dropping).",
        variable_name or "integrand",
        n_nan,
        vertical_coord,
        n_filled,
        n_remaining,
    )

    return interpolated


def convert_units(function, target, variable_name, app_logger=None):
    """
    Convert a quantified DataArray to ``target`` units with an informative error.

    Parameters
    ----------
    function : xarray.DataArray
    target : str
        e.g. ``"W/m^2"`` or ``"J/m^2"``.
    variable_name : str
    app_logger : logging.Logger, optional

    Returns
    -------
    xarray.DataArray

    Raises
    ------
    ValueError
        If the conversion is not dimensionally possible.  pint signals this
        with ``DimensionalityError``, a ``TypeError`` subclass, so the guard
        here catches ``TypeError`` as well as ``ValueError``.
    """
    log = app_logger if app_logger is not None else logging.getLogger(
        "lorenzcycletoolkit"
    )
    try:
        return function.metpy.convert_units(target)
    except _UNIT_ERRORS as exc:
        try:
            current = str(function.metpy.units)
        except Exception:  # pragma: no cover
            current = "unknown"
        message = (
            f"Unit error in {variable_name}: cannot convert from '{current}' "
            f"to '{target}'. This usually means a factor is missing from the "
            f"term's formula. Original error: {exc}"
        )
        log.error(message)
        raise ValueError(message) from exc

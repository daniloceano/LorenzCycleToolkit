#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Longitude canonicalisation helpers.

The Lorenz Energy Cycle zonal operator is an integral over one continuous
longitude arc,

.. math::
    [X] = \\frac{1}{\\lambda_e - \\lambda_w}\\int_{\\lambda_w}^{\\lambda_e} X\\,d\\lambda ,

with a positive width (Brennan and Vincent 1980, p. 963).  Longitude, however,
is cyclic, and requested domain limits may be expressed in a different
convention from the dataset (0..360 versus -180..180), or may straddle the
coordinate seam.

This module keeps that handling in one place so that a requested domain is
either (a) expressed unambiguously in the dataset's own convention, or
(b) rejected with an explicit message.  It never silently returns a different
subset from the one requested.

"""

import numpy as np


def dataset_longitude_convention(lons) -> str:
    """
    Classify the longitude convention of a coordinate array.

    Returns
    -------
    str
        ``"-180..180"`` if any longitude is negative, otherwise ``"0..360"``.
    """
    values = np.asarray(lons, dtype=float)
    return "-180..180" if np.nanmin(values) < 0 else "0..360"


def canonicalize_longitude(value: float, convention: str) -> float:
    """
    Express a longitude in the given convention.

    Parameters
    ----------
    value : float
        Longitude in degrees, any convention.
    convention : str
        ``"-180..180"`` or ``"0..360"``.

    Returns
    -------
    float
        The equivalent longitude in the requested convention.

    Examples
    --------
    >>> canonicalize_longitude(350.0, "-180..180")
    -10.0
    >>> canonicalize_longitude(-45.0, "0..360")
    315.0
    """
    if convention == "-180..180":
        return float((value + 180.0) % 360.0 - 180.0)
    if convention == "0..360":
        return float(value % 360.0)
    raise ValueError(f"Unknown longitude convention: {convention!r}")


def canonicalize_box(
    min_lon: float,
    max_lon: float,
    lons,
    app_logger=None,
    context: str = "domain",
):
    """
    Express a requested longitude interval in the dataset's convention and
    verify that it is representable as a single contiguous slice.

    Parameters
    ----------
    min_lon, max_lon : float
        Requested western and eastern limits, in degrees, any convention.
    lons : array-like
        The dataset's longitude coordinate.
    app_logger : logging.Logger, optional
    context : str
        Label used in log/error messages.

    Returns
    -------
    (float, float)
        Western and eastern limits in the dataset's convention, with
        ``west < east``.

    Raises
    ------
    ValueError
        If ``min_lon`` is greater than ``max_lon`` (swapped limits), or if the
        requested arc crosses the dataset's coordinate seam.  Wrapped domains
        are not supported by the zonal-average operator, and a silent partial
        selection would corrupt every zonal mean, so both are raised
        explicitly rather than handled implicitly.
    """
    values = np.asarray(lons, dtype=float)
    convention = dataset_longitude_convention(values)

    west = canonicalize_longitude(min_lon, convention)
    east = canonicalize_longitude(max_lon, convention)

    # A requested arc of exactly 360 degrees canonicalises to west == east.
    requested_width = (max_lon - min_lon) % 360.0
    if requested_width == 0 and max_lon != min_lon:
        raise ValueError(
            f"{context}: a full 360-degree longitude band is not supported by "
            "the limited-area zonal operator."
        )

    # A request whose own numbers run backwards is, in the overwhelming
    # majority of cases, two swapped values in the box-limits file.  The only
    # other reading -- a domain that wraps across the coordinate seam -- is not
    # supported either, so both are rejected here, with the cheap fix named
    # first.  The seam test below then handles the remaining case: limits given
    # in increasing order that still wrap once expressed in the dataset's
    # convention (e.g. 170 to 190 on a -180..180 dataset).
    if min_lon > max_lon:
        raise ValueError(
            f"{context}: min_lon ({min_lon}) is greater than max_lon "
            f"({max_lon}). If the two values are swapped, swap them back in "
            "the box-limits file. If a domain crossing the coordinate seam "
            "was intended, it is not supported: the zonal average requires "
            "one contiguous arc of increasing longitude."
        )

    if east <= west:
        raise ValueError(
            f"{context}: the requested longitude interval "
            f"[{min_lon}, {max_lon}] crosses the dataset's coordinate seam "
            f"(dataset convention {convention}, "
            f"range [{values.min():.3f}, {values.max():.3f}]). "
            "Wrapped longitude domains are not supported: the zonal average "
            "requires one contiguous arc. Re-express the data in the other "
            "longitude convention (see convert_longitude_range) so that the "
            "domain becomes contiguous, then retry."
        )

    if app_logger is not None and (west != min_lon or east != max_lon):
        app_logger.info(
            f"🌐 {context}: longitude limits [{min_lon:.3f}, {max_lon:.3f}] "
            f"canonicalised to [{west:.3f}, {east:.3f}] for dataset "
            f"convention {convention}."
        )

    return west, east


def verify_selected_domain(
    selected_lons,
    selected_lats,
    requested,
    app_logger=None,
    context: str = "domain",
):
    """
    Confirm that a sliced dataset actually covers the requested domain and
    report the requested versus realised (nearest-grid) limits.

    Parameters
    ----------
    selected_lons, selected_lats : array-like
        Coordinates of the sliced dataset.
    requested : dict
        ``{"min_lon", "max_lon", "min_lat", "max_lat"}`` in the dataset's
        convention.
    app_logger : logging.Logger, optional
    context : str

    Returns
    -------
    dict
        Realised limits and the offsets from the request.

    Raises
    ------
    ValueError
        If the slice is empty, or if it misses the requested limits by more
        than one grid spacing, which indicates the requested domain is not
        contained in the data.
    """
    lons = np.asarray(selected_lons, dtype=float)
    lats = np.asarray(selected_lats, dtype=float)

    if lons.size == 0 or lats.size == 0:
        raise ValueError(
            f"{context}: the requested domain "
            f"lon=[{requested['min_lon']}, {requested['max_lon']}], "
            f"lat=[{requested['min_lat']}, {requested['max_lat']}] "
            "selected an empty subset of the data."
        )

    dlon = float(np.min(np.diff(lons))) if lons.size > 1 else 0.0
    dlat = float(np.min(np.diff(lats))) if lats.size > 1 else 0.0

    realised = {
        "min_lon": float(lons.min()),
        "max_lon": float(lons.max()),
        "min_lat": float(lats.min()),
        "max_lat": float(lats.max()),
        "n_lon": int(lons.size),
        "n_lat": int(lats.size),
        "dlon": dlon,
        "dlat": dlat,
    }

    offsets = {
        "min_lon": realised["min_lon"] - requested["min_lon"],
        "max_lon": realised["max_lon"] - requested["max_lon"],
        "min_lat": realised["min_lat"] - requested["min_lat"],
        "max_lat": realised["max_lat"] - requested["max_lat"],
    }

    if app_logger is not None:
        app_logger.info(
            f"🗺️ {context} requested: "
            f"lon=[{requested['min_lon']:.3f}, {requested['max_lon']:.3f}], "
            f"lat=[{requested['min_lat']:.3f}, {requested['max_lat']:.3f}]"
        )
        app_logger.info(
            f"🗺️ {context} realised (nearest grid): "
            f"lon=[{realised['min_lon']:.3f}, {realised['max_lon']:.3f}] "
            f"({realised['n_lon']} pts, d={realised['dlon']:.3f}), "
            f"lat=[{realised['min_lat']:.3f}, {realised['max_lat']:.3f}] "
            f"({realised['n_lat']} pts, d={realised['dlat']:.3f})"
        )

    tol_lon = max(abs(dlon), 1e-6) * 1.001
    tol_lat = max(abs(dlat), 1e-6) * 1.001
    for key, tol in (
        ("min_lon", tol_lon),
        ("max_lon", tol_lon),
        ("min_lat", tol_lat),
        ("max_lat", tol_lat),
    ):
        if abs(offsets[key]) > tol:
            raise ValueError(
                f"{context}: the realised {key} ({realised[key]:.3f}) differs "
                f"from the requested value ({requested[key]:.3f}) by "
                f"{offsets[key]:.3f} degrees, which exceeds one grid spacing. "
                "The requested domain is not fully contained in the input data."
            )

    realised["offsets"] = offsets
    return realised

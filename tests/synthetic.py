#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Synthetic analytic fields for the Lorenz Energy Cycle regression tests.

The fields below are deliberately simple so that the primary-source
expressions can be evaluated independently, either in closed form or with a
short NumPy quadrature written directly from the published equations.  Nothing
here calls the LEC term implementations, so the expectations are not a
restatement of the code under test.

Design
------
Longitudes span exactly one period of the prescribed zonal wave, so that
    [cos(k*lambda)]      = 0
    [cos(k*lambda)^2]    = 1/2
    [sin(k*lambda)]      = 0
hold to machine precision under the (periodic) trapezoidal rule used by
``xarray.integrate``.  This makes the zonal means of the eddy products exact.

Latitude structure is linear in phi and vertical structure is linear in p, so
that d[T]/dphi, d[u]/dphi, d[u]/dp and d[v]/dp are represented exactly by
centred (and one-sided) finite differences.

Created for the audit-fixes branch.
"""

import argparse

import numpy as np
import pandas as pd
import xarray as xr
from metpy.units import units

# --- Box definition ---------------------------------------------------------
LON_MIN, LON_MAX, DLON = -60.0, -30.0, 0.5
LAT_MIN, LAT_MAX, DLAT = -50.0, -20.0, 0.5
LEVELS = np.array([20000.0, 30000.0, 50000.0, 70000.0, 85000.0, 100000.0])
NTIME = 3

# Zonal wavenumber: exactly one full period across the 30-degree box.
K = 12.0

# Field amplitudes / gradients
T0, A_T, GRAD_T = 280.0, 4.0, 15.0          # K, K, K per radian of latitude
U0, A_U, GRAD_U = 12.0, 6.0, 8.0            # m/s
V0, A_V, GRAD_V = 3.0, 5.0, 4.0             # m/s
W0, A_W = 0.02, 0.15                        # Pa/s
Z0, A_Z = 5.0e4, 800.0                      # m^2/s^2

DUDP = -1.5e-4      # m/s per Pa   (linear in p -> exact numerical derivative)
DVDP = 0.9e-4       # m/s per Pa


def _grids():
    lon = np.arange(LON_MIN, LON_MAX + 0.5 * DLON, DLON)
    lat = np.arange(LAT_MIN, LAT_MAX + 0.5 * DLAT, DLAT)
    lev = np.sort(LEVELS)
    time = pd.date_range("2020-06-01", periods=NTIME, freq="6h")
    return lon, lat, lev, time


def _fields(uniform_geopotential=False):
    lon, lat, lev, time = _grids()

    LON, LAT, LEV, _ = np.meshgrid(
        lon, lat, lev, np.arange(NTIME), indexing="ij"
    )
    lam = np.deg2rad(LON)
    phi = np.deg2rad(LAT)
    wave_c = np.cos(K * lam)
    wave_s = np.sin(K * lam)

    # Temperature: zonal-mean part linear in phi (so dT*/dphi = GRAD_T exactly
    # and dT*/dp = 0), plus a zonal wave (the eddy).
    T = T0 + 30.0 * (100000.0 - LEV) / 80000.0 + GRAD_T * phi + A_T * wave_c

    # Winds: zonal-mean parts linear in phi AND linear in p.
    u = U0 + GRAD_U * phi + DUDP * (LEV - 100000.0) + A_U * wave_c
    v = V0 + GRAD_V * phi + DVDP * (LEV - 100000.0) + A_V * wave_c

    omega = W0 + A_W * wave_c

    if uniform_geopotential:
        # Function of pressure and latitude only -> Phi' == 0 identically.
        geopt = Z0 * (100000.0 - LEV) / 80000.0 + 300.0 * phi
    else:
        geopt = Z0 * (100000.0 - LEV) / 80000.0 + A_Z * wave_s + 300.0 * phi

    return lon, lat, lev, time, T, u, v, omega, geopt


def make_dataset(uniform_geopotential=False):
    """Build the synthetic Dataset with the coordinates the toolkit expects."""
    lon, lat, lev, time, T, u, v, omega, geopt = _fields(uniform_geopotential)

    ds = xr.Dataset(
        {
            "T": (("longitude", "latitude", "level", "time"), T),
            "U": (("longitude", "latitude", "level", "time"), u),
            "V": (("longitude", "latitude", "level", "time"), v),
            "W": (("longitude", "latitude", "level", "time"), omega),
            "Z": (("longitude", "latitude", "level", "time"), geopt),
        },
        coords={"longitude": lon, "latitude": lat, "level": lev, "time": time},
    )
    ds = ds.assign_coords(
        rlats=np.deg2rad(ds["latitude"]),
        coslats=np.cos(np.deg2rad(ds["latitude"])),
        rlons=np.deg2rad(ds["longitude"]),
    )
    ds["level"] = ds["level"] * units("Pa")
    return ds


def variable_list():
    """Namelist-equivalent DataFrame for the synthetic dataset."""
    return pd.DataFrame(
        {
            "Variable": [
                "T", "Z", "W", "U", "V",
                "longitude", "latitude", "time", "level",
            ],
            "Units": [
                "K", "m**2/s**2", "Pa/s", "m/s", "m/s",
                "degrees", "degrees", "", "Pa",
            ],
        },
        index=[
            "Air Temperature", "Geopotential", "Omega Velocity",
            "Eastward Wind Component", "Northward Wind Component",
            "Longitude", "Latitude", "Time", "Vertical Level",
        ],
    )


def default_args(**overrides):
    base = dict(residuals=True, fixed=True, track=False, choose=False)
    base.update(overrides)
    return argparse.Namespace(**base)


def make_box(tmp_path, uniform_geopotential=False):
    """Construct a BoxData object over the synthetic dataset."""
    from src.utils.box_data import BoxData

    vlevels = tmp_path / "vlevels"
    vlevels.mkdir(exist_ok=True)
    ds = make_dataset(uniform_geopotential=uniform_geopotential)
    return BoxData(
        data=ds,
        variable_list_df=variable_list(),
        western_limit=LON_MIN,
        eastern_limit=LON_MAX,
        southern_limit=LAT_MIN,
        northern_limit=LAT_MAX,
        args=default_args(),
        results_subdirectory=str(tmp_path),
        results_subdirectory_vertical_levels=str(vlevels),
    )


class SilentLogger:
    """Minimal logger stand-in for the term classes."""

    def debug(self, *a, **k):
        pass

    def info(self, *a, **k):
        pass

    def warning(self, *a, **k):
        pass

    def error(self, *a, **k):
        pass

    def exception(self, *a, **k):
        pass


# --- Independent reference quadrature ---------------------------------------
# These reproduce the averaging operators of Brennan and Vincent (1980, p. 963)
# with plain NumPy, independently of src/utils/calc_averages.py.

def zonal_mean(field, rlons, axis=0):
    """[X] = (1/dlambda) int X dlambda, trapezoidal."""
    width = rlons[-1] - rlons[0]
    return np.trapezoid(field, rlons, axis=axis) / width


def area_mean(zonal_field, rlats, axis=0):
    """mean(X) = (1/(sin phi_n - sin phi_s)) int [X] cos(phi) dphi."""
    weight = np.cos(rlats)
    denom = np.sin(rlats[-1]) - np.sin(rlats[0])
    shape = [1] * zonal_field.ndim
    shape[axis] = weight.size
    return (
        np.trapezoid(zonal_field * weight.reshape(shape), rlats, axis=axis)
        / denom
    )

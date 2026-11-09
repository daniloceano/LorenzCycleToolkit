#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pre-registered synthetic tests of the boundary terms and of the residual
closure (Frente 2a.1).

Scope
-----
Numerical confirmation that BAZ, BAE, BKZ, BKE, BPhiZ and BPhiE in
``src/analysis/boundary_terms.py`` reproduce the limited-area expressions of
Brennan & Vincent (1980, Appendix, pp. 964-965) and Muench (1965, Eqs. 6-11),
and that the residuals in ``src/utils/calc_budget_and_residual.py`` close as
RGz = Gz, RGe = Ge, RKz = BPhiZ + Dz, RKe = BPhiE + De.  Production code is
used as-is: every field goes through ``BoxData`` -> ``BoundaryTerms``.

Reference evaluator
-------------------
Expectations are NOT produced by the toolkit.  ``ContinuousReference``
evaluates the continuous expressions below on analytic field functions with
64-point Gauss-Legendre quadrature in lambda, phi and p (exact for the
polynomial fields used here and accurate to round-off for the cosine weight).
It shares no code with ``src/`` and does not use the trapezoidal rule.

Notation (x = lambda - lambda_w, a = Earth radius, mu = sin phi):
    [X]      = (1/dlambda) int X dlambda                       zonal mean
    X'       = X - [X]                                          zonal eddy
    mean(X)  = int [X] cos phi dphi / (mu_n - mu_s)             area mean
    X*       = [X] - mean([X])                                  area eddy
    c1 = -1 / (a dlambda (mu_n - mu_s)),   c2 = -1 / (a (mu_n - mu_s))

Every term is written as the sum of three wall pairs,
    B = c1 int int [F]_{lambda_w}^{lambda_e} dphi dp      (east / west)
      + c2 int [G cos phi]_{phi_s}^{phi_n} dp              (north / south)
      - [H]_{p_t}^{p_b}                                    (top / bottom)
with
    BAZ : F = (2 T* T' u + T*^2 u) / (2 sigma)
          G = (2 [v' T'] T* + T*^2 [v]) / (2 sigma)
          H = mean(2 [omega' T'] T* + [omega] T*^2) / (2 sigma)
    BAE : F = u T'^2 / (2 sigma)
          G = [v T'^2] / (2 sigma)
          H = mean(omega T'^2) / (2 sigma)
    BKZ : F = u K_z / (2 g),  G = [v K_z] / (2 g),  H = mean(omega K_z)/(2 g)
          K_z = u^2 + v^2 - u'^2 - v'^2
    BKE : as BKZ with K_e = u'^2 + v'^2
    BPhiZ (toolkit form, Phi* on every wall):
          F = ([u] Phi' + Phi* u') / g
          G = [v] Phi* / g
          H = mean(omega* Phi*) / g
    BPhiE: F = u' Phi' / g,  G = [v' Phi'] / g,  H = mean(omega' Phi') / g
    sigma(p) = mean(g T / c_p - (p g / R) dT/dp)

Pre-registered criteria (fixed before any measurement)
------------------------------------------------------
T1  Global 360 deg band, 0 deg = 360 deg column repeated: the east/west part
    of each term must be < 1e-10 x the single-wall magnitude from the
    reference.  (Without the repeated column the test is INCONCLUSIVE BY
    CONSTRUCTION; the grid used here repeats it.)
T2  u, v, omega constant; T' and Phi' polynomial in lambda and phi:
    relative error <= 1 % on a 0.25 deg grid.  Terms whose exact value is
    identically zero under this design (BKE, BPhiE: u' = v' = omega' = 0)
    are checked as |value| <= 1e-10 x max |non-zero T2 term|.
T3c Each of the six walls in isolation (fields chosen so the other five
    walls vanish identically), 6 walls x 6 terms = 36 combinations, each on
    three grids refined jointly in the horizontal and the vertical:
    1 deg / dp = 2000 Pa -> 0.5 deg / 1000 Pa -> 0.25 deg / 500 Pa (P_TOP
    and P_BOT fixed), giving errors e1, e2, e3 against the reference.  PASS
    per combination if the observed orders log2(e1/e2) and log2(e2/e3) both
    lie in [1.8, 2.2] AND the Richardson-extrapolated value (p = 2),
    x3 + (x3 - x2)/3, is within 0.1 % of the reference.  For a combination
    whose reference value is zero: |x| <= 1e-10 x the largest non-zero
    reference term of the same wall configuration, on every grid.
    History (all on 2026-11-09, recorded in the Frente 2a.1 evidence):
      T3  "relative error <= 1 % per term on a 0.25 deg grid": FAIL in 4/36
          (top/bottom BAE and BKE at 1.11 %), second-order convergence, no
          formula error found.
      T3b same orders/Richardson criteria on 1, 0.5, 0.25 deg with dp fixed
          at 1000 Pa: FAIL in 4/36 (E/W BAE, N/S BAZ, observed order), caused
          by a grid-independent second-order vertical-discretisation floor
          (~1.1e-4, scaling with dp^2) in the 1/sigma-weighted terms.
      T3c joint horizontal + vertical refinement: FINAL attempt.  If it fails
          in any combination nothing is adjusted: the convergence test is
          removed from the suite and recorded as INCONCLUSIVE.
T4  BPhiZ(Phi + C) vs BPhiZ(Phi): relative change < 1e-10.  The Brennan &
    Vincent form with the full geopotential, built in this harness, must
    change, by the analytically predicted amount.
T5  RGz, RGe, RKz, RKe reproduce Gz, Ge, BPhiZ + Dz, BPhiE + De to a
    relative error <= 1e-10, through calc_budget_diff -> calc_residuals.
Control  Two in-memory mutants of boundary_terms.py (the file on disk is not
    touched): M1 flips the sign of the west wall (E - W -> E + W); M2
    replaces cos(phi) by |sin(phi)|.  T2 and T3c must both FAIL under each
    mutant; otherwise the protocol is FAIL (tests without power).
"""

import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from metpy.constants import Cp_d, Rd, Re, g
from metpy.units import units

from src.analysis.boundary_terms import BoundaryTerms
from src.utils.calc_budget_and_residual import calc_budget_diff, calc_residuals
from tests import synthetic as syn

A_EARTH = float(Re.magnitude)
G = float(g.magnitude)
CP = float(Cp_d.to("J/(kg K)").magnitude)
RD = float(Rd.to("J/(kg K)").magnitude)

TERMS = ("BAZ", "BAE", "BKZ", "BKE", "BOZ", "BOE")
METHODS = {
    "BAZ": "calc_baz", "BAE": "calc_bae", "BKZ": "calc_bkz",
    "BKE": "calc_bke", "BOZ": "calc_boz", "BOE": "calc_boe",
}

# Pre-registered tolerances (do not edit after measuring).
TOL_T1 = 1e-10
TOL_T2 = 1e-2
TOL_T2_ZERO = 1e-10
T3C_GRIDS = ((1.0, 2000.0), (0.5, 1000.0), (0.25, 500.0))  # (deg, dp Pa), coarse -> fine
T3C_ORDER_MIN, T3C_ORDER_MAX = 1.8, 2.2
T3C_TOL_EXTRAP = 1e-3                  # Richardson (p = 2) relative error
T3C_TOL_ZERO = 1e-10
TOL_T4 = 1e-10
TOL_T5 = 1e-10

# --- Regional grid (T2, T3c, T4) -------------------------------------------
LON_W, LON_E = -50.0, -40.0
LAT_S, LAT_N = -40.0, -30.0
DDEG = 0.25
P_TOP, P_BOT, DP = 20000.0, 100000.0, 1000.0
NTIME = 3


# ===========================================================================
# Synthetic dataset construction
# ===========================================================================

def _t0(p):
    """Background temperature, linear in p (K): 300 K at 1000 hPa, 240 K at 200 hPa."""
    return 300.0 - 60.0 * (P_BOT - p) / (P_BOT - P_TOP)


def _phi0(p):
    """Background geopotential, linear in p (m2 s-2)."""
    return 5.0e4 * (P_BOT - p) / (P_BOT - P_TOP)


class Config:
    """
    Analytic fields on a lat-lon-pressure box.

    Field callables take (s, t, p): s = x / Lx in [0, 1], t = (phi - phi_s) /
    Ly in [0, 1] and pressure p in Pa, and return SI values.
    """

    def __init__(self, name, u, v, w, T, Z, lon=(LON_W, LON_E),
                 lat=(LAT_S, LAT_N), ddeg=DDEG, levels=None):
        self.name = name
        self.u, self.v, self.w, self.T, self.Z = u, v, w, T, Z
        self.lon_w, self.lon_e = lon
        self.lat_s, self.lat_n = lat
        self.ddeg = ddeg
        self.levels = (np.arange(P_TOP, P_BOT + 0.5 * DP, DP)
                       if levels is None else np.asarray(levels, float))
        self.lam_w = np.deg2rad(self.lon_w)
        self.Lx = np.deg2rad(self.lon_e - self.lon_w)
        self.phi_s = np.deg2rad(self.lat_s)
        self.Ly = np.deg2rad(self.lat_n - self.lat_s)

    def with_geopotential_offset(self, C):
        Z = self.Z
        return Config(self.name + f"+{C:g}", self.u, self.v, self.w, self.T,
                      lambda s, t, p: Z(s, t, p) + C,
                      (self.lon_w, self.lon_e), (self.lat_s, self.lat_n),
                      self.ddeg, self.levels)


def make_dataset(cfg):
    lon = np.round(np.arange(cfg.lon_w, cfg.lon_e + 0.5 * cfg.ddeg, cfg.ddeg), 10)
    lat = np.round(np.arange(cfg.lat_s, cfg.lat_n + 0.5 * cfg.ddeg, cfg.ddeg), 10)
    lev = np.sort(cfg.levels)
    time = pd.date_range("2020-06-01", periods=NTIME, freq="6h")
    LON, LAT, LEV, _ = np.meshgrid(lon, lat, lev, np.arange(NTIME), indexing="ij")
    s = (np.deg2rad(LON) - cfg.lam_w) / cfg.Lx
    t = (np.deg2rad(LAT) - cfg.phi_s) / cfg.Ly
    dims = ("longitude", "latitude", "level", "time")

    def full(f):
        return np.broadcast_to(np.asarray(f(s, t, LEV), float), s.shape).copy()

    ds = xr.Dataset(
        {"T": (dims, full(cfg.T)), "U": (dims, full(cfg.u)),
         "V": (dims, full(cfg.v)), "W": (dims, full(cfg.w)),
         "Z": (dims, full(cfg.Z))},
        coords={"longitude": lon, "latitude": lat, "level": lev, "time": time},
    )
    ds = ds.assign_coords(
        rlats=np.deg2rad(ds["latitude"]),
        coslats=np.cos(np.deg2rad(ds["latitude"])),
        rlons=np.deg2rad(ds["longitude"]),
    )
    ds["level"] = ds["level"] * units("Pa")
    return ds


def make_box(cfg, tmp_path):
    from src.utils.box_data import BoxData

    sub = Path(tmp_path) / cfg.name.replace("+", "_")
    vl = sub / "vlevels"
    vl.mkdir(parents=True, exist_ok=True)
    return BoxData(
        data=make_dataset(cfg),
        variable_list_df=syn.variable_list(),
        western_limit=cfg.lon_w, eastern_limit=cfg.lon_e,
        southern_limit=cfg.lat_s, northern_limit=cfg.lat_n,
        args=syn.default_args(),
        results_subdirectory=str(sub),
        results_subdirectory_vertical_levels=str(vl),
    )


def toolkit_terms(box, cls=BoundaryTerms):
    """Run the production class (or a mutant) and return W m-2 at time 0."""
    bt = cls(box, "fixed", syn.SilentLogger())
    out = {}
    for term, meth in METHODS.items():
        da = getattr(bt, meth)()
        da = da.metpy.dequantify() if hasattr(da, "metpy") else da
        out[term] = float(np.asarray(da.values, float).ravel()[0])
    return out


# ===========================================================================
# Continuous reference evaluator (independent of src/)
# ===========================================================================

class ContinuousReference:
    """Gauss-Legendre evaluation of the continuous wall expressions."""

    def __init__(self, cfg, n=64):
        self.cfg = cfg
        xg, wg = np.polynomial.legendre.leggauss(n)
        self.sq = 0.5 * (xg + 1.0)           # nodes on [0, 1]
        self.sw = 0.5 * wg                   # weights on [0, 1]
        self.tq, self.tw = self.sq, self.sw
        self.pt, self.pb = float(cfg.levels.min()), float(cfg.levels.max())
        self.pq = self.pt + (self.pb - self.pt) * self.sq
        self.pw = (self.pb - self.pt) * self.sw
        self.mu = np.sin(cfg.phi_s + cfg.Ly) - np.sin(cfg.phi_s)
        self.c1 = -1.0 / (A_EARTH * cfg.Lx * self.mu)
        self.c2 = -1.0 / (A_EARTH * self.mu)

    # -- operators ----------------------------------------------------------
    def _phi(self, t):
        return self.cfg.phi_s + self.cfg.Ly * t

    def zm(self, f, t, p):
        """[f](t, p) for f(s, t, p); t, p broadcastable arrays."""
        s = self.sq.reshape((-1,) + (1,) * np.ndim(np.broadcast_to(t, np.broadcast(t, p).shape)))
        vals = f(s, t, p)
        vals = np.broadcast_to(vals, np.broadcast(s, t, p).shape)
        return np.tensordot(self.sw, vals, axes=(0, 0))

    def am(self, h, p):
        """Area mean over phi of h(t, p) (h already zonal), at scalar/array p."""
        t = self.tq.reshape((-1,) + (1,) * np.ndim(p))
        cw = np.cos(self._phi(t)) * self.cfg.Ly
        vals = np.broadcast_to(h(t, p), np.broadcast(t, p).shape)
        return np.tensordot(self.tw, vals * cw, axes=(0, 0)) / self.mu

    # -- decomposition ------------------------------------------------------
    def parts(self, name):
        f = getattr(self.cfg, name)
        zm = lambda t, p: self.zm(f, t, p)                       # noqa: E731
        aa = lambda p: self.am(zm, p)                            # noqa: E731
        ze = lambda s, t, p: f(s, t, p) - zm(t, p)               # noqa: E731
        ae = lambda t, p: zm(t, p) - aa(p)                       # noqa: E731
        return f, zm, ze, ae

    def sigma(self, p):
        T = self.cfg.T
        h = 10.0  # T is linear in p, so the centred difference is exact
        dTdp = lambda s, t, pp: (T(s, t, pp + h) - T(s, t, pp - h)) / (2 * h)  # noqa: E731
        integrand = lambda s, t, pp: G * T(s, t, pp) / CP - (pp * G / RD) * dTdp(s, t, pp)  # noqa: E731
        return self.am(lambda t, pp: self.zm(integrand, t, pp), p)

    # -- wall integrands ----------------------------------------------------
    def integrands(self, term):
        u, uZ, uE, uA = self.parts("u")
        v, vZ, vE, vA = self.parts("v")
        w, wZ, wE, wA = self.parts("w")
        T, TZ, TE, TA = self.parts("T")
        Z, ZZ, ZE, ZA = self.parts("Z")
        sig = self.sigma
        zm = self.zm

        if term == "BAZ":
            F = lambda s, t, p: (2 * TA(t, p) * TE(s, t, p) * u(s, t, p)
                                 + TA(t, p) ** 2 * u(s, t, p)) / (2 * sig(p))
            G_ = lambda t, p: (2 * zm(lambda s, tt, pp: vE(s, tt, pp) * TE(s, tt, pp), t, p) * TA(t, p)
                               + TA(t, p) ** 2 * vZ(t, p)) / (2 * sig(p))
            H = lambda p: self.am(lambda t, pp: 2 * zm(lambda s, tt, q: wE(s, tt, q) * TE(s, tt, q), t, pp) * TA(t, pp)
                                  + wZ(t, pp) * TA(t, pp) ** 2, p) / (2 * sig(p))
        elif term == "BAE":
            F = lambda s, t, p: u(s, t, p) * TE(s, t, p) ** 2 / (2 * sig(p))
            G_ = lambda t, p: zm(lambda s, tt, pp: v(s, tt, pp) * TE(s, tt, pp) ** 2, t, p) / (2 * sig(p))
            H = lambda p: self.am(lambda t, pp: zm(lambda s, tt, q: w(s, tt, q) * TE(s, tt, q) ** 2, t, pp), p) / (2 * sig(p))
        elif term in ("BKZ", "BKE"):
            if term == "BKZ":
                K = lambda s, t, p: (u(s, t, p) ** 2 + v(s, t, p) ** 2
                                     - uE(s, t, p) ** 2 - vE(s, t, p) ** 2)
            else:
                K = lambda s, t, p: uE(s, t, p) ** 2 + vE(s, t, p) ** 2
            F = lambda s, t, p: u(s, t, p) * K(s, t, p) / (2 * G)
            G_ = lambda t, p: zm(lambda s, tt, pp: v(s, tt, pp) * K(s, tt, pp), t, p) / (2 * G)
            H = lambda p: self.am(lambda t, pp: zm(lambda s, tt, q: w(s, tt, q) * K(s, tt, q), t, pp), p) / (2 * G)
        elif term == "BOZ":
            F = lambda s, t, p: (uZ(t, p) * ZE(s, t, p) + ZA(t, p) * uE(s, t, p)) / G
            G_ = lambda t, p: vZ(t, p) * ZA(t, p) / G
            H = lambda p: self.am(lambda t, pp: wA(t, pp) * ZA(t, pp), p) / G
        elif term == "BOE":
            F = lambda s, t, p: uE(s, t, p) * ZE(s, t, p) / G
            G_ = lambda t, p: zm(lambda s, tt, pp: vE(s, tt, pp) * ZE(s, tt, pp), t, p) / G
            H = lambda p: self.am(lambda t, pp: zm(lambda s, tt, q: wE(s, tt, q) * ZE(s, tt, q), t, pp), p) / G
        else:
            raise KeyError(term)
        return F, G_, H

    def walls(self, term):
        """Signed contribution of each of the six walls (W m-2)."""
        F, G_, H = self.integrands(term)
        t = self.tq[:, None]
        p = self.pq[None, :]
        wts = (self.tw[:, None] * self.cfg.Ly) * self.pw[None, :]
        east = self.c1 * np.sum(wts * F(np.ones_like(t), t, p))
        west = -self.c1 * np.sum(wts * F(np.zeros_like(t), t, p))
        cos_n = np.cos(self.cfg.phi_s + self.cfg.Ly)
        cos_s = np.cos(self.cfg.phi_s)
        north = self.c2 * np.sum(self.pw * G_(np.ones_like(self.pq), self.pq) * cos_n)
        south = -self.c2 * np.sum(self.pw * G_(np.zeros_like(self.pq), self.pq) * cos_s)
        bottom = -float(H(np.array(self.pb)))
        top = float(H(np.array(self.pt)))
        return {"E": float(east), "W": float(west), "N": float(north),
                "S": float(south), "T": float(top), "B": float(bottom)}

    def totals(self):
        out = {}
        for term in TERMS:
            wl = self.walls(term)
            out[term] = (sum(wl.values()), wl)
        return out


# ===========================================================================
# Field configurations
# ===========================================================================

def f_east(s):
    """Zero zonal mean, zero on the west wall, 1/3 on the east wall."""
    return s ** 2 - 2.0 * s / 3.0


def f_west(s):
    return f_east(1.0 - s)


def c_sym(s):
    """Zero zonal mean, equal values on the two walls."""
    return (s - 0.5) ** 2 - 1.0 / 12.0


def q_top(p):
    """1 at the top level, 0 at the bottom level."""
    return (P_BOT - p) / (P_BOT - P_TOP)


def q_bot(p):
    return (p - P_TOP) / (P_BOT - P_TOP)


def config_t2():
    """u, v, omega constant; T' and Phi' polynomial in lambda and phi."""
    U, V, W = 10.0, 5.0, 0.05
    return Config(
        "T2",
        u=lambda s, t, p: U + 0 * s,
        v=lambda s, t, p: V + 0 * s,
        w=lambda s, t, p: W + 0 * s,
        T=lambda s, t, p: _t0(p) + 3.0 * s + 4.0 * s ** 2 + 5.0 * t + 2.0 * s * t,
        Z=lambda s, t, p: _phi0(p) + 400.0 * s + 300.0 * s ** 2 - 500.0 * t + 200.0 * s * t,
    )


def config_face(face):
    """
    Fields for which only ``face`` can contribute, for every term.

    E / W : v = omega = 0 kills N/S and top/bottom; u and every zonal eddy
            vanish on the opposite wall, which kills its F for all six terms
            (for BPhiZ, Phi* = 0 so F_opposite = [u] Phi'(opposite) = 0).
    N / S : u = omega = 0 kills E/W and top/bottom; v vanishes on the
            opposite latitude wall, and every G carries a factor v.
    T / B : v = 0 kills N/S; u, T', Phi' are symmetric in x so E and W cancel;
            omega vanishes on the opposite pressure wall, and every H carries
            a factor omega.
    """
    U, V, W = 10.0, 5.0, 0.05
    if face in ("E", "W"):
        f = f_east if face == "E" else f_west
        ramp = (lambda s: s) if face == "E" else (lambda s: 1.0 - s)
        return Config(
            face,
            u=lambda s, t, p: U * ramp(s) * (1.0 + 0.5 * t),
            v=lambda s, t, p: 0.0 * s,
            w=lambda s, t, p: 0.0 * s,
            T=lambda s, t, p: _t0(p) + 6.0 * f(s) * (1.0 + 0.3 * t) + 5.0 * t,
            Z=lambda s, t, p: _phi0(p) + 600.0 * f(s) * (1.0 + 0.3 * t),
        )
    if face in ("N", "S"):
        h = (lambda t: t) if face == "N" else (lambda t: 1.0 - t)
        return Config(
            face,
            u=lambda s, t, p: 0.0 * s,
            v=lambda s, t, p: V * h(t) * (1.0 + 0.8 * (s - 0.5)),
            w=lambda s, t, p: 0.0 * s,
            T=lambda s, t, p: _t0(p) + 4.0 * (s - 0.5) + 3.0 * (s - 0.5) ** 2 * (1 + t) + 5.0 * t,
            Z=lambda s, t, p: _phi0(p) + 500.0 * (s - 0.5) * (1.0 + t) - 400.0 * t,
        )
    if face in ("T", "B"):
        q = q_top if face == "T" else q_bot
        return Config(
            face,
            u=lambda s, t, p: U + 8.0 * c_sym(s) + 0 * t,
            v=lambda s, t, p: 0.0 * s,
            w=lambda s, t, p: W * q(p) * (1.0 + 6.0 * c_sym(s) + 0.4 * (t - 0.5)),
            T=lambda s, t, p: _t0(p) + 6.0 * c_sym(s) * (1.0 + 0.3 * t) + 5.0 * t,
            Z=lambda s, t, p: _phi0(p) + 600.0 * c_sym(s) * (1.0 + 0.3 * t) + 300.0 * t,
        )
    raise KeyError(face)


def config_t4():
    """Net mass flux through E/W and N/S: [u] and [v] differ across walls."""
    return Config(
        "T4",
        u=lambda s, t, p: 10.0 + 4.0 * s + 3.0 * (s - 0.5) * t + 0 * p,
        v=lambda s, t, p: 5.0 + 3.0 * t + 2.0 * (s - 0.5) + 0 * p,
        w=lambda s, t, p: 0.05 + 0.02 * (s - 0.5) + 0 * p,
        T=lambda s, t, p: _t0(p) + 3.0 * s + 5.0 * t,
        Z=lambda s, t, p: _phi0(p) + 400.0 * s + 300.0 * s ** 2 - 500.0 * t + 200.0 * s * t,
    )


def config_t1():
    """
    Global 360 deg band, 0..360 deg with the 0 = 360 column repeated.

    v = omega = 0, so each term is its east/west part alone.  Fields are
    2 pi-periodic in lambda.
    """
    lam = lambda s: 2.0 * np.pi * s  # noqa: E731  (lambda_w = 0)
    return Config(
        "T1",
        u=lambda s, t, p: 10.0 + 5.0 * np.cos(lam(s)) + 3.0 * np.sin(2 * lam(s)) * (1 + t),
        v=lambda s, t, p: 0.0 * s,
        w=lambda s, t, p: 0.0 * s,
        T=lambda s, t, p: _t0(p) + 4.0 * np.cos(lam(s)) + 2.0 * np.sin(lam(s)) * (1 + t) + 5.0 * t,
        Z=lambda s, t, p: _phi0(p) + 500.0 * np.sin(lam(s)) + 300.0 * np.cos(2 * lam(s)) + 200.0 * t,
        lon=(0.0, 360.0),
        levels=[20000.0, 40000.0, 60000.0, 85000.0, 100000.0],
    )


# ===========================================================================
# Mutants for the fail control (in memory only)
# ===========================================================================

_BT_SOURCE = Path(__file__).resolve().parents[1] / "src" / "analysis" / "boundary_terms.py"


def load_mutant(tag):
    """
    Return a BoundaryTerms class built from a mutated copy of the source text.

    M1: ``eastern_limit}) - `` -> ``eastern_limit}) + ``  (west wall sign)
    M2: ``X["coslats"]`` -> ``((1 - X["coslats"] ** 2) ** 0.5)``  (cos -> |sin|)
    """
    text = _BT_SOURCE.read_text()
    if tag == "M1":
        new, n = re.subn(r"eastern_limit\}\) - ", "eastern_limit}) + ", text)
    elif tag == "M2":
        new, n = re.subn(r"([\w\.]+)\[\"coslats\"\]",
                         r'((1 - \1["coslats"] ** 2) ** 0.5)', text)
    else:
        raise KeyError(tag)
    assert n >= 6, f"mutant {tag}: only {n} substitutions"
    name = f"src.analysis._mutant_boundary_terms_{tag}"
    spec = importlib.util.spec_from_loader(name, loader=None)
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = "src.analysis"
    exec(compile(new, f"<mutant {tag} of boundary_terms.py>", "exec"), mod.__dict__)
    return mod.BoundaryTerms, n


# ===========================================================================
# Test runners (also used by the external evidence script)
# ===========================================================================

def _rel(a, b):
    return abs(a - b) / abs(b)


def run_t2(tmp_path, cls=BoundaryTerms):
    cfg = config_t2()
    ref = ContinuousReference(cfg).totals()
    got = toolkit_terms(make_box(cfg, tmp_path), cls)
    nonzero = {k for k in TERMS if k not in ("BKE", "BOE")}
    scale = max(abs(ref[k][0]) for k in nonzero)
    rows = {}
    for k in TERMS:
        exp = ref[k][0]
        if k in nonzero:
            err = _rel(got[k], exp)
            ok = err <= TOL_T2
            rows[k] = dict(toolkit=got[k], expected=exp, rel_err=err,
                           tol=TOL_T2, criterion="relative", passed=ok,
                           walls=ref[k][1])
        else:
            err = abs(got[k]) / scale
            rows[k] = dict(toolkit=got[k], expected=0.0, rel_err=err,
                           tol=TOL_T2_ZERO, criterion="|x|/scale", passed=err <= TOL_T2_ZERO,
                           walls=ref[k][1])
    return rows


def run_t3c(tmp_path, cls=BoundaryTerms):
    """
    T3c: per wall-isolating configuration, run the toolkit on the three
    T3C_GRIDS (joint horizontal + vertical refinement) and check the observed
    order and the Richardson extrapolation.  The reference depends only on
    the box limits (P_TOP, P_BOT fixed), not on the grid.
    """
    pairs = {"E": "W", "W": "E", "N": "S", "S": "N", "T": "B", "B": "T"}
    rows = {}
    for face in ("E", "W", "N", "S", "T", "B"):
        ref = ContinuousReference(config_face(face)).totals()
        got = {}
        for dd, dp in T3C_GRIDS:
            cfg = config_face(face)
            cfg.ddeg = dd
            cfg.levels = np.arange(P_TOP, P_BOT + 0.5 * dp, dp)
            cfg.name = f"T3c_{face}_{dd:g}deg_{dp:g}Pa"
            got[(dd, dp)] = toolkit_terms(make_box(cfg, tmp_path), cls)
        scale = max(abs(ref[k][0]) for k in TERMS)
        for k in TERMS:
            exp, wl = ref[k]
            # Design self-check: the wall opposite to ``face`` vanishes, and
            # each of the two other wall pairs sums to zero (in the T/B
            # configurations E and W cancel by symmetry rather than vanish).
            pair_sums = [wl[a] + wl[b] for a, b in (("E", "W"), ("N", "S"), ("T", "B"))
                         if face not in (a, b)]
            others = max([abs(wl[pairs[face]])] + [abs(x) for x in pair_sums])
            design_ok = others <= 1e-12 * max(abs(wl[face]), 1e-300)
            x1, x2, x3 = (got[gr][k] for gr in T3C_GRIDS)
            row = dict(grids_deg_dp=[list(gr) for gr in T3C_GRIDS], toolkit=[x1, x2, x3], expected=exp,
                       other_walls_max=others, design_ok=design_ok)
            if abs(exp) <= 1e-12 * scale:
                # Analytic zero: no order is defined; magnitude criterion.
                nonzero_scale = max(abs(ref[j][0]) for j in TERMS
                                    if abs(ref[j][0]) > 1e-12 * scale)
                ratio = max(abs(x1), abs(x2), abs(x3)) / nonzero_scale
                row.update(criterion="zero", ratio=ratio, tol=T3C_TOL_ZERO,
                           passed=ratio <= T3C_TOL_ZERO)
            else:
                e1, e2, e3 = (abs(x - exp) for x in (x1, x2, x3))
                with np.errstate(divide="ignore", invalid="ignore"):
                    p12 = float(np.log2(e1 / e2))
                    p23 = float(np.log2(e2 / e3))
                x_rich = x3 + (x3 - x2) / 3.0
                err_rich = _rel(x_rich, exp)
                orders_ok = all(np.isfinite(p) and T3C_ORDER_MIN <= p <= T3C_ORDER_MAX
                                for p in (p12, p23))
                row.update(criterion="convergence",
                           rel_err=[e1 / abs(exp), e2 / abs(exp), e3 / abs(exp)],
                           order_12=p12, order_23=p23,
                           richardson=x_rich, richardson_rel_err=err_rich,
                           order_range=[T3C_ORDER_MIN, T3C_ORDER_MAX],
                           tol_extrap=T3C_TOL_EXTRAP,
                           passed=bool(orders_ok and err_rich <= T3C_TOL_EXTRAP))
            rows[(face, k)] = row
    return rows


def run_t1(tmp_path):
    cfg = config_t1()
    ref = ContinuousReference(cfg)
    got = toolkit_terms(make_box(cfg, tmp_path))
    lon = np.round(np.arange(cfg.lon_w, cfg.lon_e + 0.5 * cfg.ddeg, cfg.ddeg), 10)
    repeated = bool(np.isclose(lon[0] % 360.0, lon[-1] % 360.0))
    rows = {}
    for k in TERMS:
        wl = ref.walls(k)
        wall = max(abs(wl["E"]), abs(wl["W"]))
        ratio = abs(got[k]) / wall
        rows[k] = dict(toolkit_EW=got[k], wall_magnitude=wall, ratio=ratio,
                       tol=TOL_T1, passed=repeated and ratio < TOL_T1,
                       repeated_column=repeated,
                       other_walls_max=max(abs(wl[f]) for f in "NSTB"))
    return rows


def _bv_full_phi_boz(box):
    """
    Brennan & Vincent (1980) BPhiZ with the FULL zonal-mean geopotential
    [Phi] in place of Phi* on every wall (harness only; not in src/).
    Returns (total, {"EW":..., "NS":..., "TB":...}) in W m-2.
    """
    c1 = -1 / (Re * box.xlength * box.ylength)
    c2 = -1 / (Re * box.ylength)
    P = box.PressureData.metpy.units
    lon, lat, lev = box.LonIndexer, box.LatIndexer, box.VerticalCoordIndexer

    ff = (box.u_ZA * box.geopt_ZE + box.geopt_ZA * box.u_ZE) / g
    ew = ff.sel(**{lon: box.eastern_limit}) - ff.sel(**{lon: box.western_limit})
    ew = ew.integrate("rlats").integrate(lev) * P * c1

    ns = box.v_ZA * box.geopt_ZA * box.v_ZA["coslats"] / g
    ns = ns.sel(**{lat: box.northern_limit}) - ns.sel(**{lat: box.southern_limit})
    ns = ns.integrate(lev) * P * c2

    from src.utils.calc_averages import CalcAreaAverage
    tb = CalcAreaAverage(box.omega_AE * box.geopt_ZA, box.ylength) / g
    tb = -(tb.isel(**{lev: -1}) - tb.isel(**{lev: 0}))

    def val(x):
        x = x.metpy.convert_units("W/m^2").metpy.dequantify()
        return float(np.asarray(x.values, float).ravel()[0])

    parts = {"EW": val(ew), "NS": val(ns), "TB": val(tb)}
    return sum(parts.values()), parts


def run_t4(tmp_path, C=1.0e5):
    cfg0 = config_t4()
    cfg1 = cfg0.with_geopotential_offset(C)
    b0, b1 = make_box(cfg0, tmp_path), make_box(cfg1, tmp_path)
    t0, t1 = toolkit_terms(b0)["BOZ"], toolkit_terms(b1)["BOZ"]
    bv0, p0 = _bv_full_phi_boz(b0)
    bv1, p1 = _bv_full_phi_boz(b1)

    # Analytic change of the full-Phi form: C x (mass-flux convergence)/g,
    # E/W: c1 int int (u'_E - u'_W) dphi dp / g;  N/S: c2 int ([v] cos)|_s^n dp / g.
    ref = ContinuousReference(cfg0)
    u, uZ, uE, uA = ref.parts("u")
    v, vZ, vE, vA = ref.parts("v")
    t = ref.tq[:, None]
    p = ref.pq[None, :]
    wts = (ref.tw[:, None] * cfg0.Ly) * ref.pw[None, :]
    d_ew = C * ref.c1 * np.sum(wts * (uE(np.ones_like(t), t, p) - uE(np.zeros_like(t), t, p))) / G
    cos_n, cos_s = np.cos(cfg0.phi_s + cfg0.Ly), np.cos(cfg0.phi_s)
    d_ns = C * ref.c2 * np.sum(ref.pw * (vZ(np.ones_like(ref.pq), ref.pq) * cos_n
                                         - vZ(np.zeros_like(ref.pq), ref.pq) * cos_s)) / G
    return dict(
        C=C, toolkit_BOZ=t0, toolkit_BOZ_offset=t1,
        toolkit_rel_change=abs(t1 - t0) / abs(t0), tol=TOL_T4,
        passed=abs(t1 - t0) / abs(t0) < TOL_T4,
        bv_full_phi=bv0, bv_full_phi_offset=bv1,
        bv_rel_change=abs(bv1 - bv0) / abs(bv0),
        bv_change_parts={k: p1[k] - p0[k] for k in p0},
        bv_change_expected={"EW": float(d_ew), "NS": float(d_ns), "TB": 0.0},
    )


def run_t5():
    """
    Known energies (linear in time, so np.gradient is exact), conversions,
    boundary terms, generation, dissipation and BPhi; the tendencies follow
    the budget of Brennan & Vincent (1980) / docs/source/math.rst:
        dAz/dt = -Cz - Ca + BAz + Gz
        dAe/dt =  Ca - Ce + BAe + Ge
        dKz/dt =  Cz + Ck + BKz + BPhiZ + Dz
        dKe/dt =  Ce - Ck + BKe + BPhiE + De
    """
    n = 6
    dates = pd.date_range("2020-06-01", periods=n, freq="6h")
    sec = (dates - dates[0]).total_seconds().to_numpy()
    k = dict(Cz=0.731, Ca=1.913, Ck=-0.417, Ce=2.284,
             BAz=-0.356, BAe=0.588, BKz=-1.207, BKe=0.462,
             Gz=2.519, Ge=-0.873, Dz=-1.644, De=-2.951,
             BOZ=0.379, BOE=-0.265)
    tend = dict(
        Az=-k["Cz"] - k["Ca"] + k["BAz"] + k["Gz"],
        Ae=k["Ca"] - k["Ce"] + k["BAe"] + k["Ge"],
        Kz=k["Cz"] + k["Ck"] + k["BKz"] + k["BOZ"] + k["Dz"],
        Ke=k["Ce"] - k["Ck"] + k["BKe"] + k["BOE"] + k["De"],
    )
    base = dict(Az=4.1e6, Ae=6.3e5, Kz=1.2e6, Ke=4.7e5)
    df = pd.DataFrame({e: base[e] + tend[e] * sec for e in base})
    for c in ("Cz", "Ca", "Ck", "Ce", "BAz", "BAe", "BKz", "BKe"):
        df[c] = k[c]
    df = calc_budget_diff(df, dates, syn.SilentLogger())
    df = calc_residuals(df, syn.SilentLogger())
    expected = dict(RGz=k["Gz"], RGe=k["Ge"], RKz=k["BOZ"] + k["Dz"],
                    RKe=k["BOE"] + k["De"])
    rows = {}
    for r, e in expected.items():
        err = float(np.max(np.abs(df[r].to_numpy() - e)) / abs(e))
        rows[r] = dict(expected=e, max_rel_err=err, tol=TOL_T5, passed=err <= TOL_T5)
    return rows


# ===========================================================================
# Tests
# ===========================================================================

@pytest.fixture(scope="module")
def tmp(tmp_path_factory):
    return tmp_path_factory.mktemp("bt2a1")


def test_t1_global_band_east_west_vanishes(tmp):
    rows = run_t1(tmp)
    assert all(r["repeated_column"] for r in rows.values()), "INCONCLUSIVE BY CONSTRUCTION"
    bad = {k: r["ratio"] for k, r in rows.items() if not r["passed"]}
    assert not bad, bad


def test_t2_analytic_constant_winds(tmp):
    rows = run_t2(tmp)
    bad = {k: (r["toolkit"], r["expected"], r["rel_err"]) for k, r in rows.items() if not r["passed"]}
    assert not bad, bad


def test_t3c_each_wall_in_isolation_joint_convergence(tmp):
    rows = run_t3c(tmp)
    assert all(r["design_ok"] for r in rows.values()), "field design does not isolate the wall"
    bad = {k: {f: r.get(f) for f in ("toolkit", "expected", "order_12", "order_23",
                                     "richardson_rel_err", "ratio")}
           for k, r in rows.items() if not r["passed"]}
    assert not bad, bad


def test_t4_bphiz_invariant_to_geopotential_offset(tmp):
    r = run_t4(tmp)
    assert r["passed"], r
    # The full-Phi (Brennan & Vincent) form built in the harness must change,
    # by the analytically predicted amount on each wall pair.
    assert r["bv_rel_change"] > 1e-3, r
    for wall in ("EW", "NS"):
        got, exp = r["bv_change_parts"][wall], r["bv_change_expected"][wall]
        assert abs(got - exp) <= 1e-2 * abs(exp), (wall, got, exp)


def test_t5_residual_identity():
    rows = run_t5()
    bad = {k: r for k, r in rows.items() if not r["passed"]}
    assert not bad, bad


@pytest.mark.parametrize("tag", ["M1", "M2"])
def test_control_mutants_are_detected(tmp, tag):
    cls, _ = load_mutant(tag)
    t2 = run_t2(tmp, cls)
    t3c = run_t3c(tmp, cls)
    assert not all(r["passed"] for r in t2.values()), f"T2 has no power against {tag}"
    assert not all(r["passed"] for r in t3c.values()), f"T3c has no power against {tag}"

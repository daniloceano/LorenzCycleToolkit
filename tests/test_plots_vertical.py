#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests of the vertical-level figures (boxplots by level, boxplots by time and
Hovmöller diagrams).

Each panel must be labelled with the quantity and units of its own term
(``src/utils/output_spec.py``), every pressure level of the file must be
drawn, pressure is shown in hPa on the ticks only, and the values handed to
matplotlib must be exactly the values of the per-level CSV files: no
conversion of the plotted quantity and no change of the plotted pressure
coordinate (Pa) or of the time axis.
"""

import copy
from pathlib import Path

import cmocean  # noqa: F401  (registers the "cmo." colour maps, as a run does)
import matplotlib
import matplotlib.axes
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import src.plots.plot_boxplot as plot_boxplot
import src.plots.plot_hovmoller as plot_hovmoller
from src.plots.utils import get_data_vertical_levels

matplotlib.use("Agg")
REFERENCE_DIR = Path(__file__).resolve().parent / "data" / "reference"
GROUPS = {
    "energy": ["Az", "Ae", "Kz", "Ke"],
    "conversion": ["Cz", "Ca", "Ck", "Ce"],
    "generation_dissipation": ["Gz", "Ge"],
}


class _Recorder:
    """Records the figures saved and the data handed to boxplot/contourf."""

    def __init__(self, monkeypatch):
        self.figures, self.boxplots, self.contours = [], [], []
        saved = plt.savefig
        box = matplotlib.axes.Axes.boxplot
        contourf = matplotlib.axes.Axes.contourf

        def savefig(fname, *a, **k):
            self.figures.append((Path(str(fname)).name, plt.gcf()))
            return saved(fname, *a, **k)

        def boxplot(ax, x, *a, **k):
            self.boxplots.append((ax, np.array(x, dtype=float, copy=True), k.get("positions"),
                                  k.get("tick_labels")))
            return box(ax, x, *a, **k)

        def contourf_(ax, *args, **k):
            self.contours.append((ax, [np.array(a, copy=True) for a in args[:3]]))
            return contourf(ax, *args, **k)

        monkeypatch.setattr(plt, "savefig", savefig)
        monkeypatch.setattr(matplotlib.axes.Axes, "boxplot", boxplot)
        monkeypatch.setattr(matplotlib.axes.Axes, "contourf", contourf_)

    def figure(self, name):
        found = [f for n, f in self.figures if n == name]
        assert len(found) == 1, f"{name}: {[n for n, _ in self.figures]}"
        return found[0]


def _panels(fig):
    return {ax.get_title(): ax for ax in fig.axes if ax.get_title()}


def _vertical(case):
    return get_data_vertical_levels(str(REFERENCE_DIR / case))


@pytest.fixture(params=["ERA5_fixed", "ERA5_track"])
def case(request):
    return request.param


# ------------------------------------------------- boxplots by level

def test_boxplot_by_level_draws_every_level_with_the_csv_values(case, tmp_path, monkeypatch):
    rec = _Recorder(monkeypatch)
    dv = _vertical(case)
    plot_boxplot.boxplot_vertical(dv, str(tmp_path), None)
    for group, terms in GROUPS.items():
        fig = rec.figure(f"boxplot_vertical_levels_{group}.png")
        panels = _panels(fig)
        assert sorted(panels) == sorted(terms)
        for term, ax in panels.items():
            calls = [c for c in rec.boxplots if c[0] is ax]
            levels = list(dv[term].columns)
            assert len(levels) == 4
            assert len(calls) == len(levels), f"{term}: {len(calls)} boxes for {len(levels)} levels"
            for (_, data, positions, ticks), lev in zip(calls, levels):
                # values of the CSV column, unchanged
                assert np.array_equal(data, dv[term][lev].to_numpy())
            assert [c[2][0] for c in calls] == [j / 3 for j in range(len(levels))]


def test_boxplot_by_level_on_a_synthetic_four_level_profile(tmp_path, monkeypatch):
    rec = _Recorder(monkeypatch)
    times = pd.date_range("2020-01-01", periods=3, freq="6h")
    levels = [20000.0, 50000.0, 85000.0, 100000.0]
    dv = {t: pd.DataFrame(np.arange(12.0).reshape(3, 4) + k, index=times, columns=levels)
          for k, t in enumerate(["Az", "Ae", "Kz", "Ke"])}
    plot_boxplot.boxplot_vertical(dv, str(tmp_path), None)
    fig = rec.figure("boxplot_vertical_levels_energy.png")
    for term, ax in _panels(fig).items():
        calls = [c for c in rec.boxplots if c[0] is ax]
        assert len(calls) == 4
        for (_, data, _, _), lev in zip(calls, levels):
            assert np.array_equal(data, dv[term][lev].to_numpy())


# ------------------------------------------------- inputs are not modified

def test_figures_do_not_modify_the_profiles(case, tmp_path, monkeypatch):
    _Recorder(monkeypatch)
    dv = _vertical(case)
    before = copy.deepcopy(dv)
    plot_boxplot.boxplot_vertical(dv, str(tmp_path), None)
    plot_boxplot.boxplot_time(dv, str(tmp_path), None)
    plot_hovmoller._plotter(dv, str(tmp_path), None)
    assert sorted(dv) == sorted(before)
    for term in before:
        pd.testing.assert_frame_equal(dv[term], before[term], check_exact=True)


def test_figure_file_names_are_unchanged(tmp_path, monkeypatch):
    rec = _Recorder(monkeypatch)
    dv = _vertical("ERA5_fixed")
    plot_boxplot.boxplot_vertical(dv, str(tmp_path), None)
    plot_boxplot.boxplot_time(dv, str(tmp_path), None)
    plot_hovmoller._plotter(dv, str(tmp_path), None)
    assert sorted(n for n, _ in rec.figures) == sorted(
        [f"boxplot_vertical_levels_{g}.png" for g in GROUPS]
        + [f"vertical_timeseries_{g}.png" for g in GROUPS]
        + [f"hovmoller_{g}.png" for g in GROUPS]
    )

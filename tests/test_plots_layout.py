#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Layout checks of the vertical boxplots by time and of the Hovmöller diagrams.

The figures are drawn with the Agg backend and inspected after drawing:
no two visible texts of a figure (tick labels, titles, axis and colour-bar
labels, notes) overlap, the boxes of a panel never overlap one another, and
the time axis keeps a small number of ticks. Hourly, 3-hourly, 6-hourly,
irregular and very short series are covered, besides the ERA5 references.
"""

from pathlib import Path

import cmocean  # noqa: F401  (registers the "cmo." colour maps, as a run does)
import matplotlib
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
TERMS = ["Az", "Ae", "Kz", "Ke", "Cz", "Ca", "Ck", "Ce", "Gz", "Ge"]
LEVELS = [70000.0, 85000.0, 92500.0, 100000.0]
SERIES = {
    "hourly": pd.date_range("2020-01-01", periods=5, freq="1h"),
    "3-hourly": pd.date_range("2020-01-01", periods=9, freq="3h"),
    "6-hourly-long": pd.date_range("2020-01-01", periods=40, freq="6h"),
    "irregular": pd.to_datetime(["2020-01-01 00:00", "2020-01-01 01:00", "2020-01-01 04:00",
                                 "2020-01-01 06:00", "2020-01-01 18:00"]),
    "two-times": pd.date_range("2020-01-01", periods=2, freq="6h"),
}
MAX_TIME_TICKS = 4


def _synthetic(times):
    rng = np.random.default_rng(0)
    return {t: pd.DataFrame(rng.normal(size=(len(times), len(LEVELS))), index=pd.DatetimeIndex(times), columns=LEVELS)
            for t in TERMS}


def _cases():
    cases = {name: _synthetic(times) for name, times in SERIES.items()}
    for ref in ("ERA5_fixed", "ERA5_track"):
        cases[ref] = get_data_vertical_levels(str(REFERENCE_DIR / ref))
    return cases


@pytest.fixture(params=sorted(_cases()))
def case(request):
    return request.param, _cases()[request.param]


@pytest.fixture
def saved(monkeypatch):
    """Figures handed to savefig, drawn so that text extents are final."""
    figures = []

    def savefig(fname, *args, **kwargs):
        fig = plt.gcf()
        fig.canvas.draw()
        figures.append((Path(str(fname)).name, fig))

    monkeypatch.setattr(plt, "savefig", savefig)
    return figures


def _drawn_ticklabels(axis):
    """Tick labels that are drawn: matplotlib also returns labels of ticks outside the view limits."""
    low, high = sorted(axis.get_view_interval())
    span = high - low
    return [tick.label1 for tick in axis.get_major_ticks()
            if low - 1e-9 * span <= tick.get_loc() <= high + 1e-9 * span and tick.label1.get_visible()]


def _texts(fig):
    """Drawn, non-empty text artists of a figure with their window extents."""
    renderer = fig.canvas.get_renderer()
    artists = list(fig.texts)
    for ax in fig.axes:
        artists += [ax.title, ax._left_title, ax._right_title, ax.xaxis.label, ax.yaxis.label,
                    ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()]
        artists += _drawn_ticklabels(ax.xaxis) + _drawn_ticklabels(ax.yaxis)
    found = []
    for artist in artists:
        if artist.get_visible() and artist.get_text().strip():
            bbox = artist.get_window_extent(renderer)
            if bbox.width > 0 and bbox.height > 0:
                found.append((artist.get_text(), bbox))
    return found


def _overlapping_texts(fig):
    texts = _texts(fig)
    pairs = []
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            if texts[i][1].overlaps(texts[j][1]):
                pairs.append((texts[i][0], texts[j][0]))
    return pairs


def _time_panels(fig):
    return [ax for ax in fig.axes if ax.get_title() in TERMS]


def test_boxplots_by_time_have_no_overlapping_boxes(case, saved):
    name, dv = case
    plot_boxplot.boxplot_time(dv, "unused", None)
    times = mdates.date2num(dv["Az"].index.to_pydatetime())
    for _, fig in saved:
        for ax in _time_panels(fig):
            spans = sorted((p.get_path().vertices[:, 0].min(), p.get_path().vertices[:, 0].max())
                           for p in ax.patches)
            assert len(spans) == len(times), name
            for (a0, a1), (b0, b1) in zip(spans, spans[1:]):
                assert a1 < b0, f"{name} {ax.get_title()}: boxes overlap ({a0}, {a1}) and ({b0}, {b1})"
            # every box stays centred on its own time
            centres = [(x0 + x1) / 2 for x0, x1 in spans]
            assert np.allclose(centres, times, rtol=0, atol=1e-9)


def test_box_widths_follow_the_spacing_of_the_times():
    hour = 1.0 / 24.0
    x = np.array([0.0, 1.0, 4.0, 6.0, 18.0]) * hour
    widths = plot_boxplot.box_widths(list(x))
    assert np.allclose(widths, np.array([1.0, 1.0, 2.0, 2.0, 12.0]) * 0.6 * hour)
    # half-widths of neighbours never add up to their gap
    for k in range(len(x) - 1):
        assert widths[k] / 2 + widths[k + 1] / 2 < x[k + 1] - x[k]
    assert plot_boxplot.box_widths([5.0]) == [0.6 / 24.0]


@pytest.mark.parametrize("figure", ["boxplot_time", "hovmoller"])
def test_no_overlapping_texts_and_few_time_ticks(case, saved, figure):
    name, dv = case
    if figure == "boxplot_time":
        plot_boxplot.boxplot_time(dv, "unused", None)
    else:
        plot_hovmoller._plotter(dv, "unused", None)
    assert len(saved) == 3
    for fname, fig in saved:
        pairs = _overlapping_texts(fig)
        assert not pairs, f"{name} {fname}: overlapping texts {pairs[:5]}"
        for ax in _time_panels(fig):
            ticks = [t for t in _drawn_ticklabels(ax.xaxis) if t.get_text()]
            assert 2 <= len(ticks) <= MAX_TIME_TICKS, f"{name} {fname} {ax.get_title()}: {len(ticks)} time ticks"
    plt.close("all")


def test_boxplot_by_time_with_a_single_time(saved):
    dv = _synthetic(pd.date_range("2020-01-01", periods=1, freq="6h"))
    plot_boxplot.boxplot_time(dv, "unused", None)
    for _, fig in saved:
        for ax in _time_panels(fig):
            assert len(ax.patches) == 1
        assert not _overlapping_texts(fig)

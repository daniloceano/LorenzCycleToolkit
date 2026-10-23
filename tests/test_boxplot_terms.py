#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Boxplots of the terms (boxplot_terms_*.png) are conventional boxes.

Each box is the rectangle from the first to the third quartile; the median,
the whiskers, the caps and the outliers are those of matplotlib's
boxplot_stats for the same values; colours and the order of the terms are
unchanged. Short, long, skewed and constant series are covered, besides the
ERA5 references. With notches, the notch interval of a short series reaches
beyond the quartiles and the box outline folds outward.
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib import cbook

import src.plots.plot_boxplot as plot_boxplot
from src.plots import utils

matplotlib.use("Agg")
REFERENCE_DIR = Path(__file__).resolve().parent / "data" / "reference"
TERMS = ["Az", "Ae", "Kz", "Ke", "Cz", "Ca", "Ck", "Ce"]


def _synthetic(values_for):
    return pd.DataFrame({term: values_for(i) for i, term in enumerate(TERMS)})


rng = np.random.default_rng(0)
CASES = {
    "five-times": _synthetic(lambda i: rng.normal(10.0 * (i + 1), 3.0, 5)),
    "two-times": _synthetic(lambda i: rng.normal(1.0, 1.0, 2)),
    "forty-times": _synthetic(lambda i: rng.normal(0.0, 1.0, 40)),
    "outliers": _synthetic(lambda i: np.concatenate([rng.normal(0.0, 1.0, 20), [15.0, -12.0]])),
    "constant": _synthetic(lambda i: np.full(5, 2.5)),
}
for ref in ("ERA5_fixed", "ERA5_track"):
    CASES[ref] = pd.read_csv(REFERENCE_DIR / ref / f"testdata_{ref}_results.csv", index_col=0)


@pytest.fixture
def drawn(monkeypatch):
    """(term, values, artists, colour index) of every box drawn by plot_boxplot_terms; savefig is not written."""
    boxes = []
    real = plt.boxplot

    def boxplot(values, *args, **kwargs):
        artists = real(values, *args, **kwargs)
        boxes.append((kwargs["tick_labels"][0], np.asarray(values, dtype=float), artists, kwargs["positions"][0]))
        return artists

    monkeypatch.setattr(plot_boxplot.plt, "boxplot", boxplot)
    monkeypatch.setattr(plot_boxplot.plt, "savefig", lambda *a, **k: None)
    return boxes


@pytest.mark.parametrize("case", sorted(CASES))
def test_term_boxes_are_conventional(drawn, case):
    plot_boxplot.plot_boxplot_terms(CASES[case], "unused", None)
    assert drawn
    for term, values, artists, position in drawn:
        stats = cbook.boxplot_stats(values)[0]
        q1, q3 = stats["q1"], stats["q3"]
        box = artists["boxes"][0].get_path().vertices
        scale = max(abs(q1), abs(q3), 1e-300)
        assert all(min(abs(y - q1), abs(y - q3)) <= 1e-9 * scale for y in box[:, 1]), (case, term)
        assert len({round(x, 12) for x in box[:, 0]}) <= 2, (case, term)  # straight sides: no notch
        median = artists["medians"][0].get_ydata()
        assert np.allclose(median, stats["med"], rtol=1e-12, atol=0), (case, term)
        whiskers = sorted(w.get_ydata()[1] for w in artists["whiskers"])
        assert np.allclose(whiskers, [stats["whislo"], stats["whishi"]], rtol=1e-12, atol=0), (case, term)
        caps = sorted(c.get_ydata()[0] for c in artists["caps"])
        assert np.allclose(caps, [stats["whislo"], stats["whishi"]], rtol=1e-12, atol=0), (case, term)
        assert sorted(artists["fliers"][0].get_ydata()) == sorted(stats["fliers"]), (case, term)


def test_term_boxes_keep_colours_and_order(drawn):
    """One figure per group of TERM_DETAILS, its terms in that order, box i in COLORS[i] with alpha 0.85."""
    plot_boxplot.plot_boxplot_terms(CASES["ERA5_track"], "unused", None)
    start = 0
    for group, details in utils.TERM_DETAILS.items():
        expected = [t for t in details["terms"] if t in CASES["ERA5_track"].columns]
        chunk = drawn[start:start + len(expected)]
        assert [term for term, _, _, _ in chunk] == expected, group
        for i, (term, _, artists, position) in enumerate(chunk):
            assert position == i, (group, term)
            colour = matplotlib.colors.to_rgba(utils.COLORS[i % len(utils.COLORS)], 0.85)
            assert artists["boxes"][0].get_facecolor() == pytest.approx(colour), (group, term)
        start += len(expected)
    assert start == len(drawn)

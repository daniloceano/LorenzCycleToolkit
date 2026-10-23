"""
Applicability of the cyclophaser periods to short series (E0-06, P1).

The minimum lengths (8 time steps for the high-resolution and processed
options, 3 for the low-resolution ones) were demonstrated for cyclophaser
1.7.3 with the options of plot_periods. The tests check both sides of each
limit, against the installed cyclophaser as well, the caveat recorded for
fewer than 20 steps, and that an untested version, option set or time step
raises a compatibility error instead of applying the limits silently.
These limits say when cyclophaser can run, not that its phases are valid.
"""

import math

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from src.plots.plot_periods import cyclophaser_inputs
from src.utils import cyclophaser_compat as cc

HIGH = 0.25   # ERA5-like grid spacing (degrees)
LOW = 2.5     # NCEP-R2-like grid spacing


def series(n, step_hours=1, nan_at=()):
    t = np.linspace(0.0, 1.0, n)
    zeta = list(-2e-5 - 8e-5 * np.exp(-((t - 0.5) / 0.2) ** 2))
    for i in nan_at:
        zeta[i] = float("nan")
    times = pd.date_range("2005-08-09", periods=n, freq=f"{step_hours}h")
    track = pd.DataFrame({"min_max_zeta_850": zeta, "central_lat": [-22.5] * n})
    return track, times


def inputs(n, dx, processed=False, step_hours=1, nan_at=()):
    track, times = series(n, step_hours, nan_at)
    lat = xr.DataArray(np.arange(0.0, 10 * dx, dx))
    values, branch, options, _ = cyclophaser_inputs(track, lat, processed)
    return values, times, branch, options


def run_cyclophaser(values, times, options):
    from cyclophaser import determine_periods

    return determine_periods(list(values), x=times, process_vorticity_args=options)


BOUNDARIES = [  # (dx, processed, branch, largest unavailable N, smallest applicable N)
    (HIGH, False, "high_res", 7, 8),
    (HIGH, True, "processed", 7, 8),
    (LOW, False, "low_res", 2, 3),
]


@pytest.mark.parametrize("dx, processed, branch, below, at", BOUNDARIES, ids=[b[2] for b in BOUNDARIES])
def test_minimum_length_both_sides(dx, processed, branch, below, at):
    values, times, got_branch, options = inputs(below, dx, processed)
    assert got_branch == branch
    low = cc.assess(values, times, branch, options)
    assert not low.applicable
    assert (low.n_usable, low.minimum) == (below, cc.MINIMUM_STEPS[branch])
    assert f"{below} usable time steps < minimum {at}" in low.reason
    assert "cyclophaser not called" in low.reason

    values, times, _, options = inputs(at, dx, processed)
    high = cc.assess(values, times, branch, options)
    assert high.applicable and high.reason is None
    assert cc.NOT_VALIDATED in high.notes


@pytest.mark.parametrize("dx, processed, branch, below, at", BOUNDARIES, ids=[b[2] for b in BOUNDARIES])
def test_limits_agree_with_the_installed_cyclophaser(dx, processed, branch, below, at):
    """The smallest applicable N runs; for the 8-step limits, N = 7 really fails."""
    values, times, _, options = inputs(at, dx, processed)
    assert run_cyclophaser(values, times, options) is not None
    if branch != "low_res":  # N = 2 for low_res is an approved floor, not a demonstrated failure
        values, times, _, options = inputs(below, dx, processed)
        with pytest.raises(ValueError, match="polyorder must be less than window_length"):
            run_cyclophaser(values, times, options)


@pytest.mark.parametrize("dx, processed", [(HIGH, False), (LOW, False)], ids=["high_res", "low_res"])
def test_lowpass_caveat_below_20_steps(dx, processed):
    values, times, branch, options = inputs(19, dx, processed)
    notes19 = cc.assess(values, times, branch, options).notes
    values, times, branch, options = inputs(20, dx, processed)
    notes20 = cc.assess(values, times, branch, options).notes
    assert any("N = 19 < 20" in n and "low-pass" in n for n in notes19)
    assert not any("low-pass" in n for n in notes20)
    assert cc.NOT_VALIDATED in notes19 and cc.NOT_VALIDATED in notes20


def test_no_lowpass_caveat_without_lanczos_filter():
    values, times, branch, options = inputs(19, HIGH, processed=True)
    notes = cc.assess(values, times, branch, options).notes
    assert notes == [cc.NOT_VALIDATED]


def test_only_usable_time_steps_are_counted():
    values, times, branch, options = inputs(9, HIGH, nan_at=(2, 5))
    result = cc.assess(values, times, branch, options)
    assert (result.applicable, result.n_usable, result.n_total) == (False, 7, 9)


def test_unrecognized_version_is_a_compatibility_error(monkeypatch):
    monkeypatch.setattr(cc, "installed_version", lambda: "1.8.0")
    for n in (5, 30):  # also below the minimum: the old limit is not applied silently
        values, times, branch, options = inputs(n, HIGH)
        with pytest.raises(cc.CyclophaserCompatibilityError, match="cyclophaser 1.8.0 installed"):
            cc.assess(values, times, branch, options)


def test_unrecognized_options_are_a_compatibility_error():
    values, times, branch, options = inputs(30, HIGH)
    changed = dict(options, savgol_polynomial=2)
    with pytest.raises(cc.CyclophaserCompatibilityError, match="differ from the tested ones"):
        cc.assess(values, times, branch, changed)


def test_untested_time_step_is_a_compatibility_error():
    values, times, branch, options = inputs(10, HIGH, step_hours=12)
    with pytest.raises(cc.CyclophaserCompatibilityError, match="time step"):
        cc.assess(values, times, branch, options)
    # below the minimum the failure does not depend on the step
    values, times, branch, options = inputs(5, HIGH, step_hours=24)
    assert not cc.assess(values, times, branch, options).applicable


def test_tested_options_are_those_of_plot_periods():
    """The options are pinned: changing them in plot_periods must be noticed."""
    for n in (3, 5, 8, 33, 40, 400):
        for dx, processed, branch in ((HIGH, False, "high_res"), (LOW, False, "low_res"), (HIGH, True, "processed")):
            values, times, got, options = inputs(n, dx, processed)
            assert got == branch and options == cc.tested_options(branch, n)
    assert cc.tested_options("low_res", 40) == {
        "use_filter": 6, "replace_endpoints_with_lowpass": 24, "use_smoothing": 5,
        "use_smoothing_twice": False, "savgol_polynomial": 3, "cutoff_low": 168, "cutoff_high": 12}
    assert math.isclose(cc.MINIMUM_STEPS["high_res"], 8) and cc.TESTED_VERSION == "1.7.3"

"""
Applicability of the cyclophaser life-cycle periods to a series (E0-06).

cyclophaser 1.7.3, with the options used by LorenzCycleToolkit, rejects
short series: in the high-resolution and processed-vorticity branches the
Savitzky-Golay window derived from the series length is not larger than the
polynomial order (3) below 8 time steps ("polyorder must be less than
window_length"). These minimum lengths were demonstrated for this version and
these options with time steps of 1, 3 and 6 hours (Onda 1B diagnostic):

    high_res   (grid spacing below 2 degrees)   8 time steps
    processed  (vorticity from the track file)  8 time steps
    low_res                                      3 time steps

Below the minimum the periods are not applicable: cyclophaser is not called
and the periods are reported as "unavailable". These limits are technical:
they say when cyclophaser can run, not that the phases it finds are
physically meaningful. For another version, other options or another time
step the limits are not known and a CyclophaserCompatibilityError is raised
instead of applying them silently.
"""

import math
from dataclasses import dataclass, field

import pandas as pd

TESTED_VERSION = "1.7.3"
TESTED_STEPS_HOURS = (1, 3, 6)
MINIMUM_STEPS = {"high_res": 8, "processed": 8, "low_res": 3}
# int(0.05 * N) end points are replaced by the low-pass series; below 20 that
# is 0 and cyclophaser 1.7.3 replaces the whole filtered series ([-0:]).
LOWPASS_REPLACEMENT_BELOW = 20
LANCZOS_BRANCHES = ("high_res", "low_res")

NOT_VALIDATED = (
    "technical acceptance only: cyclophaser ran on this series; the phases are "
    "not physically validated by this check"
)


class CyclophaserCompatibilityError(RuntimeError):
    """The cyclophaser version, options or time step are not the tested ones."""


def tested_options(branch, n):
    """Options of plot_periods tested with cyclophaser 1.7.3, for n time steps."""
    if branch == "high_res":
        return {"use_filter": "auto", "replace_endpoints_with_lowpass": 24, "use_smoothing": "auto",
                "use_smoothing_twice": False, "savgol_polynomial": 3, "cutoff_low": 168, "cutoff_high": 48}
    if branch == "low_res":
        return {"use_filter": n // 6 if n // 6 > 3 else 4, "replace_endpoints_with_lowpass": 24,
                "use_smoothing": n // 8 | 1 if n // 8 | 1 > 3 else 4, "use_smoothing_twice": False,
                "savgol_polynomial": 3, "cutoff_low": 168, "cutoff_high": 12}
    if branch == "processed":
        return {"use_filter": False, "use_smoothing": "auto", "use_smoothing_twice": "auto"}
    raise CyclophaserCompatibilityError(f"unknown cyclophaser options branch {branch!r}")


def installed_version():
    from importlib.metadata import version

    return version("cyclophaser")


@dataclass
class Assessment:
    """Decision for one series: apply cyclophaser, or periods unavailable."""

    applicable: bool
    n_usable: int
    n_total: int
    branch: str
    minimum: int
    reason: str = None
    notes: list = field(default_factory=list)


def assess(vorticity_data, times, branch, options):
    """
    Decide whether cyclophaser can be applied to the series.

    Args:
        vorticity_data (list): Series passed to cyclophaser.
        times (sequence): Times of the series.
        branch (str): Options branch selected by plot_periods.
        options (dict): Options selected by plot_periods.

    Returns:
        Assessment

    Raises:
        CyclophaserCompatibilityError: Version, options or time step not tested.
    """
    version = installed_version()
    if version != TESTED_VERSION:
        raise CyclophaserCompatibilityError(
            f"cyclophaser {version} installed; the applicability limits were demonstrated "
            f"for cyclophaser {TESTED_VERSION} only"
        )
    n_total = len(vorticity_data)
    expected = tested_options(branch, n_total)
    if dict(options) != expected:
        raise CyclophaserCompatibilityError(
            f"cyclophaser options of branch {branch!r} differ from the tested ones: {options} != {expected}"
        )

    times = pd.DatetimeIndex(times)
    if len(times) != n_total:
        raise CyclophaserCompatibilityError(f"{n_total} values but {len(times)} times")
    usable = [i for i, v in enumerate(vorticity_data)
              if v is not None and math.isfinite(float(v)) and not pd.isna(times[i])]
    n_usable = len(usable)
    minimum = MINIMUM_STEPS[branch]
    if n_usable < minimum:
        return Assessment(
            False, n_usable, n_total, branch, minimum,
            reason=(f"{n_usable} usable time steps < minimum {minimum} for the {branch} options "
                    f"of cyclophaser {TESTED_VERSION} (Savitzky-Golay window not larger than the "
                    f"polynomial order); cyclophaser not called"),
        )

    steps = {(times[b] - times[a]) / pd.Timedelta("1h") for a, b in zip(usable, usable[1:])}
    if len(steps) != 1 or next(iter(steps)) not in TESTED_STEPS_HOURS:
        raise CyclophaserCompatibilityError(
            f"time step(s) {sorted(steps)} h not tested; the applicability limits were demonstrated "
            f"for regular steps of {TESTED_STEPS_HOURS} h"
        )

    notes = [NOT_VALIDATED]
    if n_usable < n_total:
        notes.append(f"{n_total - n_usable} non-finite values passed to cyclophaser unchanged")
    if branch in LANCZOS_BRANCHES and n_total < LOWPASS_REPLACEMENT_BELOW:
        notes.append(
            f"N = {n_total} < {LOWPASS_REPLACEMENT_BELOW}: int(0.05 N) = 0, so cyclophaser "
            f"{TESTED_VERSION} replaces the whole band-pass filtered series with the "
            f"{options['replace_endpoints_with_lowpass']}-step low-pass series"
        )
    return Assessment(True, n_usable, n_total, branch, minimum, notes=notes)

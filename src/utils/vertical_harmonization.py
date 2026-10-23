"""
Common vertical control volume of every framework (DEC-10, option B; DEC-10-L, option B).

The fixed framework (-f), the track mode (-t) and the --choose mode decide
the pressure levels once, before any term is computed, with the same rule:

    L_n = levels of the input grid whose core fields are finite everywhere
          in the box of time step n
    L   = intersection of the L_n over the time steps processed (input order kept)

and every term at every time step is integrated over L. A level whose fields
were incomplete at one time step used to be dropped at that time step only in
the moving framework, so consecutive column integrals covered different layers
and the finite-difference tendencies and the budget residuals carried the
change of the integration domain. In the fixed framework the box is the same
at every time step. Nothing is interpolated and the input data are not modified.

Only levels of the input grid can be excluded: a level that is not in the
input is never a gap, so a deliberately sparse grid with complete data is used
as it is. A level that is in the grid but whose data are incomplete (even
entirely missing) is excluded. Excluding levels at the top or at the bottom of
the input layer only makes the layer thinner. An excluded level between the
top and the bottom of L would leave an internal gap: the toolkit does not
integrate across missing data (DEC-10-L, option B), so the run stops. Fewer
than two levels cannot be integrated over a layer and also stop the run.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import pandas as pd
import xarray as xr

# Fields that must be complete inside the box for a level to be used.
CORE_FIELDS = (
    "Air Temperature",
    "Eastward Wind Component",
    "Northward Wind Component",
    "Omega Velocity",
)

# The trapezoidal rule (and the vertical derivatives) need two levels.
MIN_LEVELS = 2

RULE = (
    "common vertical control volume (DEC-10, option B), in every framework: a pressure "
    "level of the input grid is used only if every core field is finite everywhere "
    "inside the box of every time step processed; the same levels are used at every "
    "time step and nothing is interpolated"
)

INTERNAL_GAP_POLICY = (
    "stop (DEC-10-L, option B): an excluded level between the top and the bottom of the "
    "levels used stops the run; levels that are not in the input grid are not gaps"
)

# How the messages and the metadata name each framework.
MODES = {
    "fixed": "fixed framework (-f)",
    "track": "track mode (-t)",
    "choose": "--choose mode",
}


class VerticalHarmonizationError(ValueError):
    """The levels common to all time steps cannot be integrated as a single layer."""


def core_fields(variable_list_df):
    """Namelist rows whose fields must be complete for a level to be used."""
    geopotential = (
        "Geopotential" if "Geopotential" in variable_list_df.index else "Geopotential Height"
    )
    return list(CORE_FIELDS) + [geopotential]


def complete_levels(data, variable_list_df, level_name, box):
    """
    Mask of the levels whose core fields are finite everywhere inside a box.

    Args:
        data (xr.Dataset): Fields on pressure levels (one or more times).
        variable_list_df (pd.DataFrame): Parsed namelist.
        level_name (str): Name of the vertical coordinate.
        box (dict): Latitude/longitude selector of the box (see ``nearest_box``).

    Returns:
        xr.DataArray: Boolean, along ``level_name``: True where every core
        field is finite at every point of the box and at every time of ``data``.
    """
    valid = xr.ones_like(data[level_name], dtype=bool)
    for name in core_fields(variable_list_df):
        if name not in variable_list_df.index:
            continue
        var = data[variable_list_df.loc[name]["Variable"]].sel(**box)
        other = [d for d in var.dims if d != level_name]
        valid = valid & var.notnull().all(dim=other)
    return valid


def nearest_box(data, lon_name, lat_name, west, east, south, north):
    """
    Box selector of BoxData: each requested edge is moved to the nearest grid point.

    Returns:
        dict: ``{lat_name: slice(south, north), lon_name: slice(west, east)}``
        with the edges taken from the coordinates of ``data``.
    """

    def nearest(name, value):
        return data[name].sel({name: value}, method="nearest")

    return {
        lat_name: slice(nearest(lat_name, south), nearest(lat_name, north)),
        lon_name: slice(nearest(lon_name, west), nearest(lon_name, east)),
    }


def _time_label(time):
    return pd.Timestamp(time).strftime("%Y-%m-%d %H:%M:%S")


def _pa(levels):
    return ", ".join(f"{p:.0f} Pa" for p in levels)


class VerticalHarmonization:
    """
    Levels used by every time step of a run.

    Attributes:
        mode (str): Framework of the run (a key of ``MODES``).
        level_name (str): Name of the vertical coordinate.
        available_levels (list): Levels of the input grid (Pa, input order).
        used_levels (list): Levels complete in the box of every time step.
        incomplete_at (dict): For each excluded level, the time steps whose
            box does not have complete data on it.
        complete_by_time (list): (time, levels complete at that time step).
        core_fields (list): Namelist rows checked.
    """

    def __init__(self, mode, level_name, available_levels, complete_by_time, core_fields):
        self.mode = mode
        self.level_name = level_name
        self.available_levels = list(available_levels)
        self.complete_by_time = list(complete_by_time)
        self.core_fields = list(core_fields)
        self.used_levels = [
            p for p in self.available_levels
            if all(p in levels for _, levels in self.complete_by_time)
        ]
        self.incomplete_at = {
            p: [t for t, levels in self.complete_by_time if p not in levels]
            for p in self.available_levels if p not in self.used_levels
        }

    @property
    def excluded_levels(self):
        return list(self.incomplete_at)

    def restrict(self, dataset):
        """
        Select the used levels; the object itself is returned when nothing is excluded.

        The used levels are one contiguous block of the input levels (an
        internal gap stops the run), so they are taken with a positional
        slice: a view with the memory layout of the input. A copy (selection
        by a list of levels) would change the order of the floating-point
        sums of the averages and the results in the last bits.
        """
        if not self.excluded_levels:
            return dataset
        first = self.available_levels.index(self.used_levels[0])
        selected = dataset.isel({self.level_name: slice(first, first + len(self.used_levels))})
        if [float(p) for p in selected[self.level_name].values] != self.used_levels:
            raise RuntimeError(
                f"levels {selected[self.level_name].values} selected instead of {self.used_levels}"
            )
        return selected

    def as_metadata(self):
        """Record of the harmonisation for lec_metadata.json."""
        return {
            "mode": self.mode,
            "rule": RULE,
            "internal_gap_policy": INTERNAL_GAP_POLICY,
            "minimum_levels": MIN_LEVELS,
            "time_steps": len(self.complete_by_time),
            "core_fields": self.core_fields,
            "available_levels_Pa": self.available_levels,
            "used_levels_Pa": self.used_levels,
            "excluded_levels_Pa": self.excluded_levels,
            "excluded_because_incomplete_at": [
                {"level_Pa": p, "times": [_time_label(t) for t in times]}
                for p, times in self.incomplete_at.items()
            ],
            "input_layer_Pa": {
                "top": min(self.available_levels), "bottom": max(self.available_levels),
            },
            "effective_layer_Pa": {
                "top": min(self.used_levels), "bottom": max(self.used_levels),
            },
            "complete_levels_by_time": [
                {"time": _time_label(t), "levels_Pa": levels}
                for t, levels in self.complete_by_time
            ],
        }


def harmonize_levels(steps, variable_list_df, app_logger, mode):
    """
    Decide the levels used by every time step of a run, in any framework.

    Args:
        steps (list): ``(time, data_at_time, (west, east, south, north))`` for
            every time step processed, with the box of that step (in the fixed
            framework, the same box at every step).
        variable_list_df (pd.DataFrame): Parsed namelist.
        app_logger (logging.Logger): Logger of the run.
        mode (str): Framework of the run, a key of ``MODES``.

    Returns:
        VerticalHarmonization: The levels used and why the others are excluded.

    Raises:
        VerticalHarmonizationError: If fewer than ``MIN_LEVELS`` levels are
            complete at every time step, or if an excluded level lies between
            the top and the bottom of the levels used (internal gap).
    """
    label = MODES[mode]
    names = variable_list_df["Variable"]
    level_name = names.loc["Vertical Level"]
    lon_name, lat_name = names.loc["Longitude"], names.loc["Latitude"]

    available = None
    complete_by_time = []
    for time, step, (west, east, south, north) in steps:
        levels = [float(p) for p in step[level_name].values]
        available = levels if available is None else available
        box = nearest_box(step, lon_name, lat_name, west, east, south, north)
        mask = complete_levels(step, variable_list_df, level_name, box).values
        complete_by_time.append((time, [p for p, ok in zip(levels, mask) if ok]))

    result = VerticalHarmonization(
        mode, level_name, available, complete_by_time, core_fields(variable_list_df)
    )
    _report_levels_by_time(result, label, app_logger)
    excluded = "; ".join(
        f"{p:.0f} Pa (incomplete at {', '.join(_time_label(t) for t in times)})"
        for p, times in result.incomplete_at.items()
    )

    if len(result.used_levels) < MIN_LEVELS:
        message = (
            f"Vertical harmonization ({label}): only {len(result.used_levels)} pressure "
            f"level(s) ({_pa(result.used_levels) or 'none'}) have complete data inside the "
            f"box of every time step; at least {MIN_LEVELS} are needed to integrate over a "
            f"layer (a numerical minimum). Excluded: {excluded}. No result is computed."
        )
        app_logger.error(f"❌ {message}")
        raise VerticalHarmonizationError(message)

    top, bottom = min(result.used_levels), max(result.used_levels)
    interior = [p for p in result.excluded_levels if top < p < bottom]
    if interior:
        message = (
            f"Vertical harmonization ({label}): excluding {_pa(interior)} would leave an "
            f"internal gap between {top:.0f} Pa and {bottom:.0f} Pa. The toolkit does not "
            "integrate across missing data (DEC-10-L, option B): the run stops and no "
            f"result is computed. Excluded: {excluded}. Levels that are not in the input "
            "grid are not gaps; these levels are in the grid but their data are incomplete."
        )
        app_logger.error(f"❌ {message}")
        raise VerticalHarmonizationError(message)

    if result.excluded_levels:
        app_logger.warning(
            f"⚠️ Vertical harmonization ({label}, common control volume, DEC-10 B): "
            f"{len(result.excluded_levels)} of {len(available)} pressure levels are "
            f"excluded from ALL {len(complete_by_time)} time steps because their data are "
            f"incomplete inside the box of at least one time step: {excluded}. Every term at "
            f"every time step is integrated over {len(result.used_levels)} levels, "
            f"{top:.0f} Pa (top) to {bottom:.0f} Pa (bottom): a thinner layer than the input "
            f"layer, {min(available):.0f} Pa to {max(available):.0f} Pa. Nothing is interpolated."
        )
    app_logger.info(
        f"📏 Vertical harmonization ({label}): input grid {_pa(available)}; "
        f"excluded {_pa(result.excluded_levels) or 'none'}; used {_pa(result.used_levels)} "
        f"({top:.0f} Pa to {bottom:.0f} Pa) at every time step."
    )
    return result


def _report_levels_by_time(result, label, app_logger):
    """Terminal and log: the input grid and the levels complete at each time step."""
    available, n = result.available_levels, len(result.complete_by_time)
    app_logger.info(f"📏 Vertical levels ({label}): input grid of {len(available)} levels, {_pa(available)}.")
    partial = [(t, levels) for t, levels in result.complete_by_time if levels != available]
    for time, levels in partial:
        missing = [p for p in available if p not in levels]
        app_logger.info(
            f"📏 Vertical levels ({label}) at {_time_label(time)}: complete inside the box "
            f"{_pa(levels) or 'none'}; incomplete {_pa(missing)}."
        )
    others = "the other" if partial else "all"
    if n - len(partial):
        app_logger.info(
            f"📏 Vertical levels ({label}): all {len(available)} levels complete inside the "
            f"box at {others} {n - len(partial)} time step(s)."
        )

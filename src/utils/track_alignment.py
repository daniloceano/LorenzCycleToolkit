"""
Alignment of the track times with the times of the meteorological data (track mode).

A track may be sampled more often than the data, for example a 3-hourly track
with 6-hourly ERA5 data (``--cdsapi``/``--ncar`` with ``--time-resolution 6``,
or a local file). The track is then subsampled: the original positions at the
data times are kept and the others are left out. Nothing is interpolated or
extrapolated, missing data times are never filled and the track file is not
changed.

Rules of :func:`align_track_times`:

- track and data times must be strictly increasing (no duplicates);
- when every track time is a data time nothing changes (the behaviour before
  this module, also when the data are sampled more often than the track);
- otherwise the track is subsampled only when the data are sampled less often
  than the track, both are regular inside the track period, the data step is
  a whole multiple of the track step and the data times fall on track times.
  The track times on the data grid (first data time of the track period plus
  whole data steps) are kept, and every one of them must be in the data: a
  missing one is a gap or a truncation of the data and is an error, never a
  subsampled time;
- at least :data:`MIN_TIMES` positions must remain.

Any other case is left to the checks of ``preprocessing.process_data``, which
refuse it with an error.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import json
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# The tendencies (finite differences in time) need at least two times.
MIN_TIMES = 2

# Attribute set on the dataset by preprocessing.process_data (JSON text, so that
# it survives the selections of the run) and read back by track_alignment().
ALIGNMENT_ATTR = "lectk_track_time_alignment"

ALIGNMENT_RULE = (
    "the track positions at the times of the meteorological data are kept with their original "
    "coordinates and the other positions are left out; no interpolation or extrapolation is made, "
    "missing data times are not filled and the track file is not changed"
)


class TrackAlignmentError(ValueError):
    """The track times cannot be aligned with the data times without interpolation or filling."""


@dataclass
class TrackAlignment:
    """
    Result of :func:`align_track_times`.

    Attributes:
        adapted (bool): True when the track was subsampled.
        selected (pd.DatetimeIndex): Track times used by the run.
        n_original (int): Positions in the track file.
        track_step_h (float): Track interval in hours (None if irregular or a single time).
        data_step_h (float): Data interval in hours inside the track period (None if unknown).
        excluded (pd.DatetimeIndex): Track times left out (empty when not adapted).
    """

    adapted: bool
    selected: pd.DatetimeIndex
    n_original: int
    track_step_h: float = None
    data_step_h: float = None
    excluded: pd.DatetimeIndex = field(default_factory=lambda: pd.DatetimeIndex([]))

    def record(self, trackfile=None):
        """Plain dictionary for lec_metadata.json."""
        return {
            "adapted": self.adapted,
            "rule": ALIGNMENT_RULE if self.adapted else "track times used as given (all are data times)",
            "trackfile": trackfile,
            "track_interval_h": self.track_step_h,
            "data_interval_h": self.data_step_h,
            "selected_interval_h": _single_step_h(self.selected),
            "positions_original": self.n_original,
            "positions_selected": int(len(self.selected)),
            "positions_excluded": int(len(self.excluded)),
            "selected_times": [_iso(t) for t in self.selected] if self.adapted else None,
            "spatial_interpolation": False,
            "trackfile_modified": False,
        }


def _iso(time):
    return pd.Timestamp(time).strftime("%Y-%m-%dT%H:%M:%S")


def _hours(delta):
    return float(pd.Timedelta(delta) / pd.Timedelta(hours=1))


def _steps(times):
    return sorted({_hours(d) for d in np.diff(times.values)})


def _single_step_h(times):
    steps = _steps(pd.DatetimeIndex(times)) if len(times) > 1 else []
    return steps[0] if len(steps) == 1 else None


def _times(values, what):
    times = pd.DatetimeIndex(pd.to_datetime(np.asarray(values)))
    if times.tz is not None:
        times = times.tz_localize(None)
    if times.has_duplicates:
        duplicated = [_iso(t) for t in times[times.duplicated()].unique()]
        raise TrackAlignmentError(f"{what} has duplicated times: {', '.join(duplicated)}")
    if not times.is_monotonic_increasing:
        position = int(np.argmax(np.diff(times.values) < np.timedelta64(0)))
        raise TrackAlignmentError(
            f"{what} times are not in increasing order: {_iso(times[position + 1])} comes after "
            f"{_iso(times[position])}"
        )
    return times


def _list(times, limit=6):
    shown = ", ".join(_iso(t) for t in times[:limit])
    return shown + (f" ... ({len(times)} in total)" if len(times) > limit else "")


def align_track_times(track_times, data_times, data_label="the data"):
    """
    Align the track times with the data times by exact match (see the module rules).

    Args:
        track_times (sequence): Times of the track file, as read.
        data_times (sequence): Times of the meteorological data (or the times
            a CDS/NCAR request will retrieve).
        data_label (str): How the messages name the data.

    Returns:
        TrackAlignment: ``adapted`` is False and ``selected`` holds all track
        times when no subsampling is needed or possible (the caller's checks
        then apply unchanged).

    Raises:
        TrackAlignmentError: Duplicated or unordered times, or a subsampling
            that would need interpolation, extrapolation or filling, or that
            leaves fewer than :data:`MIN_TIMES` positions.
    """
    track = _times(track_times, "The track")
    data = _times(data_times, f"The time axis of {data_label}")
    unchanged = TrackAlignment(False, track, len(track), _single_step_h(track))
    if track.isin(data).all():
        unchanged.data_step_h = _single_step_h(data[(data >= track[0]) & (data <= track[-1])])
        return unchanged

    window = data[(data >= track[0]) & (data <= track[-1])]
    track_steps = _steps(track) if len(track) > 1 else []
    # inside the track period when possible, else over the whole data axis
    data_steps = _steps(window) if len(window) > 1 else (_steps(data) if len(data) > 1 else [])
    if not track_steps or not data_steps or min(data_steps) <= min(track_steps):
        # the data are not sampled less often than the track: nothing to subsample
        return unchanged

    track_step, data_step = min(track_steps), min(data_steps)
    if len(window) < MIN_TIMES:
        raise TrackAlignmentError(
            f"Only {len(window)} time(s) of {data_label} fall inside the track period "
            f"({_iso(track[0])} to {_iso(track[-1])}); at least {MIN_TIMES} are needed."
        )
    if len(track_steps) != 1:
        raise TrackAlignmentError(
            f"The track is sampled more often than {data_label} ({data_step:g} h), but its interval is not "
            f"regular (intervals of {', '.join(f'{s:g}' for s in track_steps)} h), so it cannot be subsampled "
            f"safely. Provide a regular track or data with the times of the track."
        )
    if len(data_steps) != 1:
        gaps = [window[i + 1] for i, d in enumerate(np.diff(window.values)) if _hours(d) != data_step]
        raise TrackAlignmentError(
            f"The time axis of {data_label} is not regular inside the track period (intervals of "
            f"{', '.join(f'{s:g}' for s in data_steps)} h): times are missing before {_list(gaps)}. "
            f"Missing data times are not filled."
        )
    if data_step % track_step:
        raise TrackAlignmentError(
            f"The interval of {data_label} ({data_step:g} h) is not a whole multiple of the track interval "
            f"({track_step:g} h): the subsampled track would not keep a regular interval."
        )
    anchor = window[0]
    if _hours(anchor - track[0]) % track_step:
        raise TrackAlignmentError(
            f"The times of {data_label} ({_list(window, 3)}) do not fall on track times "
            f"({_list(track, 3)}): aligning them would need an interpolation of the track, which is not made."
        )

    on_grid = track[[(_hours(t - anchor) % data_step) == 0 for t in track]]
    missing = on_grid[~on_grid.isin(data)]
    if len(missing):
        before = missing[missing < data[0]]
        after = missing[missing > data[-1]]
        inside = missing[(missing >= data[0]) & (missing <= data[-1])]
        parts = []
        if len(before):
            parts.append(f"before the first data time {_iso(data[0])}: {_list(before)}")
        if len(after):
            parts.append(f"after the last data time {_iso(data[-1])}: {_list(after)}")
        if len(inside):
            parts.append(f"inside the data period (gap): {_list(inside)}")
        raise TrackAlignmentError(
            f"Track times on the {data_step:g}-hourly grid of {data_label} are missing from the data, "
            f"{'; '.join(parts)}. These are gaps or truncations of the data, not subsampled times; "
            f"they are not filled. Adjust the track or provide data covering these times."
        )
    if len(on_grid) < MIN_TIMES:
        raise TrackAlignmentError(
            f"Only {len(on_grid)} track time(s) coincide with the times of {data_label}; at least {MIN_TIMES} "
            f"are needed."
        )
    return TrackAlignment(True, on_grid, len(track), track_step, data_step, track[~track.isin(on_grid)])


def warn_adapted(alignment, app_logger):
    """Log, once, the WARNING that the track was subsampled (terminal and log file)."""
    if not alignment.adapted or app_logger is None:
        return
    excluded = alignment.excluded
    app_logger.warning(
        "\n" + "=" * 70 + "\n"
        "⚠️  WARNING — Temporal resolution mismatch\n"
        + "=" * 70 + "\n"
        f"Original track interval:        {alignment.track_step_h:g} h\n"
        f"Meteorological data interval:   {alignment.data_step_h:g} h\n"
        "Action: track timestamps aligned with meteorological data.\n"
        "\n"
        f"Original positions: {alignment.n_original}\n"
        f"Selected positions: {len(alignment.selected)} "
        f"({_iso(alignment.selected[0])} to {_iso(alignment.selected[-1])})\n"
        f"Excluded positions: {len(excluded)}"
        + (f" (first {_iso(excluded[0])}, last {_iso(excluded[-1])})" if len(excluded) else "") + "\n"
        "\n"
        "No spatial interpolation was performed.\n"
        "Original track file was not modified.\n"
        + "=" * 70
    )


def record_alignment(dataset, alignment, trackfile=None):
    """Attach the alignment record to the dataset (read back by :func:`track_alignment`)."""
    dataset.attrs[ALIGNMENT_ATTR] = json.dumps(alignment.record(trackfile))


def track_alignment(dataset):
    """Return the record of :func:`record_alignment`, or None if absent (fixed and --choose modes)."""
    text = dataset.attrs.get(ALIGNMENT_ATTR)
    return None if text is None else json.loads(text)


def planned_request_times(track_times, time_step):
    """
    Times that a ``--cdsapi``/``--ncar`` request retrieves for a track (before any download).

    Same hours as the requests (see ``ncar_gdex.daily_hours``, the rule of the CDS API route).
    """
    from .ncar_gdex import daily_hours

    return pd.DatetimeIndex([
        pd.Timestamp(day) + pd.Timedelta(hours=hour)
        for day, hours in daily_hours(track_times, time_step)
        for hour in hours
    ])

"""
Alignment of the track times with the times of the meteorological data (track mode).

A track sampled more often than the data (e.g. 3-hourly track, 6-hourly data)
is subsampled at the data times with its original positions; nothing is
interpolated, filled or written back to the track file. Genuine gaps,
truncations and ambiguous configurations stay errors. The end-to-end tests
run the same 6-hourly data through the three sources (local file, --ncar,
--cdsapi) with a 3-hourly track that crosses midnight.
"""

import hashlib
import json
import logging
import shutil
import sys

import numpy as np
import pandas as pd
import pytest

import lorenzcycletoolkit
from src.utils import ncar_gdex as g
from src.utils import preprocessing
from src.utils.track_alignment import (MIN_TIMES, TrackAlignmentError, align_track_times,
                                       planned_request_times)
from tests.ncss_emulator import FakeGdex

LOGGER = logging.getLogger("test_track_alignment")


def _range(start, end, hours):
    return pd.date_range(start, end, freq=f"{hours}h")


# ------------------------------------------------------------------ rules


def test_3h_track_with_6h_data_keeps_the_positions_at_the_data_times():
    track = _range("2018-08-09 00:00", "2018-08-09 21:00", 3)
    data = _range("2018-08-09 00:00", "2018-08-09 18:00", 6)
    alignment = align_track_times(track, data)
    assert alignment.adapted
    assert list(alignment.selected) == list(data)
    assert list(alignment.excluded) == list(_range("2018-08-09 03:00", "2018-08-09 21:00", 6))
    record = alignment.record("inputs/track")
    assert (record["track_interval_h"], record["data_interval_h"], record["selected_interval_h"]) == (3.0, 6.0, 6.0)
    assert (record["positions_original"], record["positions_selected"], record["positions_excluded"]) == (8, 4, 4)
    assert record["spatial_interpolation"] is False and record["trackfile_modified"] is False


def test_3h_track_with_12h_data():
    track = _range("2018-08-09 00:00", "2018-08-11 00:00", 3)
    data = _range("2018-08-08 00:00", "2018-08-12 00:00", 12)
    alignment = align_track_times(track, data)
    assert alignment.adapted
    assert list(alignment.selected) == list(_range("2018-08-09 00:00", "2018-08-11 00:00", 12))
    assert len(alignment.excluded) == len(track) - 5


def test_6h_track_with_3h_data_keeps_every_track_position():
    track = _range("2018-08-09 00:00", "2018-08-10 00:00", 6)
    data = _range("2018-08-08 21:00", "2018-08-10 03:00", 3)
    alignment = align_track_times(track, data)
    assert not alignment.adapted and list(alignment.selected) == list(track)


def test_identical_times_are_used_as_given():
    times = _range("2018-08-09 00:00", "2018-08-09 18:00", 6)
    alignment = align_track_times(times, times)
    assert not alignment.adapted and list(alignment.selected) == list(times)
    assert alignment.record()["selected_times"] is None


def test_crossing_midnight_and_month_end():
    track = _range("2018-08-31 15:00", "2018-09-01 09:00", 3)
    data = _range("2018-08-31 12:00", "2018-09-01 12:00", 6)
    alignment = align_track_times(track, data)
    assert [t.strftime("%m-%d %H") for t in alignment.selected] == ["08-31 18", "09-01 00", "09-01 06"]


def test_data_times_off_the_track_times_need_interpolation_and_fail():
    # 3-hourly track at 01, 04, 07 ... UTC; 6-hourly data at 00, 06, 12 ... UTC: no common time
    track = _range("2018-08-09 01:00", "2018-08-10 01:00", 3)
    data = _range("2018-08-09 00:00", "2018-08-10 06:00", 6)
    with pytest.raises(TrackAlignmentError, match="do not fall on track times.*interpolation"):
        align_track_times(track, data)


def test_too_few_common_times_fail():
    track = _range("2018-08-09 00:00", "2018-08-09 06:00", 3)  # 00, 03, 06
    data = _range("2018-08-09 00:00", "2018-08-10 00:00", 12)  # only 00 inside the track period
    with pytest.raises(TrackAlignmentError, match=f"at least {MIN_TIMES} are needed"):
        align_track_times(track, data)


def test_an_internal_gap_of_the_data_is_an_error_not_a_subsampled_time():
    track = _range("2018-08-09 00:00", "2018-08-10 00:00", 3)
    data = _range("2018-08-09 00:00", "2018-08-10 00:00", 6).delete(2)  # 12 UTC missing
    with pytest.raises(TrackAlignmentError, match="not regular inside the track period.*not filled"):
        align_track_times(track, data)


def test_truncated_data_is_an_error():
    track = _range("2018-08-09 00:00", "2018-08-10 21:00", 3)
    data = _range("2018-08-09 00:00", "2018-08-10 06:00", 6)  # ends before the track
    with pytest.raises(TrackAlignmentError, match=r"after the last data time 2018-08-10T06:00:00: 2018-08-10T12:00:00, "
                                                  r"2018-08-10T18:00:00\. These are gaps or truncations"):
        align_track_times(track, data)


def test_an_irregular_track_is_not_subsampled():
    track = _range("2018-08-09 00:00", "2018-08-10 00:00", 3).delete(3)
    data = _range("2018-08-09 00:00", "2018-08-10 00:00", 6)
    with pytest.raises(TrackAlignmentError, match="interval is not regular"):
        align_track_times(track, data)


def test_an_irregular_track_with_every_time_in_the_data_fails():
    # 00, 06, 18, 24 UTC: all in the 6-hourly data, but 12 UTC is missing from the track
    track = pd.DatetimeIndex(["2018-08-09 00:00", "2018-08-09 06:00", "2018-08-09 18:00", "2018-08-10 00:00"])
    data = _range("2018-08-09 00:00", "2018-08-10 00:00", 6)
    with pytest.raises(TrackAlignmentError, match=r"interval is not regular \(intervals of 6, 12 h\)"):
        align_track_times(track, data)


def test_an_irregular_track_fails_also_with_data_sampled_more_often():
    # 19 UTC is not a data time, so no subsampling applies; the track gap must still be reported
    track = pd.DatetimeIndex(["2018-08-09 00:00", "2018-08-09 06:00", "2018-08-09 18:00", "2018-08-09 19:00"])
    data = _range("2018-08-08 21:00", "2018-08-09 21:00", 3)
    with pytest.raises(TrackAlignmentError, match="interval is not regular"):
        align_track_times(track, data)


def test_a_regular_track_with_every_time_in_the_data_is_unchanged():
    track = _range("2018-08-09 00:00", "2018-08-10 00:00", 6)
    data = _range("2018-08-08 00:00", "2018-08-11 00:00", 6)
    alignment = align_track_times(track, data)
    assert not alignment.adapted and list(alignment.selected) == list(track)
    assert alignment.track_step_h == 6


def test_a_track_with_fewer_than_two_positions_fails():
    data = _range("2018-08-09 00:00", "2018-08-10 00:00", 6)
    with pytest.raises(TrackAlignmentError, match=f"at least {MIN_TIMES} are needed"):
        align_track_times(data[:1], data)


def test_a_data_step_that_is_not_a_multiple_of_the_track_step_fails():
    track = _range("2018-08-09 00:00", "2018-08-10 00:00", 4)
    data = _range("2018-08-09 00:00", "2018-08-10 00:00", 6)
    with pytest.raises(TrackAlignmentError, match="not a whole multiple"):
        align_track_times(track, data)


@pytest.mark.parametrize("which", ["track", "data"])
def test_duplicated_or_unordered_times_fail(which):
    good = _range("2018-08-09 00:00", "2018-08-09 18:00", 6)
    duplicated = good.append(good[-1:])
    unordered = good[[0, 2, 1, 3]]
    for bad, message in ((duplicated, "duplicated times"), (unordered, "not in increasing order")):
        args = (bad, good) if which == "track" else (good, bad)
        with pytest.raises(TrackAlignmentError, match=message):
            align_track_times(*args)


def test_planned_request_times_follow_the_cds_rule():
    track = _range("2018-08-09 15:00", "2018-08-10 09:00", 3)
    planned = planned_request_times(track, 6)
    assert list(planned) == list(_range("2018-08-09 12:00", "2018-08-10 06:00", 6))


# ------------------------------------------------------------------ the three data sources

TRACK_TIMES = _range("2018-08-09 12:00", "2018-08-10 12:00", 3)  # 9 positions, crosses midnight
TRACK = "time;Lat;Lon\n" + "".join(
    f"{t:%Y-%m-%d-%H%M};{-23.0 + 0.125 * i:.3f};{-45.5 + 0.125 * i:.3f}\n" for i, t in enumerate(TRACK_TIMES)
)
SELECTED = _range("2018-08-09 12:00", "2018-08-10 12:00", 6)


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _workdir(lectk_workdir):
    workdir = lectk_workdir(namelist="namelist_ERA5-cdsapi")
    shutil.copyfile(workdir / "inputs" / "namelist", workdir / "inputs" / "namelist_ERA5-cdsapi")
    (workdir / "inputs" / "track").write_text(TRACK)
    return workdir


def _run(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["lorenzcycletoolkit.py", *argv])
    return lorenzcycletoolkit.main()


def _outputs(workdir, stem):
    results = workdir / "LEC_Results" / f"{stem}_track"
    frame = pd.read_csv(results / f"{stem}_track_results.csv", index_col=0)
    metadata = json.loads((results / "lec_metadata.json").read_text())
    trackfile = pd.read_csv(results / f"{stem}_track_trackfile", sep=";")
    log = next(results.glob("log.*")).read_text()
    return frame, metadata, trackfile, log


def test_the_three_sources_align_the_3h_track_with_6h_data_in_the_same_way(lectk_workdir, monkeypatch):
    workdir = _workdir(lectk_workdir)
    track_sha = _sha(workdir / "inputs" / "track")

    # --ncar: the request is every 6 h; before the correction the run failed after the download
    fake = FakeGdex()
    monkeypatch.setattr(g, "default_opener", lambda: fake.opener)
    assert _run(monkeypatch, ["ncar_case.nc", "-r", "-t", "--ncar", "--time-resolution", "6"]) == 0
    assert len(fake.calls) == 10  # two days, five fields
    data = workdir / "ncar_case.nc"

    # local file: the same 6-hourly data
    shutil.copyfile(data, workdir / "local_case.nc")
    assert _run(monkeypatch, ["local_case.nc", "-r", "-t"]) == 0

    # --cdsapi: the acquisition is replaced by a copy of the same data (offline)
    requested = []

    def fake_cds(args, track, app_logger):
        requested.append((args.time_resolution, len(track)))
        shutil.copyfile(data, args.infile)
        return args.infile

    monkeypatch.setattr(preprocessing, "get_cdsapi_data", fake_cds)
    assert _run(monkeypatch, ["cds_case.nc", "-r", "-t", "--cdsapi", "--time-resolution", "6"]) == 0
    assert requested == [(6, 9)]

    outputs = {stem: _outputs(workdir, stem) for stem in ("ncar_case", "local_case", "cds_case")}
    reference = outputs["ncar_case"]
    original = pd.read_csv(workdir / "inputs" / "track", sep=";", index_col="time", parse_dates=["time"],
                           date_format="%Y-%m-%d-%H%M")
    for stem, (frame, metadata, trackfile, log) in outputs.items():
        assert list(pd.to_datetime(frame.index)) == list(SELECTED), stem
        assert np.isfinite(frame.to_numpy(dtype=float)).all(), stem
        pd.testing.assert_frame_equal(frame, reference[0], check_exact=True)
        alignment = metadata["track_time_alignment"]
        assert alignment == reference[1]["track_time_alignment"], stem
        assert alignment["adapted"] is True
        assert (alignment["track_interval_h"], alignment["data_interval_h"]) == (3.0, 6.0)
        assert (alignment["positions_original"], alignment["positions_selected"],
                alignment["positions_excluded"]) == (9, 5, 4)
        assert alignment["selected_times"] == [f"{t:%Y-%m-%dT%H:%M:%S}" for t in SELECTED]
        assert metadata["time"]["n"] == 5 and metadata["time"]["distinct_steps_s"] == [21600]
        # one warning per run, in the log file, with the counts of this run
        assert log.count("WARNING — Temporal resolution mismatch") == 1, stem
        for line in ("Original track interval:        3 h", "Meteorological data interval:   6 h",
                     "Original positions: 9", "Selected positions: 5", "Excluded positions: 4",
                     "No spatial interpolation was performed.", "Original track file was not modified."):
            assert line in log, (stem, line)
        # the positions kept are the original coordinates
        times = pd.to_datetime(trackfile["time"], format="%Y-%m-%d-%H%M")
        assert list(times) == list(SELECTED)
        np.testing.assert_array_equal(trackfile["Lat"].to_numpy(), original.loc[SELECTED, "Lat"].to_numpy())
        np.testing.assert_array_equal(trackfile["Lon"].to_numpy(), original.loc[SELECTED, "Lon"].to_numpy())
    assert "after the download the track will be aligned" in outputs["ncar_case"][3]
    assert _sha(workdir / "inputs" / "track") == track_sha


def test_compatible_times_record_no_adaptation_and_no_warning(lectk_workdir, monkeypatch):
    workdir = _workdir(lectk_workdir)
    (workdir / "inputs" / "track").write_text(
        "time;Lat;Lon\n" + "".join(TRACK.splitlines(keepends=True)[1::2]))  # the 6-hourly positions
    fake = FakeGdex()
    monkeypatch.setattr(g, "default_opener", lambda: fake.opener)
    assert _run(monkeypatch, ["ncar_case.nc", "-r", "-t", "--ncar", "--time-resolution", "6"]) == 0
    frame, metadata, trackfile, log = _outputs(workdir, "ncar_case")
    assert list(pd.to_datetime(frame.index)) == list(SELECTED)
    assert metadata["track_time_alignment"]["adapted"] is False
    assert metadata["track_time_alignment"]["positions_selected"] == 5
    assert "Temporal resolution mismatch" not in log


def test_ncar_refuses_before_any_request_a_track_the_request_cannot_serve(lectk_workdir, monkeypatch):
    workdir = _workdir(lectk_workdir)
    shifted = TRACK.replace("-1200;", "-1300;").replace("-1500;", "-1600;").replace("-1800;", "-1900;") \
        .replace("-2100;", "-2200;").replace("-0000;", "-0100;").replace("-0300;", "-0400;") \
        .replace("-0600;", "-0700;").replace("-0900;", "-1000;")  # 13, 16, 19 ... UTC
    (workdir / "inputs" / "track").write_text(shifted)
    fake = FakeGdex()
    monkeypatch.setattr(g, "default_opener", lambda: fake.opener)
    with pytest.raises(TrackAlignmentError, match="do not fall on track times"):
        _run(monkeypatch, ["ncar_case.nc", "-r", "-t", "--ncar", "--time-resolution", "6"])
    assert fake.calls == [] and not (workdir / "ncar_case.nc").exists()


def test_cdsapi_refuses_before_the_download_a_track_the_request_cannot_serve(lectk_workdir, monkeypatch):
    workdir = _workdir(lectk_workdir)
    (workdir / "inputs" / "track").write_text(TRACK.replace("2018-08-09-2100", "2018-08-09-2000"))  # irregular
    monkeypatch.setattr(preprocessing, "get_cdsapi_data", lambda *a: pytest.fail("no download"))
    with pytest.raises(TrackAlignmentError, match="interval is not regular"):
        _run(monkeypatch, ["cds_case.nc", "-r", "-t", "--cdsapi", "--time-resolution", "6"])


def test_a_local_file_with_a_gap_is_refused(lectk_workdir, monkeypatch):
    workdir = _workdir(lectk_workdir)
    fake = FakeGdex()
    monkeypatch.setattr(g, "default_opener", lambda: fake.opener)
    assert _run(monkeypatch, ["ncar_case.nc", "-r", "-t", "--ncar", "--time-resolution", "6"]) == 0
    import xarray as xr
    with xr.open_dataset(workdir / "ncar_case.nc") as dataset:
        dataset.drop_isel(valid_time=2).to_netcdf(workdir / "gap_case.nc")  # 2018-08-10 00 UTC removed
    with pytest.raises(TrackAlignmentError, match="not regular inside the track period"):
        _run(monkeypatch, ["gap_case.nc", "-r", "-t"])

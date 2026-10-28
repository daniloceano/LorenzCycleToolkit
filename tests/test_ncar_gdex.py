"""
Tests of the optional NCAR GDEX route (src/utils/ncar_gdex.py, --ncar).

The service is emulated in-process (tests/ncss_emulator.py), so no test opens
a connection, except the one that checks that an unreachable service fails
cleanly: the session guard blocks that attempt and the test expects it.
"""

import argparse
import contextlib
import datetime
import email.message
import email.utils
import io
import json
import logging
import math
import os
import random
import shutil
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import lorenzcycletoolkit
from src.utils import ncar_gdex as g
from tests.ncss_emulator import FakeClock, FakeGdex, FakeWire, field_values

try:
    import fcntl
except ImportError:  # Windows
    fcntl = None

REPO_ROOT = Path(__file__).resolve().parents[1]
LOGGER = logging.getLogger("test_ncar_gdex")
SMALL = (-20.0, -50.0, -22.0, -48.0)
TWO_DAYS = [(datetime.date(2018, 8, 9), (0, 12)), (datetime.date(2018, 8, 10), (0, 12))]
E1_AREA = (-17.5, -50, -27.5, -40)
E1_DAYS = [(datetime.date(2018, 8, 9), (0, 6, 12, 18)), (datetime.date(2018, 8, 10), (0, 6, 12, 18))]

# The ten E1 requests of the P2-A pilot (2026-10-26), as sent to the service.
P2A_BASE = "https://tds.gdex.ucar.edu/thredds/ncss/grid/files/d633000/e5.oper.an.pl/201808/"
P2A_E1_QUERY = ("?var={v}&north=-17.5&south=-27.5&west=-50&east=-40&horizStride=1&timeStride=6&accept=netcdf3"
                "&time_start=2018-08-{d}T00:00:00Z&time_end=2018-08-{d}T18:00:00Z")
P2A_E1_FILES = {"T": "128_130_t.ll025sc", "Z": "128_129_z.ll025sc", "U": "128_131_u.ll025uv",
                "V": "128_132_v.ll025uv", "W": "128_135_w.ll025sc"}
P2A_E1_URLS = [
    f"{P2A_BASE}e5.oper.an.pl.{code}.201808{d}00_201808{d}23.nc" + P2A_E1_QUERY.format(v=v, d=d)
    for v, code in P2A_E1_FILES.items() for d in ("09", "10")
]


@pytest.fixture(autouse=True)
def no_waits(monkeypatch):
    """Retries without waiting."""
    monkeypatch.setattr(g, "DEFAULT_BACKOFF", 0.0)


def _track(times, lats, lons):
    return pd.DataFrame({"Lat": lats, "Lon": lons}, index=pd.DatetimeIndex(times, name="time"))


def _names(namelist="namelist_ERA5"):
    return g.output_names(str(REPO_ROOT / "inputs" / namelist))


def _retrieve(tmp_path, fake, area=SMALL, days=TWO_DAYS, name="out.nc", **options):
    plan = g.plan_retrieval(area, days)
    path = tmp_path / name
    record = g.retrieve(plan, str(path), _names(), logger=LOGGER, opener=fake.opener, **options)
    return plan, path, record


def _expected(plan, field):
    return field_values(field, list(plan.times), plan.latitudes, plan.longitudes)


# ------------------------------------------------- same request as --cdsapi
def _cds_requests(track, time_resolution, tmp_path):
    """Requests that get_cdsapi_data sends for a track (client and file handling stubbed)."""
    from src.utils import tools

    captured = []

    class Client:
        def __init__(self, *args, **kwargs):
            pass

        def retrieve(self, name, request, target):
            captured.append(request)
            Path(target).write_bytes(b"")

    tmp_path.mkdir(parents=True, exist_ok=True)
    infile = tmp_path / "cds.nc"
    infile.write_bytes(b"")
    args = argparse.Namespace(infile=str(infile), time_resolution=time_resolution)
    with patch("cdsapi.Client", Client), patch("xarray.open_dataset"), patch("xarray.concat"):
        tools.get_cdsapi_data(args, track, LOGGER)
    return captured


def _assert_same_request(track, step, tmp_path):
    cds = _cds_requests(track, step, tmp_path)
    days = g.daily_hours(track.index, step)
    assert [list(g.area_from_track(track))] * len(cds) == [r["area"] for r in cds]
    cds_days = [(datetime.date(int(r["year"]), int(r["month"]), int(r["day"])),
                 tuple(int(t[:2]) for t in r["time"])) for r in cds]
    assert cds_days == days


@pytest.mark.parametrize("start, end, freq, step", [
    ("1992-06-06 18:00", "1992-06-08 15:00", "3h", 3),   # example of the CDS API documentation
    ("2018-08-09 00:00", "2018-08-10 18:00", "6h", 6),
    ("2018-08-09 00:00", "2018-08-10 18:00", "6h", 3),
    ("2020-02-28 21:00", "2020-03-01 03:00", "1h", 1),   # leap day
    ("2019-12-31 12:00", "2020-01-01 06:00", "6h", 12),  # year change
    ("2005-08-09 00:00", "2005-08-09 04:00", "1h", 1),   # one day
    ("2005-08-09 05:00", "2005-08-09 05:00", "1h", 6),   # one time
    ("2010-03-01 00:00", "2010-03-03 00:00", "6h", 24),  # ends at midnight
])
def test_same_area_and_hours_as_cdsapi(start, end, freq, step, tmp_path):
    times = pd.date_range(start, end, freq=freq)
    track = _track(times, np.linspace(-35.3, -27.1, len(times)), np.linspace(-48.6, -38.2, len(times)))
    _assert_same_request(track, step, tmp_path)


def test_same_area_and_hours_as_cdsapi_random_tracks(tmp_path):
    rng = random.Random(20261026)
    for case in range(60):
        start = datetime.datetime(1979, 1, 1) + datetime.timedelta(hours=rng.randrange(0, 46 * 365 * 24))
        times = [start + datetime.timedelta(hours=h) for h in range(0, rng.randrange(1, 240), rng.choice([1, 3, 6]))]
        lats = [rng.uniform(-60, -15) for _ in times]
        lons = [rng.uniform(-70, -10) for _ in times]
        _assert_same_request(_track(times, lats, lons), rng.choice([1, 2, 3, 4, 6, 8, 12, 24]), tmp_path / str(case))


def test_requests_carry_the_parameters_of_the_p2a_pilot():
    plan = g.plan_retrieval(E1_AREA, E1_DAYS)

    def parsed(url):
        parts = urllib.parse.urlsplit(url)
        return parts.scheme, parts.netloc, parts.path, dict(urllib.parse.parse_qsl(parts.query))

    assert sorted(parsed(r.url) for r in plan.requests) == sorted(parsed(u) for u in P2A_E1_URLS)
    assert all(r.shape == (4, 37, 41, 41) for r in plan.requests)


# ---------------------------------------------------------------- longitudes
@pytest.mark.parametrize("west, east, requests, output", [
    (-50, -40, [(-50, -40)], (310.0, 320.0)),        # west of Greenwich: negative limits, GDEX longitudes
    (20, 35, [(20, 35)], (20.0, 35.0)),              # east of Greenwich
    (310, 320, [(-50, -40)], (310.0, 320.0)),        # 0..360 limits
    (-20, 20, [(-20, -0.25), (0, 20)], (-20.0, 20.0)),  # crosses Greenwich: two requests, -180..180 output
    (350, 370, [(-10, -0.25), (0, 10)], (-10.0, 10.0)),
    (150, 180, [(150, 180)], (150.0, 180.0)),
])
def test_longitude_segments(west, east, requests, output):
    segments = g.longitude_segments(west, east)
    assert [(s.west, s.east) for s in segments] == requests
    lons = [lon for s in segments for lon in s.output]
    assert (lons[0], lons[-1]) == output
    assert np.all(np.diff(lons) == 0.25)
    assert all(0 <= lon < 360 for s in segments for lon in s.response)


@pytest.mark.parametrize("west, east", [(170, 190), (-190, -170), (175, 185), (-40, -50), (10, 10), (0, 360),
                                        (-50.1, -40)])
def test_longitude_intervals_refused(west, east):
    with pytest.raises(ValueError):
        g.longitude_segments(west, east)


def test_latitudes_clipped_at_the_pole():
    plan = g.plan_retrieval((95, -50, 80, -40), TWO_DAYS)
    assert plan.latitudes[0] == 90.0 and plan.latitudes[-1] == 80.0


# --------------------------------------------------------------------- names
def test_output_names_follow_the_namelist():
    assert _names("namelist_ERA5") == {"T": "T", "Z": "Z", "W": "W", "U": "U", "V": "V", "time": "time",
                                       "level": "level", "latitude": "latitude", "longitude": "longitude"}
    default = _names("namelist.default")
    assert default == _names("namelist_ERA5-cdsapi")
    assert (default["T"], default["time"], default["level"]) == ("t", "valid_time", "pressure_level")


@pytest.mark.parametrize("old, new, message", [
    ("Geopotential;geopotential;Z;m**2/s**2", "Geopotential Height;geopotential_height;Z;m", "geopotential height"),
    ("Omega Velocity;omega;W;Pa/s", "Omega Velocity;omega;W;hPa/s", "wrong units"),
    ("Air Temperature;air_temperature;T;K", "Air Temperature;air_temperature;T;degC", "wrong units"),
])
def test_namelists_that_cannot_describe_gdex_fields_are_refused(old, new, message, tmp_path):
    text = (REPO_ROOT / "inputs" / "namelist_ERA5").read_text().replace(old, new)
    namelist = tmp_path / "namelist"
    namelist.write_text(text)
    with pytest.raises(ValueError):
        g.output_names(str(namelist))


# ---------------------------------------------------------------- retrieval
def test_retrieval_copies_the_values_and_records_provenance(tmp_path):
    fake = FakeGdex()
    plan, path, record = _retrieve(tmp_path, fake)
    assert len(fake.calls) == len(plan.requests) == 10
    assert not (tmp_path / "out.nc.parts").exists()
    with xr.open_dataset(path) as ds:
        for field in "TZWUV":
            assert ds[field].dtype == np.float32
            assert np.array_equal(ds[field].transpose("time", "level", "latitude", "longitude").values,
                                  _expected(plan, field))
        assert ds["level"].attrs["units"] == "hPa"
        assert list(ds.time.values) == list(np.array(plan.times, dtype="datetime64[ns]"))
    saved = json.loads((tmp_path / "out.nc.provenance.json").read_text())
    assert saved == json.loads(json.dumps(record))
    assert saved["output"]["sha256"] == g._sha256(str(path))
    assert [r["url"] for r in saved["responses"]] == [r.url for r in plan.requests]
    assert saved["totals"] == {"requests": 10, "downloaded_responses": 10, "reused_responses": 0,
                               "downloaded_bytes": sum(r["bytes"] for r in saved["responses"]), "retries": 0,
                               "received_bytes_counted": sum(r["bytes"] + r["http"]["header_bytes"]
                                                             for r in saved["responses"]),
                               "received_bytes_limit": g.DEFAULT_MAX_TOTAL_BYTES}


def test_retrieval_across_greenwich(tmp_path):
    fake = FakeGdex()
    plan, path, _ = _retrieve(tmp_path, fake, area=(-30.0, -2.0, -32.0, 3.0))
    assert len(fake.calls) == 20
    with xr.open_dataset(path) as ds:
        assert float(ds.longitude[0]) == -2.0 and float(ds.longitude[-1]) == 3.0
        assert np.all(np.diff(ds.longitude.values) == 0.25)
        assert np.array_equal(ds["U"].values, _expected(plan, "U"))


def test_float64_copy_is_exact(tmp_path):
    plan, path, record = _retrieve(tmp_path, FakeGdex(), float64=True)
    with xr.open_dataset(path) as ds:
        assert ds["T"].dtype == np.float64
        assert np.array_equal(ds["T"].values, _expected(plan, "T").astype(np.float64))
    assert "float64" in record["output"]["values"]


def test_uncompressed_output_has_the_same_values(tmp_path):
    plan, path, record = _retrieve(tmp_path, FakeGdex(), compress=False)
    with xr.open_dataset(path) as ds:
        assert np.array_equal(ds["Z"].values, _expected(plan, "Z"))
    assert record["output"]["format"] == "NETCDF3_64BIT_OFFSET"


def test_concurrent_requests_give_the_same_file_values(tmp_path):
    fake = FakeGdex(delay=lambda n_times: 0.05)
    plan, path, _ = _retrieve(tmp_path, fake, workers=3)
    assert len({thread for _, _, thread in fake.calls}) > 1
    with xr.open_dataset(path) as ds:
        assert np.array_equal(ds["V"].values, _expected(plan, "V"))


# -------------------------------------------------------------------- faults
@pytest.mark.parametrize("fault", ["http503", "http500", "http502", "http504", "http429", "timeout", "refused",
                                   "reset", "read-timeout", "truncate", "garbage"])
def test_transient_failures_are_retried(fault, tmp_path):
    fake = FakeGdex(faults={"T_20180809": [fault]})
    plan, path, record = _retrieve(tmp_path, fake)
    assert path.exists()
    assert len(fake.calls) == len(plan.requests) + 1
    response = next(r for r in record["responses"] if r["key"] == "T_20180809_0")
    assert response["attempts"] == 2 and record["totals"]["retries"] == 1


def test_a_too_slow_response_is_retried(tmp_path, monkeypatch):
    # an attempt may last its timeout plus the largest response at MIN_TRANSFER_RATE; this one trickles
    monkeypatch.setattr(g, "MIN_TRANSFER_RATE", 1e15)
    fake = FakeGdex()
    slow = []

    def opener(request, timeout=None):
        response = fake.opener(request, timeout)
        if "var=T" in request.full_url and "20180809" in request.full_url and not slow:
            slow.append(request.full_url)
            original = response.read
            response.read = lambda size=-1: time.sleep(0.02) or original(min(size, 64))
        return response

    plan = g.plan_retrieval(SMALL, TWO_DAYS)
    record = g.retrieve(plan, str(tmp_path / "out.nc"), _names(), logger=LOGGER, opener=opener, workers=1,
                        timeout=0.1)
    response = next(r for r in record["responses"] if r["key"] == "T_20180809_0")
    assert response["attempts"] == 2 and slow
    assert any(a["result"].startswith("response too slow") for a in record["attempts"])


def test_retry_after_longer_than_the_wait_stops_without_retrying(tmp_path):
    fake = FakeGdex()
    calls = []

    def opener(request, timeout=None):
        calls.append(request.full_url)
        headers = email.message.Message()
        headers["Retry-After"] = str(int(g.MAX_WAIT) + 300)
        raise urllib.error.HTTPError(request.full_url, 503, "busy", headers, io.BytesIO(b"busy"))

    with pytest.raises(g.NcarTransferError, match="asked to wait 600 s"):
        g.retrieve(g.plan_retrieval(SMALL, TWO_DAYS), str(tmp_path / "out.nc"), _names(), logger=LOGGER,
                   opener=opener, workers=1)
    assert len(calls) == 1 and fake.calls == []


@pytest.mark.parametrize("value, low, high", [("120", 120, 120), ("0", 0, 0), ("soon", None, None)])
def test_retry_after_in_seconds(value, low, high):
    headers = email.message.Message()
    headers["Retry-After"] = value
    wait = g._retry_after(headers)
    assert (wait is None) if low is None else low <= wait <= high


def test_retry_after_as_an_http_date():
    headers = email.message.Message()
    later = datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(seconds=120)
    headers["Retry-After"] = email.utils.format_datetime(later, usegmt=True)
    assert 110 <= g._retry_after(headers) <= 121
    headers.replace_header("Retry-After", "Wed, 21 Oct 2015 07:28:00 GMT")
    assert g._retry_after(headers) == 0.0
    assert g._retry_after(email.message.Message()) is None and g._retry_after(None) is None


def test_a_response_without_length_is_refused_when_its_bytes_cannot_be_counted(tmp_path):
    # an opener that does not count what its connections receive cannot see the chunk framing
    fake = FakeGdex(faults={"Z_20180810": ["no-length"] * 4})
    with pytest.raises(g.NcarContentError, match="no Content-Length.*cannot be counted"):
        _retrieve(tmp_path, fake, workers=1)
    assert sum(1 for _, fault, _ in fake.calls if fault == "no-length") == 1


def test_persistent_failure_stops_without_a_file_and_keeps_the_verified_responses(tmp_path):
    fake = FakeGdex(faults={"V_20180810": ["http500"] * 4})
    with pytest.raises(g.NcarTransferError, match="failed 4 times"):
        _retrieve(tmp_path, fake)
    assert not (tmp_path / "out.nc").exists()
    assert not (tmp_path / "out.nc.part").exists()
    assert not (tmp_path / "out.nc.provenance.json").exists()
    kept = sorted(p.name for p in (tmp_path / "out.nc.parts").glob("*.nc"))
    assert "V_20180810_0.nc" not in kept and not list((tmp_path / "out.nc.parts").glob("*.part"))


@pytest.mark.parametrize("fault, error, text", [
    ("http404", g.NcarRequestError, "--cdsapi"),
    ("http400", g.NcarRequestError, "refused the request"),
    ("http403", g.NcarRequestError, "refused access"),
    ("http302", g.NcarRequestError, "redirect"),
    ("html", g.NcarContentError, "instead of NetCDF"),
    ("wrong-grid", g.NcarContentError, "latitude values differ"),
    ("missing", g.NcarContentError, "missing values"),
    ("nan", g.NcarContentError, "missing values"),
    ("wrong-var", g.NcarContentError, "variable T missing"),
    ("oversize", g.NcarContentError, "more than"),
])
def test_wrong_answers_stop_at_once(fault, error, text, tmp_path):
    fake = FakeGdex(faults={"T_20180809": [fault] * 4})
    with pytest.raises(error, match=text):
        _retrieve(tmp_path, fake, workers=1)
    assert sum(1 for url, f, _ in fake.calls if f == fault) == 1
    assert not (tmp_path / "out.nc").exists()


def test_dates_not_yet_published_fail_on_the_first_request(tmp_path):
    fake = FakeGdex(last_day=datetime.date(2018, 8, 9))
    with pytest.raises(g.NcarRequestError, match="not be published yet"):
        _retrieve(tmp_path, fake, workers=1)
    assert len(fake.calls) == 1 and "2018081000" in fake.calls[0][0]
    assert not (tmp_path / "out.nc.parts").exists()


def test_limits_are_checked_before_any_request(tmp_path, monkeypatch):
    fake = FakeGdex()
    with pytest.raises(g.NcarRetrievalError, match="more than the limit"):
        _retrieve(tmp_path, fake, max_total_bytes=1000)
    monkeypatch.setattr(shutil, "disk_usage", lambda path: shutil._ntuple_diskusage(10, 9, 1))
    with pytest.raises(g.NcarRetrievalError, match="not enough free disk space"):
        _retrieve(tmp_path, fake)
    assert fake.calls == []


def _reservation(area=SMALL, days=TWO_DAYS):
    sizes = {r.max_bytes for r in g.plan_retrieval(area, days).requests}
    assert len(sizes) == 1
    return sizes.pop()


def test_received_bytes_are_counted_and_recorded(tmp_path):
    fake = FakeGdex()
    budget = g.ByteBudget(10 ** 9)
    plan, path, record = _retrieve(tmp_path, fake, budget=budget)
    responses = {r["key"]: r for r in record["responses"]}
    expected = sum(r["bytes"] + r["http"]["header_bytes"] for r in responses.values())
    assert budget.used == expected == record["totals"]["received_bytes_counted"]
    assert budget.reserved == 0 and len(record["attempts"]) == len(plan.requests)
    for attempt in record["attempts"]:
        assert attempt["result"] == "ok" and attempt["charged_bytes"] < plan.requests[0].max_bytes
        assert attempt["body_bytes"] == responses[attempt["key"]]["bytes"]


def test_failed_attempts_count_against_the_byte_limit(tmp_path):
    reservation = _reservation()
    fake = FakeGdex(faults={"T_20180810": ["http503"] * 4})
    budget = g.ByteBudget(2 * reservation)
    with pytest.raises(g.NcarTransferError, match="not requested"):
        _retrieve(tmp_path, fake, workers=1, budget=budget)
    # two failed attempts are charged their whole reservation; the third is not made
    assert len(fake.calls) == 2 and budget.used == 2 * reservation and budget.reserved == 0
    assert [a["result"].startswith("HTTP 503") for a in budget.attempts] == [True, True]


def test_the_default_limit_counts_retries(tmp_path):
    reservation = _reservation()
    fake = FakeGdex(faults={"T_20180809": ["truncate"] * 4})
    with pytest.raises(g.NcarTransferError, match="not requested"):
        _retrieve(tmp_path, fake, days=TWO_DAYS[:1], workers=1, max_total_bytes=3 * reservation, attempts=4)
    assert len(fake.calls) == 3


@pytest.mark.parametrize("fault", ["timeout", "refused", "http503", "reset", "read-timeout", "truncate", "garbage"])
def test_each_kind_of_attempt_is_charged_as_announced(fault, tmp_path):
    # with an opener that cannot count (see the FakeWire tests for the toolkit's own), a failed attempt is
    # charged its whole reservation and a complete response its status line, headers and body
    reservation = _reservation(days=TWO_DAYS[:1])
    budget = g.ByteBudget(10 ** 9)
    fake = FakeGdex(faults={"T_20180809": [fault]})
    plan, path, record = _retrieve(tmp_path, fake, days=TWO_DAYS[:1], workers=1, budget=budget)
    assert path.exists() and budget.reserved == 0
    assert budget.used == sum(a["charged_bytes"] for a in budget.attempts) == record["totals"]["received_bytes_counted"]
    failed, retried = [a for a in budget.attempts if a["key"] == "T_20180809_0"]
    assert failed["result"] != "ok" and failed["charged_bytes"] == reservation
    assert failed.get("header_bytes", 0) + failed.get("body_bytes", 0) <= reservation
    assert retried["result"] == "ok" and retried["charged_bytes"] == retried["header_bytes"] + retried["body_bytes"]
    assert retried["charged_bytes"] < reservation and len(fake.calls) == len(plan.requests) + 1


def test_a_limit_that_cannot_hold_every_request_is_announced(tmp_path, caplog):
    budget = g.ByteBudget(2 * _reservation())
    with caplog.at_level(logging.WARNING, logger=LOGGER.name), pytest.raises(g.NcarTransferError, match="not requested"):
        _retrieve(tmp_path, FakeGdex(), workers=1, budget=budget)
    assert "leaves room for 2 of the 10 request(s)" in caplog.text


def test_a_retrieval_stops_before_a_request_the_limit_cannot_hold_and_resumes(tmp_path):
    # a limit of one reservation lets exactly one request through
    budget = g.ByteBudget(_reservation())
    with pytest.raises(g.NcarTransferError, match="not requested"):
        _retrieve(tmp_path, FakeGdex(), workers=1, budget=budget)
    assert len(budget.attempts) == 1 and budget.used <= budget.limit
    assert len(list((tmp_path / "out.nc.parts").glob("*.nc"))) == 1
    # a larger limit that keeps the bytes already received resumes from the kept response
    budget = g.ByteBudget(10 ** 9, used=budget.used)
    fake = FakeGdex()
    plan, path, record = _retrieve(tmp_path, fake, workers=1, budget=budget)
    assert len(fake.calls) == len(plan.requests) - 1 and record["totals"]["reused_responses"] == 1
    with xr.open_dataset(path) as ds:
        assert np.array_equal(ds["T"].values, _expected(plan, "T"))


def test_one_limit_serves_several_retrievals(tmp_path):
    budget = g.ByteBudget(10 ** 9)
    _retrieve(tmp_path, FakeGdex(), budget=budget)
    budget.limit = budget.used + _reservation() - 1
    fake = FakeGdex()
    with pytest.raises(g.NcarTransferError, match="not requested"):
        _retrieve(tmp_path, fake, name="again.nc", budget=budget)
    assert fake.calls == [] and not (tmp_path / "again.nc.parts").exists()


def test_requests_in_progress_reserve_their_bytes(tmp_path):
    budget = g.ByteBudget(2 * _reservation())
    fake = FakeGdex(delay=lambda n: 0.3)
    with pytest.raises(g.NcarTransferError, match="not requested"):
        _retrieve(tmp_path, fake, workers=4, budget=budget)
    assert len(fake.calls) <= 2 and budget.used <= budget.limit and budget.reserved == 0


# ------------------------------------------------- bytes counted on the wire
WIRE_BASE = "http://gdex.test/thredds/ncss/grid/files/d633000/e5.oper.an.pl"


def _wire_retrieve(tmp_path, monkeypatch, wire, days=TWO_DAYS[:1], name="out.nc", **options):
    """Retrieval through the toolkit's own opener, whose connections are FakeWire sockets."""
    monkeypatch.setattr(g, "_open_socket", wire)
    plan = g.plan_retrieval(SMALL, days, base_url=WIRE_BASE)
    path = tmp_path / name
    record = g.retrieve(plan, str(path), _names(), logger=LOGGER, workers=1, **options)
    return plan, path, record


def _first_wire_request():
    return g.plan_retrieval(SMALL, TWO_DAYS[:1], base_url=WIRE_BASE).requests[0]


def _limit(request, chunked=False):
    """Most bytes a connection may receive for a request (one more detects a longer answer)."""
    return request.max_bytes + (g._framing_allowance(request) if chunked else 0)


def test_the_toolkit_opener_counts_every_byte_its_connections_receive(tmp_path, monkeypatch):
    wire = FakeWire()
    budget = g.ByteBudget(10 ** 9)
    plan, path, record = _wire_retrieve(tmp_path, monkeypatch, wire, budget=budget)
    assert path.exists() and len(wire.sent) == len(plan.requests)
    assert [a["received_bytes"] for a in budget.attempts] == wire.sent
    assert all(a["charged_bytes"] == a["received_bytes"] == a["header_bytes"] + a["body_bytes"]
               for a in budget.attempts)
    assert budget.used == sum(wire.sent) == record["totals"]["received_bytes_counted"] and budget.reserved == 0
    with xr.open_dataset(path) as ds:
        assert np.array_equal(ds["T"].values, _expected(plan, "T"))


def test_a_chunked_response_is_counted_with_its_framing(tmp_path, monkeypatch):
    wire = FakeWire(mode=lambda key, n: "chunked:1000" if key == "T_20180809" else "length")
    budget = g.ByteBudget(10 ** 9)
    plan, path, record = _wire_retrieve(tmp_path, monkeypatch, wire, budget=budget)
    chunked = next(a for a in budget.attempts if a["key"] == "T_20180809_0")
    framing = chunked["received_bytes"] - chunked["header_bytes"] - chunked["body_bytes"]
    assert chunked["result"] == "ok" and framing > 0 and chunked["charged_bytes"] == chunked["received_bytes"]
    assert chunked["reserved_bytes"] == _limit(plan.requests[0], chunked=True)
    assert budget.used == sum(wire.sent) and budget.reserved == 0
    with xr.open_dataset(path) as ds:
        assert np.array_equal(ds["T"].values, _expected(plan, "T"))


@pytest.mark.parametrize("mode, text", [
    ("extensions", "chunk framing and body"),   # one-byte chunks with 4000-byte extensions: framing explodes
    ("endless", "larger than"),                 # chunked body that never ends
    ("headers", "headers alone"),               # 600 kB of headers
])
def test_a_misbehaving_answer_stops_at_the_reservation_and_is_charged_what_arrived(mode, text, tmp_path,
                                                                                   monkeypatch):
    wire = FakeWire(mode=lambda key, n: mode)
    budget = g.ByteBudget(10 ** 9)
    with pytest.raises(g.NcarContentError, match=text):
        _wire_retrieve(tmp_path, monkeypatch, wire, budget=budget)
    request = _first_wire_request()
    (attempt,) = budget.attempts
    # the connection never read more than the reservation (+1 byte); nothing is retried
    assert attempt["received_bytes"] == wire.sent[0] <= _limit(request, chunked=mode != "headers") + 1
    assert len(wire.sent) == 1
    # the answer was abandoned: what the network may still have delivered is charged too
    assert attempt["receive_window_bytes"] == wire.window
    assert attempt["charged_bytes"] == attempt["received_bytes"] + wire.window == budget.used


def test_an_interrupted_answer_is_charged_its_receive_window(tmp_path, monkeypatch):
    wire = FakeWire(mode=lambda key, n: "stall" if n == 1 else "length", window=12345)
    budget = g.ByteBudget(10 ** 9)
    plan, path, record = _wire_retrieve(tmp_path, monkeypatch, wire, budget=budget, attempts=2)
    failed, retried = [a for a in budget.attempts if a["key"] == "T_20180809_0"]
    assert "timed out" in failed["result"] and failed["received_bytes"] == wire.sent[0]
    assert failed["charged_bytes"] == failed["received_bytes"] + 12345
    assert retried["result"] == "ok" and "receive_window_bytes" not in retried
    assert budget.used == sum(wire.sent) + 12345 * sum(1 for m in wire.modes if m == "stall")


def test_a_counted_error_answer_is_charged_what_arrived(tmp_path, monkeypatch):
    wire = FakeWire(mode=lambda key, n: "http503" if n == 1 else "length")
    budget = g.ByteBudget(10 ** 9)
    _wire_retrieve(tmp_path, monkeypatch, wire, budget=budget, attempts=2)
    failed = budget.attempts[0]
    assert failed["result"].startswith("HTTP 503") and failed["received_bytes"] == wire.sent[0] > 3000
    assert failed["charged_bytes"] == failed["received_bytes"] + wire.window


def test_bytes_past_the_limit_are_announced_and_stop_the_retrieval(tmp_path, monkeypatch, caplog):
    # the receive window of an interrupted answer may take the total past the limit: said, never silent
    request = _first_wire_request()
    wire = FakeWire(mode=lambda key, n: "stall", window=10 * request.max_bytes)
    budget = g.ByteBudget(request.max_bytes + 1000)
    with caplog.at_level(logging.WARNING, logger=LOGGER.name), pytest.raises(g.NcarTransferError, match="not requested"):
        _wire_retrieve(tmp_path, monkeypatch, wire, budget=budget, attempts=3)
    assert len(wire.sent) == 1 and budget.used > budget.limit
    assert f"pass the limit by {budget.used - budget.limit}" in caplog.text


def test_no_room_for_the_chunk_framing_stops_before_the_body(tmp_path, monkeypatch):
    request = _first_wire_request()
    wire = FakeWire(mode=lambda key, n: "chunked:8192")
    budget = g.ByteBudget(request.max_bytes)
    with pytest.raises(g.NcarTransferError, match="response not read, because its chunk framing"):
        _wire_retrieve(tmp_path, monkeypatch, wire, budget=budget)
    (attempt,) = budget.attempts
    assert len(wire.sent) == 1 and attempt["received_bytes"] == wire.sent[0] < request.max_bytes
    assert attempt["reserved_bytes"] == request.max_bytes and budget.reserved == 0


def test_the_counting_socket_counts_and_limits_every_read():
    wire = FakeWire()
    raw = wire(("gdex.test", 80))
    raw.sendall(b"GET /missing HTTP/1.1\r\nHost: gdex.test\r\n\r\n")
    counted = g._CountingSocket(raw, g._Received(50))
    assert len(counted.recv(30)) == 30
    # the next read is cut to what the limit leaves (+1 byte), and that byte is detected
    with pytest.raises(g._ResponseTooLarge):
        counted.recv(1000)
    assert counted.lectk_received.bytes == 51 == wire.sent[0]
    with pytest.raises(g._ResponseTooLarge):
        counted.recv(1000)
    assert wire.sent[0] == 51
    # without a limit, every read is counted
    free = g._CountingSocket(wire(("gdex.test", 80)), g._Received(None))
    free._sock.sendall(b"GET /missing HTTP/1.1\r\nHost: gdex.test\r\n\r\n")
    data = free.recv(10 ** 6)
    assert free.lectk_received.bytes == len(data) == wire.sent[1] > 51
    free.close()
    assert free.lectk_received.window == wire.window and free._sock.closed


class _AdapterTime:
    """The adapter's ``time`` module, with a simulated monotonic clock."""

    def __init__(self, clock):
        self.monotonic = clock.monotonic

    def __getattr__(self, name):
        return getattr(time, name)


@pytest.mark.parametrize("mode", ["trickle", "trickle:chunked:1000"])
def test_the_time_limit_of_an_attempt_holds_inside_one_slow_read(mode, tmp_path, monkeypatch):
    # Each socket read returns 64 bytes and takes 1 s of a simulated clock: no read comes near the timeout
    # without data, but the attempt passes its time limit long before the answer ends. One read of the
    # HTTP library gathers many socket reads, so the limit must hold inside it, not only between reads.
    clock = FakeClock()
    monkeypatch.setattr(g, "time", _AdapterTime(clock))
    wire = FakeWire(mode=lambda key, n: mode if n == 1 else "length", clock=clock, trickle=(64, 1.0))
    request = _first_wire_request()
    timeout = 10.0
    allowed = timeout + request.max_bytes / g.MIN_TRANSFER_RATE
    budget = g.ByteBudget(10 ** 9)
    plan, path, record = _wire_retrieve(tmp_path, monkeypatch, wire, budget=budget, attempts=2, timeout=timeout)
    slow, retried = [a for a in budget.attempts if a["key"] == request.key]
    assert slow["result"].startswith("response too slow") and retried["result"] == "ok" and path.exists()
    # stopped when the time was over (at most one socket read later), far from the end of the answer
    assert slow["seconds"] <= allowed + 1.0
    assert slow["received_bytes"] == wire.sent[0] <= 64 * (math.ceil(allowed) + 1) < wire.sizes[0] // 4
    # each wait of the socket was bounded by what was left of the attempt, so a silent server cannot hold it either
    assert wire.sockets[0].timeouts and max(wire.sockets[0].timeouts) <= allowed
    # the abandoned answer is charged what arrived plus the receive buffer
    assert slow["charged_bytes"] == slow["received_bytes"] + wire.window


# ---------------------------------------------------------- one at a time
needs_flock = pytest.mark.skipif(fcntl is None, reason="flock is not available")


@needs_flock
def test_a_second_retrieval_of_the_same_file_stops_at_once(tmp_path):
    lock = tmp_path / "out.nc.lock"
    handle = os.open(lock, os.O_CREAT | os.O_RDWR)
    fcntl.flock(handle, fcntl.LOCK_EX)
    os.write(handle, b"4242 on elsewhere")
    fake = FakeGdex()
    try:
        with pytest.raises(g.NcarRetrievalError, match="another retrieval is writing .*process 4242 on elsewhere"):
            _retrieve(tmp_path, fake)
    finally:
        os.close(handle)
    assert fake.calls == [] and not (tmp_path / "out.nc.parts").exists() and lock.exists()


@needs_flock
def test_the_lock_of_a_process_that_died_does_not_block(tmp_path):
    lock = tmp_path / "out.nc.lock"
    code = ("import fcntl, os, sys, time; h = os.open(sys.argv[1], os.O_CREAT | os.O_RDWR); "
            "fcntl.flock(h, fcntl.LOCK_EX); print('locked', flush=True); time.sleep(60)")
    child = subprocess.Popen([sys.executable, "-c", code, str(lock)], stdout=subprocess.PIPE, text=True)
    try:
        assert child.stdout.readline().strip() == "locked"
        with pytest.raises(g.NcarRetrievalError, match="another retrieval is writing"):
            _retrieve(tmp_path, FakeGdex())
    finally:
        child.kill()
        child.wait()
        child.stdout.close()
    # the system released the lock of the dead process; its lock file is reused, then removed
    assert lock.exists()
    plan, path, _ = _retrieve(tmp_path, FakeGdex())
    assert path.exists() and not lock.exists()


def test_two_retrievals_started_together_write_one_file(tmp_path):
    fake = FakeGdex(delay=lambda n: 0.2)
    start, outcomes = threading.Barrier(2), []

    def run(index):
        start.wait()
        time.sleep(0.05 * index)
        try:
            _retrieve(tmp_path, fake, days=TWO_DAYS[:1], workers=1)
            outcomes.append("ok")
        except g.NcarRetrievalError as error:
            outcomes.append(str(error))

    threads = [threading.Thread(target=run, args=(i,)) for i in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert sorted(outcomes)[-1] == "ok" and "another retrieval is writing" in sorted(outcomes)[0]
    assert len(fake.calls) == 5 and not (tmp_path / "out.nc.lock").exists()


def test_the_lock_is_released_when_a_retrieval_fails(tmp_path):
    with pytest.raises(g.NcarTransferError):
        _retrieve(tmp_path, FakeGdex(faults={"V_20180810": ["http500"] * 4}), workers=1)
    assert not (tmp_path / "out.nc.lock").exists()
    plan, path, record = _retrieve(tmp_path, FakeGdex(), workers=1)
    assert path.exists() and record["totals"]["reused_responses"] >= 1


@needs_flock
@pytest.mark.parametrize("recreated", [True, False])
def test_a_lock_file_replaced_while_locking_is_locked_again(recreated, tmp_path, monkeypatch):
    # a retrieval that ends removes its lock file; one that locked the old file must lock the one there now
    real, calls = fcntl.flock, []
    lock = tmp_path / "out.nc.lock"

    def flock(handle, operation):
        real(handle, operation)
        calls.append(operation)
        if len(calls) == 1:
            lock.unlink()
            if recreated:
                lock.write_text("")

    monkeypatch.setattr(fcntl, "flock", flock)
    plan, path, _ = _retrieve(tmp_path, FakeGdex())
    assert path.exists() and len(calls) == 2 and not lock.exists()


def test_without_flock_the_lock_is_the_exclusive_creation_of_the_lock_file(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "fcntl", None)  # as on Windows
    lock = tmp_path / "out.nc.lock"
    plan, path, _ = _retrieve(tmp_path, FakeGdex(), days=TWO_DAYS[:1])
    assert path.exists() and not lock.exists()
    # a lock file left behind (by a crash, or a retrieval in progress) stops the next retrieval and says how to go on
    lock.write_text("")
    fake = FakeGdex()
    with pytest.raises(g.NcarRetrievalError, match="remove .*out.nc.lock if no retrieval is running"):
        _retrieve(tmp_path, fake, name="out.nc")
    assert fake.calls == [] and lock.exists()


def _track_and_args(tmp_path, name="case.nc"):
    trackfile = tmp_path / "track"
    trackfile.write_text(TRACK)
    track = g._read_track(str(trackfile), LOGGER)
    args = argparse.Namespace(infile=str(tmp_path / name), trackfile=str(trackfile), time_resolution=6)
    return track, args


@needs_flock
def test_runs_that_read_an_existing_file_do_not_need_the_lock(tmp_path, monkeypatch):
    track, args = _track_and_args(tmp_path)
    namelist = str(REPO_ROOT / "inputs" / "namelist_ERA5")
    fake = FakeGdex()
    assert g.get_ncar_data(args, track, LOGGER, namelist=namelist, opener=fake.opener) == args.infile
    handle = os.open(f"{args.infile}.lock", os.O_CREAT | os.O_RDWR)
    fcntl.flock(handle, fcntl.LOCK_EX)
    monkeypatch.setattr(g, "default_opener", lambda: pytest.fail("no request to reuse a file"))
    try:
        # reading the existing file works while another run holds the lock ...
        assert g.get_ncar_data(args, track, LOGGER, namelist=namelist) == args.infile
        # ... retrieving a missing file does not
        other = argparse.Namespace(**dict(vars(args), infile=str(tmp_path / "other.nc")))
        handle_other = os.open(f"{other.infile}.lock", os.O_CREAT | os.O_RDWR)
        fcntl.flock(handle_other, fcntl.LOCK_EX)
        with pytest.raises(g.NcarRetrievalError, match="another retrieval is writing"):
            g.get_ncar_data(other, track, LOGGER, namelist=namelist)
        os.close(handle_other)
    finally:
        os.close(handle)


def test_a_file_written_by_another_run_while_waiting_for_the_lock_is_reused(tmp_path, monkeypatch):
    track, args = _track_and_args(tmp_path)
    namelist = str(REPO_ROOT / "inputs" / "namelist_ERA5")
    ready = tmp_path / "ready"
    ready.mkdir()
    done, _ = _track_and_args(ready)
    g.get_ncar_data(argparse.Namespace(**dict(vars(args), infile=str(ready / "case.nc"))), done, LOGGER,
                    namelist=namelist, opener=FakeGdex().opener)
    real = g._exclusive

    @contextlib.contextmanager
    def other_run_ends_first(path):
        # the other run moves its verified file in just before this run gets the lock
        for suffix in ("", ".provenance.json"):
            shutil.copy(ready / f"case.nc{suffix}", f"{path}{suffix}")
        with real(path):
            yield

    monkeypatch.setattr(g, "_exclusive", other_run_ends_first)
    monkeypatch.setattr(g, "default_opener", lambda: pytest.fail("no request: the file of the other run is reused"))
    assert g.get_ncar_data(args, track, LOGGER, namelist=namelist) == args.infile


def test_a_missing_output_directory_is_reported(tmp_path):
    fake = FakeGdex()
    with pytest.raises(g.NcarRetrievalError, match="does not exist"):
        _retrieve(tmp_path, fake, name="missing/out.nc")
    assert fake.calls == []


def test_an_existing_file_is_never_overwritten(tmp_path):
    (tmp_path / "out.nc").write_bytes(b"user data")
    with pytest.raises(g.NcarRetrievalError, match="never overwritten"):
        _retrieve(tmp_path, FakeGdex())
    assert (tmp_path / "out.nc").read_bytes() == b"user data"


# ------------------------------------------------------------------- resume
def test_an_interrupted_retrieval_resumes_from_the_verified_responses(tmp_path):
    with pytest.raises(g.NcarTransferError):
        _retrieve(tmp_path, FakeGdex(faults={"V_20180810": ["http503"] * 4}), workers=1)
    parts = tmp_path / "out.nc.parts"
    stored = sorted(p.name for p in parts.glob("*.nc"))
    # a damaged response is not reused
    damaged = parts / stored[0]
    damaged.write_bytes(damaged.read_bytes()[:-10] + b"0123456789")

    assert "V_20180810_0.nc" not in stored and len(stored) >= 1
    fake = FakeGdex()
    plan, path, record = _retrieve(tmp_path, fake, workers=1)
    # only the responses not kept, plus the damaged one, are requested again
    assert len(fake.calls) == len(plan.requests) - len(stored) + 1
    assert record["totals"]["reused_responses"] == len(stored) - 1
    with xr.open_dataset(path) as ds:
        for field in "TZWUV":
            assert np.array_equal(ds[field].values, _expected(plan, field))
    assert not parts.exists()


def test_reused_responses_keep_when_and_by_which_code_they_were_received(tmp_path):
    with pytest.raises(g.NcarTransferError):
        _retrieve(tmp_path, FakeGdex(faults={"V_20180810": ["http503"] * 4}), workers=1)
    kept = {p.name[:-8]: json.loads(p.read_text()) for p in (tmp_path / "out.nc.parts").glob("*.nc.json")}
    assert kept and all(r["received_utc"] and "lectk_commit" in r for r in kept.values())
    plan, path, record = _retrieve(tmp_path, FakeGdex(), workers=1)
    for response in record["responses"]:
        if response["reused"]:
            assert response["received_utc"] == kept[response["key"]]["received_utc"]
        else:
            assert response["received_utc"] >= max(r["received_utc"] for r in kept.values())
        assert response["lectk_commit"] == g._lectk_commit()


# -------------------------------------------------------------------- reuse
def test_existing_file_is_reused_only_when_verified(tmp_path):
    plan, path, _ = _retrieve(tmp_path, FakeGdex())
    names = _names()
    assert g.check_existing(str(path), plan, names)["output"]["sha256"] == g._sha256(str(path))
    # a shorter request is covered
    assert g.check_existing(str(path), g.plan_retrieval(SMALL, TWO_DAYS[:1]), names)
    # a request with other times, area or names is not
    for other, other_names, text in (
        (g.plan_retrieval(SMALL, [(datetime.date(2018, 8, 11), (0,))]), names, "times"),
        (g.plan_retrieval((-19.0, -50.0, -22.0, -48.0), TWO_DAYS), names, "latitudes"),
        (plan, _names("namelist.default"), "names"),
    ):
        with pytest.raises(g.NcarReuseError, match=text):
            g.check_existing(str(path), other, other_names)
    # a changed file is not
    with open(path, "ab") as handle:
        handle.write(b"\0")
    with pytest.raises(g.NcarReuseError, match="changed after it was written"):
        g.check_existing(str(path), plan, names)
    # nor a file without provenance record
    os.remove(f"{path}.provenance.json")
    with pytest.raises(g.NcarReuseError, match="no provenance record"):
        g.check_existing(str(path), plan, names)


@pytest.mark.parametrize("written, asked", [(True, False), (False, True)])
def test_a_file_of_another_numeric_type_is_not_reused(written, asked, tmp_path):
    # the toolkit computes in the type of its input: a float64 file gives slightly different results
    plan, path, _ = _retrieve(tmp_path, FakeGdex(), float64=written)
    names = _names()
    assert g.check_existing(str(path), plan, names, float64=written)
    with pytest.raises(g.NcarReuseError, match="float64" if written else "float32, but this run asks for float64"):
        g.check_existing(str(path), plan, names, float64=asked)
    # records written before the dtype field are read from their description of the values
    record_path = f"{path}.provenance.json"
    record = json.loads(Path(record_path).read_text())
    del record["output"]["dtype"]
    Path(record_path).write_text(json.dumps(record))
    assert g.check_existing(str(path), plan, names, float64=written)
    with pytest.raises(g.NcarReuseError):
        g.check_existing(str(path), plan, names, float64=asked)


# ---------------------------------------------------------- command line
def test_cdsapi_and_ncar_are_exclusive(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["lorenzcycletoolkit.py", "x.nc", "-t", "--cdsapi", "--ncar"])
    with pytest.raises(SystemExit) as exit_info:
        lorenzcycletoolkit.main()
    assert exit_info.value.code == 2
    assert "alternative data sources" in capsys.readouterr().err


def test_ncar_flag_defaults_to_false():
    args = lorenzcycletoolkit.create_arg_parser().parse_args(["x.nc", "-f"])
    assert args.ncar is False and args.cdsapi is False


TRACK = "time;Lat;Lon\n2018-08-09-0000;-23.00;-45.50\n2018-08-09-0600;-22.75;-45.25\n" \
        "2018-08-09-1200;-22.50;-45.00\n2018-08-09-1800;-22.25;-44.75\n"
BOX = "min_lon;-48\nmax_lon;-42\nmin_lat;-25.5\nmax_lat;-19.5\n"


@pytest.mark.parametrize("mode", ["-f", "-t"])
def test_ncar_run_end_to_end(mode, lectk_workdir, monkeypatch):
    workdir = lectk_workdir(namelist="namelist_ERA5")
    (workdir / "inputs" / "track").write_text(TRACK)
    (workdir / "inputs" / "box_limits").write_text(BOX)
    fake = FakeGdex()
    monkeypatch.setattr(g, "default_opener", lambda: fake.opener)
    argv = ["lorenzcycletoolkit.py", "ncar_case.nc", "-r", mode, "--ncar", "--time-resolution", "6"]
    monkeypatch.setattr(sys, "argv", argv)
    assert lorenzcycletoolkit.main() == 0
    assert len(fake.calls) == 5  # one day, five fields
    method = "fixed" if mode == "-f" else "track"
    results = workdir / "LEC_Results" / f"ncar_case_{method}"
    status = json.loads((results / "run_status.json").read_text())
    metadata = json.loads((results / "lec_metadata.json").read_text())
    provenance = json.loads((workdir / "ncar_case.nc.provenance.json").read_text())
    assert status["status"] == "success"
    assert metadata["input"]["sha256"] == provenance["output"]["sha256"]
    assert provenance["request"]["area"] == {"north": -7.0, "west": -61.0, "south": -38.0, "east": -29.0}
    frame = pd.read_csv(results / f"ncar_case_{method}_results.csv", index_col=0)
    assert len(frame) == 4 and np.isfinite(frame.to_numpy(dtype=float)).all()

    # second run: the verified file is reused, nothing is requested
    monkeypatch.setattr(sys, "argv", argv + ["--overwrite"])
    assert lorenzcycletoolkit.main() == 0
    assert len(fake.calls) == 5


def test_ncar_run_refuses_an_unverified_existing_file(lectk_workdir, monkeypatch):
    workdir = lectk_workdir(namelist="namelist_ERA5")
    (workdir / "inputs" / "track").write_text(TRACK)
    (workdir / "ncar_case.nc").write_bytes(b"not provenance-checked")
    fake = FakeGdex()
    monkeypatch.setattr(g, "default_opener", lambda: fake.opener)
    monkeypatch.setattr(sys, "argv", ["lorenzcycletoolkit.py", "ncar_case.nc", "-r", "-t", "--ncar"])
    with pytest.raises(g.NcarReuseError, match="no provenance record"):
        lorenzcycletoolkit.main()
    assert fake.calls == []
    assert (workdir / "ncar_case.nc").read_bytes() == b"not provenance-checked"


@pytest.mark.expect_network_attempts
def test_unreachable_service_fails_cleanly(tmp_path, monkeypatch, network_attempts):
    monkeypatch.setattr(g, "DEFAULT_ATTEMPTS", 2)
    plan = g.plan_retrieval(SMALL, TWO_DAYS[:1])
    with pytest.raises(g.NcarTransferError, match="failed 2 times"):
        g.retrieve(plan, str(tmp_path / "out.nc"), _names(), logger=LOGGER)
    assert len(network_attempts) >= 1
    assert not (tmp_path / "out.nc").exists()


def test_dry_run_lists_the_requests_without_retrieving(tmp_path, monkeypatch, capsys):
    track = tmp_path / "track"
    track.write_text(TRACK)
    monkeypatch.setattr(g, "default_opener", lambda: pytest.fail("no request in a dry run"))
    code = g.main([str(tmp_path / "out.nc"), "--trackfile", str(track), "--time-resolution", "6",
                   "--namelist", str(REPO_ROOT / "inputs" / "namelist_ERA5"), "--dry-run"])
    lines = capsys.readouterr().out.splitlines()
    assert code == 0
    assert len([line for line in lines if line.startswith("https://")]) == 5
    assert not (tmp_path / "out.nc").exists()


def test_standalone_command_retrieves_then_reuses_and_keeps_the_numeric_type(tmp_path, monkeypatch, capsys):
    track = tmp_path / "track"
    track.write_text(TRACK)
    fake = FakeGdex()
    monkeypatch.setattr(g, "default_opener", lambda: fake.opener)
    out = tmp_path / "out.nc"
    command = [str(out), "--trackfile", str(track), "--time-resolution", "6",
               "--namelist", str(REPO_ROOT / "inputs" / "namelist_ERA5"), "--float64"]
    assert g.main(command) == 0 and len(fake.calls) == 5
    record = json.loads(Path(f"{out}.provenance.json").read_text())
    assert record["output"]["dtype"] == "float64" and record["request"]["trackfile_sha256"] == g._sha256(str(track))
    # the same command reuses the file without a request; the other numeric type is refused
    assert g.main(command) == 0 and len(fake.calls) == 5
    assert g.main(command[:-1]) == 1 and len(fake.calls) == 5
    assert "holds the values as float64, but this run asks for float32" in capsys.readouterr().err


def test_standalone_command_reports_errors_without_a_traceback(tmp_path, capsys):
    track = tmp_path / "track"
    track.write_text(TRACK)
    out = tmp_path / "out.nc"
    out.write_bytes(b"not provenance-checked")
    code = g.main([str(out), "--trackfile", str(track), "--namelist", str(REPO_ROOT / "inputs" / "namelist_ERA5")])
    assert code == 1
    assert "no provenance record" in capsys.readouterr().err
    assert out.read_bytes() == b"not provenance-checked"

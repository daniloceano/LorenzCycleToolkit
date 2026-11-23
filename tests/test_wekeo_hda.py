"""
Offline tests of the WEkEO route (src/utils/wekeo_hda.py, --wekeo).

No network and no hda package are needed: the hda client is simulated. The
CDS and WEkEO routes are fed from the same synthetic "ERA5" (values that
depend only on time, level, latitude and longitude), so that the file of the
WEkEO route can be compared with the file the CDS route
(src/utils/tools.py:get_cdsapi_data, with a simulated cdsapi) writes for the
same track.
"""

import argparse
import io
import json
import logging
import os
import subprocess
import sys
import types
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from src.utils import wekeo_hda as w
from src.utils.ncar_gdex import area_from_track, daily_hours

LOGGER = logging.getLogger("test_wekeo")
STEP = 5.0  # coarse synthetic grid (both routes alike); the toolkit grid is 0.25 degree


# ------------------------------------------------------------------ synthetic ERA5
def truth(times, north, west, south, east):
    """Dataset shaped like a CDS/WEkEO NetCDF (names, dims, order, attributes)."""
    times = pd.DatetimeIndex(times)
    levels = np.array(sorted((float(x) for x in w.PRESSURE_LEVELS), reverse=True))
    lat = np.arange(north, south - 1e-9, -STEP)
    lon = np.arange(west, east + 1e-9, STEP)
    hours = ((times - pd.Timestamp("1979-01-01")) / pd.Timedelta(hours=1)).values
    shape = (len(times), len(levels), len(lat), len(lon))
    base = (hours[:, None, None, None] * 1e-3 + levels[None, :, None, None] * 1e-2
            + lat[None, None, :, None] * 1e-1 + lon[None, None, None, :] * 1e-4)
    units = {"u": "m s**-1", "v": "m s**-1", "t": "K", "w": "Pa s**-1", "z": "m**2 s**-2"}
    data = {name: (("valid_time", "pressure_level", "latitude", "longitude"),
                   (base * (k + 1)).astype("float32"), {"units": units[name], "GRIB_paramId": k})
            for k, name in enumerate(w.SHORT_NAMES)}
    assert all(v[1].shape == shape for v in data.values())
    return xr.Dataset(data, coords={
        "number": 0, "valid_time": times, "pressure_level": ("pressure_level", levels, {"units": "hPa"}),
        "latitude": ("latitude", lat, {"units": "degrees_north"}),
        "longitude": ("longitude", lon, {"units": "degrees_east"}),
        "expver": ("valid_time", np.array(["0001"] * len(times)))}, attrs={"Conventions": "CF-1.7"})


def netcdf_bytes(ds, tmp_path, name):
    path = tmp_path / name
    ds.to_netcdf(path)
    return path.read_bytes()


# ------------------------------------------------------------------ simulated hda
class FakeResponse:
    def __init__(self, body, with_length=True, status=200):
        self.body, self.status_code = body, status
        self.headers = {"Content-Length": str(len(body))} if with_length else {}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def raise_for_status(self):
        if self.status_code >= 400:
            import requests
            raise requests.exceptions.HTTPError(str(self.status_code))

    def iter_content(self, size):
        stream = io.BytesIO(self.body)
        while block := stream.read(size):
            yield block


class FakeClient:
    """Answers like hda.Client: search -> one result, order -> id, HEAD 200, GET the NetCDF."""

    timeout = 5

    def __init__(self, tmp_path, with_length=True, pending=0):
        self.tmp_path, self.with_length, self.pending = tmp_path, with_length, pending
        self.calls, self.orders = [], {}
        self.session = types.SimpleNamespace(get=self._get)

    def search(self, query):
        self.calls.append(("search", query))
        return types.SimpleNamespace(results=[{"id": f"item-{query['month'][0]}",
                                               "properties": {"location": "loc"}}], query=query)

    def post(self, order_query, action):
        self.calls.append(("post", order_query))
        order_id = f"order-{len(self.orders) + 1}"
        query = [q for kind, q in self.calls if kind == "search"][-1]
        self.orders[order_id] = query
        return {"download_id": order_id}

    def head(self, action, order_id):
        self.calls.append(("head", order_id))
        if self.pending:
            self.pending -= 1
            return types.SimpleNamespace(status_code=202)
        return types.SimpleNamespace(status_code=200)

    def full_url(self, *parts):
        return "fake://" + "/".join(parts)

    def _get(self, url, stream, timeout):
        order_id = url.rsplit("/", 1)[-1]
        self.calls.append(("get", order_id))
        q = self.orders[order_id]
        y, m = int(q["year"][0]), int(q["month"][0])
        times = [pd.Timestamp(y, m, int(d), int(t[:2])) for d in q["day"] for t in q["time"]]
        west, south, east, north = q["bbox"]
        body = netcdf_bytes(truth(times, north, west, south, east), self.tmp_path, f"{order_id}.nc")
        return FakeResponse(body, with_length=self.with_length)


# ------------------------------------------------------------------ fixtures
def make_track(start, end, lat=(-30, -25), lon=(-50, -45)):
    index = pd.date_range(start, end, freq="6h")
    return pd.DataFrame({"Lat": np.linspace(*lat, len(index)), "Lon": np.linspace(*lon, len(index))}, index=index)


@pytest.fixture
def track():
    # crosses a month (Jan 30 12 UTC to Feb 1 06 UTC): two monthly requests, extra hours on Jan 30
    return make_track("2020-01-30 12:00", "2020-02-01 06:00")


def args_for(path, max_gib=None):
    return argparse.Namespace(infile=str(path), time_resolution=6, max_gib=max_gib, cds_max_seconds=None)


# ------------------------------------------------------------------ request
def test_plan_monthly_with_cds_area_hours_and_levels(track):
    plan = w.plan_requests(track, 6)
    assert plan["area"] == area_from_track(track) == (-10, -65, -45, -30)
    assert [r["month"] for r in plan["requests"]] == ["2020-01", "2020-02"]
    jan, feb = (r["query"] for r in plan["requests"])
    assert jan["bbox"] == [-65, -45, -30, -10]  # [W, S, E, N]
    assert jan["day"] == ["30", "31"] and jan["time"] == ["00:00", "06:00", "12:00", "18:00"]
    assert feb["day"] == ["01"] and feb["time"] == ["00:00", "06:00"]
    assert jan["pressure_level"] == list(w.PRESSURE_LEVELS) and len(w.PRESSURE_LEVELS) == 37
    assert jan["variable"] == list(w.VARIABLES)
    assert (jan["data_format"], jan["download_format"]) == ("netcdf", "unarchived")
    expected = [pd.Timestamp(d) + pd.Timedelta(hours=h) for d, hs in daily_hours(track.index, 6) for h in hs]
    assert plan["times"] == expected and plan["times"][0] == pd.Timestamp("2020-01-30 12:00")


def test_plan_greenwich_as_cds(track):
    """Across Greenwich the longitudes are passed as the CDS route passes them (west < 0 < east)."""
    plan = w.plan_requests(make_track("2020-03-01", "2020-03-02", lon=(-5, 5)), 6)
    north, west, south, east = plan["area"]
    assert (west, east) == (-20, 20)
    assert plan["requests"][0]["query"]["bbox"] == [-20, south, 20, north]


def test_same_fields_and_levels_as_cds_route(track, tmp_path):
    requests = []

    class Cds:
        def __init__(self, **kw):
            pass

        def retrieve(self, dataset, request, target):
            requests.append(request)
            day = pd.Timestamp(f"{request['year']}-{request['month']}-{request['day']}")
            times = [day + pd.Timedelta(hours=int(t[:2])) for t in request["time"]]
            north, west, south, east = request["area"]
            truth(times, north, west, south, east).to_netcdf(target)

    from src.utils.tools import get_cdsapi_data
    with patch("cdsapi.Client", Cds):
        get_cdsapi_data(args_for(tmp_path / "cds.nc"), track, LOGGER)
    plan = w.plan_requests(track, 6)
    assert {tuple(r["area"]) for r in requests} == {plan["area"]}
    assert {tuple(r["pressure_level"]) for r in requests} == {w.PRESSURE_LEVELS}
    assert {tuple(r["variable"]) for r in requests} == {w.VARIABLES}
    cds_times = [pd.Timestamp(f"{r['year']}-{r['month']}-{r['day']}") + pd.Timedelta(hours=int(t[:2]))
                 for r in requests for t in r["time"]]
    assert cds_times == plan["times"]


# ------------------------------------------------------------------ output contract
def test_file_equivalent_to_cds_route(track, tmp_path):
    """Same synthetic ERA5 through both routes: identical datasets (names, units, dims, order, longitudes)."""
    class Cds:
        def __init__(self, **kw):
            pass

        def retrieve(self, dataset, request, target):
            day = pd.Timestamp(f"{request['year']}-{request['month']}-{request['day']}")
            times = [day + pd.Timedelta(hours=int(t[:2])) for t in request["time"]]
            north, west, south, east = request["area"]
            truth(times, north, west, south, east).to_netcdf(target)

    from src.utils.tools import get_cdsapi_data
    with patch("cdsapi.Client", Cds):
        cds = get_cdsapi_data(args_for(tmp_path / "cds.nc"), track, LOGGER)
    client = FakeClient(tmp_path)
    out = w.get_wekeo_data(args_for(tmp_path / "wekeo.nc"), track, LOGGER, client=client, sleep=lambda s: None)
    with xr.open_dataset(cds) as a, xr.open_dataset(out) as b:
        xr.testing.assert_identical(a, b)
        assert b["t"].dims == ("valid_time", "pressure_level", "latitude", "longitude")
        assert b.longitude.values.min() == -65  # same longitude convention as the request
        assert len(b.valid_time) == 8  # the extra hours of the monthly request were dropped
    prov = json.loads((tmp_path / "wekeo.nc.provenance.json").read_text())
    assert prov["route"] == "wekeo" and prov["dataset"] == w.WEKEO_DATASET
    assert [r["order_id"] for r in prov["requests"]] == ["order-1", "order-2"]
    assert all(len(r["sha256"]) == 64 for r in prov["requests"]) and len(prov["output"]["sha256"]) == 64
    assert prov["bytes_received_this_run"] == sum(r["bytes"] for r in prov["requests"])
    assert "hda_version" in prov
    assert not os.path.exists(w.ledger_path(out)) and not os.path.exists(w.parts_dir(out))


def test_pending_order_is_polled(track, tmp_path):
    waits = []
    client = FakeClient(tmp_path, pending=2)
    w.get_wekeo_data(args_for(tmp_path / "x.nc"), track, LOGGER, client=client, sleep=waits.append)
    assert waits == [w.POLL_FIRST, w.POLL_FIRST * 1.5]


# ------------------------------------------------------------------ package and credential
def test_missing_package(monkeypatch):
    real = w.importlib.import_module

    def fake_import(name, *a):
        if name == "hda":
            raise ImportError("No module named 'hda'")
        return real(name, *a)

    monkeypatch.setattr(w.importlib, "import_module", fake_import)
    with pytest.raises(w.WekeoSetupError, match=r"lorenzcycletoolkit\[wekeo\]"):
        w.make_client()


def test_missing_credential():
    class ConfigurationError(Exception):
        pass

    def configuration():
        raise ConfigurationError("Missing or incomplete configuration")

    fake = types.SimpleNamespace(Configuration=configuration, Client=None)
    with pytest.raises(w.WekeoSetupError) as error:
        w.make_client(fake)
    text = str(error.value)
    assert "~/.hdarc" in text and "HDA_USER" in text and w.REGISTER_URL in text


def test_client_has_limited_retries_and_quiet_token_log():
    seen = {}
    fake = types.SimpleNamespace(Configuration=lambda: "cfg", Client=lambda config, **kw: seen.update(kw) or "c")
    assert w.make_client(fake) == "c"
    assert seen["retry_max"] == 3 and seen["progress"] is False
    assert logging.getLogger("hda").level == logging.WARNING  # hda logs its token at DEBUG


# ------------------------------------------------------------------ volume
def test_estimate_above_limit_requests_nothing(track, tmp_path):
    client = FakeClient(tmp_path)
    with pytest.raises(w.WekeoVolumeError, match="--max-gib"):
        w.get_wekeo_data(args_for(tmp_path / "x.nc", max_gib=0.001), track, LOGGER, client=client)
    assert client.calls == [] and not os.path.exists(tmp_path / "x.nc")


def test_response_past_limit_is_not_read(track, tmp_path):
    plan = w.plan_requests(track, 6)
    budget = w._Budget(1000)
    client = FakeClient(tmp_path, with_length=False)
    client.orders["o"] = plan["requests"][0]["query"]
    with pytest.raises(w.WekeoVolumeError, match="passed the byte limit"):
        w.download_order(client, "o", tmp_path / "part", budget)
    assert budget.used <= budget.limit + w.CHUNK
    client.with_length = True
    with pytest.raises(w.WekeoVolumeError, match="not read"):
        w.download_order(client, "o", tmp_path / "part", w._Budget(1000))


def test_disk_space_check(track, tmp_path, monkeypatch):
    monkeypatch.setattr(w.shutil, "disk_usage", lambda d: types.SimpleNamespace(free=10))
    with pytest.raises(w.WekeoVolumeError, match="disk"):
        w.check_volume(w.plan_requests(track, 6), 10 ** 12, str(tmp_path))


# ------------------------------------------------------------------ orders (no orphan)
def test_unknown_submission_blocks(track, tmp_path):
    path = tmp_path / "x.nc"
    plan = w.plan_requests(track, 6)
    w._save_ledger(str(path), {"2020-01": {"query": plan["requests"][0]["query"], "status": "submitting"}})
    client = FakeClient(tmp_path)
    with pytest.raises(w.WekeoOrderError, match="unknown"):
        w.get_wekeo_data(args_for(path), track, LOGGER, client=client)
    assert client.calls == []


def test_recorded_order_is_resumed_not_resubmitted(track, tmp_path):
    path = tmp_path / "x.nc"
    plan = w.plan_requests(track, 6)
    jan = plan["requests"][0]["query"]
    w._save_ledger(str(path), {"2020-01": {"query": jan, "status": "submitted", "order_id": "old-1"}})
    client = FakeClient(tmp_path)
    client.orders["old-1"] = jan
    w.get_wekeo_data(args_for(path), track, LOGGER, client=client, sleep=lambda s: None)
    searched = [q["month"] for kind, q in client.calls if kind == "search"]
    assert searched == [["02"]]  # only February was ordered; January resumed order old-1
    assert ("get", "old-1") in client.calls


def test_ledger_written_before_order(track, tmp_path):
    path = tmp_path / "x.nc"
    seen = []

    class Recording(FakeClient):
        def post(self, order_query, action):
            seen.append(json.loads(open(w.ledger_path(str(path))).read())["orders"])
            return super().post(order_query, action)

    w.get_wekeo_data(args_for(path), track, LOGGER, client=Recording(tmp_path), sleep=lambda s: None)
    assert seen[0]["2020-01"]["status"] == "submitting"


def test_failed_order_status(tmp_path):
    client = types.SimpleNamespace(head=lambda a, i: types.SimpleNamespace(status_code=500))
    with pytest.raises(w.WekeoOrderError, match="failed"):
        w.wait_order(client, "o", sleep=lambda s: None)


def test_order_wait_has_a_deadline():
    clock = iter(range(0, 10 ** 6, 100))
    client = types.SimpleNamespace(head=lambda a, i: types.SimpleNamespace(status_code=202))
    with pytest.raises(w.WekeoOrderError, match="resumes it"):
        w.wait_order(client, "o", max_seconds=300, sleep=lambda s: None, clock=lambda: next(clock))


def test_response_not_netcdf(tmp_path):
    bad = tmp_path / "bad"
    bad.write_bytes(b"<html>error</html>")
    with pytest.raises(w.WekeoContentError):
        w.netcdf_members(str(bad), str(tmp_path))


# ------------------------------------------------------------------ command line
LAUNCHER = os.path.join(os.path.dirname(__file__), "run_cli_offline.py")
SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lorenzcycletoolkit.py")


@pytest.mark.parametrize("extra, message", [
    (["--wekeo", "--cdsapi"], "alternative data sources"),
    (["--wekeo", "--ncar"], "alternative data sources"),
    (["--max-gib", "5"], "--max-gib needs --ncar or --wekeo"),
])
def test_cli_rejects_combinations(tmp_path, extra, message):
    result = subprocess.run([sys.executable, LAUNCHER, SCRIPT, str(tmp_path / "x.nc"), "-t", *extra],
                            capture_output=True, text=True, timeout=300)
    assert result.returncode == 2 and message in result.stderr

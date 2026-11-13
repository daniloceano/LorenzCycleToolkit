"""
Tests of the experimental GDEX access methods (src/utils/gdex_sources.py).

Get a Subset (API and file downloads) and the kerchunk byte ranges are
emulated in-process; no test opens a connection. Each method must give the
file that NCSS gives for the same plan (tests/ncss_emulator.py), value by value
and with the same coordinate metadata.
"""

import datetime
import email.message
import functools
import io
import json
import logging
import os
import urllib.error
import urllib.parse
import urllib.request
import zlib

import netCDF4
import numpy as np
import pytest
import xarray as xr

from src.utils import gdex_sources as gs
from src.utils import ncar_gdex as g
from tests.ncss_emulator import FakeGdex, FakeResponse, field_values

LOGGER = logging.getLogger("test_gdex_sources")
POLAR = (90.0, -50.0, 88.0, -48.0)          # kerchunk test chunks hold only the rows 90N..80N
POLAR_GREENWICH = (90.0, -2.0, 88.0, 3.0)
DAYS = [(datetime.date(2018, 8, 9), (12, 18)), (datetime.date(2018, 8, 10), (0, 6))]
TOKEN = "s3cr3t-token-value"


@pytest.fixture(autouse=True)
def no_waits(monkeypatch):
    monkeypatch.setattr(g, "DEFAULT_BACKOFF", 0.0)
    monkeypatch.setattr(gs, "SUBSET_POLL", (0.0, 0.0, 1.0))


def _names():
    return g.output_names(os.path.join(os.path.dirname(__file__), "..", "inputs", "namelist_ERA5"))


def _ncss_file(tmp_path, plan):
    path = tmp_path / "ncss.nc"
    g.retrieve(plan, str(path), _names(), logger=LOGGER, opener=FakeGdex().opener)
    return path


def _assert_same_as_ncss(tmp_path, plan, path):
    reference = _ncss_file(tmp_path, plan)
    with xr.open_dataset(path) as a, xr.open_dataset(reference) as b:
        for name in b.variables:
            assert np.array_equal(a[name].values, b[name].values), name
            if name in b.coords:
                assert a[name].attrs == b[name].attrs, name
        for name in b.data_vars:
            assert a[name].dtype == np.float32


# ---------------------------------------------------------------- kerchunk
HOURS = {datetime.datetime.combine(day, datetime.time(h)) for day, hours in DAYS for h in hours}


@functools.lru_cache(maxsize=1)
def _kerchunk_store():
    """GDEX files (byte strings) and kerchunk references of the plan's hours, rows 90N..80N only."""
    rows = FakeKerchunk.ROWS
    lats = 90.0 - 0.25 * np.arange(rows)
    lons = 0.25 * np.arange(1440)
    files, refs = {}, {}
    for field in g.FIELDS:
        entries = {f"{field.name}/.zarray": json.dumps({
            "chunks": [1, 37, rows, 1440], "compressor": {"id": "zlib", "level": 1}, "dtype": "<f4",
            "filters": [{"elementsize": 4, "id": "shuffle"}], "order": "C", "zarr_format": 2,
            "shape": [87648, 37, rows, 1440], "fill_value": 9.999e20}),
            f"{field.name}/.zattrs": json.dumps({"units": field.units, "long_name": field.row,
                                                 "_ARRAY_DIMENSIONS": ["time", "level", "latitude", "longitude"],
                                                 "missing_value": 9.999e20})}
        for day in sorted({m.date() for m in HOURS}):
            url = (f"https://data.rda.ucar.edu/d633000/e5.oper.an.pl/{day:%Y%m}/e5.oper.an.pl.{field.code}."
                   f"{field.grid}.{day:%Y%m%d}00_{day:%Y%m%d}23.nc")
            body = bytearray(b"HDF-header" * 10)
            for moment in sorted(m for m in HOURS if m.date() == day):
                values = field_values(field.name, [moment], lats, lons)[0]
                raw = np.ascontiguousarray(values, dtype="<f4").view(np.uint8).reshape(-1, 4).T.tobytes()
                chunk = zlib.compress(raw, 1)
                index = (moment - datetime.datetime(2010, 1, 1)) // datetime.timedelta(hours=1)
                entries[f"{field.name}/{index}.0.0.0"] = [url, len(body), len(chunk)]
                body += chunk
            files[url.replace("data.rda.ucar.edu", "data.gdex.ucar.edu")] = bytes(body)
        refs[f"{gs.KERCHUNK_NAMES[field.name]}.2010-2019-remote.json"] = entries
    return files, refs


class FakeKerchunk:
    """GDEX NetCDF-4 files as byte ranges, and their kerchunk references written to ``refs_dir``."""

    ROWS = 41

    def __init__(self, refs_dir, faults=None):
        self.calls, self.faults = [], dict(faults or {})
        self.files, self.refs = _kerchunk_store()
        for name, entries in self.refs.items():
            (refs_dir / name).write_text(json.dumps({"version": 1, "refs": entries}))

    def opener(self, request, timeout=None):
        url = request.full_url
        self.calls.append(url)
        fault = self.faults.pop(url, None) if url in self.faults else None
        if fault == "timeout":
            raise TimeoutError("timed out (emulated)")
        start, end = (int(x) for x in request.headers["Range"].split("=")[1].split("-"))
        body = self.files[url][start:end + 1]
        return FakeResponse(body, {"Content-Length": len(body), "Content-Type": "application/octet-stream"},
                            status=206)


def _kerchunk(tmp_path, area=POLAR, faults=None, **options):
    refs = tmp_path / "refs"
    refs.mkdir(exist_ok=True)
    fake = FakeKerchunk(refs, faults)
    plan = g.plan_retrieval(area, DAYS)
    source = gs.KerchunkSource(str(refs), opener=fake.opener)
    path = tmp_path / "kerchunk.nc"
    record = g.retrieve(plan, str(path), _names(), logger=LOGGER, source=source, **options)
    return plan, path, record, fake


def test_kerchunk_gives_the_ncss_file(tmp_path):
    plan, path, record, fake = _kerchunk(tmp_path)
    _assert_same_as_ncss(tmp_path, plan, path)
    assert len(fake.calls) == 5 * len(plan.times)          # one range per field and hour
    assert all("data.gdex.ucar.edu" in url for url in fake.calls)
    assert record["method"]["method"] == "kerchunk" and record["method"]["chunks"] == 20
    assert set(record["method"]["reference_files"]) == set(fake.refs)
    assert record["totals"]["received_bytes_counted"] >= record["method"]["chunk_bytes"]
    assert record["source"]["service"] == gs.KerchunkSource.service
    assert record["responses"][0]["http"]["chunks"][0]["status"] == 206


def _assert_contiguous_longitudes(plan, path):
    """Longitudes increasing, without duplicates or gaps, and those of the plan."""
    with xr.open_dataset(path) as data:
        lons = data["longitude"].values
    assert np.array_equal(lons, np.asarray(plan.longitudes))
    assert np.allclose(np.diff(lons), g.GRID_STEP) and len(np.unique(lons)) == len(lons)


def test_kerchunk_across_greenwich(tmp_path):
    plan, path, _, _ = _kerchunk(tmp_path, area=POLAR_GREENWICH)
    _assert_same_as_ncss(tmp_path, plan, path)


def test_kerchunk_across_greenwich_reads_each_chunk_once(tmp_path):
    # D-026 defect 3: with two longitude segments each piece read the same global chunk again
    plan, path, record, fake = _kerchunk(tmp_path, area=POLAR_GREENWICH)
    assert len({r.segment_index for r in plan.requests}) == 2
    assert len(fake.calls) == 5 * len(plan.times) == record["method"]["chunks"]
    assert record["totals"]["received_bytes_counted"] <= gs.KerchunkSource(
        str(tmp_path / "refs"), opener=fake.opener).estimated_bytes(plan)
    # the second segment's pieces say which piece read their chunks
    shared = [r["http"] for r in record["responses"] if r["http"].get("chunks_read_with")]
    assert len(shared) == 5 * len(DAYS)
    _assert_same_as_ncss(tmp_path, plan, path)
    _assert_contiguous_longitudes(plan, path)


def test_kerchunk_across_greenwich_fits_a_limit_counting_each_chunk_once(tmp_path):
    # the Frente 2b harness allowed 1.5x the estimate (each chunk once); the duplicated reads exhausted it
    refs = tmp_path / "refs"
    refs.mkdir()
    fake = FakeKerchunk(refs)
    plan = g.plan_retrieval(POLAR_GREENWICH, DAYS)
    source = gs.KerchunkSource(str(refs), opener=fake.opener)
    estimate = source.estimated_bytes(plan)
    chunk_bytes = sum(size for _, _, size, _ in source._chunks.values())
    limit = estimate - g.HEADER_ALLOWANCE * len(source._chunks) + chunk_bytes // 2   # < 2x the chunks
    g.retrieve(plan, str(tmp_path / "k.nc"), _names(), logger=LOGGER, source=source, max_total_bytes=limit)
    _assert_same_as_ncss(tmp_path, plan, tmp_path / "k.nc")


def test_kerchunk_resumes_one_segment_reading_its_chunks_once(tmp_path):
    plan, path, _, fake = _kerchunk(tmp_path, area=POLAR_GREENWICH, keep_parts=True)
    parts = tmp_path / "kerchunk.nc.parts"
    path.unlink()
    (tmp_path / "kerchunk.nc.provenance.json").unlink()
    for request in plan.requests:
        if request.segment_index == 1:
            (parts / f"{request.key}.nc").unlink()
    fake.calls.clear()
    source = gs.KerchunkSource(str(tmp_path / "refs"), opener=fake.opener)
    g.retrieve(plan, str(path), _names(), logger=LOGGER, source=source)
    assert len(fake.calls) == 5 * len(plan.times)
    _assert_same_as_ncss(tmp_path, plan, path)


def test_kerchunk_retries_a_failed_range(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    url = next(iter(FakeKerchunk(refs).files))
    plan, path, record, fake = _kerchunk(tmp_path, faults={url: "timeout"})
    assert record["totals"]["retries"] == 1
    _assert_same_as_ncss(tmp_path, plan, path)


def test_kerchunk_limit_is_checked_before_any_range(tmp_path):
    with pytest.raises(g.NcarRetrievalError, match="more than the limit"):
        _kerchunk(tmp_path, max_total_bytes=1024)
    assert not (tmp_path / "kerchunk.nc").exists()


def test_kerchunk_refuses_a_reference_to_another_day(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    fake = FakeKerchunk(refs)
    name = "temperature.2010-2019-remote.json"
    content = json.loads((refs / name).read_text())
    index = (datetime.datetime(2018, 8, 9, 12) - datetime.datetime(2010, 1, 1)) // datetime.timedelta(hours=1)
    content["refs"][f"T/{index}.0.0.0"][0] = content["refs"][f"T/{index}.0.0.0"][0].replace("20180809", "20180810")
    (refs / name).write_text(json.dumps(content))
    source = gs.KerchunkSource(str(refs), opener=fake.opener)
    with pytest.raises(g.NcarContentError, match="not to the T file"):
        g.retrieve(g.plan_retrieval(POLAR, DAYS), str(tmp_path / "k.nc"), _names(), logger=LOGGER, source=source)
    assert fake.calls == []


# ------------------------------------------------------------ Get a Subset
class FakeSubset:
    """GDEX Get a Subset API and file server."""

    def __init__(self, listed=None, lose_submit=False, failing_file=None, listing_fails_after_purge=False,
                 redirect_without_slash=False, redirect_to=None, crossing_fault=True):
        self.listed = list(listed or [])
        # GDEX since 2026-11-10: status/<id>, get_req_files/<id> and purge/<id> answer 301 to the same path with
        # a trailing slash (relative Location, query kept, as Django's APPEND_SLASH does)
        self.redirect_without_slash = redirect_without_slash
        self.redirect_to = redirect_to                  # callable(url) -> Location of every API answer (tests)
        # What Frente 2b observed (mechanism inferred, D-025a): for a box crossing Greenwich in the GDEX
        # 0..360 grid, a file larger than the size listed for it arrives. Emulated with a padded body.
        self.crossing_fault = crossing_fault
        self.crossing = set()
        self.controls = []
        self.lose_submit = lose_submit
        self.failing_file = failing_file
        self.listing_fails_after_purge = listing_fails_after_purge
        self.purged = False
        self.calls, self.files, self.next_id = [], {}, 960001
        self.tokens = set()

    def _json(self, content, status=200):
        body = json.dumps(content).encode()
        return FakeResponse(body, {"Content-Length": len(body), "Content-Type": "application/json"}, status=status)

    def _make_files(self, rid, control):
        start = datetime.datetime.strptime(control["date"].split("/to/")[0], "%Y%m%d%H%M")
        end = datetime.datetime.strptime(control["date"].split("/to/")[1], "%Y%m%d%H%M")
        north, south, west, east = (float(control[k]) for k in ("nlat", "slat", "wlon", "elon"))
        lats = north - 0.25 * np.arange(int(round((north - south) / 0.25)) + 1)
        lons = np.array(sorted((west + 0.25 * i) % 360.0 for i in range(int(round((east - west) / 0.25)) + 1)))
        hours = []
        moment = start
        while moment <= end:
            hours.append(moment)
            moment += datetime.timedelta(hours=1)
        for field in g.FIELDS:
            for day in sorted({h.date() for h in hours}):
                moments = [h for h in hours if h.date() == day]
                name = f"{rid}.{field.name}.e5.oper.an.pl.{field.code}.{field.grid}.{day:%Y%m%d}00_{day:%Y%m%d}23.nc"
                buffer = netCDF4.Dataset("inmemory.nc", "w", memory=1024, format="NETCDF4")
                for dim, size in (("time", len(moments)), ("latitude", len(lats)), ("level", 37),
                                  ("longitude", len(lons))):
                    buffer.createDimension(dim, size)
                for dim, dtype, data, units in (
                        ("latitude", "f8", lats, "degrees_north"), ("level", "f8", g.PRESSURE_LEVELS_HPA, "hPa"),
                        ("longitude", "f8", lons, "degrees_east")):
                    v = buffer.createVariable(dim, dtype, (dim,))
                    v.units = units
                    v[:] = np.asarray(data)
                t = buffer.createVariable("time", "i4", ("time",))
                t.units, t.calendar = "hours since 1900-01-01", "gregorian"
                t[:] = g._hours_since_origin(moments)
                v = buffer.createVariable(field.name, "f4", ("time", "level", "latitude", "longitude"),
                                          fill_value=np.float32(9.999e20))
                v.units, v.long_name, v.missing_value = field.units, field.row, np.float32(9.999e20)
                v[:] = field_values(field.name, moments, lats, lons)
                self.files[f"https://data.gdex.ucar.edu/download/{rid}/{name}"] = bytes(buffer.close())

    def opener(self, request, timeout=None):
        url = request.full_url
        parsed = urllib.parse.urlparse(url)
        query = urllib.parse.parse_qs(parsed.query)
        self.calls.append((request.get_method(), parsed.path))
        if url in self.files:
            if url == self.failing_file:
                raise ConnectionResetError(54, "reset (emulated)")
            body = self.files[url]
            if url.split("/download/", 1)[1].split("/", 1)[0] in self.crossing:
                body += b"\0" * (g.HEADER_ALLOWANCE + 1024)
            return FakeResponse(body, {"Content-Length": len(body), "Content-Type": "application/x-netcdf"})
        if parsed.netloc != "gdex.ucar.edu":
            raise AssertionError(f"request sent to another host: {parsed.scheme}://{parsed.netloc}")
        self.tokens.update(query.get("token", []))
        endpoint = parsed.path.split("/api/", 1)[1]
        target = None if self.redirect_to is None else self.redirect_to(url)
        if target is not None:
            raise self._redirect(url, target)
        if (self.redirect_without_slash and not endpoint.endswith("/")
                and endpoint.split("/", 1)[0] in ("status", "get_req_files", "purge")):
            location = parsed.path + "/" + (f"?{parsed.query}" if parsed.query else "")
            raise self._redirect(url, location)
        if endpoint == "status/":
            if self.purged and self.listing_fails_after_purge:
                raise urllib.error.HTTPError(url, 503, "unavailable", {}, io.BytesIO(b"busy " + TOKEN.encode()))
            return self._json({"data": self.listed} if self.listed else
                              {"data": None, "error_messages": ["No requests for given email"]})
        if endpoint == "submit/":
            control = json.loads(request.data)
            self.controls.append(control)
            rid = str(self.next_id)
            self.next_id += 1
            if self.crossing_fault and float(control["wlon"]) % 360.0 > float(control["elon"]) % 360.0:
                self.crossing.add(rid)
            parameters = ",".join(f"8!{control['dataset']}:{p}" for p in control["param"].split("/"))
            rinfo = (f"dsnum={control['dataset']};startdate=" + "{0}-{1}-{2} {3}:{4}".format(
                control["date"][0:4], control["date"][4:6], control["date"][6:8], control["date"][8:10],
                control["date"][10:12]) + ";enddate=" + "{0}-{1}-{2} {3}:{4}".format(
                *(control["date"].split("/to/")[1][i:j] for i, j in ((0, 4), (4, 6), (6, 8), (8, 10), (10, 12))))
                + f";parameters={parameters};product=1;level=;nlat={control['nlat']};slat={control['slat']};"
                  f"wlon={control['wlon']};elon={control['elon']};tindex={control['groupindex']};ofmt=netCDF")
            self.listed.append({"request_index": int(rid), "status": "Completed", "rinfo": rinfo})
            self._make_files(rid, control)
            if self.lose_submit:
                raise TimeoutError("submission answer lost (emulated)")
            return self._json({"data": {"request_id": rid}})
        rid = endpoint.rstrip("/").rsplit("/", 1)[-1]
        if endpoint.startswith("status/"):
            status = next((r["status"] for r in self.listed if str(r["request_index"]) == rid), None)
            return self._json({"data": {"status": status}})
        if endpoint.startswith("get_req_files/"):
            files = [{"web_path": u, "size": len(b)} for u, b in self.files.items() if f"/{rid}/" in u]
            return self._json({"data": {"web_files": files}})
        if endpoint.startswith("purge/"):
            self.purged = True
            for r in self.listed:
                if str(r["request_index"]) == rid:
                    r["status"] = "Set for Purge"
            return self._json({"data": {"purge_successful": "true"}})
        raise urllib.error.HTTPError(url, 404, "not found", {}, io.BytesIO(b""))

    @staticmethod
    def _redirect(url, location):
        headers = email.message.Message()
        headers["Location"] = location
        headers["Content-Length"] = "0"
        return urllib.error.HTTPError(url, 301, "Moved Permanently", headers, io.BytesIO(b""))

    def submits(self):
        return sum(1 for method, path in self.calls if path.endswith("/submit/"))


def _subset(tmp_path, fake, area=POLAR, name="subset.nc", **source_options):
    plan = g.plan_retrieval(area, DAYS)
    source = gs.SubsetSource(token=TOKEN, opener=fake.opener, api="https://gdex.ucar.edu/api/", **source_options)
    path = tmp_path / name
    record = g.retrieve(plan, str(path), _names(), logger=LOGGER, source=source)
    return plan, path, record


def _listed_like(fake_area_plan, status):
    """A listed request identical to the plan's Get a Subset request."""
    fake = FakeSubset()
    control = gs.subset_control(fake_area_plan)
    request = urllib.request.Request(f"https://gdex.ucar.edu/api/submit/?token={TOKEN}",
                                     data=json.dumps(control).encode(), method="POST")
    fake.opener(request)
    entry = dict(fake.listed[0], status=status)
    return fake, entry


def test_subset_gives_the_ncss_file_without_leaking_the_token(tmp_path, caplog):
    fake = FakeSubset()
    with caplog.at_level(logging.DEBUG):
        plan, path, record = _subset(tmp_path, fake)
    _assert_same_as_ncss(tmp_path, plan, path)
    assert fake.submits() == 1 and fake.tokens == {TOKEN}
    assert record["method"]["created_here"] is True and record["method"]["purge"] == {"requested": False}
    assert TOKEN not in caplog.text
    assert TOKEN not in (tmp_path / "subset.nc.provenance.json").read_text()
    assert not (tmp_path / "subset.nc.parts").exists()


def test_subset_across_greenwich(tmp_path):
    plan, path, _ = _subset(tmp_path, FakeSubset(), area=POLAR_GREENWICH)
    _assert_same_as_ncss(tmp_path, plan, path)


# --- D-026 defect 1: HTTP 301 of the GDEX API
def test_subset_follows_the_trailing_slash_redirect_of_the_api(tmp_path, caplog):
    fake = FakeSubset(redirect_without_slash=True)
    with caplog.at_level(logging.DEBUG):
        plan, path, record = _subset(tmp_path, fake, purge=True)
    paths = [p for _, p in fake.calls]
    assert "/api/status/960001" in paths and "/api/status/960001/" in paths
    assert "/api/get_req_files/960001" in paths and "/api/get_req_files/960001/" in paths
    assert [p for m, p in fake.calls if m == "DELETE"] == ["/api/purge/960001", "/api/purge/960001/"]
    assert record["method"]["purge"]["answer"] == {"purge_successful": "true"}
    assert fake.tokens == {TOKEN} and TOKEN not in caplog.text
    assert TOKEN not in (tmp_path / "subset.nc.provenance.json").read_text()
    _assert_same_as_ncss(tmp_path, plan, path)


def test_subset_follows_an_absolute_redirect_to_the_same_host(tmp_path):
    def location(url):
        parsed = urllib.parse.urlparse(url)
        return None if parsed.path.endswith("/") else f"https://gdex.ucar.edu{parsed.path}/?{parsed.query}"

    fake = FakeSubset(redirect_to=location)
    plan, path, _ = _subset(tmp_path, fake)
    assert "/api/status/960001" in [p for _, p in fake.calls] and fake.tokens == {TOKEN}
    _assert_same_as_ncss(tmp_path, plan, path)


@pytest.mark.parametrize("location", [
    "https://evil.example.org/api/status/960001/?token={token}",
    "http://gdex.ucar.edu/api/status/960001/?token={token}",          # downgrade to plain HTTP
    "https://gdex.ucar.edu.evil.example.org/api/status/960001/",
    "//evil.example.org/api/status/960001/",
])
def test_subset_never_follows_a_redirect_to_another_host_or_scheme(tmp_path, caplog, location):
    fake = FakeSubset(redirect_to=lambda url: location.format(token=TOKEN))
    with caplog.at_level(logging.DEBUG), pytest.raises(g.NcarRequestError, match="not followed") as error:
        _subset(tmp_path, fake)
    assert TOKEN not in str(error.value) and TOKEN not in caplog.text
    assert all(path.startswith("/api/") for _, path in fake.calls)   # the fake fails on any other host
    assert fake.submits() == 0


def test_subset_stops_after_too_many_redirects(tmp_path):
    fake = FakeSubset(redirect_to=lambda url: urllib.parse.urlparse(url).path + "x/")
    with pytest.raises(g.NcarRequestError, match="redirect"):
        _subset(tmp_path, fake)
    assert len(fake.calls) == 1 + gs.SUBSET_MAX_REDIRECTS


# --- D-026 defect 2: boxes crossing Greenwich
def test_subset_never_submits_a_box_crossing_greenwich(tmp_path):
    fake = FakeSubset()
    plan, path, record = _subset(tmp_path, fake, area=POLAR_GREENWICH, purge=True)
    assert fake.submits() == 2
    boxes = [(float(c["wlon"]), float(c["elon"])) for c in fake.controls]
    assert boxes == [(-2.0, -0.25), (0.0, 3.0)]
    assert all(w % 360 <= e % 360 for w, e in boxes)
    # one request per segment, both created and purged by this retrieval
    method = record["method"]
    assert [r["request_id"] for r in method["requests"]] == ["960001", "960002"]
    assert all(r["created_here"] for r in method["requests"])
    assert [p for m, p in fake.calls if m == "DELETE"] == ["/api/purge/960001", "/api/purge/960002"]
    assert [p["listed_after"] for p in method["purge"]["requests"]] == ["Set for Purge"] * 2
    _assert_same_as_ncss(tmp_path, plan, path)
    _assert_contiguous_longitudes(plan, path)


def test_subset_controls():
    plan = g.plan_retrieval(POLAR, DAYS)
    assert gs.subset_controls(plan) == [gs.subset_control(plan)]          # one box: the request is unchanged
    assert gs.subset_control(plan) == {
        "dataset": "d633000", "groupindex": "26", "date": "201808091200/to/201808100600",
        "param": "T/Z/W/U/V", "product": "Analysis", "oformat": "netCDF",
        "nlat": "90", "slat": "88", "wlon": "-50", "elon": "-48"}
    crossing = g.plan_retrieval(POLAR_GREENWICH, DAYS)
    west, east = gs.subset_controls(crossing)
    assert (west["wlon"], west["elon"], east["wlon"], east["elon"]) == ("-2", "-0.25", "0", "3")
    assert {k: v for k, v in west.items() if k not in ("wlon", "elon")} == \
        {k: v for k, v in gs.subset_control(crossing).items() if k not in ("wlon", "elon")}


def test_subset_across_greenwich_needs_room_for_both_requests(tmp_path):
    fake = FakeSubset(listed=[{"request_index": i, "status": "Completed", "rinfo": "dsnum=x"} for i in range(6)])
    with pytest.raises(g.NcarRequestError, match="8 open requests"):
        _subset(tmp_path, fake, area=POLAR_GREENWICH)
    assert fake.submits() == 0


def test_subset_across_greenwich_resumes_with_the_same_requests(tmp_path):
    fake = FakeSubset()
    plan = g.plan_retrieval(POLAR_GREENWICH, DAYS)
    original = fake._make_files

    def make_and_break(rid, control):
        original(rid, control)
        if rid == "960002":
            fake.failing_file = list(fake.files)[-1]

    fake._make_files = make_and_break
    with pytest.raises(g.NcarTransferError):
        g.retrieve(plan, str(tmp_path / "s.nc"), _names(), logger=LOGGER,
                   source=gs.SubsetSource(token=TOKEN, opener=fake.opener))
    parts = tmp_path / "s.nc.parts"
    assert json.loads((parts / gs.SUBSET_STATE).read_text())["request_id"] == "960001"
    assert json.loads((parts / "subset_request.segment1.json").read_text())["request_id"] == "960002"
    fake.failing_file = None
    g.retrieve(plan, str(tmp_path / "s.nc"), _names(), logger=LOGGER,
               source=gs.SubsetSource(token=TOKEN, opener=fake.opener))
    assert fake.submits() == 2
    _assert_same_as_ncss(tmp_path, plan, tmp_path / "s.nc")


def test_subset_without_greenwich_keeps_one_request_and_its_record(tmp_path):
    fake = FakeSubset()
    plan, path, record = _subset(tmp_path, fake)
    assert fake.controls == [gs.subset_control(plan)] and fake.submits() == 1
    assert "requests" not in record["method"] and record["method"]["request_id"] == "960001"
    assert set(record["method"]) >= {"request_id", "created_here", "reused_existing_request", "submit_answer_lost",
                                     "files", "listed_bytes", "remote_status"}


def test_subset_purges_only_its_own_request_and_records_the_listing(tmp_path):
    fake = FakeSubset()
    _, _, record = _subset(tmp_path, fake, purge=True)
    purge = record["method"]["purge"]
    assert [p for m, p in fake.calls if m == "DELETE"] == ["/api/purge/960001"]
    assert purge["listed_after"] == "Set for Purge" and purge["removed_from_listing"] is False
    saved = json.loads((tmp_path / "subset.nc.provenance.json").read_text())
    assert saved["method"]["purge"] == purge


def test_an_identical_completed_request_is_reused_and_never_purged(tmp_path):
    plan = g.plan_retrieval(POLAR, DAYS)
    fake, entry = _listed_like(plan, "Completed")
    fake.calls.clear()
    plan, path, record = _subset(tmp_path, fake, purge=True)
    assert fake.submits() == 0
    assert record["method"]["reused_existing_request"] is True
    assert record["method"]["purge"]["done"] is False
    assert not any(m == "DELETE" for m, _ in fake.calls)
    _assert_same_as_ncss(tmp_path, plan, path)


def test_an_identical_request_waiting_for_purge_stops_without_submitting(tmp_path):
    plan = g.plan_retrieval(POLAR, DAYS)
    fake, entry = _listed_like(plan, "Set for Purge")
    fake.listed = [entry]
    fake.calls.clear()
    with pytest.raises(g.NcarRequestError, match="deduplication"):
        _subset(tmp_path, fake)
    assert fake.submits() == 0


def test_nothing_is_submitted_with_too_many_listed_requests(tmp_path):
    fake = FakeSubset(listed=[{"request_index": i, "status": "Completed", "rinfo": "dsnum=x"} for i in range(7)])
    with pytest.raises(g.NcarRequestError, match="8 open requests"):
        _subset(tmp_path, fake)
    assert fake.submits() == 0


def test_a_lost_submission_answer_adopts_the_new_request_without_resubmitting_or_purging(tmp_path):
    fake = FakeSubset(lose_submit=True)
    plan, path, record = _subset(tmp_path, fake, purge=True)
    assert fake.submits() == 1
    method = record["method"]
    assert method["request_id"] == "960001" and method["submit_answer_lost"] is True
    assert method["created_here"] is False and method["reused_existing_request"] is False
    assert method["purge"]["done"] is False and "origin not proven" in method["purge"]["reason"]
    assert not any(m == "DELETE" for m, _ in fake.calls)
    _assert_same_as_ncss(tmp_path, plan, path)


def test_resuming_with_every_piece_kept_makes_no_api_call_before_the_purge(tmp_path):
    fake = FakeSubset()
    plan = g.plan_retrieval(POLAR, DAYS)
    path = tmp_path / "s.nc"
    g.retrieve(plan, str(path), _names(), logger=LOGGER, keep_parts=True,
               source=gs.SubsetSource(token=TOKEN, opener=fake.opener))
    path.unlink()
    (tmp_path / "s.nc.provenance.json").unlink()
    fake.calls.clear()
    record = g.retrieve(plan, str(path), _names(), logger=LOGGER, keep_parts=True,
                        source=gs.SubsetSource(token=TOKEN, opener=fake.opener, purge=True))
    assert [m for m, _ in fake.calls] == ["DELETE", "GET"]          # purge, then the listing that follows
    assert record["method"]["request_id"] == "960001" and record["method"]["created_here"] is True
    assert record["method"]["purge"]["listed_after"] == "Set for Purge"
    # without the request state, nothing is purged and nothing fails
    path.unlink()
    (tmp_path / "s.nc.provenance.json").unlink()
    (tmp_path / "s.nc.parts" / gs.SUBSET_STATE).unlink()
    fake.calls.clear()
    record = g.retrieve(plan, str(path), _names(), logger=LOGGER,
                        source=gs.SubsetSource(token=TOKEN, opener=fake.opener, purge=True))
    assert fake.calls == [] and record["method"]["purge"]["done"] is False
    _assert_same_as_ncss(tmp_path, plan, path)


def test_a_failed_listing_after_the_purge_is_recorded_and_the_file_kept(tmp_path, caplog):
    fake = FakeSubset(listing_fails_after_purge=True)
    with caplog.at_level(logging.DEBUG):
        plan, path, record = _subset(tmp_path, fake, purge=True)
    purge = record["method"]["purge"]
    assert purge["answer"] == {"purge_successful": "true"} and "HTTP 503" in purge["error"]
    assert TOKEN not in purge["error"] and TOKEN not in caplog.text
    assert "removed_from_listing" not in purge
    _assert_same_as_ncss(tmp_path, plan, path)


def test_an_interrupted_download_resumes_with_the_same_request(tmp_path):
    fake = FakeSubset()
    plan = g.plan_retrieval(POLAR, DAYS)
    first = gs.SubsetSource(token=TOKEN, opener=fake.opener)
    fake.failing_file = "pending"
    original = fake._make_files

    def make_and_break(rid, control):
        original(rid, control)
        fake.failing_file = list(fake.files)[-1]

    fake._make_files = make_and_break
    with pytest.raises(g.NcarTransferError):
        g.retrieve(plan, str(tmp_path / "s.nc"), _names(), logger=LOGGER, source=first)
    state = json.loads((tmp_path / "s.nc.parts" / gs.SUBSET_STATE).read_text())
    assert state["request_id"] == "960001" and len(state["files"]) == len(fake.files) - 1
    fake.failing_file = None
    downloads_before = sum(1 for _, p in fake.calls if "/download/" in p)
    second = gs.SubsetSource(token=TOKEN, opener=fake.opener)
    g.retrieve(plan, str(tmp_path / "s.nc"), _names(), logger=LOGGER, source=second)
    assert fake.submits() == 1
    assert sum(1 for _, p in fake.calls if "/download/" in p) - downloads_before == 1
    _assert_same_as_ncss(tmp_path, plan, tmp_path / "s.nc")


def test_the_token_file_must_be_private(tmp_path):
    path = tmp_path / "token"
    path.write_text(TOKEN + "\n")
    path.chmod(0o644)
    with pytest.raises(g.NcarRequestError, match="chmod 600") as error:
        gs.read_token(environ={}, path=str(path))
    assert TOKEN not in str(error.value)
    path.chmod(0o600)
    assert gs.read_token(environ={}, path=str(path)) == TOKEN
    assert gs.read_token(environ={gs.TOKEN_ENV: " env-token "}, path=str(path)) == "env-token"
    with pytest.raises(g.NcarRequestError, match="GDEX_API_TOKEN"):
        gs.read_token(environ={}, path=str(tmp_path / "missing"))


def test_rinfo_matching():
    plan = g.plan_retrieval(POLAR, DAYS)
    control = gs.subset_control(plan)
    assert control["date"] == "201808091200/to/201808100600"
    # the format GDEX lists (D-Bench records)
    rinfo = ("dsnum=d633000;startdate=2018-08-09 12:00;enddate=2018-08-10 06:00;"
             "parameters=8!d633000:T,8!d633000:Z,8!d633000:U,8!d633000:V,8!d633000:W;product=1;level=;"
             "nlat=90.0;slat=88.0;wlon=-50.0;elon=-48.0;tindex=26;ofmt=netCDF")
    assert gs.rinfo_matches(rinfo, control)
    assert not gs.rinfo_matches(rinfo.replace("slat=88", "slat=87"), control)
    assert not gs.rinfo_matches(rinfo.replace(",8!d633000:W", ""), control)        # other fields
    assert not gs.rinfo_matches(rinfo.replace("level=;", "level=850;"), control)   # other levels
    assert not gs.rinfo_matches(rinfo.replace("ofmt=netCDF", "ofmt=grib"), control)
    assert not gs.rinfo_matches("garbage", control)


def test_token_redaction_covers_its_url_encoded_forms():
    token = "a+b/c=d"
    text = f"https://gdex.ucar.edu/api/status/?token={urllib.parse.quote_plus(token)} {token} " \
           f"{urllib.parse.quote(token, safe='')}"
    assert token not in gs._redact(text, token) and "%2B" not in gs._redact(text, token)

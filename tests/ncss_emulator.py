"""
In-process emulation of the GDEX THREDDS NetCDF Subset Service (test helper).

``FakeGdex`` answers NCSS URLs of the GDEX ERA5 pressure-level files the way
the real service answered the P2-A pilot of 2026-10-26: NetCDF-3 classic, a
float32 field with dimensions (time, level, latitude, longitude), latitudes
from north to south, longitudes in 0..359.75 whatever the sign of the
requested limits, and the attributes of the real responses. The values come
from smooth analytic functions of time, level, latitude and longitude, so any
area and period can be requested and the expected values recomputed exactly
with :func:`field_values`.

Faults can be scripted per field and day (HTTP errors, timeouts, connection
resets, truncated or oversized bodies, web pages, wrong grids, missing
values). The opener is called in-process, so the network guard of the test
suite is not involved.
"""

import datetime
import email.message
import io
import re
import threading
import time
import urllib.error
import urllib.parse

import numpy as np
import scipy.io

LEVELS = (
    1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 125.0, 150.0,
    175.0, 200.0, 225.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0,
    650.0, 700.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0,
    975.0, 1000.0,
)
FILL = np.float32(9.999e20)
ORIGIN = datetime.datetime(1900, 1, 1)
# name -> (file code, grid tag, units, long name, ECMWF parameter)
FIELDS = {
    "T": ("128_130_t", "ll025sc", "K", "Temperature", 130),
    "Z": ("128_129_z", "ll025sc", "m**2 s**-2", "Geopotential", 129),
    "W": ("128_135_w", "ll025sc", "Pa s**-1", "Vertical velocity", 135),
    "U": ("128_131_u", "ll025uv", "m s**-1", "U component of wind", 131),
    "V": ("128_132_v", "ll025uv", "m s**-1", "V component of wind", 132),
}
PATH = re.compile(
    r"/(?P<month>\d{6})/e5\.oper\.an\.pl\.(?P<code>128_\d{3}_[a-z])\.(?P<grid>ll025(?:sc|uv))\."
    r"(?P<day>\d{8})00_(?P=day)23\.nc$"
)


def _standard_atmosphere(p):
    """Temperature (K) and geopotential height (m) of a simple standard atmosphere."""
    p = np.asarray(p, dtype=np.float64)
    temperature = np.maximum(216.65, 288.15 * (p / 1013.25) ** 0.190263)
    height = np.where(
        p >= 226.32,
        44330.8 * (1.0 - (p / 1013.25) ** 0.190263),
        11000.0 + 6341.6 * np.log(226.32 / p),
    )
    return temperature, height


def field_values(name, moments, latitudes, longitudes):
    """Values (float32, time x level x lat x lon) of the analytic field ``name``."""
    hours = np.array([(m - ORIGIN).total_seconds() / 3600.0 for m in moments])[:, None, None, None]
    p = np.asarray(LEVELS)[None, :, None, None]
    lat = np.deg2rad(np.asarray(latitudes, dtype=np.float64))[None, None, :, None]
    lon = np.deg2rad(np.asarray(longitudes, dtype=np.float64) % 360.0)[None, None, None, :]
    phase = 2.0 * np.pi * (hours % 24.0) / 24.0
    temperature, height = _standard_atmosphere(p)
    if name == "T":
        values = temperature + 8.0 * np.cos(lat) + 1.5 * np.sin(2.0 * lon) * np.cos(lat) + 0.5 * np.sin(phase)
    elif name == "Z":
        values = 9.80665 * height + 400.0 * np.cos(lat) + 60.0 * np.sin(2.0 * lon + phase) * np.cos(lat)
    elif name == "U":
        values = 15.0 * np.cos(lat) * (1.0 - p / 1100.0) + 3.0 * np.sin(lon + phase)
    elif name == "V":
        values = 4.0 * np.sin(2.0 * lon) * np.cos(lat) + 0.0 * p + 0.2 * np.cos(phase)
    elif name == "W":
        values = 0.05 * np.sin(lon + phase) * np.cos(lat) * (p / 1000.0)
    else:
        raise KeyError(name)
    shape = (len(moments), len(LEVELS), len(latitudes), len(longitudes))
    return np.broadcast_to(values, shape).astype(np.float32)


def grid_latitudes(north, south):
    count = int(round((north - south) / 0.25)) + 1
    return np.array([north - 0.25 * i for i in range(count)])


def grid_longitudes(west, east):
    """GDEX longitudes (0..359.75) between two limits, as NCSS selects them."""
    start = west % 360.0
    stop = start + (east - west)
    count = int(round((stop - start) / 0.25)) + 1
    return np.array([start + 0.25 * i for i in range(count)])


def netcdf3_response(name, moments, latitudes, longitudes, values=None, variable_name=None):
    """
    Bytes of a NetCDF-3 classic file shaped like a GDEX NCSS response.

    Written with scipy's pure-Python NetCDF-3 writer: the netCDF-C library is
    not thread-safe, and the emulator answers requests from several threads
    while the retrieval reads responses with netCDF-C.
    """
    code, grid, unit, long_name, parameter = FIELDS[name]
    values = field_values(name, moments, latitudes, longitudes) if values is None else values
    buffer = io.BytesIO()
    nc = scipy.io.netcdf_file(buffer, "w", version=1)
    nc.createDimension("level", len(LEVELS))
    nc.createDimension("latitude", len(latitudes))
    nc.createDimension("time", len(moments))
    nc.createDimension("longitude", len(longitudes))
    coordinates = (
        ("level", "d", LEVELS, {"long_name": "pressure level", "short_name": "plev", "units": "hPa",
                                "alternate_units": "millibar"}),
        ("latitude", "d", latitudes, {"long_name": "latitude", "short_name": "lat", "units": "degrees_north",
                                      "qualifier": "Gaussian", "standard_name": "latitude"}),
        ("time", "i", [int((m - ORIGIN).total_seconds() // 3600) for m in moments],
         {"long_name": "time", "units": "hours since 1900-01-01 00:00:00", "calendar": "gregorian"}),
        ("longitude", "d", longitudes, {"long_name": "longitude", "short_name": "lon", "units": "degrees_east",
                                        "standard_name": "longitude"}),
    )
    for coordinate, kind, data, attributes in coordinates:
        variable = nc.createVariable(coordinate, kind, (coordinate,))
        variable[:] = np.asarray(data)
        for key, value in attributes.items():
            setattr(variable, key, value)
    variable = nc.createVariable(variable_name or name, "f", ("time", "level", "latitude", "longitude"))
    variable[:] = np.asarray(values, dtype=np.float32)
    attributes = {
        "long_name": long_name, "short_name": name.lower(), "units": unit,
        "original_format": "WMO GRIB 1 with ECMWF local table", "ecmwf_local_table": np.int32(128),
        "ecmwf_parameter": np.int32(parameter), "_FillValue": FILL, "missing_value": FILL,
        "rda_dataset": "ds633.0", "rda_dataset_doi": "DOI: 10.5065/BH6N-5N20",
        "coordinates": "time level latitude longitude ",
    }
    for key, value in attributes.items():
        setattr(variable, key, value)
    nc.DATA_SOURCE = "ECMWF: https://cds.climate.copernicus.eu, Copernicus Climate Data Store"
    nc.Conventions = "CF-1.6"
    nc.emulator = "tests/ncss_emulator.py"
    nc.flush()
    body = buffer.getvalue()
    nc.close()
    return body


class FakeResponse(io.BytesIO):
    """Response object with the attributes the retrieval uses (status, headers, read)."""

    def __init__(self, body, headers, status=200, cut_at=None, fail_at=None):
        super().__init__(body)
        self.status = status
        self.headers = email.message.Message()
        for key, value in headers.items():
            self.headers[key] = str(value)
        self._cut_at = cut_at
        self._fail_at = fail_at
        self._fail_with = ConnectionResetError(54, "Connection reset by peer (emulated)")

    def read(self, size=-1):
        position = self.tell()
        if self._fail_at is not None and position >= self._fail_at:
            raise self._fail_with
        if self._cut_at is not None:
            if position >= self._cut_at:
                return b""
            size = self._cut_at - position if size is None or size < 0 else min(size, self._cut_at - position)
        if self._fail_at is not None:
            size = self._fail_at - position if size is None or size < 0 else min(size, self._fail_at - position)
        return super().read(size)


class FakeGdex:
    """
    NCSS emulator. Call :meth:`opener` like ``urllib.request.urlopen``.

    Args:
        faults (dict): "<field>_<YYYYMMDD>" -> list of faults, one used per
            request of that field and day (``None`` = normal answer).
        last_day (datetime.date): Days after it answer HTTP 404 (not published).
        delay (callable): Seconds to wait before answering a request, given the
            number of times requested (simulated server time).
    """

    def __init__(self, faults=None, last_day=None, delay=None):
        self.faults = {key: list(value) for key, value in (faults or {}).items()}
        self.last_day = last_day
        self.delay = delay
        self.calls = []
        self.lock = threading.Lock()

    # ------------------------------------------------------------------ core
    def respond(self, url):
        """Return (status, headers, body, fault) for an NCSS URL."""
        parsed = urllib.parse.urlsplit(url)
        query = dict(urllib.parse.parse_qsl(parsed.query))
        match = PATH.search(parsed.path)
        if match is None:
            return 404, {"Content-Type": "text/plain"}, b"Not Found", None
        name = query.get("var")
        if name not in FIELDS or FIELDS[name][:2] != (match["code"], match["grid"]):
            return 400, {"Content-Type": "text/plain"}, f"variable {name} not in this dataset".encode(), None
        day = datetime.datetime.strptime(match["day"], "%Y%m%d").date()
        key = f"{name}_{day:%Y%m%d}"
        with self.lock:
            fault = self.faults[key].pop(0) if self.faults.get(key) else None
            self.calls.append((url, fault, threading.current_thread().name))
        if self.last_day is not None and day > self.last_day:
            return 404, {"Content-Type": "text/plain"}, b"Not Found", fault
        if query.get("accept") != "netcdf3":
            return 400, {"Content-Type": "text/plain"}, b"accept netcdf3 only", fault
        north, south = float(query["north"]), float(query["south"])
        west, east = float(query["west"]), float(query["east"])
        if (west % 360.0) + (east - west) >= 360.0:
            return 400, {"Content-Type": "text/plain"}, b"seam crossing not emulated", fault
        start = datetime.datetime.strptime(query["time_start"], "%Y-%m-%dT%H:%M:%SZ")
        end = datetime.datetime.strptime(query["time_end"], "%Y-%m-%dT%H:%M:%SZ")
        stride = int(query.get("timeStride", 1))
        if start.date() != day or end.date() != day or end < start:
            return 400, {"Content-Type": "text/plain"}, b"time range outside the file", fault
        moments = []
        moment = start
        while moment <= end:
            moments.append(moment)
            moment += datetime.timedelta(hours=stride)
        latitudes = grid_latitudes(north, south)
        longitudes = grid_longitudes(west, east)
        if self.delay is not None:
            time.sleep(self.delay(len(moments)))
        values = None
        variable_name = None
        if fault == "wrong-grid":
            latitudes = latitudes[:-1]
        elif fault in ("missing", "nan"):
            values = field_values(name, moments, latitudes, longitudes).copy()
            values[0, 0, 0, 0] = FILL if fault == "missing" else np.nan
        elif fault == "wrong-var":
            variable_name = "Q"
        body = netcdf3_response(name, moments, latitudes, longitudes, values, variable_name)
        headers = {
            "Content-Type": "application/x-netcdf",
            "Content-Length": len(body),
            "Date": "Mon, 26 Oct 2026 21:05:32 GMT",
            "Content-Disposition": f"attachment; filename={parsed.path.rsplit('/', 1)[-1]}",
        }
        return 200, headers, body, fault

    # ---------------------------------------------------------------- opener
    def opener(self, request, timeout=None):
        url = request.full_url if hasattr(request, "full_url") else request
        status, headers, body, fault = self.respond(url)
        if fault == "timeout":
            raise TimeoutError("timed out (emulated)")
        if fault == "refused":
            raise urllib.error.URLError(ConnectionRefusedError(61, "Connection refused (emulated)"))
        for code in (500, 502, 503, 504, 429, 404, 400, 403, 302):
            if fault == f"http{code}":
                status, body = code, f"emulated HTTP {code}".encode()
                headers = {"Content-Type": "text/plain", "Retry-After": "0"}
                if code == 302:
                    headers["Location"] = "https://example.invalid/login"
        if status != 200:
            message = email.message.Message()
            for key, value in headers.items():
                message[key] = str(value)
            raise urllib.error.HTTPError(url, status, f"emulated {status}", message, io.BytesIO(body))
        if fault == "html":
            page = b"<html><body>Scheduled maintenance</body></html>"
            return FakeResponse(page, {"Content-Type": "text/html", "Content-Length": len(page)})
        if fault == "garbage":
            junk = b"\x00" * len(body)
            return FakeResponse(junk, headers)
        if fault == "oversize":
            headers = dict(headers, **{"Content-Length": 10 * len(body) + 10 ** 7})
            return FakeResponse(body, headers)
        if fault == "truncate":
            return FakeResponse(body, headers, cut_at=len(body) // 2)
        if fault in ("reset", "read-timeout"):
            response = FakeResponse(body, headers, fail_at=len(body) // 3)
            response._fail_with = (ConnectionResetError(54, "Connection reset by peer (emulated)")
                                   if fault == "reset" else TimeoutError("read timed out (emulated)"))
            return response
        if fault == "no-length":
            headers = {k: v for k, v in headers.items() if k != "Content-Length"}
            return FakeResponse(body, headers)
        return FakeResponse(body, headers)

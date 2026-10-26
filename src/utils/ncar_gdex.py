# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    ncar_gdex.py                                       :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2026/10/26 23:00:00 by daniloceano       #+#    #+#              #
#    Updated: 2026/10/26 23:00:00 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
ERA5 pressure-level data from the NCAR Geoscience Data Exchange (optional).

Alternative to the CDS API route (:func:`src.utils.tools.get_cdsapi_data`)
that reads the same ERA5 fields from the NCAR GDEX copy of ERA5 (dataset
d633000, DOI 10.5065/BH6N-5N20) through its THREDDS NetCDF Subset Service
(NCSS). The request is the one the CDS route makes for the same track and
``--time-resolution``: the five fields, the 37 pressure levels, the track
extent with a 15-degree buffer rounded outwards to whole degrees, and the same
hours of each day. NCSS cuts the area on the server and answers on the same
connection, without a queue or an account. GDEX keeps one file per field and
day, so each request is one field for one day (two requests when the area
crosses the Greenwich meridian, where the 0-359.75 longitudes of the GDEX grid
wrap).

Every response is checked before use: HTTP status, number of bytes, NetCDF
signature, variable, units, coordinates (times, levels, latitudes and
longitudes exactly those requested) and absence of missing values. The values
are copied unchanged (float32, no interpolation or unit conversion) into one
NetCDF file whose variable and coordinate names are those of the namelist used
by the run. The file gets its final name only when complete and checked, next
to a provenance record ``<file>.provenance.json`` with the URL, size and
SHA-256 of every response and of the file. An existing file is reused only if
its provenance record matches its SHA-256 and covers the request; otherwise
the program stops and says why.

Responses are kept in ``<file>.parts/`` until the file is complete, so an
interrupted retrieval resumes from the verified responses already received.

GDEX publishes ERA5 some months after the CDS: dates that are not available
yet stop the retrieval with a message that suggests ``--cdsapi``.

Run on its own (from the repository root) to retrieve a file without running
the energetics, or to see the requests first with ``--dry-run``::

    python -m src.utils.ncar_gdex OUTFILE.nc --trackfile inputs/track --dry-run

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astronomia, Geofísica e Ciências Atmosféricas
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import argparse
import contextlib
import datetime
import hashlib
import http.client
import json
import logging
import math
import os
import platform
import shutil
import ssl
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from metpy.units import units

from .input_files import resolve_input_file

NCSS_BASE_URL = "https://tds.gdex.ucar.edu/thredds/ncss/grid/files/d633000/e5.oper.an.pl"
SOURCE = {
    "provider": "NCAR Geoscience Data Exchange (GDEX)",
    "dataset": "d633000",
    "title": "ERA5 Reanalysis (0.25 Degree Latitude-Longitude Grid)",
    "doi": "10.5065/BH6N-5N20",
    "license": "CC-BY-4.0",
    "service": "THREDDS NetCDF Subset Service (NCSS)",
}
PROVENANCE_SCHEMA = "lectk-era5-retrieval/1"
USER_AGENT = "LorenzCycleToolkit-ncar-gdex/1 (+https://github.com/daniloceano/LorenzCycleToolkit)"

# Same buffer as get_cdsapi_data.
BUFFER_DEGREES = 15
GRID_STEP = 0.25
PRESSURE_LEVELS_HPA = (
    1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 125.0, 150.0,
    175.0, 200.0, 225.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0,
    650.0, 700.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0,
    975.0, 1000.0,
)
TIME_UNITS = "hours since 1900-01-01 00:00:00"
TIME_ORIGIN = datetime.datetime(1900, 1, 1)
DIMENSIONS = ("time", "level", "latitude", "longitude")

# A NetCDF-3 response holds the float32 values plus coordinates and header.
HEADER_ALLOWANCE = 64 * 1024
TYPICAL_HEADER_BYTES = 8 * 1024
NETCDF_SIGNATURES = (b"CDF\x01", b"CDF\x02", b"CDF\x05", b"\x89HDF")

DEFAULT_WORKERS = 1
MAX_WORKERS = 4
DEFAULT_TIMEOUT = 300.0
DEFAULT_ATTEMPTS = 4
DEFAULT_BACKOFF = 10.0
MAX_WAIT = 300.0
DEFAULT_MAX_TOTAL_BYTES = 8 * 1024 ** 3
RETRY_STATUS = {429, 500, 502, 503, 504}
# The netCDF-C library is not thread-safe: responses are checked one at a time.
_NETCDF_LOCK = threading.Lock()
NETWORK_ERRORS = (urllib.error.URLError, TimeoutError, ConnectionError, ssl.SSLError,
                  http.client.HTTPException)


@dataclass(frozen=True)
class SourceField:
    """One field of the GDEX pressure-level files."""

    row: str
    name: str
    code: str
    grid: str
    units: str


FIELDS = (
    SourceField("Air Temperature", "T", "128_130_t", "ll025sc", "K"),
    SourceField("Geopotential", "Z", "128_129_z", "ll025sc", "m**2 s**-2"),
    SourceField("Omega Velocity", "W", "128_135_w", "ll025sc", "Pa s**-1"),
    SourceField("Eastward Wind Component", "U", "128_131_u", "ll025uv", "m s**-1"),
    SourceField("Northward Wind Component", "V", "128_132_v", "ll025uv", "m s**-1"),
)
COORDINATE_ROWS = (
    ("time", "Time"),
    ("level", "Vertical Level"),
    ("latitude", "Latitude"),
    ("longitude", "Longitude"),
)


class NcarRetrievalError(RuntimeError):
    """The GDEX retrieval could not produce a verified file."""


class NcarRequestError(NcarRetrievalError):
    """The service refused a request (not retried)."""


class NcarTransferError(NcarRetrievalError):
    """A network or server failure persisted after the retries."""


class NcarContentError(NcarRetrievalError):
    """A response is not the data that was requested."""


class NcarReuseError(NcarRetrievalError):
    """An existing file cannot be verified as the data of this request."""


class _Retryable(Exception):
    """Failure of one attempt that may succeed if repeated."""

    def __init__(self, message, retry_after=None):
        super().__init__(message)
        self.retry_after = retry_after


class _Cancelled(Exception):
    """The retrieval was stopped while a response was being received."""


def _describe(error):
    """Text of a failed attempt for the user (no internal class names)."""
    return str(error) if isinstance(error, _Retryable) else f"{type(error).__name__}: {error}"


# --------------------------------------------------------------------- plan
@dataclass(frozen=True)
class LongitudeSegment:
    """One NCSS longitude interval and its place in the output file."""

    west: float
    east: float
    response: tuple
    output: tuple
    offset: int


@dataclass(frozen=True)
class PieceRequest:
    """One NCSS request: one field, one day, one longitude segment."""

    field: SourceField
    day: datetime.date
    hours: tuple
    segment: LongitudeSegment
    segment_index: int
    latitudes: tuple
    time_offset: int
    url: str

    @property
    def shape(self):
        return (len(self.hours), len(PRESSURE_LEVELS_HPA), len(self.latitudes), len(self.segment.response))

    @property
    def n_values(self):
        return int(np.prod(self.shape))

    @property
    def max_bytes(self):
        return 4 * self.n_values + HEADER_ALLOWANCE

    @property
    def key(self):
        return f"{self.field.name}_{self.day:%Y%m%d}_{self.segment_index}"


@dataclass(frozen=True)
class RetrievalPlan:
    """Every request of a retrieval and the coordinates of the output file."""

    north: float
    west: float
    south: float
    east: float
    times: tuple
    latitudes: tuple
    longitudes: tuple
    requests: tuple
    base_url: str

    @property
    def estimated_bytes(self):
        return sum(4 * r.n_values + TYPICAL_HEADER_BYTES for r in self.requests)

    @property
    def output_bytes(self):
        return 4 * len(FIELDS) * len(self.times) * len(PRESSURE_LEVELS_HPA) * len(self.latitudes) * len(self.longitudes)


def area_from_track(track):
    """
    Area of the CDS API request for a track, as [north, west, south, east].

    Same rule as :func:`src.utils.tools.get_cdsapi_data`: the extent of the
    track positions with a 15-degree buffer, rounded outwards to whole degrees.

    Args:
        track (pd.DataFrame): Track with 'Lat' and 'Lon' columns.

    Returns:
        tuple: (north, west, south, east) in degrees (integers).
    """
    min_lat, max_lat = float(track["Lat"].min()), float(track["Lat"].max())
    min_lon, max_lon = float(track["Lon"].min()), float(track["Lon"].max())
    return (
        math.ceil(max_lat + BUFFER_DEGREES),
        math.floor(min_lon - BUFFER_DEGREES),
        math.floor(min_lat - BUFFER_DEGREES),
        math.ceil(max_lon + BUFFER_DEGREES),
    )


def daily_hours(times, time_step):
    """
    Hours of each day that the CDS API route requests for a track.

    Same rule as :func:`src.utils.tools.get_cdsapi_data`: every day from the
    first to the last track date; the first day starts at the first track hour
    rounded down to a multiple of ``time_step``, the last day ends with the
    last track hour (included), and the days in between are complete.

    Args:
        times (sequence): Track times.
        time_step (int): Time resolution in hours (``--time-resolution``).

    Returns:
        list: (datetime.date, tuple of hours) for each day with hours.

    Raises:
        ValueError: If ``time_step`` is not a positive whole number of hours.
    """
    if int(time_step) != time_step or time_step < 1:
        raise ValueError(f"time resolution must be a positive whole number of hours, not {time_step!r}")
    time_step = int(time_step)
    index = pd.DatetimeIndex(times)
    first, last = index.min(), index.max()
    days = pd.date_range(start=first.date(), end=last.date(), freq="D")
    plan = []
    for day in days:
        if day.date() == first.date():
            start = (first.hour // time_step) * time_step
            if day.date() == last.date():
                end = ((last.hour // time_step) + 1) * time_step
                if end <= start:
                    end = start + time_step
                end = min(end, 24)
            else:
                end = 24
        elif day.date() == last.date():
            start = 0
            end = ((last.hour // time_step) + 1) * time_step
            if end == 0:
                end = time_step
            end = min(end, 24)
        else:
            start, end = 0, 24
        hours = tuple(range(start, end, time_step))
        if hours:
            plan.append((day.date(), hours))
    return plan


def _points(start, stop):
    """Grid values from start to stop (both included), GRID_STEP apart."""
    step = GRID_STEP if stop >= start else -GRID_STEP
    count = int(round((stop - start) / step)) + 1
    return tuple(start + step * i for i in range(count))


def _on_grid(value):
    return math.isclose(value / GRID_STEP, round(value / GRID_STEP), abs_tol=1e-9)


def longitude_segments(west, east):
    """
    Split a longitude interval into the NCSS requests of the GDEX grid.

    The GDEX grid runs from 0 to 359.75 degrees east. An interval inside it is
    one request; an interval that crosses the Greenwich meridian is two, and
    its output longitudes run from -180 to 180 so that they stay increasing.
    Intervals entirely west of Greenwich are requested with negative
    longitudes, as in the P2-A pilot. Intervals crossing 180 degrees are
    refused: the toolkit works in -180..180 and needs one contiguous arc.

    Args:
        west (float): Western limit in degrees (any convention).
        east (float): Eastern limit in degrees, greater than ``west``.

    Returns:
        tuple: One or two :class:`LongitudeSegment`.

    Raises:
        ValueError: For an empty, reversed, full-circle or 180-crossing
            interval, or limits that are not on the 0.25-degree grid.
    """
    width = east - west
    if not 0 < width < 360:
        raise ValueError(
            f"longitude interval [{west}, {east}] must have the western limit first and "
            "span less than 360 degrees"
        )
    if not (_on_grid(west) and _on_grid(east)):
        raise ValueError(f"longitude limits [{west}, {east}] are not on the 0.25-degree grid")
    if ((west + 180.0) % 360.0) - 180.0 + width > 180.0:
        raise ValueError(
            f"longitude interval [{west}, {east}] crosses 180 degrees: the toolkit works in "
            "-180..180 and needs one contiguous arc, so this area cannot be processed"
        )
    start = west % 360.0
    stop = start + width
    if stop < 360.0:
        response = _points(start, stop)
        if start >= 180.0:
            return (LongitudeSegment(start - 360.0, stop - 360.0, response, response, 0),)
        return (LongitudeSegment(start, stop, response, response, 0),)
    west_part = _points(start, 360.0 - GRID_STEP)
    east_part = _points(0.0, stop - 360.0)
    return (
        LongitudeSegment(start - 360.0, -GRID_STEP, west_part, tuple(lon - 360.0 for lon in west_part), 0),
        LongitudeSegment(0.0, stop - 360.0, east_part, east_part, len(west_part)),
    )


def _number(value):
    return f"{value + 0.0:g}"


def request_url(base_url, field, day, hours, north, south, west, east):
    """
    NCSS URL of one field, one day and one area.

    Args:
        base_url (str): NCSS address of the GDEX pressure-level files.
        field (SourceField): Field to request.
        day (datetime.date): Day (one GDEX file).
        hours (tuple): Hours of the day, equally spaced.
        north, south, west, east (float): Area, on the 0.25-degree grid.

    Returns:
        str: The URL.
    """
    stride = hours[1] - hours[0] if len(hours) > 1 else 1
    if any(b - a != stride for a, b in zip(hours, hours[1:])) or stride < 1:
        raise ValueError(f"hours {hours} are not equally spaced")
    stamp = f"{day:%Y%m%d}"
    path = f"{base_url}/{day:%Y%m}/e5.oper.an.pl.{field.code}.{field.grid}.{stamp}00_{stamp}23.nc"
    return (
        f"{path}?var={field.name}&north={_number(north)}&south={_number(south)}"
        f"&west={_number(west)}&east={_number(east)}&horizStride=1"
        f"&time_start={day:%Y-%m-%d}T{hours[0]:02d}:00:00Z"
        f"&time_end={day:%Y-%m-%d}T{hours[-1]:02d}:00:00Z"
        f"&timeStride={stride}&accept=netcdf3"
    )


def plan_retrieval(area, days_hours, base_url=NCSS_BASE_URL):
    """
    Plan the NCSS requests for an area and the hours of each day.

    Args:
        area (sequence): (north, west, south, east) in degrees, on the
            0.25-degree grid; latitudes beyond the poles are clipped.
        days_hours (list): (datetime.date, hours) pairs, as returned by
            :func:`daily_hours`.
        base_url (str): NCSS address of the GDEX pressure-level files.

    Returns:
        RetrievalPlan: The requests, in day, field and segment order.

    Raises:
        ValueError: For an empty or invalid area or time selection.
    """
    north, west, south, east = (float(v) for v in area)
    north, south = min(north, 90.0), max(south, -90.0)
    if not (_on_grid(north) and _on_grid(south)):
        raise ValueError(f"latitude limits [{south}, {north}] are not on the 0.25-degree grid")
    if north <= south:
        raise ValueError(f"the northern limit ({north}) must be north of the southern limit ({south})")
    if not days_hours:
        raise ValueError("no times to retrieve")
    segments = longitude_segments(west, east)
    latitudes = _points(north, south)
    times, requests = [], []
    for day, hours in days_hours:
        hours = tuple(int(h) for h in hours)
        if not hours or hours[0] < 0 or hours[-1] > 23:
            raise ValueError(f"hours {hours} of {day} are outside 0-23")
        for field in FIELDS:
            for index, segment in enumerate(segments):
                url = request_url(base_url, field, day, hours, north, south, segment.west, segment.east)
                requests.append(PieceRequest(field, day, hours, segment, index, latitudes, len(times), url))
        times.extend(_moments(day, hours))
    if len(set(times)) != len(times) or times != sorted(times):
        raise ValueError("the times to retrieve must be increasing and unique")
    longitudes = tuple(lon for segment in segments for lon in segment.output)
    return RetrievalPlan(north, west, south, east, tuple(times), latitudes, longitudes, tuple(requests), base_url)


# ------------------------------------------------------------------- names
def _same_units(declared, expected):
    try:
        magnitude = units.Quantity(1.0, str(declared)).to(expected).magnitude
    except Exception:
        return False
    return math.isclose(magnitude, 1.0, rel_tol=1e-12)


def output_names(namelist):
    """
    Variable and coordinate names of the output file, from a namelist.

    The output file is written with the names of the namelist that the run
    reads, so that the data are found under the names the run expects.

    Args:
        namelist (str): Path to the namelist (``;``-separated, as
            ``inputs/namelist``).

    Returns:
        dict: GDEX name (T, Z, W, U, V, time, level, latitude, longitude) ->
        output name.

    Raises:
        ValueError: If the namelist lacks a field or coordinate, describes
            geopotential height instead of geopotential, gives units other
            than those of the GDEX fields, or repeats a name.
    """
    table = pd.read_csv(namelist, sep=";", index_col=0, header=0)
    names = {}
    for field in FIELDS:
        if field.row not in table.index:
            if field.row == "Geopotential" and "Geopotential Height" in table.index:
                raise ValueError(
                    f"{namelist} describes geopotential height, but GDEX provides geopotential "
                    f"({field.units}); use a namelist with a 'Geopotential' row in m**2/s**2 "
                    "(e.g. inputs/namelist_ERA5)"
                )
            raise ValueError(f"{namelist} has no '{field.row}' row")
        declared = table.loc[field.row]["Units"]
        if not _same_units(declared, field.units):
            raise ValueError(
                f"{namelist} gives '{declared}' for '{field.row}', but the GDEX field {field.name} "
                f"is in {field.units}; the values would be read with the wrong units"
            )
        names[field.name] = str(table.loc[field.row]["Variable"]).strip()
    for coordinate, row in COORDINATE_ROWS:
        if row not in table.index:
            raise ValueError(f"{namelist} has no '{row}' row")
        names[coordinate] = str(table.loc[row]["Variable"]).strip()
    if len(set(names.values())) != len(names):
        raise ValueError(f"{namelist} gives the same name to two variables or coordinates: {names}")
    return names


# ------------------------------------------------------------------ fetch
@dataclass
class PieceResult:
    """A verified response stored in the parts directory."""

    request: PieceRequest
    path: str
    sha256: str
    size: int
    attempts: int
    elapsed: float
    http: dict
    reused: bool


class _RefuseRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def default_opener():
    """URL opener that does not follow redirects (a redirect is reported, not followed)."""
    return urllib.request.build_opener(_RefuseRedirect).open


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _error_text(error):
    try:
        text = error.read(2048).decode("utf-8", "replace")
    except Exception:
        return ""
    return " ".join(text.split())[:500]


def _retry_after(headers):
    try:
        return max(0.0, float(headers.get("Retry-After")))
    except (TypeError, ValueError, AttributeError):
        return None


def _refusal(code, reason, body, url):
    if code == 404:
        return (
            f"GDEX has no file for this request (HTTP 404): the date may not be published yet "
            f"(GDEX follows the CDS by some months; use --cdsapi for recent dates) or the file "
            f"was renamed. URL: {url}"
        )
    if code in (301, 302, 303, 307, 308):
        return f"the service redirected the request (HTTP {code}); redirects are not followed. URL: {url}"
    if code in (401, 403):
        return f"the service refused access (HTTP {code} {reason}); GDEX NCSS needs no account. URL: {url}"
    return f"the service refused the request (HTTP {code} {reason}): {body or 'no message'}. URL: {url}"


def _download(request, part, timeout, opener, cancel):
    """One attempt: stream the response to ``part`` and return (sha256, size, http)."""
    message = urllib.request.Request(request.url, headers={"User-Agent": USER_AGENT})
    try:
        response = opener(message, timeout=timeout)
    except urllib.error.HTTPError as error:
        body = _error_text(error)
        if error.code in RETRY_STATUS:
            raise _Retryable(f"HTTP {error.code} {error.reason}: {body or 'no message'}",
                             _retry_after(error.headers))
        raise NcarRequestError(_refusal(error.code, error.reason, body, request.url))
    digest, size = hashlib.sha256(), 0
    with response:
        headers = response.headers
        content_type = headers.get("Content-Type") or ""
        length = headers.get("Content-Length")
        try:
            expected = int(length) if length is not None else None
        except ValueError:
            expected = None
        if expected is not None and expected > request.max_bytes:
            raise NcarContentError(
                f"{request.key}: the service announced {expected} bytes, more than the "
                f"{request.max_bytes} bytes this request can hold"
            )
        with open(part, "wb") as handle:
            while True:
                if cancel.is_set():
                    raise _Cancelled()
                block = response.read(1 << 20)
                if not block:
                    break
                size += len(block)
                if size > request.max_bytes:
                    raise NcarContentError(
                        f"{request.key}: response larger than the {request.max_bytes} bytes this request can hold"
                    )
                digest.update(block)
                try:
                    handle.write(block)
                except OSError as error:
                    raise NcarRetrievalError(f"cannot write {part}: {error}") from error
        if expected is not None and size != expected:
            raise _Retryable(f"incomplete response: {size} of {expected} bytes")
        if not content_type.startswith("application/x-netcdf"):
            with open(part, "rb") as handle:
                start = handle.read(300).decode("utf-8", "replace")
            raise NcarContentError(
                f"{request.key}: the service answered '{content_type}' instead of NetCDF "
                f"(GDEX may be under maintenance): {' '.join(start.split())[:200]}"
            )
        status = getattr(response, "status", None)
        http = {
            "status": status,
            "content_type": content_type,
            "content_length": expected,
            "date": headers.get("Date"),
            "content_disposition": headers.get("Content-Disposition"),
        }
    with open(part, "rb") as handle:
        signature = handle.read(4)
    if signature not in NETCDF_SIGNATURES:
        raise _Retryable(f"response does not start with a NetCDF signature ({signature!r})")
    return digest.hexdigest(), size, http


def _hours_since_origin(moments):
    """Values of the GDEX time coordinate (whole hours since 1900-01-01) of datetimes."""
    return np.array([(m - TIME_ORIGIN) // datetime.timedelta(hours=1) for m in moments], dtype=np.int64)


def _moments(day, hours):
    return [datetime.datetime.combine(day, datetime.time(hour)) for hour in hours]


def check_piece(path, request):
    """
    Check that a stored response holds exactly the data of its request.

    Args:
        path (str): NetCDF file of the response.
        request (PieceRequest): The request.

    Raises:
        NcarContentError: If the variable, units, shape, coordinates or values
            differ from the request, or values are missing.
        _Retryable: If the file cannot be read as NetCDF.
    """
    with _NETCDF_LOCK:
        _check_piece(path, request)


def _check_piece(path, request):
    try:
        dataset = xr.open_dataset(path, decode_cf=False)
    except Exception as error:
        raise _Retryable(f"response is not a readable NetCDF file ({type(error).__name__}: {error})")
    with dataset:
        name = request.field.name
        if name not in dataset.variables:
            raise NcarContentError(f"{request.key}: variable {name} missing (found {sorted(dataset.data_vars)})")
        variable = dataset[name]
        problems = []
        if variable.dims != DIMENSIONS:
            problems.append(f"dimensions {variable.dims} instead of {DIMENSIONS}")
        if variable.shape != request.shape:
            problems.append(f"shape {variable.shape} instead of {request.shape}")
        if variable.dtype != np.float32:
            problems.append(f"type {variable.dtype} instead of float32")
        if variable.attrs.get("units") != request.field.units:
            problems.append(f"units {variable.attrs.get('units')!r} instead of {request.field.units!r}")
        checks = (
            ("level", PRESSURE_LEVELS_HPA, "hPa"),
            ("latitude", request.latitudes, "degrees_north"),
            ("longitude", request.segment.response, "degrees_east"),
        )
        for coordinate, expected, unit in checks:
            if coordinate not in dataset.variables:
                problems.append(f"coordinate {coordinate} missing")
                continue
            values = np.asarray(dataset[coordinate].values, dtype=np.float64)
            if values.shape != (len(expected),) or not np.array_equal(values, np.asarray(expected)):
                problems.append(f"{coordinate} values differ from the request")
            if dataset[coordinate].attrs.get("units") != unit:
                problems.append(f"{coordinate} units {dataset[coordinate].attrs.get('units')!r} instead of {unit!r}")
        if "time" not in dataset.variables:
            problems.append("coordinate time missing")
        else:
            times = dataset["time"]
            if times.attrs.get("units") != TIME_UNITS:
                problems.append(f"time units {times.attrs.get('units')!r} instead of {TIME_UNITS!r}")
            if times.attrs.get("calendar") not in ("gregorian", "standard"):
                problems.append(f"time calendar {times.attrs.get('calendar')!r}")
            expected = _hours_since_origin(_moments(request.day, request.hours))
            if times.shape != expected.shape or not np.array_equal(np.asarray(times.values, dtype=np.int64), expected):
                problems.append("times differ from the request")
        if not problems:
            values = variable.values
            missing = ~np.isfinite(values)
            for attribute in ("_FillValue", "missing_value"):
                if attribute in variable.attrs:
                    missing |= values == np.float32(variable.attrs[attribute])
            if missing.any():
                problems.append(f"{int(missing.sum())} missing values")
        if problems:
            raise NcarContentError(f"{request.key}: " + "; ".join(problems))


def _write_json(path, content):
    temporary = f"{path}.tmp"
    with open(temporary, "w", encoding="utf-8") as handle:
        json.dump(content, handle, indent=2, sort_keys=False)
        handle.write("\n")
    os.replace(temporary, path)


def _stored_piece(request, parts_dir):
    """Return the PieceResult of a verified response kept from an earlier attempt, or None."""
    path = os.path.join(parts_dir, f"{request.key}.nc")
    record_path = f"{path}.json"
    if not (os.path.isfile(path) and os.path.isfile(record_path)):
        return None
    try:
        with open(record_path, encoding="utf-8") as handle:
            record = json.load(handle)
        if record.get("url") != request.url or record.get("sha256") != _sha256(path):
            return None
        check_piece(path, request)
    except (OSError, ValueError, NcarContentError, _Retryable):
        return None
    return PieceResult(request, path, record["sha256"], record["bytes"], 0, 0.0, record.get("http", {}), True)


def _wait(seconds, cancel):
    if cancel.wait(seconds):
        raise _Cancelled()


def fetch_piece(request, parts_dir, opener, timeout, attempts, backoff, cancel, logger):
    """
    Obtain one verified response, from the parts directory or from NCSS.

    Retries network failures, incomplete responses and HTTP 429/5xx with
    increasing waits (or the server's Retry-After); a refused request (other
    4xx) or a response that is not the requested data stops at once.

    Returns:
        PieceResult: The verified response.

    Raises:
        NcarRequestError, NcarTransferError, NcarContentError; _Cancelled if
        the retrieval was stopped (no request is made after a failure).
    """
    if cancel.is_set():
        raise _Cancelled()
    try:
        return _fetch_piece(request, parts_dir, opener, timeout, attempts, backoff, cancel, logger)
    except _Cancelled:
        raise
    except BaseException:
        cancel.set()
        raise


def _fetch_piece(request, parts_dir, opener, timeout, attempts, backoff, cancel, logger):
    stored = _stored_piece(request, parts_dir)
    if stored is not None:
        logger.debug(f"♻️ {request.key}: verified response kept from an earlier attempt")
        return stored
    path = os.path.join(parts_dir, f"{request.key}.nc")
    part = f"{path}.part"
    for attempt in range(1, attempts + 1):
        started = time.monotonic()
        try:
            sha256, size, http = _download(request, part, timeout, opener, cancel)
            check_piece(part, request)
        except (_Retryable, *NETWORK_ERRORS) as error:
            if os.path.exists(part):
                os.remove(part)
            if attempt == attempts:
                raise NcarTransferError(
                    f"{request.key}: failed {attempts} times; last error: {_describe(error)}. URL: {request.url}"
                ) from error
            delay = min(MAX_WAIT, max(backoff * 3 ** (attempt - 1), getattr(error, "retry_after", None) or 0.0))
            logger.warning(
                f"⚠️ {request.key}: {_describe(error)}; retrying in {delay:.0f} s (attempt {attempt + 1}/{attempts})"
            )
            _wait(delay, cancel)
            continue
        except BaseException:
            if os.path.exists(part):
                os.remove(part)
            raise
        elapsed = time.monotonic() - started
        os.replace(part, path)
        _write_json(f"{path}.json", {"url": request.url, "sha256": sha256, "bytes": size, "http": http})
        return PieceResult(request, path, sha256, size, attempt, elapsed, http, False)
    raise AssertionError("unreachable")


# --------------------------------------------------------------- assemble
def _copy_attributes(source, target, skip=("_FillValue",)):
    for attribute in source.ncattrs():
        if attribute not in skip:
            target.setncattr(attribute, source.getncattr(attribute))


def write_output(plan, pieces, path, names, compress=True, float64=False):
    """
    Write the verified responses into one NetCDF file, values unchanged.

    The file is written piece by piece (one response in memory at a time).

    Args:
        plan (RetrievalPlan): The plan the responses belong to.
        pieces (list): PieceResult of every request of the plan.
        path (str): File to write.
        names (dict): Output names, from :func:`output_names`.
        compress (bool): NetCDF-4 with lossless zlib compression (level 1,
            shuffle, one time step per chunk); otherwise NetCDF-3 64-bit offset.
        float64 (bool): Store the float32 values as float64 (exact). The
            toolkit computes in the type of its input: results from float32
            and float64 copies of the same values differ by about 1e-7
            (relative), and the regression references of the test suite were
            computed from float64 input.
    """
    by_key = {piece.request.key: piece for piece in pieces}
    if set(by_key) != {request.key for request in plan.requests}:
        raise NcarRetrievalError("the responses do not match the plan")
    first = {field.name: next(p for p in pieces if p.request.field.name == field.name) for field in FIELDS}
    n_lat, n_lon = len(plan.latitudes), len(plan.longitudes)
    file_format = "NETCDF4" if compress else "NETCDF3_64BIT_OFFSET"
    value_type = np.float64 if float64 else np.float32
    with netCDF4.Dataset(path, "w", format=file_format) as out:
        out.setncatts({
            "Conventions": "CF-1.6",
            "title": "ERA5 pressure-level fields for the Lorenz Energy Cycle Toolkit",
            "source": (
                f"ERA5 reanalysis (ECMWF, Copernicus Climate Change Service) from {SOURCE['provider']} "
                f"dataset {SOURCE['dataset']} (DOI {SOURCE['doi']}), {SOURCE['service']}"
            ),
            "history": (
                f"{datetime.datetime.now(datetime.timezone.utc):%Y-%m-%dT%H:%M:%SZ}: "
                f"{len(pieces)} NCSS responses copied unchanged by LorenzCycleToolkit src/utils/ncar_gdex.py"
            ),
            "lectk_provenance": f"{os.path.basename(path).removesuffix('.part')}.provenance.json",
        })
        out.createDimension(names["time"], len(plan.times))
        out.createDimension(names["level"], len(PRESSURE_LEVELS_HPA))
        out.createDimension(names["latitude"], n_lat)
        out.createDimension(names["longitude"], n_lon)
        with netCDF4.Dataset(first["T"].path) as source:
            for coordinate, dtype, values in (
                ("time", "i4", _hours_since_origin(plan.times)),
                ("level", "f8", PRESSURE_LEVELS_HPA),
                ("latitude", "f8", plan.latitudes),
                ("longitude", "f8", plan.longitudes),
            ):
                variable = out.createVariable(names[coordinate], dtype, (names[coordinate],))
                _copy_attributes(source.variables[coordinate], variable)
                variable[:] = np.asarray(values)
        dimensions = tuple(names[d] for d in DIMENSIONS)
        options = {"zlib": True, "complevel": 1, "shuffle": True,
                   "chunksizes": (1, len(PRESSURE_LEVELS_HPA), n_lat, n_lon)} if compress else {}
        for field in FIELDS:
            with netCDF4.Dataset(first[field.name].path) as source:
                original = source.variables[field.name]
                fill = original.getncattr("_FillValue") if "_FillValue" in original.ncattrs() else None
                fill = None if fill is None else value_type(fill)
                variable = out.createVariable(names[field.name], value_type, dimensions, fill_value=fill, **options)
                _copy_attributes(original, variable, skip=("_FillValue", "coordinates", "missing_value"))
                if "missing_value" in original.ncattrs():
                    variable.setncattr("missing_value", value_type(original.getncattr("missing_value")))
                variable.setncattr("coordinates", " ".join(dimensions))
        for piece in pieces:
            request = piece.request
            with netCDF4.Dataset(piece.path) as source:
                source.set_auto_maskandscale(False)
                values = source.variables[request.field.name][:]
            target = out.variables[names[request.field.name]]
            target.set_auto_maskandscale(False)
            start = request.segment.offset
            target[request.time_offset:request.time_offset + len(request.hours), :, :,
                   start:start + len(request.segment.response)] = values.astype(value_type)


def check_output(path, plan, names):
    """
    Check the written file against the plan: coordinates and no missing values.

    Raises:
        NcarRetrievalError: If anything differs.
    """
    problems = []
    with xr.open_dataset(path) as dataset:
        expected = {
            names["level"]: np.asarray(PRESSURE_LEVELS_HPA),
            names["latitude"]: np.asarray(plan.latitudes),
            names["longitude"]: np.asarray(plan.longitudes),
        }
        for coordinate, values in expected.items():
            if not np.array_equal(np.asarray(dataset[coordinate].values, dtype=np.float64), values):
                problems.append(f"{coordinate} differs from the plan")
        times = np.asarray(dataset[names["time"]].values, dtype="datetime64[ns]")
        if not np.array_equal(times, np.asarray(plan.times, dtype="datetime64[ns]")):
            problems.append("times differ from the plan")
        for field in FIELDS:
            variable = dataset[names[field.name]]
            for index in range(variable.shape[0]):
                block = variable.isel({names["time"]: index}).values
                if not np.isfinite(block).all():
                    problems.append(f"{names[field.name]} has missing values at time index {index}")
                    break
    if problems:
        raise NcarRetrievalError(f"the written file {path} failed its check: " + "; ".join(problems))


# -------------------------------------------------------------- provenance
def _lectk_commit():
    directory = os.path.dirname(os.path.abspath(__file__))
    try:
        result = subprocess.run(["git", "-C", directory, "rev-parse", "HEAD"], capture_output=True,
                                text=True, timeout=10, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    commit = result.stdout.strip()
    return commit if result.returncode == 0 and len(commit) == 40 else None


def _software():
    from .run_status import _version

    return {
        "lectk_version": _version(),
        "lectk_commit": _lectk_commit(),
        "module": "src/utils/ncar_gdex.py",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "xarray": xr.__version__,
        "netCDF4": netCDF4.__version__,
        "netcdf_c": netCDF4.__netcdf4libversion__,
    }


def _iso(moment):
    return f"{moment:%Y-%m-%dT%H:%M:%SZ}"


def provenance_path(path):
    return f"{path}.provenance.json"


def check_existing(path, plan, names):
    """
    Decide whether an existing file can be reused for a plan.

    The file is reused only if its provenance record exists, names this
    module's source, matches the file's SHA-256, was written with the same
    variable names and covers every time, level, latitude and longitude of the
    plan.

    Returns:
        dict: The provenance record.

    Raises:
        NcarReuseError: Saying why the file cannot be reused.
    """
    record_path = provenance_path(path)
    advice = (f"Run without --ncar to use {path} as a local input file, or move it away "
              f"(with {os.path.basename(record_path)}) to retrieve the data again.")
    if not os.path.isfile(record_path):
        raise NcarReuseError(f"{path} exists but has no provenance record ({record_path}), "
                             f"so it cannot be verified as the data of this request. {advice}")
    try:
        with open(record_path, encoding="utf-8") as handle:
            record = json.load(handle)
        output = record["output"]
        recorded_sha = output["sha256"]
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise NcarReuseError(f"the provenance record {record_path} cannot be read ({error}). {advice}")
    if record.get("schema") != PROVENANCE_SCHEMA or record.get("source", {}).get("dataset") != SOURCE["dataset"]:
        raise NcarReuseError(f"{record_path} is not a GDEX provenance record of this toolkit. {advice}")
    if _sha256(path) != recorded_sha:
        raise NcarReuseError(f"{path} changed after it was written (SHA-256 differs from {record_path}). {advice}")
    if output.get("names") != names:
        raise NcarReuseError(
            f"{path} was written with the names {output.get('names')}, but the namelist of this run "
            f"gives {names}. {advice}"
        )
    missing = []
    if not {_iso(t) for t in plan.times} <= set(output.get("times_utc", [])):
        missing.append("times")
    if not set(plan.latitudes) <= set(output.get("latitudes", [])):
        missing.append("latitudes")
    if not {lon % 360.0 for lon in plan.longitudes} <= {lon % 360.0 for lon in output.get("longitudes", [])}:
        missing.append("longitudes")
    if list(output.get("levels_hPa", [])) != list(PRESSURE_LEVELS_HPA):
        missing.append("pressure levels")
    if missing:
        raise NcarReuseError(f"{path} does not cover the {', '.join(missing)} of this request. {advice}")
    return record


# ---------------------------------------------------------------- retrieve
def retrieve(plan, path, names, logger=None, workers=None, opener=None, timeout=None, attempts=None,
             backoff=None, compress=True, float64=False, keep_parts=False, max_total_bytes=None,
             request_info=None, namelist=None):
    """
    Retrieve every request of a plan and write the verified file ``path``.

    Args:
        plan (RetrievalPlan): Requests to make.
        path (str): Output NetCDF file; must not exist.
        names (dict): Output names, from :func:`output_names`.
        logger (logging.Logger): Logger (defaults to the toolkit logger).
        workers (int): Requests in flight at the same time (1 to 4;
            default DEFAULT_WORKERS).
        opener (callable): ``opener(urllib.request.Request, timeout=...)``
            returning a response; defaults to :func:`default_opener`.
        timeout (float): Seconds without data before an attempt fails
            (default DEFAULT_TIMEOUT).
        attempts (int): Attempts per request (default DEFAULT_ATTEMPTS).
        backoff (float): First wait between attempts, in seconds, then
            three times longer each time (default DEFAULT_BACKOFF).
        compress (bool): See :func:`write_output`.
        float64 (bool): See :func:`write_output`.
        keep_parts (bool): Keep the responses after the file is written.
        max_total_bytes (int): Largest retrieval accepted, in bytes
            (default DEFAULT_MAX_TOTAL_BYTES).
        request_info (dict): How the request was derived, for the record.
        namelist (str): Namelist the names came from, for the record.

    Returns:
        dict: The provenance record (also written next to the file).

    Raises:
        NcarRetrievalError: If the file could not be produced; the verified
            responses stay in ``<path>.parts/`` for the next attempt.
    """
    logger = logger or logging.getLogger("lorenzcycletoolkit")
    opener = opener or default_opener()
    workers = DEFAULT_WORKERS if workers is None else workers
    timeout = DEFAULT_TIMEOUT if timeout is None else timeout
    attempts = DEFAULT_ATTEMPTS if attempts is None else attempts
    backoff = DEFAULT_BACKOFF if backoff is None else backoff
    max_total_bytes = DEFAULT_MAX_TOTAL_BYTES if max_total_bytes is None else max_total_bytes
    if not 1 <= int(workers) <= MAX_WORKERS:
        raise ValueError(f"workers must be between 1 and {MAX_WORKERS}, not {workers}")
    if os.path.exists(path):
        raise NcarRetrievalError(f"{path} already exists; it is never overwritten")
    estimate = plan.estimated_bytes
    if estimate > max_total_bytes:
        raise NcarRetrievalError(
            f"the retrieval would transfer about {estimate / 1024 ** 3:.2f} GiB, more than the limit of "
            f"{max_total_bytes / 1024 ** 3:.2f} GiB. Shorten the period, use a larger --time-resolution, "
            f"or raise the limit (python -m src.utils.ncar_gdex ... --max-gib)."
        )
    directory = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(directory):
        raise NcarRetrievalError(f"the directory of {path} does not exist")
    parts_dir = f"{path}.parts"
    free = shutil.disk_usage(directory).free
    needed = estimate + plan.output_bytes
    if free < needed:
        raise NcarRetrievalError(
            f"not enough free disk space in {directory}: about {needed / 1024 ** 3:.2f} GiB needed "
            f"(responses and file), {free / 1024 ** 3:.2f} GiB free"
        )
    os.makedirs(parts_dir, exist_ok=True)
    days = sorted({r.day for r in plan.requests})
    logger.info(
        f"📥 NCAR GDEX: {len(plan.requests)} request(s) for {len(days)} day(s) "
        f"({days[0]} to {days[-1]}), {len(plan.times)} time(s), {len(PRESSURE_LEVELS_HPA)} levels, "
        f"area N{plan.north:g} W{plan.west:g} S{plan.south:g} E{plan.east:g} "
        f"({len(plan.latitudes)}x{len(plan.longitudes)} points), about {estimate / 1024 ** 2:.1f} MiB"
    )
    started_utc = datetime.datetime.now(datetime.timezone.utc)
    started = time.monotonic()
    cancel = threading.Event()
    results = {}
    # The last day goes first: a period that ends after the dates GDEX has
    # published fails on the first request instead of after the other days.
    ordered = sorted(plan.requests, key=lambda request: request.day != days[-1])
    executor = ThreadPoolExecutor(max_workers=int(workers), thread_name_prefix="ncar-gdex")
    futures = []
    try:
        futures = [executor.submit(fetch_piece, request, parts_dir, opener, timeout, attempts, backoff,
                                   cancel, logger) for request in ordered]
        for future in as_completed(futures):
            if future.exception() is not None:
                cancel.set()
                break
            piece = future.result()
            results[piece.request.key] = piece
            note = "kept from an earlier attempt" if piece.reused else f"{piece.elapsed:.1f} s"
            logger.info(f"  ⬇️ [{len(results)}/{len(plan.requests)}] {piece.request.key}: "
                        f"{piece.size / 1024 ** 2:.2f} MiB ({note})")
    except BaseException:
        cancel.set()
        executor.shutdown(wait=True, cancel_futures=True)
        logger.error(f"❌ NCAR GDEX retrieval interrupted; {len(results)} verified response(s) kept in {parts_dir}")
        raise
    executor.shutdown(wait=True, cancel_futures=True)
    finished = [f for f in futures if f.done() and not f.cancelled()]
    results.update({f.result().request.key: f.result() for f in finished if f.exception() is None})
    errors = [f.exception() for f in finished if f.exception() is not None]
    if errors or len(results) != len(plan.requests):
        real = [error for error in errors if not isinstance(error, _Cancelled)]
        error = real[0] if real else NcarRetrievalError("retrieval stopped before every response was received")
        logger.error(f"❌ NCAR GDEX retrieval failed: {error}")
        if results:
            logger.error(f"   {len(results)} verified response(s) kept in {parts_dir}; running again resumes from them.")
        else:
            with contextlib.suppress(OSError):
                os.rmdir(parts_dir)
        if isinstance(error, NcarRetrievalError):
            raise error
        raise NcarRetrievalError(f"{type(error).__name__}: {error}") from error
    transfer = time.monotonic() - started

    pieces = [results[r.key] for r in plan.requests]
    part_path = f"{path}.part"
    assembly_started = time.monotonic()
    try:
        write_output(plan, pieces, part_path, names, compress=compress, float64=float64)
        check_output(part_path, plan, names)
    except BaseException:
        if os.path.exists(part_path):
            os.remove(part_path)
        raise
    assembly = time.monotonic() - assembly_started
    sha256 = _sha256(part_path)
    downloaded = [p for p in pieces if not p.reused]
    record = {
        "schema": PROVENANCE_SCHEMA,
        "source": dict(SOURCE, base_url=plan.base_url),
        "request": dict(request_info or {}, **{
            "area": {"north": plan.north, "west": plan.west, "south": plan.south, "east": plan.east},
            "times_utc": [_iso(t) for t in plan.times],
            "levels_hPa": list(PRESSURE_LEVELS_HPA),
            "fields": [{"row": f.row, "source_name": f.name, "units": f.units, "output_name": names[f.name]}
                       for f in FIELDS],
        }),
        "responses": [
            {"key": p.request.key, "url": p.request.url, "bytes": p.size, "sha256": p.sha256,
             "attempts": p.attempts, "seconds": round(p.elapsed, 3), "reused": p.reused, "http": p.http}
            for p in pieces
        ],
        "output": {
            "file": os.path.basename(path),
            "sha256": sha256,
            "size_bytes": os.path.getsize(part_path),
            "format": "NETCDF4 (zlib level 1, shuffle)" if compress else "NETCDF3_64BIT_OFFSET",
            "values": ("float32 values of the responses stored as float64 (exact conversion)" if float64
                       else "float32 copied unchanged from the responses"),
            "names": names,
            "namelist": namelist,
            "namelist_sha256": _sha256(namelist) if namelist and os.path.isfile(namelist) else None,
            "times_utc": [_iso(t) for t in plan.times],
            "levels_hPa": list(PRESSURE_LEVELS_HPA),
            "latitudes": list(plan.latitudes),
            "longitudes": list(plan.longitudes),
        },
        "software": _software(),
        "timing": {
            "started_utc": _iso(started_utc),
            "finished_utc": _iso(datetime.datetime.now(datetime.timezone.utc)),
            "transfer_seconds": round(transfer, 3),
            "assembly_seconds": round(assembly, 3),
            "workers": int(workers),
        },
        "totals": {
            "requests": len(pieces),
            "downloaded_responses": len(downloaded),
            "reused_responses": len(pieces) - len(downloaded),
            "downloaded_bytes": sum(p.size for p in downloaded),
            "retries": sum(max(0, p.attempts - 1) for p in downloaded),
        },
    }
    _write_json(f"{provenance_path(path)}.tmp", record)
    os.replace(part_path, path)
    os.replace(f"{provenance_path(path)}.tmp", provenance_path(path))
    if not keep_parts:
        shutil.rmtree(parts_dir, ignore_errors=True)
    logger.info(
        f"✅ NCAR GDEX data ready: {path} ({record['output']['size_bytes'] / 1024 ** 2:.1f} MiB, "
        f"{len(downloaded)} response(s) in {transfer:.1f} s, written in {assembly:.1f} s); "
        f"provenance in {provenance_path(path)}"
    )
    return record


def get_ncar_data(args, track, app_logger, namelist="inputs/namelist", **options):
    """
    Retrieve the ERA5 data of a track from NCAR GDEX into ``args.infile``.

    Counterpart of :func:`src.utils.tools.get_cdsapi_data` with the same
    area, levels and hours. If ``args.infile`` exists it is reused only when
    its provenance record verifies it (see :func:`check_existing`).

    Args:
        args (argparse.Namespace): Command-line arguments (``infile``,
            ``trackfile``, ``time_resolution``).
        track (pd.DataFrame): Track with 'Lat' and 'Lon' columns and a
            datetime index.
        app_logger (logging.Logger): Logger for the application.
        namelist (str): Namelist of the run; the output names come from it.
        **options: Passed to :func:`retrieve`.

    Returns:
        str: Path to the verified NetCDF file.
    """
    namelist = resolve_input_file(namelist)
    try:
        names = output_names(namelist)
        plan = plan_retrieval(area_from_track(track), daily_hours(track.index, args.time_resolution))
    except ValueError as error:
        app_logger.error(f"❌ NCAR GDEX request cannot be made: {error}")
        raise
    if os.path.exists(args.infile):
        try:
            check_existing(args.infile, plan, names)
        except NcarReuseError as error:
            app_logger.error(f"❌ {error}")
            raise
        app_logger.info(f"✅ NCAR GDEX data already exists and matches its provenance record, "
                        f"skipping download: {args.infile}")
        return args.infile
    trackfile = getattr(args, "trackfile", None)
    request_info = {
        "rule": "area, levels and hours of the CDS API route (get_cdsapi_data) for the track",
        "time_resolution_hours": int(args.time_resolution),
        "track_extent": {"min_lat": float(track["Lat"].min()), "max_lat": float(track["Lat"].max()),
                         "min_lon": float(track["Lon"].min()), "max_lon": float(track["Lon"].max())},
        "trackfile": trackfile,
        "trackfile_sha256": _sha256(trackfile) if trackfile and os.path.isfile(trackfile) else None,
    }
    retrieve(plan, args.infile, names, logger=app_logger, request_info=request_info, namelist=namelist, **options)
    return args.infile


# --------------------------------------------------------------- command
def _read_track(trackfile, logger):
    from .validation import validate_track_file

    delimiter, _ = validate_track_file(trackfile, logger)
    return pd.read_csv(trackfile, delimiter=delimiter, index_col="time", parse_dates=["time"],
                       date_format="%Y-%m-%d-%H%M")


def main(argv=None):
    """Retrieve a file (or show the requests) without running the energetics."""
    parser = argparse.ArgumentParser(
        prog="python -m src.utils.ncar_gdex",
        description=("Retrieve the ERA5 data of a track from NCAR GDEX (d633000) with the area, levels and "
                     "hours that --cdsapi would request, and write a verified NetCDF file with a provenance "
                     "record."),
    )
    parser.add_argument("outfile", help="NetCDF file to write (must not exist).")
    parser.add_argument("--trackfile", default="inputs/track", help="Track file (default: inputs/track).")
    parser.add_argument("--time-resolution", type=int, default=3, help="Hours between times (default: 3).")
    parser.add_argument("--namelist", default="inputs/namelist",
                        help="Namelist whose names the file uses (default: inputs/namelist or its .default).")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help=f"Requests in flight at the same time, 1 to {MAX_WORKERS} (default: {DEFAULT_WORKERS}).")
    parser.add_argument("--max-gib", type=float, default=DEFAULT_MAX_TOTAL_BYTES / 1024 ** 3,
                        help="Largest transfer accepted, in GiB (default: %(default)g).")
    parser.add_argument("--no-compress", action="store_true", help="Write NetCDF-3 instead of compressed NetCDF-4.")
    parser.add_argument("--float64", action="store_true",
                        help="Store the float32 values as float64 (exact); the toolkit then computes in float64.")
    parser.add_argument("--keep-parts", action="store_true", help="Keep the responses after writing the file.")
    parser.add_argument("--dry-run", action="store_true", help="Show the requests and sizes; retrieve nothing.")
    parser.add_argument("-v", "--verbosity", action="store_true", help="Debug messages.")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbosity else logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    logger = logging.getLogger("lorenzcycletoolkit.ncar_gdex")
    track = _read_track(args.trackfile, logger)
    namelist = resolve_input_file(args.namelist)
    names = output_names(namelist)
    plan = plan_retrieval(area_from_track(track), daily_hours(track.index, args.time_resolution))
    if args.dry_run:
        print(f"area N{plan.north:g} W{plan.west:g} S{plan.south:g} E{plan.east:g}: "
              f"{len(plan.latitudes)} x {len(plan.longitudes)} points, {len(plan.times)} times, "
              f"{len(PRESSURE_LEVELS_HPA)} levels")
        print(f"{len(plan.requests)} requests, about {plan.estimated_bytes / 1024 ** 2:.1f} MiB to transfer, "
              f"{plan.output_bytes / 1024 ** 2:.1f} MiB of values; names from {namelist}: {names}")
        for request in plan.requests:
            print(request.url)
        return 0
    if os.path.exists(args.outfile):
        check_existing(args.outfile, plan, names)
        logger.info(f"✅ {args.outfile} already exists and matches its provenance record.")
        return 0
    request_info = {
        "rule": "area, levels and hours of the CDS API route (get_cdsapi_data) for the track",
        "time_resolution_hours": args.time_resolution,
        "trackfile": args.trackfile,
        "trackfile_sha256": _sha256(args.trackfile),
    }
    try:
        retrieve(plan, args.outfile, names, logger=logger, workers=args.workers, compress=not args.no_compress,
                 float64=args.float64, keep_parts=args.keep_parts, max_total_bytes=int(args.max_gib * 1024 ** 3),
                 request_info=request_info, namelist=namelist)
    except NcarRetrievalError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""
ERA5 pressure-level data from the DestinE Earth Data Hub (Zarr), optional.

Alternative to the CDS API route (:func:`src.utils.tools.get_cdsapi_data`) that
reads the copy of ERA5 pressure-level data that the Destination Earth (DestinE)
Earth Data Hub serves as a Zarr store over HTTPS. There is no request queue:
the toolkit reads the chunks it needs, and only those, each one once.

What the product is, and what that means for the cycle (read before using it):

- **its own pressure levels** (19 in the 2026-11 measurements, 1000 to 1 hPa),
  used as they are, without interpolation; the toolkit itself drops the levels
  above 10 hPa, as for every source. Terms of the Lorenz Energy Cycle can differ
  from routes with 37 levels.
- **values stored with reduced precision** (bit rounding of the float32
  mantissa). The v3 store keeps 13 bits (T within 0.008 K to 0.016 K, z within
  16 m2 s-2 of the original); the v2 store keeps 10 bits (T within 0.0625 K or
  0.125 K, z within 64 m2 s-2), and with it, in the 2b benchmark (phase 2), 18
  of 240 terms of the cycle differed between 1 % and 2.7 % from the CDS
  (boundary, generation and residual terms); energy contents and conversions
  did not change.
- **hourly**: any ``--time-resolution`` of whole hours is served.
- **requests:** every chunk is one HTTP request, and the account has a monthly
  request quota; the number of requests is estimated and limited before any
  chunk is read (``--max-requests``) and counted.
- **volume:** a chunk of the v3 store holds 1440 hours (60 days) of one level over
  a 60 x 60 tile (15 x 15 degrees); one of the v2 store, 120 hours over 256 x 256.
  What is transferred is many times what the file keeps (14 to 20 times with the
  v2 store in the 2b benchmark). The transfer is counted and limited by
  ``--max-gib``.

The request is the one of the CDS route for the same track and
``--time-resolution``: the five fields, the track extent with a 15-degree buffer
rounded outwards to whole degrees (:func:`src.utils.ncar_gdex.area_from_track`),
the same hours of each day (:func:`src.utils.ncar_gdex.daily_hours`) and the
longitudes as the CDS route returns them (``west`` to ``east``, increasing; a box
that crosses 180 degrees is refused). The file follows the input contract of the
CDS route: variables ``u v t w z`` (float32), dimensions ``valid_time``,
``pressure_level`` (hPa, 1000 first), ``latitude`` (north to south) and
``longitude``, CDS units, the values unchanged; ``inputs/namelist_ERA5-cdsapi`` is
used.

Two stores (``--destine-store``, decision D-043 of 2026-11-25, stage 4d):

- ``v3`` (default): the Zarr v3 store ``era5-pressure-levels-v0.zarr``, to which
  the Earth Data Hub moved ERA5 on 2026-07-06 and which it updates monthly;
- ``v2``: the Zarr v2 store of the 2b benchmark
  (``reanalysis-era5-pressure-levels-v0.zarr``), no longer updated since
  2026-07-06 (it ends in mid-2026) and possibly withdrawn without notice; kept
  while it is served, because its smaller chunks need fewer requests and bytes
  for short cases. Its values are rounded more coarsely (10 bits, not 13): the
  two stores do not hold the same values.

A time beyond the end of the store is refused before any chunk is read.

Optional dependencies: ``pip install "lorenzcycletoolkit[destine]"`` (zarr 3,
numcodecs, fsspec, aiohttp; the same as ``[arco]``). Credential: the personal
Earth Data Hub API key ("token"), read from the ``~/.netrc`` entry of
``data.earthdatahub.destine.eu`` (its ``password``, as the Earth Data Hub
documentation recommends; ``NETRC`` names another file) or from the toolkit's own
variable ``DESTINE_EDH_TOKEN``, and sent as HTTP basic authentication (user
``edh``). It is never logged or recorded. The monthly quota of the account is
500 000 requests (Earth Data Hub documentation, 2026-11); usage is shown on the
account page.

Protections:

- **volume:** before reading any chunk, the chunks that the selection touches are
  counted from the store metadata, whole chunks, uncompressed, plus a margin for
  the codec (an upper bound); if that passes the limit (``--max-gib``, default
  20 GiB) or the free disk space, the run stops saying the ``--max-gib`` needed.
  Afterwards every byte received is counted and no chunk is requested once the
  limit is reached;
- **requests:** the requests (metadata already read plus one per chunk touched)
  are shown before any chunk is read; above ``--max-requests`` (default 10000, 2 % of the monthly quota of 500 000)
  the run stops saying the value needed. Every attempt is counted against the
  limit, and the HTTP requests really made are recorded;
- **limited tries:** a chunk is requested at most 3 times (authorisation and
  quota errors are not retried), with at most 4 reads in parallel;
- **no orphan order:** there is no remote queue, so no request outlives the run.

An existing destination file is reused without a check, as with ``--cdsapi``.
The record ``<file>.provenance.json`` has the route, store, versions, commit,
box, times, levels, warnings, bytes, requests and SHA-256.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astronomia, Geofísica e Ciências Atmosféricas
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import collections.abc
import concurrent.futures
import datetime
import importlib
import importlib.metadata
import math
import netrc
import os
import platform
import shutil
import threading
import time

import numpy as np
import pandas as pd
import xarray as xr

from . import arco_zarr
from .ncar_gdex import DEFAULT_MAX_TOTAL_BYTES, GRID_STEP, _lectk_commit, area_from_track, longitude_segments
from .track_alignment import planned_request_times

STORE_HOST = "data.earthdatahub.destine.eu"
# Stores of the Earth Data Hub (zarr format, URL); v3 is the default (D-043, 2026-11-25)
STORES = {
    "v3": {"url": f"https://{STORE_HOST}/era5/era5-pressure-levels-v0.zarr", "zarr_format": 3, "updated": True},
    "v2": {"url": f"https://{STORE_HOST}/era5/reanalysis-era5-pressure-levels-v0.zarr", "zarr_format": 2,
           "updated": False},
}
DEFAULT_STORE = "v3"
STORE_URL = STORES[DEFAULT_STORE]["url"]
PROVENANCE_SCHEMA = "lectk-era5-destine/1"
SHORT_NAMES = arco_zarr.SHORT_NAMES
UNITS = arco_zarr.UNITS
OUTPUT_DIMENSIONS = arco_zarr.OUTPUT_DIMENSIONS
# Levels of the product in the 2026-11 measurements (LZV, 2b); the levels used are those of the store.
EXPECTED_LEVELS = (1000.0, 925.0, 850.0, 700.0, 600.0, 500.0, 400.0, 300.0, 250.0, 200.0, 150.0, 100.0, 70.0, 50.0,
                   30.0, 20.0, 10.0, 5.0, 1.0)
LEVELS_WARNING = ("product with its own pressure levels ({n}, {top:g} to {bottom:g} hPa; routes with 37 levels "
                  "differ), used without interpolation")
PRECISION_WARNING = ("values stored with reduced precision: in the comparison of Front 2b, 18 of 240 terms of the "
                     "cycle differed between 1 % and 2.7 % from the CDS (boundary, generation and residual terms); "
                     "energy contents and conversions did not change")  # the v2 store (10 bits of mantissa)
PRECISION_WARNING_V3 = ("values stored with reduced precision (13 bits of mantissa: T within 0.016 K, z within "
                        "16 m2 s-2); with the coarser v2 store (10 bits), 18 of 240 terms of the cycle differed "
                        "between 1 % and 2.7 % from the CDS in the comparison of Front 2b (boundary, generation and "
                        "residual terms), energy contents and conversions did not change")
PRECISION_WARNINGS = {"v3": PRECISION_WARNING_V3, "v2": PRECISION_WARNING}
QUOTA_NOTE = ("every chunk is one HTTP request and counts against the monthly request quota of the Earth Data Hub "
              "account")
ACCESS_HINT = 'pip install "lorenzcycletoolkit[destine]" (zarr>=3, numcodecs, fsspec, aiohttp)'
TOKEN_ENV = "DESTINE_EDH_TOKEN"
TOKEN_URL = "https://earthdatahub.destine.eu/quota-api-keys"
DEFAULT_MAX_REQUESTS = 10000
WORKERS = 4
ATTEMPTS = 3
BACKOFF = 2.0
HTTP_TIMEOUT = 300


class DestineError(RuntimeError):
    """Base class of the errors of the DestinE route."""


class DestineSetupError(DestineError):
    """A package or the token is missing, or the token is refused."""


class DestineRequestError(DestineError):
    """The request cannot be served by the product (resolution, area, times)."""


class DestineVolumeError(DestineError):
    """The transfer would pass the byte limit, the request limit or the free disk space."""


class DestineTransferError(DestineError):
    """A chunk could not be read."""


class DestineContentError(DestineError):
    """The store or the assembled file is not what was requested."""


# ARCO errors raised by the shared selection helpers, as DestinE errors (same message)
_FROM_ARCO = {
    arco_zarr.ArcoSetupError: DestineSetupError,
    arco_zarr.ArcoRequestError: DestineRequestError,
    arco_zarr.ArcoVolumeError: DestineVolumeError,
    arco_zarr.ArcoTransferError: DestineTransferError,
    arco_zarr.ArcoContentError: DestineContentError,
}


def _as_destine(error):
    return _FROM_ARCO.get(type(error), DestineError)(str(error))


# -------------------------------------------------------------------- plan
def check_time_resolution(time_step):
    """Refuse a time resolution that is not a whole number of hours (the product is hourly), before any read."""
    if int(time_step) != time_step or time_step < 1:
        raise DestineRequestError(f"--time-resolution {time_step} is not a whole number of hours. Nothing was read.")


def plan_request(track, time_step):
    """
    What the DestinE route reads for a track (no network); the levels come from the store.

    Args:
        track (pd.DataFrame): Track with 'Lat' and 'Lon' columns and a datetime index.
        time_step (int): ``--time-resolution`` in hours.

    Returns:
        dict: ``area`` (north, west, south, east), ``times`` (the times of the CDS route),
        ``latitudes`` (north to south) and ``longitudes`` (west to east, as returned to the toolkit).

    Raises:
        DestineRequestError: If the time step is not whole hours or the box cannot be served.
    """
    check_time_resolution(time_step)
    north, west, south, east = area_from_track(track)
    try:
        longitude_segments(west, east)  # refuses a box that crosses 180 degrees or leaves the grid
    except ValueError as error:
        raise DestineRequestError(str(error)) from None
    latitudes = north - GRID_STEP * np.arange(int(round((north - south) / GRID_STEP)) + 1)
    longitudes = west + GRID_STEP * np.arange(int(round((east - west) / GRID_STEP)) + 1)
    return {"area": (north, west, south, east), "times": planned_request_times(track.index, time_step),
            "latitudes": latitudes, "longitudes": longitudes}


def levels_warning(levels):
    levels = [float(x) for x in levels]
    return LEVELS_WARNING.format(n=len(levels), top=max(levels), bottom=min(levels))


# -------------------------------------------------------------------- credential and packages
def read_token(environ=None, netrc_file=None):
    """
    The Earth Data Hub token, from ``DESTINE_EDH_TOKEN`` or from the ``~/.netrc`` entry of the store host.

    The value is returned to the caller that builds the authentication; it is not logged here or
    anywhere else, and error messages never contain it.

    Raises:
        DestineSetupError: If no token is found.
    """
    environ = os.environ if environ is None else environ
    token = (environ.get(TOKEN_ENV) or "").strip()
    if not token:  # the .netrc entry has no login line in the Earth Data Hub example
        path = netrc_file or environ.get("NETRC") or os.path.join(os.path.expanduser("~"), ".netrc")
        try:
            entry = netrc.netrc(path).authenticators(STORE_HOST)
        except (OSError, netrc.NetrcParseError):
            entry = None
        token = (entry[2] or "").strip() if entry else ""
    if not token:
        raise DestineSetupError(
            f"Earth Data Hub token not found: --destine needs the personal token of a DestinE account. Set "
            f"{TOKEN_ENV}, or add to ~/.netrc the lines 'machine {STORE_HOST}' and '  password <your API key>' (file "
            f"readable only by you). "
            f"Account and token: {TOKEN_URL}"
        )
    return token


def load_packages():
    """Import the optional packages of the route, or explain how to install them."""
    modules = {}
    for name in ("zarr", "fsspec", "aiohttp"):
        try:
            modules[name] = importlib.import_module(name)
        except ImportError:
            raise DestineSetupError(f"--destine needs the optional package {name}: {ACCESS_HINT}") from None
    major = int(str(modules["zarr"].__version__).split(".")[0])
    if major < 3:
        raise DestineSetupError(f"--destine needs zarr 3 and found zarr {modules['zarr'].__version__}: {ACCESS_HINT}")
    return modules


def package_versions():
    out = {}
    for name in ("zarr", "numcodecs", "fsspec", "aiohttp"):
        try:
            out[name] = importlib.metadata.version(name)
        except Exception:  # noqa: BLE001
            out[name] = None
    return out


# -------------------------------------------------------------------- counted store
class Budget:
    """Bytes received and requests made, counted against two limits (thread safe)."""

    def __init__(self, limit, max_requests):
        self.limit, self.max_requests = int(limit), int(max_requests)
        self.used, self.reads, self.attempts, self.retries = 0, 0, 0, 0
        self._lock = threading.Lock()

    def start(self, key):
        """Count one request before it is made; refuse it if a limit is reached."""
        with self._lock:
            if self.used >= self.limit:
                raise DestineVolumeError(f"the byte limit ({self.limit}) was reached; {key} not requested.")
            if self.attempts >= self.max_requests:
                raise DestineVolumeError(
                    f"the request limit ({self.max_requests}) was reached; {key} not requested.")
            self.attempts += 1

    def add(self, size):
        with self._lock:
            self.used += size
            self.reads += 1
            return self.used > self.limit


class HttpCounter:
    """HTTP requests really made by the client (aiohttp trace), with their status codes."""

    def __init__(self):
        self.requests, self.status, self.exceptions = 0, {}, 0
        self._lock = threading.Lock()

    def end(self, status):
        with self._lock:
            self.requests += 1
            self.status[str(status)] = self.status.get(str(status), 0) + 1

    def failed(self):
        with self._lock:
            self.requests += 1
            self.exceptions += 1

    def summary(self):
        with self._lock:
            return {"requests": self.requests, "status": dict(self.status), "exceptions": self.exceptions}


_is_metadata = arco_zarr._is_metadata  # Zarr v2 (.zarray, .zattrs, .zgroup, .zmetadata) and v3 (zarr.json)


class CountingStore(collections.abc.MutableMapping):
    """
    Read-only mapping of zarr keys that counts every request and byte, retries a few times and stops at the limits.

    It wraps any mapping of keys to bytes (an fsspec mapper in production, a dict in tests), and zarr reads it
    through :func:`src.utils.arco_zarr.open_store_dataset`. Each key read is one request of the remote store. A
    missing chunk is recorded in ``missing`` and refused (DestineContentError): zarr would fill it with NaN. At
    most ``parallel`` keys are read at a time.
    """

    def __init__(self, inner, budget, attempts=ATTEMPTS, backoff=BACKOFF, sleep=time.sleep, parallel=WORKERS):
        self.inner, self.budget = inner, budget
        self.attempts, self.backoff, self.sleep = attempts, backoff, sleep
        self.in_flight = arco_zarr.InFlight(parallel)
        self.missing = []
        self.errors = []
        self.keys_read = collections.Counter()

    def __getitem__(self, key):
        last = None
        for attempt in range(1, self.attempts + 1):
            self.budget.start(key)
            try:
                with self.in_flight:
                    value = self.inner[key]
            except KeyError:
                if _is_metadata(key):
                    raise
                self.missing.append(key)
                raise DestineContentError(f"the chunk {key} is absent from the store") from None
            except Exception as exc:  # noqa: BLE001 - network errors of several libraries
                status = getattr(exc, "status", None)
                if status in (401, 403) or isinstance(exc, PermissionError):
                    raise DestineSetupError(
                        f"the Earth Data Hub refused the token (HTTP {status or 'permission'}): check the token "
                        f"({TOKEN_ENV} or ~/.netrc) and the account ({TOKEN_URL})."
                    ) from None
                if status == 429:
                    raise DestineVolumeError(
                        "the Earth Data Hub answered HTTP 429 (too many requests): the request quota of the "
                        "account may be exhausted; nothing more was requested."
                    ) from None
                last = exc
                self.errors.append(f"{key}: {type(exc).__name__}")
                with self.budget._lock:
                    self.budget.retries += 1
                if attempt < self.attempts:
                    self.sleep(self.backoff * attempt)
                continue
            self.keys_read[key] += 1
            if self.budget.add(len(value)):
                raise DestineVolumeError(f"the byte limit ({self.budget.limit}) was passed reading {key}; stopped.")
            return value
        raise DestineTransferError(f"{key} not read after {self.attempts} attempts: {type(last).__name__}")

    def __contains__(self, key):
        # zarr tests ``key in store`` before reading a chunk; on a remote store that is one more request per
        # chunk. A chunk is declared present here and read next: if it is absent, reading it raises
        # DestineContentError (the route never accepts a missing chunk, which zarr would fill with NaN).
        if _is_metadata(key):
            return key in self.inner
        return True

    def __iter__(self):
        return iter(self.inner)

    def __len__(self):
        return len(self.inner)

    def __setitem__(self, key, value):
        raise DestineError("the DestinE store is read-only")

    def __delitem__(self, key):
        raise DestineError("the DestinE store is read-only")


def open_remote_mapping(token, counter, url=STORE_URL, packages=None):
    """fsspec mapper of the store with the token as basic authentication and a trace that counts HTTP requests."""
    packages = load_packages() if packages is None else packages
    aiohttp = packages["aiohttp"]

    async def on_end(session, context, params):
        counter.end(params.response.status)

    async def on_exception(session, context, params):
        counter.failed()

    trace = aiohttp.TraceConfig()
    trace.on_request_end.append(on_end)
    trace.on_request_exception.append(on_exception)
    return packages["fsspec"].get_mapper(url, client_kwargs={
        "auth": aiohttp.BasicAuth("edh", token), "trace_configs": [trace],
        "timeout": aiohttp.ClientTimeout(total=HTTP_TIMEOUT)})


# -------------------------------------------------------------------- selection and estimates
def store_levels(ds):
    """Pressure levels of the store in hPa, 1000 first (any unit: Pa values are converted)."""
    name = arco_zarr._dimension(ds, "level")
    values = np.asarray(ds[name].values, dtype=float)
    hpa = values / 100.0 if values.max() > 2000.0 else values
    return np.sort(hpa)[::-1], name, values.max() > 2000.0


def select(ds, plan, store_name=DEFAULT_STORE):
    """Positions in the store of the planned times, the store levels, latitudes and longitudes."""
    levels, name, in_pa = store_levels(ds)
    plan["levels"] = levels
    plan_for_store = dict(plan, levels=levels * 100.0 if in_pa else levels)
    try:
        return arco_zarr.select(ds, plan_for_store)
    except arco_zarr.ArcoError as error:
        if "time(s) are not in the store" in str(error):
            if STORES[store_name]["updated"]:
                raise DestineRequestError(
                    f"{error}. The Earth Data Hub store is updated monthly; use another data source for the "
                    f"latest dates. No chunk was read.") from None
            raise DestineRequestError(
                f"{error}. The Earth Data Hub Zarr v2 store (--destine-store v2) is no longer updated (since "
                f"2026-07-06); use --destine-store v3 or another data source for later dates. No chunk was "
                f"read.") from None
        raise _as_destine(error) from None


def chunks_touched(ds, selection):
    """Number of chunks the selection touches, by store variable (each one read once)."""
    positions = selection["positions"]
    out = {}
    for name in selection["fields"].values():
        variable = ds[name]
        chunks = variable.encoding.get("chunks")
        if not chunks:
            raise DestineContentError(f"{name}: the store does not say its chunk shape")
        if variable.encoding.get("shards"):
            raise DestineContentError(f"{name}: the store is sharded, which this route does not support")
        count = 1
        for dim, chunk in zip(variable.dims, chunks):
            count *= len(np.unique(positions[dim] // chunk)) if dim in positions else math.ceil(variable.sizes[dim] / chunk)
        out[name] = count
    return out


def codec_margin(raw_chunk_bytes):
    """Most bytes a compressor can add to a chunk that does not compress (Blosc and zstd: well under this)."""
    return raw_chunk_bytes // 64 + 4096


def estimate(ds, selection, budget):
    """
    Upper bound of the bytes, and number of requests, before any chunk is read.

    Bytes: what was already received (metadata and coordinates) plus every chunk touched, whole, uncompressed,
    with a codec margin. Requests: those already made plus one per chunk touched (retries not included).
    """
    try:
        raw_total, raw_by_variable = arco_zarr.estimate_bytes(ds, selection)
    except arco_zarr.ArcoError as error:
        raise _as_destine(error) from None
    chunks = chunks_touched(ds, selection)
    by_variable = {}
    for name, raw in raw_by_variable.items():
        per_chunk = raw // chunks[name] if chunks[name] else 0
        by_variable[name] = raw + chunks[name] * codec_margin(per_chunk)
    return {
        "bytes": budget.used + sum(by_variable.values()),
        "bytes_by_variable": by_variable,
        "bytes_raw_chunks": raw_total,
        "chunks": sum(chunks.values()),
        "chunks_by_variable": chunks,
        "requests": budget.attempts + sum(chunks.values()),
    }


def output_bytes(plan):
    """Size of the float32 values of the output file."""
    return 4 * len(SHORT_NAMES) * len(plan["times"]) * len(plan["levels"]) * len(plan["latitudes"]) * len(plan["longitudes"])


def check_limits(est, plan, limit, max_requests, directory):
    """Refuse before reading chunks if the estimate passes a limit or the output does not fit the disk."""
    if est["requests"] > max_requests:
        raise DestineVolumeError(
            f"this request needs about {est['requests']} HTTP requests to the Earth Data Hub ({est['chunks']} chunks "
            f"plus the metadata already read), above the limit of {max_requests}; no chunk was read. Use "
            f"--max-requests {est['requests']} or more to accept it (the account has a monthly request quota)."
        )
    if est["bytes"] > limit:
        needed = math.ceil(est["bytes"] / 1024 ** 3 * 10) / 10
        raise DestineVolumeError(
            f"the DestinE chunks of this request may weigh up to {est['bytes'] / 1024 ** 3:.2f} GiB ({est['bytes']} "
            f"bytes, whole chunks, uncompressed upper bound), above the limit of {limit / 1024 ** 3:.2f} GiB; no "
            f"chunk was read. Use --max-gib {needed:g} or more to accept it."
        )
    needed_disk = 2 * output_bytes(plan)
    free = shutil.disk_usage(directory).free
    if needed_disk > free:
        raise DestineVolumeError(
            f"the output file may need up to {needed_disk / 1024 ** 3:.2f} GiB of disk and {directory} has "
            f"{free / 1024 ** 3:.2f} GiB free; no chunk was read."
        )


# -------------------------------------------------------------------- reading and writing
def read_level(ds, selection, store_name, level_position):
    """
    Values of one field at one level for every planned time, latitude and longitude (float32).

    One orthogonal selection on the store array: zarr reads each chunk it touches once, whatever the
    number of times in the store (no dask graph is built).
    """
    names, positions = selection["names"], selection["positions"]
    variable = ds[store_name]
    index = {names["time"]: positions[names["time"]], names["level"]: level_position,
             names["latitude"]: positions[names["latitude"]], names["longitude"]: positions[names["longitude"]]}
    values = variable.isel({d: index[d] for d in variable.dims}).transpose(
        names["time"], names["latitude"], names["longitude"]).values
    return np.asarray(values, dtype="float32")


def write_output(ds, plan, selection, path, notes, workers=WORKERS, store_name=DEFAULT_STORE):
    """
    Read the selection and write the NetCDF file through a temporary name.

    The (field, level) blocks are read ``workers`` at a time; the file is written only from the calling
    thread (the netCDF-C library is not thread safe).
    """
    import netCDF4

    tmp = f"{path}.destine-tmp.nc"
    names = selection["names"]
    level_positions = selection["positions"][names["level"]]
    times = np.asarray(pd.DatetimeIndex(plan["times"]).values).astype("datetime64[s]").astype("int64")
    attrs = {
        "Conventions": "CF-1.7",
        "title": "ERA5 pressure-level fields for the Lorenz Energy Cycle Toolkit",
        "source": ("ERA5 reanalysis (ECMWF, Copernicus Climate Change Service), DestinE Earth Data Hub Zarr "
                   f"{store_name} store ({len(plan['levels'])} pressure levels, hourly, values stored with reduced "
                   f"precision)"),
        "history": (f"{datetime.datetime.now(datetime.timezone.utc):%Y-%m-%dT%H:%M:%SZ}: selected from the store "
                    f"by LorenzCycleToolkit src/utils/destine_edh.py, values unchanged"),
        "lectk_provenance_note": " | ".join(notes),
    }
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as out:
        out.setncatts(attrs)
        sizes = (len(times), len(plan["levels"]), len(plan["latitudes"]), len(plan["longitudes"]))
        for dim, size in zip(OUTPUT_DIMENSIONS, sizes):
            out.createDimension(dim, size)
        coordinates = {
            "valid_time": out.createVariable("valid_time", "i8", ("valid_time",)),
            "pressure_level": out.createVariable("pressure_level", "f8", ("pressure_level",)),
            "latitude": out.createVariable("latitude", "f8", ("latitude",)),
            "longitude": out.createVariable("longitude", "f8", ("longitude",)),
        }
        coordinates["valid_time"].setncatts({"units": "seconds since 1970-01-01", "calendar": "proleptic_gregorian"})
        coordinates["valid_time"][:] = times
        for name, units, values in (("pressure_level", "hPa", plan["levels"]),
                                    ("latitude", "degrees_north", plan["latitudes"]),
                                    ("longitude", "degrees_east", plan["longitudes"])):
            coordinates[name].setncattr("units", units)
            coordinates[name][:] = np.asarray(values, dtype="f8")
        variables = {}
        for short in SHORT_NAMES:
            variables[short] = out.createVariable(short, "f4", OUTPUT_DIMENSIONS, fill_value=np.float32(np.nan))
            variables[short].setncattr("units", UNITS[short])
        jobs = [(short, store_name, j) for short, store_name in selection["fields"].items()
                for j in range(len(level_positions))]
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(read_level, ds, selection, store_name, int(level_positions[j])): (short, j)
                       for short, store_name, j in jobs}
            try:
                for future in concurrent.futures.as_completed(futures):
                    short, j = futures[future]
                    variables[short][:, j, :, :] = future.result()
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
    os.replace(tmp, path)


def check_output(path, plan):
    """Check the written file against the plan: coordinates, and no missing value."""
    problems = []
    with xr.open_dataset(path) as ds:
        for name, expected in (("pressure_level", plan["levels"]), ("latitude", plan["latitudes"]),
                               ("longitude", plan["longitudes"])):
            if not np.array_equal(np.asarray(ds[name].values, dtype=float), np.asarray(expected, dtype=float)):
                problems.append(f"{name} differs from the plan")
        if not np.array_equal(ds["valid_time"].values.astype("datetime64[ns]"),
                              pd.DatetimeIndex(plan["times"]).values.astype("datetime64[ns]")):
            problems.append("times differ from the plan")
        for short in SHORT_NAMES:
            if short not in ds:
                problems.append(f"{short} missing")
                continue
            for i in range(ds[short].sizes["valid_time"]):
                if not np.isfinite(ds[short].isel(valid_time=i).values).all():
                    problems.append(f"{short} has missing values at time index {i}")
                    break
    if problems:
        raise DestineContentError(f"the written file {path} failed its check: " + "; ".join(problems))


def provenance_path(path):
    return f"{path}.provenance.json"


# -------------------------------------------------------------------- entry point
def get_destine_data(args, track, app_logger, mapping=None, sleep=time.sleep):
    """
    Retrieve the ERA5 file of a track from the DestinE Earth Data Hub Zarr store into ``args.infile``.

    Args:
        args (argparse.Namespace): ``infile``, ``time_resolution``, ``max_gib`` (None: 20 GiB),
            ``max_requests`` (None: 10000) and ``destine_store`` (None: v3).
        track (pd.DataFrame): Track with 'Lat', 'Lon' and a datetime index.
        app_logger (logging.Logger): Logger.
        mapping: Mapping of zarr keys to bytes standing for the remote store (tests pass a local one);
            default: the Earth Data Hub store, opened with the token.
        sleep (callable): wait between tries (tests pass a stub).

    Returns:
        str: ``args.infile``.
    """
    path = args.infile
    max_gib = getattr(args, "max_gib", None)
    limit = DEFAULT_MAX_TOTAL_BYTES if max_gib is None else int(max_gib * 1024 ** 3)
    max_requests = getattr(args, "max_requests", None)
    max_requests = DEFAULT_MAX_REQUESTS if max_requests is None else int(max_requests)
    store_name = getattr(args, "destine_store", None) or DEFAULT_STORE
    if store_name not in STORES:
        raise DestineRequestError(f"unknown Earth Data Hub store {store_name!r}: use one of {sorted(STORES)}")
    store_url, precision_warning = STORES[store_name]["url"], PRECISION_WARNINGS[store_name]
    plan = plan_request(track, args.time_resolution)
    north, west, south, east = plan["area"]
    app_logger.warning(f"⚠️ DestinE: {levels_warning(EXPECTED_LEVELS)} (levels confirmed from the store below)")
    app_logger.warning(f"⚠️ DestinE: {precision_warning}")
    app_logger.warning(f"⚠️ DestinE: {QUOTA_NOTE}")
    app_logger.info(f"🌐 DestinE: store {store_name}, area N{north} W{west} S{south} E{east}, {len(plan['times'])} "
                    f"times (limits {limit / 1024 ** 3:.2f} GiB and {max_requests} requests)")
    directory = os.path.dirname(os.path.abspath(path))
    started = datetime.datetime.now(datetime.timezone.utc)
    remote = mapping is None
    counter = HttpCounter()
    if remote:
        packages = load_packages()
        mapping = open_remote_mapping(read_token(), counter, url=store_url, packages=packages)
    else:
        load_packages()
    budget = Budget(limit, max_requests)
    store = CountingStore(mapping, budget, sleep=sleep)

    t0 = time.monotonic()
    try:
        ds = arco_zarr.open_store_dataset(store, STORES[store_name]["zarr_format"])  # no dask: see read_level
        selection = select(ds, plan, store_name)
        warning_levels = levels_warning(plan["levels"])
        app_logger.warning(f"⚠️ DestinE: {warning_levels}")
        est = estimate(ds, selection, budget)
        app_logger.info(f"📦 DestinE: {est['chunks']} chunks, about {est['requests']} HTTP requests, up to "
                        f"{est['bytes'] / 1024 ** 3:.2f} GiB ({output_bytes(plan) / 1024 ** 2:.1f} MiB of useful "
                        f"values)")
        check_limits(est, plan, limit, max_requests, directory)
        metadata = {"bytes": budget.used, "requests": budget.attempts}
        write_output(ds, plan, selection, path, [warning_levels, precision_warning], workers=WORKERS,
                     store_name=store_name)
    except DestineError:
        _discard(path)
        raise
    except Exception as exc:  # noqa: BLE001
        _discard(path)
        cause = exc
        while cause is not None and not isinstance(cause, (DestineError, arco_zarr.ArcoError)):
            cause = cause.__cause__ or cause.__context__
        if isinstance(cause, DestineError):
            raise cause from None
        if cause is not None:
            raise _as_destine(cause) from None
        raise DestineTransferError(f"reading the DestinE store failed: {type(exc).__name__}: {exc}") from None
    transfer_seconds = time.monotonic() - t0
    if store.missing:
        _discard(path, final=True)
        raise DestineContentError(f"{len(store.missing)} chunk(s) are absent from the store, first {store.missing[0]}")
    try:
        check_output(path, plan)
    except DestineError:
        _discard(path, final=True)
        raise
    fields = set(selection["fields"].values())
    repeated = sum(1 for key, n in store.keys_read.items()
                   if n > 1 and key.split("/")[0] in fields and not _is_metadata(key))  # data chunks only
    provenance = {
        "schema": PROVENANCE_SCHEMA,
        "route": "destine",
        "store": store_url if remote else "local store (test)",
        "store_version": store_name,
        "zarr_format": STORES[store_name]["zarr_format"],
        "warnings": [warning_levels, precision_warning, QUOTA_NOTE],
        "packages": package_versions(),
        "lectk_commit": _lectk_commit(),
        "python": platform.python_version(),
        "xarray": xr.__version__,
        "area_nwse": list(plan["area"]),
        "times": [f"{t:%Y-%m-%dT%H:%MZ}" for t in plan["times"]],
        "time_resolution_h": args.time_resolution,
        "levels_hpa": [float(x) for x in plan["levels"]],
        "limit_bytes": limit,
        "limit_requests": max_requests,
        "estimate_bytes": est["bytes"],
        "estimate_bytes_by_variable": est["bytes_by_variable"],
        "estimate_raw_chunk_bytes": est["bytes_raw_chunks"],
        "estimate_chunks": est["chunks"],
        "estimate_requests": est["requests"],
        "bytes_received": budget.used,
        "bytes_metadata": metadata["bytes"],
        "bytes_counted": "bodies of the chunk and metadata responses (HTTP headers and TLS not counted)",
        "requests_store": budget.attempts,
        "requests_metadata": metadata["requests"],
        "requests_http": counter.summary() if remote else None,
        "reads": budget.reads,
        "retries": budget.retries,
        "chunks_read_more_than_once": repeated,
        "workers": WORKERS,
        "max_parallel_reads": store.in_flight.peak,
        "transfer_seconds": round(transfer_seconds, 3),
        "started_utc": started.isoformat(),
        "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "output": {"path": os.path.abspath(path), "bytes": os.path.getsize(path), "sha256": arco_zarr._sha256(path)},
    }
    arco_zarr._write_json(provenance_path(path), provenance)
    app_logger.info(f"🎉 DestinE file ready: {path} ({provenance['output']['bytes'] / 1024 ** 2:.2f} MB, "
                    f"{budget.used / 1024 ** 2:.1f} MB received in {budget.attempts} requests, sha256 "
                    f"{provenance['output']['sha256'][:12]}…); provenance in {provenance_path(path)}")
    return path


def _discard(path, final=False):
    """Remove the temporary file; with ``final``, also the output file (never leave a file that fails its check)."""
    for name in (f"{path}.destine-tmp.nc", path if final else None):
        if name and os.path.exists(name):
            os.remove(name)

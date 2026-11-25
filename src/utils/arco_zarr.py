"""
ERA5 pressure-level data from the CDS ARCO Zarr store (time-chunked), optional.

Alternative to the CDS API route (:func:`src.utils.tools.get_cdsapi_data`) that
reads the *Analysis Ready, Cloud Optimised* copy of ERA5 pressure-level data that
the Copernicus Climate Data Store (CDS) serves as a Zarr store over HTTPS. There
is no request queue: the toolkit reads the chunks it needs, and only those.

What the product is, and what that means for the cycle (read before using it):

- **13 pressure levels**, 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150,
  100 and 50 hPa; the toolkit does not interpolate or complete them. Terms of the
  Lorenz Energy Cycle can differ from routes with 37 levels (in the 2026-11 LZV
  case, Ca was 5 % lower on average).
- **Every 6 hours** (00, 06, 12, 18 UTC): a ``--time-resolution`` that is not a
  multiple of 6 is refused before anything is read.
- **Values** equal to the CDS NetCDF for ``u``, ``v`` and ``w`` and within
  rounding (``t`` up to 2.4e-4 K, ``z`` up to 0.0625 m2 s-2) for ``t`` and ``z``
  (2b benchmark, phase 2); the cycle changes by at most 0.04 %.
- **Volume:** the store is chunked by time, a chunk holding whole latitude
  bands; what is transferred is several times what the file keeps (7.5 times
  the useful volume, 2.7 GB, in the LZV case). The transfer is counted
  and limited by ``--max-gib``.

The request is the one of the CDS route for the same track and
``--time-resolution``: the five fields, the track extent with a 15-degree buffer
rounded outwards to whole degrees (:func:`src.utils.ncar_gdex.area_from_track`),
the same hours of each day (:func:`src.utils.ncar_gdex.daily_hours`) and the
longitudes as the CDS route returns them (``west`` to ``east``, increasing, so
that a box that crosses Greenwich is one contiguous arc; a box that crosses 180
degrees is refused). The file follows the input contract of the CDS route:
variables ``u v t w z`` (float32), dimensions ``valid_time``, ``pressure_level``,
``latitude`` (north to south) and ``longitude``, CDS units, the values unchanged;
``inputs/namelist_ERA5-cdsapi`` is used.

Optional dependencies: ``pip install "lorenzcycletoolkit[arco]"`` (zarr 3,
numcodecs, fsspec, aiohttp; the store is in Zarr format 2, which zarr 3 reads).
Credential: the CDS API key, sent as a bearer token, read from ``CDSAPI_KEY``
or from the file named by ``CDSAPI_RC`` (default ``~/.cdsapirc``, line
``key: ...``). It is never logged or recorded. The store
URL is that of the ECMWF documentation (ERA5 ARCO user guide) and of the 2b
benchmark; the CDS terms of use of ERA5 apply.

Protections:

- **volume:** before reading any chunk, the chunks that the selection touches are
  counted from the store metadata, whole chunks, uncompressed (an upper bound);
  if that passes the limit (``--max-gib``, default 20 GiB, as for ``--ncar``) or
  the free disk space, the run stops saying the ``--max-gib`` needed. Afterwards
  every byte received is counted (headers and TLS excluded) and no chunk is
  requested once the limit is reached;
- **limited tries:** a chunk is requested at most 3 times (authorisation errors
  are not retried), with at most 4 reads in parallel;
- **no orphan order:** there is no remote queue, so no request outlives the run.

An existing destination file is reused without a check, as with ``--cdsapi``.
The record ``<file>.provenance.json`` has the route, store, versions, commit,
box, times, levels, bytes and SHA-256.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astronomia, Geofísica e Ciências Atmosféricas
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import asyncio
import collections.abc
import datetime
import hashlib
import importlib
import importlib.metadata
import json
import math
import os
import platform
import shutil
import threading
import time

import numpy as np
import pandas as pd
import xarray as xr

from .ncar_gdex import DEFAULT_MAX_TOTAL_BYTES, GRID_STEP, _lectk_commit, area_from_track, longitude_segments
from .track_alignment import planned_request_times

STORE_URL = ("https://arco.datastores.ecmwf.int/cadl-arco-time-048/arco/"
             "reanalysis_era5_pressure_levels/pl/timeChunked.zarr")
STORE_ZARR_FORMAT = 2  # the ARCO store is Zarr v2 (``.zmetadata``); zarr-python 3 reads it
PROVENANCE_SCHEMA = "lectk-era5-arco/1"
# Output short name -> names the store may give the variable (same fields as src/utils/tools.py:
# get_cdsapi_data). The real store (observed 2026-11-24) uses the short names; the user guide lists the long ones.
STORE_VARIABLES = {
    "u": ("u", "u_component_of_wind"),
    "v": ("v", "v_component_of_wind"),
    "t": ("t", "temperature"),
    "w": ("w", "vertical_velocity"),
    "z": ("z", "geopotential"),
}
SHORT_NAMES = tuple(STORE_VARIABLES)
# CDS NetCDF units (the contract of inputs/namelist_ERA5-cdsapi)
UNITS = {"u": "m s**-1", "v": "m s**-1", "t": "K", "w": "Pa s**-1", "z": "m**2 s**-2"}
PRESSURE_LEVELS = (1000.0, 925.0, 850.0, 700.0, 600.0, 500.0, 400.0, 300.0, 250.0, 200.0, 150.0, 100.0, 50.0)
TIME_STEP_HOURS = 6
# Names the store may use, in order of preference (the real store calls the pressure dimension "pressureLevel", the user guide "pressureLevels")
DIMENSION_CANDIDATES = {
    "time": ("time", "valid_time"),
    "level": ("pressureLevel", "pressureLevels", "pressure_level", "level", "isobaricInhPa"),
    "latitude": ("latitude", "lat"),
    "longitude": ("longitude", "lon"),
}
OUTPUT_DIMENSIONS = ("valid_time", "pressure_level", "latitude", "longitude")
LEVELS_WARNING = ("product with 13 levels (1000 to 50 hPa) every 6 h; terms of the cycle may differ from routes "
                  "with 37 levels (in the LZV case, Ca -5 % on average)")
ACCESS_HINT = 'pip install "lorenzcycletoolkit[arco]" (zarr>=3, numcodecs, fsspec, aiohttp)'
CDS_PROFILE_URL = "https://cds.climate.copernicus.eu/"
WORKERS = 4
ATTEMPTS = 3
BACKOFF = 2.0
HTTP_TIMEOUT = 300
GRID_TOLERANCE = 1e-6
CHUNK_READ = 1024 ** 2


class ArcoError(RuntimeError):
    """Base class of the errors of the ARCO route."""


class ArcoSetupError(ArcoError):
    """A package or the CDS credential is missing, or the credential is refused."""


class ArcoRequestError(ArcoError):
    """The request cannot be served by the product (resolution, area, times)."""


class ArcoVolumeError(ArcoError):
    """The transfer would pass the byte limit or the free disk space."""


class ArcoTransferError(ArcoError):
    """A chunk could not be read."""


class ArcoContentError(ArcoError):
    """The store or the assembled file is not what was requested."""


# -------------------------------------------------------------------- plan
def plan_request(track, time_step):
    """
    What the ARCO route reads for a track (no network).

    Args:
        track (pd.DataFrame): Track with 'Lat' and 'Lon' columns and a datetime index.
        time_step (int): ``--time-resolution`` in hours.

    Returns:
        dict: ``area`` (north, west, south, east), ``times`` (pandas DatetimeIndex, the times of the
        CDS route), ``latitudes`` (north to south), ``longitudes`` (west to east, as returned to the
        toolkit) and ``levels``.

    Raises:
        ArcoRequestError: If ``time_step`` is not a multiple of 6 h, or the box cannot be served.
    """
    check_time_resolution(time_step)
    north, west, south, east = area_from_track(track)
    try:
        longitude_segments(west, east)  # refuses a box that crosses 180 degrees or leaves the grid
    except ValueError as error:
        raise ArcoRequestError(str(error)) from None
    latitudes = north - GRID_STEP * np.arange(int(round((north - south) / GRID_STEP)) + 1)
    longitudes = west + GRID_STEP * np.arange(int(round((east - west) / GRID_STEP)) + 1)
    return {"area": (north, west, south, east), "times": planned_request_times(track.index, time_step),
            "latitudes": latitudes, "longitudes": longitudes, "levels": np.asarray(PRESSURE_LEVELS)}


def check_time_resolution(time_step):
    """Refuse a time resolution that the 6-hourly product cannot serve, before any read."""
    if int(time_step) != time_step or time_step < TIME_STEP_HOURS or time_step % TIME_STEP_HOURS:
        raise ArcoRequestError(
            f"--arco serves the 6-hourly ERA5 ARCO product (00, 06, 12 and 18 UTC): --time-resolution "
            f"{time_step} is not a multiple of {TIME_STEP_HOURS} h. Use 6 (or 12, 18, 24), or another "
            f"data source (--cdsapi, --wekeo or --ncar) for 3 h. Nothing was read."
        )


# -------------------------------------------------------------------- credential and packages
def read_token(environ=None, home=None):
    """
    The CDS API key, from ``CDSAPI_KEY`` or from the file named by ``CDSAPI_RC`` (default ``~/.cdsapirc``).

    The value is returned to the caller that builds the request header; it is not logged here
    or anywhere else, and error messages never contain it.

    Raises:
        ArcoSetupError: If no key is found.
    """
    environ = os.environ if environ is None else environ
    key = (environ.get("CDSAPI_KEY") or "").strip()
    if not key:
        rc = environ.get("CDSAPI_RC") or os.path.join(home or os.path.expanduser("~"), ".cdsapirc")
        try:
            with open(rc, encoding="utf-8") as f:
                for line in f:
                    name, _, value = line.partition(":")
                    if name.strip() == "key" and value.strip():
                        key = value.strip()
                        break
        except OSError:
            key = ""
    if not key:
        raise ArcoSetupError(
            "CDS credential not found: --arco reads the same CDS API key as --cdsapi. Create ~/.cdsapirc with the "
            "line 'key: <your key>' (or point CDSAPI_RC to such a file, or set CDSAPI_KEY). A free account: "
            f"{CDS_PROFILE_URL}"
        )
    return key


def load_packages():
    """Import the optional packages of the route, or explain how to install them."""
    modules = {}
    for name in ("zarr", "fsspec", "aiohttp"):
        try:
            modules[name] = importlib.import_module(name)
        except ImportError:
            raise ArcoSetupError(f"--arco needs the optional package {name}: {ACCESS_HINT}") from None
    major = int(str(modules["zarr"].__version__).split(".")[0])
    if major < 3:
        raise ArcoSetupError(f"--arco needs zarr 3 and found zarr {modules['zarr'].__version__}: {ACCESS_HINT}")
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
    """Bytes received from the store, counted against a limit (thread safe)."""

    def __init__(self, limit):
        self.limit, self.used, self.reads, self.retries = int(limit), 0, 0, 0
        self._lock = threading.Lock()

    def room(self):
        with self._lock:
            return self.limit - self.used

    def add(self, size):
        with self._lock:
            self.used += size
            self.reads += 1
            return self.used > self.limit


def _is_metadata(key):
    """Metadata keys of Zarr v2 (``.zarray``, ``.zattrs``, ``.zgroup``, ``.zmetadata``) and v3 (``zarr.json``)."""
    name = os.path.basename(key)
    return name.startswith(".z") or name == "zarr.json"


class InFlight:
    """At most ``limit`` reads of the remote store at a time, whoever asks (zarr 3 reads chunks concurrently)."""

    def __init__(self, limit=WORKERS):
        self.limit, self.now, self.peak = int(limit), 0, 0
        self._slots = threading.BoundedSemaphore(self.limit)
        self._lock = threading.Lock()

    def __enter__(self):
        self._slots.acquire()
        with self._lock:
            self.now += 1
            self.peak = max(self.peak, self.now)

    def __exit__(self, *exc):
        with self._lock:
            self.now -= 1
        self._slots.release()


class CountingStore(collections.abc.MutableMapping):
    """
    Read-only mapping of zarr keys that counts every byte it returns, retries a few times and stops at the limit.

    It wraps any mapping of keys to bytes (an fsspec mapper in production, a dict in tests), and zarr reads it
    through :func:`open_store_dataset`. A missing chunk is recorded in ``missing`` and refused
    (ArcoContentError): zarr would fill it with NaN. At most ``parallel`` keys are read at a time.
    """

    def __init__(self, inner, budget, attempts=ATTEMPTS, backoff=BACKOFF, sleep=time.sleep, parallel=WORKERS):
        self.inner, self.budget = inner, budget
        self.attempts, self.backoff, self.sleep = attempts, backoff, sleep
        self.in_flight = InFlight(parallel)
        self.missing = []
        self.errors = []

    def __getitem__(self, key):
        last = None
        for attempt in range(1, self.attempts + 1):
            if self.budget.room() <= 0:
                raise ArcoVolumeError(
                    f"the byte limit ({self.budget.limit}) was reached; {key} not requested.")
            try:
                with self.in_flight:
                    value = self.inner[key]
            except KeyError:
                if _is_metadata(key):
                    raise
                self.missing.append(key)
                raise ArcoContentError(f"the chunk {key} is absent from the store") from None
            except Exception as exc:  # noqa: BLE001 - network errors of several libraries
                status = getattr(exc, "status", None)
                if status in (401, 403) or isinstance(exc, PermissionError):
                    raise ArcoSetupError(
                        f"the CDS ARCO store refused the credential (HTTP {status or 'permission'}): check the CDS "
                        f"API key and that your CDS account has accepted the ERA5 licence ({CDS_PROFILE_URL})."
                    ) from None
                last = exc
                self.errors.append(f"{key}: {type(exc).__name__}")
                self.budget.retries += 1
                if attempt < self.attempts:
                    self.sleep(self.backoff * attempt)
                continue
            if self.budget.add(len(value)):
                raise ArcoVolumeError(f"the byte limit ({self.budget.limit}) was passed reading {key}; stopped.")
            return value
        raise ArcoTransferError(f"{key} not read after {self.attempts} attempts: {type(last).__name__}")

    def __contains__(self, key):
        # zarr tests ``key in store`` before reading a chunk; on a remote store that is one more request per
        # chunk. A chunk is declared present here and read next: if it is absent, reading it raises
        # ArcoContentError (the route never accepts a missing chunk, which zarr would fill with NaN).
        if _is_metadata(key):
            return key in self.inner
        return True

    def __iter__(self):
        return iter(self.inner)

    def __len__(self):
        return len(self.inner)

    def __setitem__(self, key, value):
        raise ArcoError("the ARCO store is read-only")

    def __delitem__(self, key):
        raise ArcoError("the ARCO store is read-only")


_STORE_CLASS = None


def _store_class():
    """
    zarr 3 store class that reads a counted mapping (built once, on first use: zarr is an optional package).

    zarr 3 asks for keys asynchronously; each key is read from the mapping in a worker thread, so the counting,
    the limits, the retries and the cap on parallel reads of :class:`CountingStore` apply unchanged.
    """
    global _STORE_CLASS
    if _STORE_CLASS is None:
        from zarr.abc.store import Store

        class MappingStore(Store):
            supports_writes = False
            supports_deletes = False
            supports_partial_writes = False
            supports_listing = False

            def __init__(self, mapping):
                super().__init__(read_only=True)
                self.mapping = mapping

            def __eq__(self, other):
                return self is other

            async def get(self, key, prototype, byte_range=None):
                if byte_range is not None:  # only sharded stores read parts of a key; they are refused
                    raise ArcoContentError(f"partial read of {key} requested: sharded stores are not supported")
                value = await asyncio.to_thread(self.mapping.get, key)
                return None if value is None else prototype.buffer.from_bytes(value)

            async def get_partial_values(self, prototype, key_ranges):
                return [await self.get(key, prototype, byte_range) for key, byte_range in key_ranges]

            async def exists(self, key):
                return await asyncio.to_thread(self.mapping.__contains__, key)

            async def set(self, key, value):
                raise ArcoError("the store is read-only")

            async def delete(self, key):
                raise ArcoError("the store is read-only")

            async def list(self):
                return
                yield

            async def list_prefix(self, prefix):
                return
                yield

            async def list_dir(self, prefix):
                return
                yield

        _STORE_CLASS = MappingStore
    return _STORE_CLASS


def open_store_dataset(store, zarr_format):
    """
    The store as a lazy xarray Dataset, from its consolidated metadata, with no dask array (``chunks=None``).

    Args:
        store (CountingStore): counted mapping of the store keys.
        zarr_format (int): 2 or 3 (given, so that no key of the other format is requested).
    """
    import warnings

    with warnings.catch_warnings():  # zarr 3 notes that consolidated metadata is not yet in the v3 specification
        warnings.filterwarnings("ignore", message=".*[Cc]onsolidated metadata.*")
        return xr.open_zarr(_store_class()(store), consolidated=True, chunks=None, zarr_format=zarr_format)


def open_remote_mapping(token, url=STORE_URL, packages=None):
    """fsspec mapper of the store with the CDS key as bearer token (the header is not kept anywhere else)."""
    packages = load_packages() if packages is None else packages
    aiohttp = packages["aiohttp"]
    return packages["fsspec"].get_mapper(
        url, headers={"Authorization": f"Bearer {token}"},
        client_kwargs={"timeout": aiohttp.ClientTimeout(total=HTTP_TIMEOUT)})


# -------------------------------------------------------------------- selection
def _dimension(ds, key):
    for name in DIMENSION_CANDIDATES[key]:
        if name in ds.dims or name in ds.coords:
            return name
    raise ArcoContentError(f"the store has no {key} dimension (looked for {DIMENSION_CANDIDATES[key]}; "
                           f"found {list(ds.dims)})")


def _positions(index, wanted, what):
    """Positions of ``wanted`` values in a store coordinate (exact on the grid); error if any is absent."""
    found = pd.Index(index).get_indexer(wanted, method="nearest", tolerance=GRID_TOLERANCE)
    absent = np.flatnonzero(found < 0)
    if len(absent):
        raise ArcoContentError(f"{len(absent)} {what} value(s) of the request are not in the store, first "
                               f"{np.asarray(wanted)[absent[0]]} (store range {np.min(index)} to {np.max(index)})")
    return found


def store_variable_names(ds):
    """Name that the store gives each of the five fields (short name -> store name)."""
    names, missing = {}, []
    for short, candidates in STORE_VARIABLES.items():
        found = next((c for c in candidates if c in ds.data_vars), None)
        if found is None:
            missing.append(short)
        else:
            names[short] = found
    if missing:
        raise ArcoContentError(f"the store lacks the field(s) {missing} (found {list(ds.data_vars)})")
    return names


def select(ds, plan):
    """
    Positions in the store of the planned times, levels, latitudes and longitudes.

    Returns:
        dict: ``names`` (dimension names of the store) and ``positions`` (integer arrays by dimension name).
    """
    names = {key: _dimension(ds, key) for key in DIMENSION_CANDIDATES}
    fields = store_variable_names(ds)
    time_index = pd.DatetimeIndex(ds[names["time"]].values)
    wanted_times = pd.DatetimeIndex(plan["times"])
    present = wanted_times.isin(time_index)
    if not present.all():
        raise ArcoContentError(f"{int((~present).sum())} requested time(s) are not in the store, first "
                               f"{wanted_times[~present][0]} (store {time_index.min()} to {time_index.max()})")
    lon = np.asarray(ds[names["longitude"]].values, dtype=float)
    west_to_store = plan["longitudes"] % 360.0 if lon.max() > 180.0 else ((plan["longitudes"] + 180.0) % 360.0) - 180.0
    positions = {
        names["time"]: time_index.get_indexer(wanted_times),
        names["level"]: _positions(ds[names["level"]].values, plan["levels"], "pressure level"),
        names["latitude"]: _positions(ds[names["latitude"]].values, plan["latitudes"], "latitude"),
        names["longitude"]: _positions(lon, west_to_store, "longitude"),
    }
    return {"names": names, "fields": fields, "positions": positions}


def estimate_bytes(ds, selection):
    """
    Most bytes the selection can transfer: the chunks it touches, whole and uncompressed.

    Returns:
        tuple: (total bytes, dict of bytes by store variable).
    """
    positions = selection["positions"]
    by_variable = {}
    for name in selection["fields"].values():
        variable = ds[name]
        chunks = variable.encoding.get("chunks")
        if not chunks:
            raise ArcoContentError(f"{name}: the store does not say its chunk shape")
        if variable.encoding.get("shards"):
            raise ArcoContentError(f"{name}: the store is sharded, which this route does not support")
        item = np.dtype(variable.encoding.get("dtype", variable.dtype)).itemsize
        count, volume = 1, 1
        for dim, chunk in zip(variable.dims, chunks):
            count *= len(np.unique(positions[dim] // chunk)) if dim in positions else math.ceil(variable.sizes[dim] / chunk)
            volume *= chunk
        by_variable[name] = count * volume * item
    return sum(by_variable.values()), by_variable


def output_bytes(plan):
    """Size of the float32 values of the output file."""
    return 4 * len(SHORT_NAMES) * len(plan["times"]) * len(plan["levels"]) * len(plan["latitudes"]) * len(plan["longitudes"])


def check_volume(estimate, plan, limit, directory):
    """Refuse before reading chunks if the estimate passes ``limit`` or the output does not fit the disk."""
    if estimate > limit:
        needed = math.ceil(estimate / 1024 ** 3 * 10) / 10
        raise ArcoVolumeError(
            f"the ARCO chunks of this request may weigh up to {estimate / 1024 ** 3:.2f} GiB ({estimate} bytes, "
            f"whole chunks, uncompressed upper bound), above the limit of {limit / 1024 ** 3:.2f} GiB; no chunk "
            f"was read. Use --max-gib {needed:g} or more to accept it."
        )
    needed_disk = 2 * output_bytes(plan)
    free = shutil.disk_usage(directory).free
    if needed_disk > free:
        raise ArcoVolumeError(
            f"the output file may need up to {needed_disk / 1024 ** 3:.2f} GiB of disk and {directory} has "
            f"{free / 1024 ** 3:.2f} GiB free; no chunk was read."
        )


# -------------------------------------------------------------------- output
def build_output(ds, plan, selection):
    """Lazy dataset with the CDS contract: names, units, dimension order and coordinates of the plan."""
    names, positions = selection["names"], selection["positions"]
    out = {}
    for short, store_name in selection["fields"].items():
        variable = ds[store_name].isel({d: positions[d] for d in names.values() if d in ds[store_name].dims})
        variable = variable.transpose(names["time"], names["level"], names["latitude"], names["longitude"])
        # Dask only now, after the selection: one task per time step and level, each reading its chunks once.
        # (Chunking the whole store first gave ~1.65 million blocks per field and a 50-minute, 20 GB run.)
        variable = variable.chunk({names["time"]: 1, names["level"]: 1, names["latitude"]: -1,
                                   names["longitude"]: -1}).astype("float32")
        out[short] = (OUTPUT_DIMENSIONS, variable.data, {"units": UNITS[short]})
    result = xr.Dataset(out, coords={
        "valid_time": ("valid_time", pd.DatetimeIndex(plan["times"])),
        "pressure_level": ("pressure_level", plan["levels"], {"units": "hPa"}),
        "latitude": ("latitude", plan["latitudes"], {"units": "degrees_north"}),
        "longitude": ("longitude", plan["longitudes"], {"units": "degrees_east"}),
    }, attrs={
        "Conventions": "CF-1.7",
        "title": "ERA5 pressure-level fields for the Lorenz Energy Cycle Toolkit",
        "source": ("ERA5 reanalysis (ECMWF, Copernicus Climate Change Service), CDS ARCO Zarr store "
                   "(time-chunked, 13 pressure levels, 6-hourly)"),
        "history": (f"{datetime.datetime.now(datetime.timezone.utc):%Y-%m-%dT%H:%M:%SZ}: selected from the store "
                    f"by LorenzCycleToolkit src/utils/arco_zarr.py, values unchanged"),
        "lectk_provenance_note": LEVELS_WARNING,
    })
    return result


def write_output(result, path, workers=WORKERS):
    """
    Read the selection and write the NetCDF file through a temporary name.

    Chunks are read ``workers`` at a time, one time step at a time, and the file is written only from the
    calling thread (the netCDF-C library is not thread safe: a failure while other threads still wrote
    into an open file crashed the interpreter in testing).
    """
    import dask
    import netCDF4

    tmp = f"{path}.arco-tmp.nc"
    times = np.asarray(result["valid_time"].values).astype("datetime64[s]").astype("int64")
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as out:
        out.setncatts(result.attrs)
        for dim, size in zip(OUTPUT_DIMENSIONS, (len(times),) + tuple(result.sizes[d] for d in OUTPUT_DIMENSIONS[1:])):
            out.createDimension(dim, size)
        coordinates = {
            "valid_time": out.createVariable("valid_time", "i8", ("valid_time",)),
            "pressure_level": out.createVariable("pressure_level", "f8", ("pressure_level",)),
            "latitude": out.createVariable("latitude", "f8", ("latitude",)),
            "longitude": out.createVariable("longitude", "f8", ("longitude",)),
        }
        coordinates["valid_time"].setncatts({"units": "seconds since 1970-01-01", "calendar": "proleptic_gregorian"})
        coordinates["valid_time"][:] = times
        for name, units in (("pressure_level", "hPa"), ("latitude", "degrees_north"), ("longitude", "degrees_east")):
            coordinates[name].setncattr("units", units)
            coordinates[name][:] = np.asarray(result[name].values, dtype="f8")
        variables = {}
        for short in SHORT_NAMES:
            variables[short] = out.createVariable(short, "f4", OUTPUT_DIMENSIONS, fill_value=np.float32(np.nan))
            variables[short].setncattr("units", UNITS[short])
        with dask.config.set(scheduler="threads", num_workers=workers):
            for i in range(len(times)):
                step = result.isel(valid_time=i).compute()
                for short in SHORT_NAMES:
                    variables[short][i] = step[short].values
    os.replace(tmp, path)


def check_output(path, plan):
    """Check the written file against the plan: coordinates, and no missing value (chunk by time)."""
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
        raise ArcoContentError(f"the written file {path} failed its check: " + "; ".join(problems))


def provenance_path(path):
    return f"{path}.provenance.json"


def _write_json(path, content):
    tmp = f"{path}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(content, f, indent=1, default=str)
    os.replace(tmp, path)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(CHUNK_READ), b""):
            digest.update(block)
    return digest.hexdigest()


# -------------------------------------------------------------------- entry point
def get_arco_data(args, track, app_logger, mapping=None, sleep=time.sleep):
    """
    Retrieve the ERA5 file of a track from the CDS ARCO Zarr store into ``args.infile``.

    Args:
        args (argparse.Namespace): ``infile``, ``time_resolution`` and ``max_gib`` (None: 20 GiB).
        track (pd.DataFrame): Track with 'Lat', 'Lon' and a datetime index.
        app_logger (logging.Logger): Logger.
        mapping: Mapping of zarr keys to bytes standing for the remote store (tests pass a local one);
            default: the CDS store, opened with the credential.
        sleep (callable): wait between tries (tests pass a stub).

    Returns:
        str: ``args.infile``.
    """
    path = args.infile
    max_gib = getattr(args, "max_gib", None)
    limit = DEFAULT_MAX_TOTAL_BYTES if max_gib is None else int(max_gib * 1024 ** 3)
    plan = plan_request(track, args.time_resolution)
    north, west, south, east = plan["area"]
    app_logger.warning(f"⚠️ ARCO: {LEVELS_WARNING}")
    app_logger.info(f"🌐 ARCO: area N{north} W{west} S{south} E{east}, {len(plan['times'])} times, "
                    f"{len(plan['levels'])} levels (limit {limit / 1024 ** 3:.2f} GiB)")
    directory = os.path.dirname(os.path.abspath(path))
    started = datetime.datetime.now(datetime.timezone.utc)
    remote = mapping is None
    if remote:
        packages = load_packages()
        mapping = open_remote_mapping(read_token(), packages=packages)
    else:
        load_packages()
    budget = Budget(limit)
    store = CountingStore(mapping, budget, sleep=sleep)

    t0 = time.monotonic()
    try:
        ds = open_store_dataset(store, STORE_ZARR_FORMAT)  # no dask yet: see build_output
        selection = select(ds, plan)
        estimate, by_variable = estimate_bytes(ds, selection)
        app_logger.info(f"📦 ARCO: the chunks touched weigh up to {estimate / 1024 ** 3:.2f} GiB "
                        f"({output_bytes(plan) / 1024 ** 2:.1f} MiB of useful values)")
        check_volume(estimate, plan, limit, directory)
        metadata_bytes = budget.used
        result = build_output(ds, plan, selection)
        write_output(result, path, workers=WORKERS)
    except ArcoError:
        _discard(path)
        raise
    except Exception as exc:  # noqa: BLE001
        _discard(path)
        cause = exc
        while cause is not None and not isinstance(cause, ArcoError):
            cause = cause.__cause__ or cause.__context__
        if cause is not None:
            raise cause from None
        raise ArcoTransferError(f"reading the ARCO store failed: {type(exc).__name__}: {exc}") from None
    transfer_seconds = time.monotonic() - t0
    if store.missing:
        _discard(path, final=True)
        raise ArcoContentError(f"{len(store.missing)} chunk(s) are absent from the store, first {store.missing[0]}")
    try:
        check_output(path, plan)
    except ArcoError:
        _discard(path, final=True)
        raise
    provenance = {
        "schema": PROVENANCE_SCHEMA,
        "route": "arco",
        "store": STORE_URL if remote else "local store (test)",
        "warning": LEVELS_WARNING,
        "packages": package_versions(),
        "lectk_commit": _lectk_commit(),
        "python": platform.python_version(),
        "xarray": xr.__version__,
        "area_nwse": list(plan["area"]),
        "times": [f"{t:%Y-%m-%dT%H:%MZ}" for t in plan["times"]],
        "time_resolution_h": args.time_resolution,
        "levels_hpa": [float(x) for x in plan["levels"]],
        "limit_bytes": limit,
        "estimate_bytes": estimate,
        "estimate_by_variable": by_variable,
        "bytes_received": budget.used,
        "bytes_metadata": metadata_bytes,
        "bytes_counted": "bodies of the chunk and metadata responses (HTTP headers and TLS not counted)",
        "reads": budget.reads,
        "retries": budget.retries,
        "workers": WORKERS,
        "max_parallel_reads": store.in_flight.peak,
        "transfer_seconds": round(transfer_seconds, 3),
        "started_utc": started.isoformat(),
        "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "output": {"path": os.path.abspath(path), "bytes": os.path.getsize(path), "sha256": _sha256(path)},
    }
    _write_json(provenance_path(path), provenance)
    app_logger.info(f"🎉 ARCO file ready: {path} ({provenance['output']['bytes'] / 1024 ** 2:.2f} MB, "
                    f"{budget.used / 1024 ** 2:.1f} MB received, sha256 {provenance['output']['sha256'][:12]}…); "
                    f"provenance in {provenance_path(path)}")
    return path


def _discard(path, final=False):
    """Remove the temporary file; with ``final``, also the output file (never leave a file that fails its check)."""
    for name in (f"{path}.arco-tmp.nc", path if final else None):
        if name and os.path.exists(name):
            os.remove(name)

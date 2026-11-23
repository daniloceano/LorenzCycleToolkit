"""
ERA5 pressure-level data from WEkEO (Harmonised Data Access, HDA), optional.

Alternative to the CDS API route (:func:`src.utils.tools.get_cdsapi_data`)
that asks WEkEO for the same ERA5 product (dataset
``EO:ECMWF:DAT:REANALYSIS_ERA5_PRESSURE_LEVELS``). The request is the one the
CDS route makes for the same track and ``--time-resolution``: the five fields,
the 37 pressure levels, the track extent with a 15-degree buffer rounded
outwards to whole degrees (:func:`src.utils.ncar_gdex.area_from_track`), with
the longitudes passed as the CDS route passes them, and the same hours of each
day (:func:`src.utils.ncar_gdex.daily_hours`). WEkEO is asked once per month of
the track (as in the 2b benchmark); each monthly request covers its days times
its hours, and the extra times at the ends are dropped when the file is
assembled, so the file has exactly the times of the CDS route.

The file follows the input contract of the CDS route: the variable, coordinate
and dimension names, units, dimension order and longitude convention of the
NetCDF that WEkEO delivers (the same as the CDS; ``inputs/namelist_ERA5-cdsapi``
is used), the values unchanged, concatenated along ``valid_time``.

Optional dependency: the ``hda`` client (``pip install
"lorenzcycletoolkit[wekeo]"``). Credential: a WEkEO (or EUMETSAT/CMEMS)
account in ``~/.hdarc`` (lines ``user:`` and ``password:``) or in the
environment variables ``HDA_USER`` and ``HDA_PASSWORD``. The toolkit never
reads, logs or records the credential; the ``hda`` client reads it.

Protections:

- **volume:** before any request, the largest possible size of the responses
  (float32, uncompressed) is compared with the limit (``--max-gib``, default
  20 GiB, as for ``--ncar``) and the free disk space; every byte received,
  including failed attempts, is counted, and a response is not read past the
  limit;
- **limited tries:** the ``hda`` client retries each HTTP call a few times
  (not its default 500); a download is tried at most twice, on the same order;
- **no orphan order:** before each order is submitted, a ledger
  ``<file>.wekeo-orders.json`` records it; its id is added once WEkEO answers.
  A later run with the same destination resumes the orders it recorded instead
  of submitting new ones. An order whose submission was recorded but whose id
  never arrived has an unknown remote state: the route stops until the user
  removes the ledger (same policy as D-033 for the CDS). The ledger and the
  monthly files are deleted once the file is complete; the provenance record
  ``<file>.provenance.json`` keeps the orders, bytes and SHA-256 sums.

An existing destination file is reused without a check, as with ``--cdsapi``.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astronomia, Geofísica e Ciências Atmosféricas
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import datetime
import hashlib
import importlib
import importlib.metadata
import json
import logging
import math
import os
import platform
import shutil
import time
import zipfile

import pandas as pd
import xarray as xr

from .ncar_gdex import DEFAULT_MAX_TOTAL_BYTES, GRID_STEP, _lectk_commit, area_from_track, daily_hours

WEKEO_DATASET = "EO:ECMWF:DAT:REANALYSIS_ERA5_PRESSURE_LEVELS"
PROVENANCE_SCHEMA = "lectk-era5-wekeo/1"
# Same fields and levels as src/utils/tools.py:get_cdsapi_data
VARIABLES = ("u_component_of_wind", "v_component_of_wind", "temperature", "vertical_velocity", "geopotential")
SHORT_NAMES = ("u", "v", "t", "w", "z")
PRESSURE_LEVELS = (
    "1", "2", "3", "5", "7", "10", "20", "30", "50", "70",
    "100", "125", "150", "175", "200", "225", "250", "300",
    "350", "400", "450", "500", "550", "600", "650", "700",
    "750", "775", "800", "825", "850", "875", "900", "925",
    "950", "975", "1000",
)
TIME_DIMENSION = "valid_time"
REGISTER_URL = "https://www.wekeo.eu/register"
INSTALL_HINT = 'pip install "lorenzcycletoolkit[wekeo]" (or pip install hda)'
# hda.Client options: a few HTTP retries (the client's default is 500), no progress bar, one download thread
CLIENT_OPTIONS = {"timeout": 120, "retry_max": 3, "sleep_max": 60, "progress": False, "max_workers": 1}
ORDER_MAX_SECONDS = 3600      # longest wait for one order to be ready
POLL_FIRST, POLL_MAX = 5.0, 60.0
DOWNLOAD_ATTEMPTS = 2         # tries of the download of one ready order (never a new order)
HEADER_ALLOWANCE = 1024 ** 2  # per response, on top of the float32 values
CHUNK = 1024 ** 2
NETCDF_SIGNATURES = (b"CDF\x01", b"CDF\x02", b"CDF\x05", b"\x89HDF")
ZIP_SIGNATURE = b"PK\x03\x04"


class WekeoError(RuntimeError):
    """Base class of the errors of the WEkEO route."""


class WekeoSetupError(WekeoError):
    """The hda package or the WEkEO credential is missing."""


class WekeoVolumeError(WekeoError):
    """The request would pass the byte limit or the free disk space."""


class WekeoOrderError(WekeoError):
    """An order failed, is still pending, or has an unknown remote state."""


class WekeoContentError(WekeoError):
    """A response or the assembled file is not what was requested."""


# -------------------------------------------------------------------- plan
def _grid_count(first, last):
    return int(round(abs(last - first) / GRID_STEP)) + 1


def plan_requests(track, time_step):
    """
    Monthly WEkEO requests for a track (no network).

    Args:
        track (pd.DataFrame): Track with 'Lat' and 'Lon' columns and a datetime index.
        time_step (int): ``--time-resolution`` in hours.

    Returns:
        dict: ``area`` (north, west, south, east), ``times`` (the times of the
        CDS route, as pandas Timestamps), ``requests`` (list of dicts with
        ``month``, ``query``, ``times`` and ``max_bytes``) and ``max_bytes``.
    """
    north, west, south, east = area_from_track(track)
    days_hours = daily_hours(track.index, time_step)
    months = {}
    for day, hours in days_hours:
        months.setdefault((day.year, day.month), []).append((day, hours))
    points = _grid_count(north, south) * _grid_count(west, east)
    requests = []
    for (year, month), entries in sorted(months.items()):
        days = sorted({f"{day.day:02d}" for day, _ in entries})
        hours = sorted({hour for _, hs in entries for hour in hs})
        query = {
            "dataset_id": WEKEO_DATASET,
            "product_type": ["reanalysis"],
            "variable": list(VARIABLES),
            "year": [f"{year}"],
            "month": [f"{month:02d}"],
            "day": days,
            "time": [f"{hour:02d}:00" for hour in hours],
            "pressure_level": list(PRESSURE_LEVELS),
            "data_format": "netcdf",
            "download_format": "unarchived",
            "bbox": [west, south, east, north],
        }
        n_times = len(days) * len(hours)
        requests.append({
            "month": f"{year}-{month:02d}",
            "query": query,
            "times": [pd.Timestamp(day) + pd.Timedelta(hours=h) for day, hs in entries for h in hs],
            "max_bytes": 4 * len(VARIABLES) * n_times * len(PRESSURE_LEVELS) * points + HEADER_ALLOWANCE,
        })
    times = [t for request in requests for t in request["times"]]
    return {"area": (north, west, south, east), "times": times, "requests": requests,
            "max_bytes": sum(r["max_bytes"] for r in requests)}


def check_volume(plan, limit, directory):
    """
    Refuse the plan before any request if it may pass ``limit`` bytes or the free disk space.

    The disk must hold the monthly responses and the assembled file (twice the estimate).
    """
    estimate = plan["max_bytes"]
    if estimate > limit:
        needed = math.ceil(estimate / 1024 ** 3 * 10) / 10
        raise WekeoVolumeError(
            f"the WEkEO request may take up to {estimate / 1024 ** 3:.2f} GiB ({estimate} bytes, float32 "
            f"uncompressed upper bound), above the limit of {limit / 1024 ** 3:.2f} GiB; nothing was requested. "
            f"Use --max-gib {needed:g} or more to accept it."
        )
    free = shutil.disk_usage(directory).free
    if 2 * estimate > free:
        raise WekeoVolumeError(
            f"the WEkEO request may need up to {2 * estimate / 1024 ** 3:.2f} GiB of disk (responses and the "
            f"assembled file) and {directory} has {free / 1024 ** 3:.2f} GiB free; nothing was requested."
        )


# -------------------------------------------------------------------- client
def load_hda():
    """Import the optional hda package, or explain how to install it."""
    try:
        return importlib.import_module("hda")
    except ImportError:
        raise WekeoSetupError(f"--wekeo needs the optional hda package: {INSTALL_HINT}") from None


def make_client(hda=None):
    """
    hda.Client with limited retries. The credential is read by hda only
    (~/.hdarc, or HDA_USER and HDA_PASSWORD); its values are never seen here.
    """
    hda = load_hda() if hda is None else hda
    # The hda client logs its access token at DEBUG level: keep its logger above that
    logging.getLogger("hda").setLevel(logging.WARNING)
    try:
        config = hda.Configuration()
    except Exception as exc:  # hda.api.ConfigurationError; its text is not shown (may not be ours to show)
        if type(exc).__name__ != "ConfigurationError":
            raise
        raise WekeoSetupError(
            "WEkEO credential not found: create ~/.hdarc with the two lines 'user:<username>' and "
            "'password:<password>' (or set HDA_USER and HDA_PASSWORD). An account is free: "
            f"{REGISTER_URL}"
        ) from None
    return hda.Client(config=config, **CLIENT_OPTIONS)


def hda_version():
    try:
        return importlib.metadata.version("hda")
    except Exception:  # noqa: BLE001
        return None


# -------------------------------------------------------------------- ledger
def ledger_path(path):
    return f"{path}.wekeo-orders.json"


def parts_dir(path):
    return f"{path}.wekeo-parts"


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
        for block in iter(lambda: f.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def read_ledger(path, plan):
    """
    Orders recorded by an earlier run for the same destination, checked against ``plan``.

    Raises:
        WekeoOrderError: if the ledger is for another request, or records a
            submission without an order id (unknown remote state).
    """
    ledger = ledger_path(path)
    if not os.path.exists(ledger):
        return {}
    with open(ledger, encoding="utf-8") as f:
        entries = json.load(f).get("orders", {})
    queries = {r["month"]: r["query"] for r in plan["requests"]}
    for month, entry in entries.items():
        if queries.get(month) != entry.get("query"):
            raise WekeoOrderError(
                f"{ledger} records orders of another request ({month}); nothing was requested. Use another "
                f"destination file, or remove the ledger after checking those orders on WEkEO."
            )
        if entry.get("status") == "submitting":
            raise WekeoOrderError(
                f"{ledger}: the order of {month} was being submitted when an earlier run stopped, and WEkEO's "
                f"answer (its order id) was not recorded, so its remote state is unknown. Nothing was requested. "
                f"Check your WEkEO account, then remove the ledger to submit again."
            )
    return entries


def _save_ledger(path, entries):
    _write_json(ledger_path(path), {"schema": "lectk-wekeo-orders/1", "dataset": WEKEO_DATASET, "orders": entries})


# -------------------------------------------------------------------- order and download
class _Budget:
    def __init__(self, limit):
        self.limit, self.used = int(limit), 0

    def remaining(self):
        return self.limit - self.used


def submit_order(client, query):
    """Search the request and return (order query, order id): the order is submitted here."""
    results = client.search(query)
    items = list(getattr(results, "results", results))
    if len(items) != 1:
        raise WekeoContentError(f"WEkEO found {len(items)} results for the request, expected 1")
    item = items[0]
    return {"dataset_id": WEKEO_DATASET, "product_id": item["id"], "location": item["properties"]["location"]}


def wait_order(client, order_id, max_seconds=ORDER_MAX_SECONDS, sleep=time.sleep, clock=time.monotonic):
    """Wait until the order is ready (HTTP 200), up to ``max_seconds``; return the seconds waited."""
    start, delay = clock(), POLL_FIRST
    while True:
        status = client.head("dataaccess/download", order_id).status_code
        if status == 200:
            return clock() - start
        if status != 202:
            raise WekeoOrderError(f"WEkEO order {order_id} failed (HTTP {status})")
        if clock() - start + delay > max_seconds:
            raise WekeoOrderError(
                f"WEkEO order {order_id} not ready after {max_seconds:g} s; it is recorded in the ledger and "
                f"a new run with the same destination resumes it (no new order)."
            )
        sleep(delay)
        delay = min(delay * 1.5, POLL_MAX)


def download_order(client, order_id, target, budget, attempts=DOWNLOAD_ATTEMPTS, logger=None):
    """
    Download a ready order to ``target``, counting every byte against ``budget``.

    Returns:
        int: bytes of the file.
    """
    import requests

    url = client.full_url("dataaccess/download", order_id)
    last = None
    for attempt in range(1, attempts + 1):
        received = 0
        try:
            with client.session.get(url, stream=True, timeout=client.timeout) as response:
                response.raise_for_status()
                length = response.headers.get("Content-Length")
                if length is not None and int(length) > budget.remaining():
                    raise WekeoVolumeError(
                        f"order {order_id}: the response has {int(length)} bytes and only {budget.remaining()} "
                        f"of the {budget.limit} allowed remain; not read."
                    )
                with open(target, "wb") as f:
                    for block in response.iter_content(CHUNK):
                        received += len(block)
                        if received > budget.remaining():
                            raise WekeoVolumeError(
                                f"order {order_id}: the response passed the byte limit ({budget.limit}); stopped."
                            )
                        f.write(block)
                if length is not None and received != int(length):
                    raise requests.exceptions.ConnectionError(
                        f"incomplete response: {received} of {length} bytes")
            budget.used += received
            return received
        except requests.exceptions.RequestException as exc:
            budget.used += received
            last = exc
            if logger is not None:
                logger.warning(f"⚠️ WEkEO download of order {order_id}, attempt {attempt}/{attempts} failed: "
                               f"{type(exc).__name__}")
        except WekeoVolumeError:
            budget.used += received
            raise
    raise WekeoOrderError(f"download of WEkEO order {order_id} failed after {attempts} attempts: {last}")


def netcdf_members(path, directory):
    """NetCDF files of a response: the file itself, or the .nc members of a zip."""
    with open(path, "rb") as f:
        head = f.read(8)
    if head.startswith(NETCDF_SIGNATURES):
        return [path]
    if head.startswith(ZIP_SIGNATURE):
        out = []
        with zipfile.ZipFile(path) as archive:
            for name in archive.namelist():
                if name.endswith(".nc") and os.path.basename(name) == name:
                    archive.extract(name, directory)
                    out.append(os.path.join(directory, name))
        if out:
            return sorted(out)
    raise WekeoContentError(f"{path}: the WEkEO response is not NetCDF")


# -------------------------------------------------------------------- assembly
def assemble(files, plan, path):
    """
    Concatenate the monthly files along valid_time, keep the times of the CDS
    route and write ``path`` (through a temporary name). Values are not changed.
    """
    datasets = [xr.open_dataset(f) for f in files]
    try:
        for f, ds in zip(files, datasets):
            missing = [v for v in SHORT_NAMES if v not in ds.data_vars]
            if missing or TIME_DIMENSION not in ds.dims:
                raise WekeoContentError(f"{f}: missing {missing or [TIME_DIMENSION]}")
        combined = xr.concat(datasets, dim=TIME_DIMENSION) if len(datasets) > 1 else datasets[0]
        wanted = pd.DatetimeIndex(plan["times"])
        have = pd.DatetimeIndex(combined[TIME_DIMENSION].values)
        absent = wanted.difference(have)
        if len(absent):
            raise WekeoContentError(f"{len(absent)} requested time(s) missing from the WEkEO files, first {absent[0]}")
        combined = combined.sel({TIME_DIMENSION: wanted})
        _check_grid(combined, plan)
        tmp = f"{path}.wekeo-tmp.nc"
        combined.to_netcdf(tmp)
        combined.close()
    finally:
        for ds in datasets:
            ds.close()
    os.replace(tmp, path)


def _check_grid(ds, plan):
    north, west, south, east = plan["area"]
    levels = sorted(float(x) for x in ds["pressure_level"].values)
    problems = []
    if levels != sorted(float(x) for x in PRESSURE_LEVELS):
        problems.append(f"pressure levels {levels}")
    lat, lon = ds["latitude"].values, ds["longitude"].values
    if not (math.isclose(lat.max(), north) and math.isclose(lat.min(), south)):
        problems.append(f"latitudes {lat.min()}..{lat.max()} instead of {south}..{north}")
    if not (math.isclose(lon.min(), west) and math.isclose(lon.max(), east)):
        problems.append(f"longitudes {lon.min()}..{lon.max()} instead of {west}..{east}")
    if problems:
        raise WekeoContentError("the WEkEO files do not match the request: " + "; ".join(problems))


# -------------------------------------------------------------------- entry point
def get_wekeo_data(args, track, app_logger, client=None, sleep=time.sleep):
    """
    Retrieve the ERA5 file of a track from WEkEO into ``args.infile``.

    Args:
        args (argparse.Namespace): ``infile``, ``time_resolution`` and ``max_gib`` (None: 20 GiB).
        track (pd.DataFrame): Track with 'Lat', 'Lon' and a datetime index.
        app_logger (logging.Logger): Logger.
        client: hda.Client (tests pass a simulated one); default: :func:`make_client`.
        sleep (callable): wait between status checks (tests pass a stub).

    Returns:
        str: ``args.infile``.
    """
    path = args.infile
    max_gib = getattr(args, "max_gib", None)
    limit = DEFAULT_MAX_TOTAL_BYTES if max_gib is None else int(max_gib * 1024 ** 3)
    plan = plan_requests(track, args.time_resolution)
    directory = os.path.dirname(os.path.abspath(path))
    north, west, south, east = plan["area"]
    app_logger.info(f"🌐 WEkEO: {len(plan['requests'])} monthly request(s), area N{north} W{west} S{south} E{east}, "
                    f"{len(plan['times'])} times, up to {plan['max_bytes'] / 1024 ** 3:.2f} GiB "
                    f"(limit {limit / 1024 ** 3:.2f} GiB)")
    check_volume(plan, limit, directory)
    entries = read_ledger(path, plan)
    if client is None:
        client = make_client()

    started = datetime.datetime.now(datetime.timezone.utc)
    budget = _Budget(limit)
    parts = parts_dir(path)
    os.makedirs(parts, exist_ok=True)
    files, records = [], []
    for request in plan["requests"]:
        month = request["month"]
        entry = entries.get(month)
        target = os.path.join(parts, f"wekeo_{month}.download")
        if entry and entry.get("status") == "downloaded" and os.path.exists(target) \
                and _sha256(target) == entry.get("sha256"):
            app_logger.info(f"  ♻️ {month}: response already downloaded and verified (ledger)")
        else:
            if entry and entry.get("order_id"):
                order_id = entry["order_id"]
                app_logger.info(f"  ↩️ {month}: resuming WEkEO order {order_id} recorded in the ledger")
            else:
                order_query = submit_order(client, request["query"])
                entry = {"query": request["query"], "status": "submitting", "order_query": order_query,
                         "submitted_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
                entries[month] = entry
                _save_ledger(path, entries)  # recorded before the order exists
                order_id = client.post(order_query, "dataaccess/download")["download_id"]
                entry.update(status="submitted", order_id=order_id)
                _save_ledger(path, entries)
                app_logger.info(f"  📨 {month}: WEkEO order {order_id} submitted")
            waited = wait_order(client, order_id, sleep=sleep)
            t0 = time.monotonic()
            size = download_order(client, order_id, target, budget, logger=app_logger)
            entry.update(status="downloaded", bytes=size, sha256=_sha256(target), wait_seconds=round(waited, 3),
                         transfer_seconds=round(time.monotonic() - t0, 3))
            _save_ledger(path, entries)
            app_logger.info(f"  ✓ {month}: {size / 1024 ** 2:.2f} MB after {waited:.0f} s in the WEkEO queue")
        files += netcdf_members(target, parts)
        records.append({"month": month, "query": request["query"], "order_id": entry.get("order_id"),
                        "bytes": entry.get("bytes"), "sha256": entry.get("sha256"),
                        "wait_seconds": entry.get("wait_seconds"), "transfer_seconds": entry.get("transfer_seconds")})

    t0 = time.monotonic()
    assemble(files, plan, path)
    provenance = {
        "schema": PROVENANCE_SCHEMA,
        "route": "wekeo",
        "dataset": WEKEO_DATASET,
        "hda_version": hda_version(),
        "lectk_commit": _lectk_commit(),
        "python": platform.python_version(),
        "xarray": xr.__version__,
        "area_nwse": list(plan["area"]),
        "times": [f"{t:%Y-%m-%dT%H:%MZ}" for t in plan["times"]],
        "time_resolution_h": args.time_resolution,
        "limit_bytes": limit,
        "bytes_received_this_run": budget.used,
        "requests": records,
        "assembly_seconds": round(time.monotonic() - t0, 3),
        "started_utc": started.isoformat(),
        "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "output": {"path": os.path.abspath(path), "bytes": os.path.getsize(path), "sha256": _sha256(path)},
    }
    _write_json(provenance_path(path), provenance)
    os.remove(ledger_path(path))
    shutil.rmtree(parts, ignore_errors=True)
    app_logger.info(f"🎉 WEkEO file ready: {path} ({provenance['output']['bytes'] / 1024 ** 2:.2f} MB, "
                    f"sha256 {provenance['output']['sha256'][:12]}…); provenance in {provenance_path(path)}")
    return path

"""
Offline tests of the DestinE Earth Data Hub route (src/utils/destine_edh.py, --destine).

No network is used. The remote store is a local zarr store (a dict of keys to
bytes) holding a synthetic ERA5 shaped like the Earth Data Hub stores (19 levels
in a dimension ``isobaricInhPa``, hourly ``valid_time``, 0.25-degree grid,
latitudes north to south, longitudes 0..360, chunks of time x 1 level x a tile),
in Zarr format 3 as the default store (v3) or 2 as the store of the 2b benchmark
(v2), whose values depend only on variable, time, level, latitude and longitude.
The tests that write or read zarr need zarr 3 (``pip install
"lorenzcycletoolkit[destine]"``) and are skipped without it; the others
(resolution, credential, counting, retries, command line) do not.
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import types

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from src.utils import destine_edh as d
from src.utils.ncar_gdex import area_from_track
from src.utils.track_alignment import planned_request_times

pytestmark = pytest.mark.filterwarnings("ignore:Consolidated metadata is currently not part:UserWarning")  # test stores (v3)
LOGGER = logging.getLogger("test_destine")
LEVELS = d.EXPECTED_LEVELS  # 1000 ... 1 hPa, as the store
UNITS = {"u": "m s**-1", "v": "m s**-1", "t": "K", "w": "Pa s**-1", "z": "m**2 s**-2"}


def args_for(path, time_resolution=6, max_gib=None, max_requests=None, destine_store=None):
    return argparse.Namespace(infile=str(path), time_resolution=time_resolution, max_gib=max_gib,
                              max_requests=max_requests, destine_store=destine_store)


FORMATS = {"v3": 3, "v2": 2}


def make_track(start, end, lat=(-30, -28), lon=(-50, -48), freq="6h"):
    index = pd.date_range(start, end, freq=freq)
    return pd.DataFrame({"Lat": np.linspace(*lat, len(index)), "Lon": np.linspace(*lon, len(index))}, index=index)


@pytest.fixture
def track():
    return make_track("2020-01-30 00:00", "2020-01-30 12:00")  # 3 times; area N-13 W-65 S-45 E-33


def field(short, times, levels, lat, lon):
    """Synthetic value: depends on the variable, time, level, latitude and longitude (-180..180)."""
    k = d.SHORT_NAMES.index(short) + 1
    hours = ((pd.DatetimeIndex(times) - pd.Timestamp("2020-01-01")) / pd.Timedelta(hours=1)).values
    lon = ((np.asarray(lon) + 180.0) % 360.0) - 180.0
    return (k * (hours[:, None, None, None] * 1e-2 + np.asarray(levels, float)[None, :, None, None] * 1e-3
                 + np.asarray(lat)[None, None, :, None] * 1e-1 + lon[None, None, None, :] * 1e-3)).astype("float32")


def memory_store(mapping):
    """A zarr 3 MemoryStore holding the keys and bytes of ``mapping``."""
    import zarr
    from zarr.core.buffer import default_buffer_prototype
    return zarr.storage.MemoryStore(
        {key: default_buffer_prototype().buffer.from_bytes(value) for key, value in mapping.items()})


def key_bytes(memory):
    return {key: value.to_bytes() for key, value in memory._store_dict.items()}


def build_store(lat, lon_signed, times, convention="360", levels=LEVELS, chunk_time=120, chunk_lat=32, chunk_lon=32,
                level_name="isobaricInhPa", version="v3", shards=None):
    """
    Synthetic DestinE store as a dict of keys to bytes: consolidated metadata, chunks of time x 1 level x a tile,
    Zarr format 3 (``version`` "v3", the default store) or 2 ("v2"); ``shards``: shard shape (v3 only).
    """
    pytest.importorskip("zarr")
    lon_signed = np.asarray(lon_signed, float)
    lon_store = lon_signed % 360.0 if convention == "360" else lon_signed
    order = np.argsort(lon_store)
    lon_store, lon_signed = lon_store[order], lon_signed[order]
    dims = ("valid_time", level_name, "latitude", "longitude")
    data, encoding = {}, {}
    for short in d.SHORT_NAMES:
        data[short] = (dims, field(short, times, levels, lat, lon_signed), {"units": UNITS[short]})
        encoding[short] = {"chunks": (chunk_time, 1, chunk_lat, chunk_lon)}
        if shards:
            encoding[short]["shards"] = shards
    ds = xr.Dataset(data, coords={"valid_time": pd.DatetimeIndex(times), level_name: np.asarray(levels, float),
                                  "latitude": np.asarray(lat, float), "longitude": lon_store})
    import zarr
    memory = zarr.storage.MemoryStore()
    ds.to_zarr(memory, mode="w", consolidated=True, encoding=encoding, zarr_format=FORMATS[version])
    return key_bytes(memory)


def lengthen_time_axis(store, periods):
    """
    The same store with ``periods`` hourly steps: the time coordinate is extended and the data arrays are resized
    in their metadata only (the chunks beyond the first time chunk are never written nor read).
    """
    import zarr
    memory = memory_store(store)
    group = zarr.open_group(memory, mode="r+", use_consolidated=False)
    time = group["valid_time"]
    values = time[:]
    step = values[1] - values[0]
    time.resize((periods,))
    time[:] = values[0] + step * np.arange(periods)
    for short in d.SHORT_NAMES:
        group[short].resize((periods,) + group[short].shape[1:])
    zarr.consolidate_metadata(memory)
    return key_bytes(memory)


class Recorder(dict):
    """Store mapping that remembers the keys read."""

    def __init__(self, *args):
        super().__init__(*args)
        self.reads = []

    def __getitem__(self, key):
        self.reads.append(key)
        return super().__getitem__(key)


def chunk_reads(reads):
    return [k for k in reads if k.split("/")[0] in d.SHORT_NAMES and not d._is_metadata(k)]


LAT = np.arange(-12.0, -46.0 - 1e-9, -0.25)  # covers the area -13..-45, north to south (kept small: memory)
LON = np.arange(-66.0, -32.0 + 1e-9, 0.25)
HOURS = pd.date_range("2020-01-30 00:00", periods=13, freq="1h")  # the 3 times of the track, hourly


@pytest.fixture
def store():
    return Recorder(build_store(LAT, LON, HOURS))


# ------------------------------------------------------------------ request
def test_plan_as_cds_route(track):
    plan = d.plan_request(track, 6)
    assert plan["area"] == area_from_track(track) == (-13, -65, -45, -33)
    assert list(plan["times"]) == list(planned_request_times(track.index, 6))
    assert plan["latitudes"][0] == -13 and plan["latitudes"][-1] == -45 and np.all(np.diff(plan["latitudes"]) < 0)
    assert plan["longitudes"][0] == -65 and plan["longitudes"][-1] == -33 and np.all(np.diff(plan["longitudes"]) > 0)


def test_plan_greenwich_is_one_arc_west_to_east():
    plan = d.plan_request(make_track("2020-03-01", "2020-03-02", lon=(-5, 5)), 6)
    assert (plan["area"][1], plan["area"][3]) == (-20, 20)
    assert plan["longitudes"][0] == -20 and plan["longitudes"][-1] == 20 and np.all(np.diff(plan["longitudes"]) > 0)


def test_plan_refuses_180_crossing():
    with pytest.raises(d.DestineRequestError, match="crosses 180"):
        d.plan_request(make_track("2020-03-01", "2020-03-02", lon=(170, 178)), 6)


@pytest.mark.parametrize("step", [0, -3, 1.5])
def test_time_resolution_not_whole_hours_refused(step):
    with pytest.raises(d.DestineRequestError, match="whole number of hours"):
        d.check_time_resolution(step)


@pytest.mark.parametrize("step", [1, 3, 6])
def test_hourly_product_serves_any_whole_hour_step(step, store, tmp_path):
    track = make_track("2020-01-30 00:00", "2020-01-30 12:00", freq=f"{step}h")
    out = d.get_destine_data(args_for(tmp_path / f"h{step}.nc", step), track, LOGGER, mapping=store)
    with xr.open_dataset(out) as ds:
        assert list(ds.valid_time.values) == list(planned_request_times(track.index, step).values)


# ------------------------------------------------------------------ output contract, values, levels
@pytest.mark.parametrize("version", ["v3", "v2"])
@pytest.mark.parametrize("convention", ["360", "signed"])
def test_file_follows_the_cds_contract_with_the_values_of_the_store(convention, version, track, tmp_path):
    mapping = build_store(LAT, LON, HOURS, convention=convention, version=version)
    out = d.get_destine_data(args_for(tmp_path / "destine.nc", destine_store=version), track, LOGGER,
                             mapping=mapping)
    plan = d.plan_request(track, 6)
    with xr.open_dataset(out) as ds:
        assert set(ds.data_vars) == set("uvtwz")
        for short in "uvtwz":
            v = ds[short]
            assert v.dims == ("valid_time", "pressure_level", "latitude", "longitude") and v.dtype == np.float32
            assert v.attrs["units"] == UNITS[short]
            expected = field(short, plan["times"], LEVELS, plan["latitudes"], plan["longitudes"])
            np.testing.assert_array_equal(v.values, expected)  # values unchanged
        assert list(ds.pressure_level.values) == list(LEVELS)  # the 19 levels of the product, nothing interpolated
        np.testing.assert_array_equal(ds.latitude.values, plan["latitudes"])  # north to south
        np.testing.assert_array_equal(ds.longitude.values, plan["longitudes"])  # west to east, -180..180 as CDS
        assert ds.valid_time.encoding["units"] == "seconds since 1970-01-01"


def test_greenwich_with_the_store_in_0_360(tmp_path):
    track = make_track("2020-03-01 00:00", "2020-03-01 12:00", lat=(-10, -8), lon=(-5, 5))
    lat = np.arange(8.0, -26.0 - 1e-9, -0.25)
    lon = np.arange(-25.0, 25.0 + 1e-9, 0.25)  # stored as 0..25 and 335..359.75
    mapping = build_store(lat, lon, pd.date_range("2020-03-01 00:00", periods=13, freq="1h"))
    out = d.get_destine_data(args_for(tmp_path / "g.nc"), track, LOGGER, mapping=mapping)
    plan = d.plan_request(track, 6)
    with xr.open_dataset(out) as ds:
        assert ds.longitude.values[0] == -20 and ds.longitude.values[-1] == 20 and 0.0 in ds.longitude.values
        np.testing.assert_array_equal(ds["t"].values, field("t", plan["times"], LEVELS, plan["latitudes"],
                                                            plan["longitudes"]))


def test_levels_of_the_store_in_pa_and_in_another_number(track, tmp_path):
    levels = (1000.0, 850.0, 500.0, 250.0, 100.0)
    mapping = build_store(LAT, LON, HOURS,
                          levels=[x * 100 for x in levels], level_name="level")
    out = d.get_destine_data(args_for(tmp_path / "pa.nc"), track, LOGGER, mapping=mapping)
    with xr.open_dataset(out) as ds:
        assert list(ds.pressure_level.values) == list(levels)
    prov = json.loads((tmp_path / "pa.nc.provenance.json").read_text())
    assert prov["levels_hpa"] == list(levels) and "(5, 1000 to 100 hPa" in prov["warnings"][0]


def test_namelist_of_the_cds_contract_is_enough(store, track, tmp_path):
    from src.utils.validation import validate_namelist_file, validate_variable_match
    from src.utils import preprocessing
    out = d.get_destine_data(args_for(tmp_path / "n.nc"), track, LOGGER, mapping=store)
    namelist = validate_namelist_file(preprocessing.ERA5_NAMELIST, LOGGER)
    with xr.open_dataset(out) as ds:
        validate_variable_match(namelist, ds, preprocessing.ERA5_NAMELIST, LOGGER)  # raises if a name is absent
    ns = argparse.Namespace(cdsapi=False, wekeo=False, arco=False, destine=True, ncar=False)
    assert preprocessing.namelist_for_run(ns) == preprocessing.ERA5_NAMELIST


def test_warnings_in_the_log_and_in_the_provenance(store, track, tmp_path, caplog, monkeypatch):
    monkeypatch.setenv(d.TOKEN_ENV, "SECRET-TOKEN-VALUE")
    with caplog.at_level(logging.INFO, logger="test_destine"):
        out = d.get_destine_data(args_for(tmp_path / "p.nc"), track, LOGGER, mapping=store)
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any("19, 1000 to 1 hPa" in w and "37 levels" in w for w in warnings)  # (a) number of levels
    assert any("13 bits" in w and "18 of 240 terms" in w and "between 1 % and 2.7 %" in w
               and "conversions did not change" in w for w in warnings)  # (b) reduced precision
    text = (tmp_path / "p.nc.provenance.json").read_text()
    prov = json.loads(text)
    assert prov["route"] == "destine" and len(prov["levels_hpa"]) == 19
    assert prov["store_version"] == "v3" and prov["zarr_format"] == 3
    assert prov["warnings"][0] == d.levels_warning(LEVELS) and prov["warnings"][1] == d.PRECISION_WARNING_V3
    assert prov["area_nwse"] == [-13, -65, -45, -33] and len(prov["times"]) == 3
    assert prov["lectk_commit"] is None or len(prov["lectk_commit"]) == 40
    assert set(prov["packages"]) == {"zarr", "numcodecs", "fsspec", "aiohttp"}
    assert len(prov["output"]["sha256"]) == 64 and prov["output"]["bytes"] == os.path.getsize(out)
    with xr.open_dataset(out) as ds:
        assert "18 of 240 terms" in ds.attrs["lectk_provenance_note"]
    assert "SECRET-TOKEN-VALUE" not in text and "SECRET-TOKEN-VALUE" not in caplog.text


# ------------------------------------------------------------------ chunks read once (the graph does not grow)
@pytest.mark.parametrize("version, chunk_time", [("v3", 1440), ("v2", 120)])
def test_each_chunk_is_read_once_whatever_the_length_of_the_store(version, chunk_time, track, tmp_path):
    """No dask graph over the store: the chunks touched are read once each, for a short and a long store."""
    counts = []
    for periods in (13, 24 * 400):  # 4 days and 400 days of hours
        mapping = build_store(LAT, LON, HOURS, version=version, chunk_time=chunk_time)
        if periods > 13:
            mapping = lengthen_time_axis(mapping, periods)  # metadata only: no data held in memory
        mapping = Recorder(mapping)
        out = d.get_destine_data(args_for(tmp_path / f"once{periods}.nc", destine_store=version), track, LOGGER,
                                 mapping=mapping)
        reads = chunk_reads(mapping.reads)
        prov = json.loads((tmp_path / f"once{periods}.nc.provenance.json").read_text())
        assert len(reads) == len(set(reads)) == prov["estimate_chunks"]  # each chunk once, all foreseen
        assert prov["chunks_read_more_than_once"] == 0
        counts.append(len(reads))
        assert os.path.exists(out)
    assert counts[0] == counts[1]  # the length of the store does not change what is read


def test_no_dask_array_is_built(store, track):
    pytest.importorskip("zarr")
    ds = d.arco_zarr.open_store_dataset(d.CountingStore(store, d.Budget(10 ** 12, 10 ** 6)), 3)
    assert all(v.chunks is None for v in ds.data_vars.values())
    selection = d.select(ds, d.plan_request(track, 6))
    block = d.read_level(ds, selection, "t", 0)
    assert isinstance(block, np.ndarray) and block.shape == (3, 129, 129)


# ------------------------------------------------------------------ estimates and limits
def test_estimate_is_an_upper_bound_and_counts_requests(store, track, tmp_path):
    out = d.get_destine_data(args_for(tmp_path / "e.nc"), track, LOGGER, mapping=store)
    prov = json.loads((tmp_path / "e.nc.provenance.json").read_text())
    assert prov["bytes_received"] <= prov["estimate_bytes"]
    assert prov["estimate_raw_chunk_bytes"] < prov["estimate_bytes"] - prov["bytes_metadata"]  # codec margin
    assert prov["requests_store"] == prov["estimate_requests"] == len(store.reads)  # no retry: exact
    assert prov["requests_store"] == prov["requests_metadata"] + prov["estimate_chunks"]
    # chunks: 1 time chunk x 19 levels x latitude tiles x longitude tiles, 5 fields
    assert prov["estimate_chunks"] % (5 * 19) == 0 and os.path.exists(out)
    assert prov["estimate_bytes"] > d.output_bytes(dict(d.plan_request(track, 6), levels=LEVELS))


def test_requests_above_the_limit_read_no_chunk(store, track, tmp_path):
    out = tmp_path / "r.nc"
    with pytest.raises(d.DestineVolumeError, match=r"--max-requests \d+ or more") as error:
        d.get_destine_data(args_for(out, max_requests=20), track, LOGGER, mapping=store)
    assert "no chunk was read" in str(error.value) and "monthly request quota" in str(error.value)
    assert chunk_reads(store.reads) == [] and not out.exists()
    needed = int(str(error.value).split("--max-requests ")[1].split()[0])
    assert os.path.exists(d.get_destine_data(args_for(out, max_requests=needed), track, LOGGER, mapping=store))


def test_bytes_above_the_limit_read_no_chunk(store, track, tmp_path):
    out = tmp_path / "v.nc"
    with pytest.raises(d.DestineVolumeError, match=r"--max-gib [0-9.]+ or more") as error:
        d.get_destine_data(args_for(out, max_gib=0.001), track, LOGGER, mapping=store)
    assert "no chunk was read" in str(error.value)
    assert chunk_reads(store.reads) == [] and not out.exists()
    needed = float(str(error.value).split("--max-gib ")[1].split()[0])
    assert os.path.exists(d.get_destine_data(args_for(out, max_gib=needed), track, LOGGER, mapping=store))


def test_disk_space_checked_before_reading(store, track, tmp_path, monkeypatch):
    monkeypatch.setattr(d.shutil, "disk_usage", lambda directory: types.SimpleNamespace(free=1000))
    with pytest.raises(d.DestineVolumeError, match="of disk"):
        d.get_destine_data(args_for(tmp_path / "d.nc"), track, LOGGER, mapping=store)
    assert chunk_reads(store.reads) == []


def test_real_count_stops_at_the_byte_limit():
    inner = Recorder({"a": b"x" * 600, "b": b"y" * 600, "c": b"z" * 600})
    store = d.CountingStore(inner, d.Budget(1000, 100))
    assert len(store["a"]) == 600
    with pytest.raises(d.DestineVolumeError, match="passed reading b"):
        store["b"]
    with pytest.raises(d.DestineVolumeError, match="not requested"):
        store["c"]
    assert inner.reads == ["a", "b"] and store.budget.used == 1200


def test_real_count_stops_at_the_request_limit():
    inner = Recorder({"a": b"x", "b": b"y", "c": b"z"})
    store = d.CountingStore(inner, d.Budget(10 ** 6, 2))
    store["a"], store["b"]
    with pytest.raises(d.DestineVolumeError, match="request limit"):
        store["c"]
    assert inner.reads == ["a", "b"] and store.budget.attempts == 2


def test_failed_attempts_are_limited_and_count_as_requests():
    calls = []

    class Flaky(dict):
        def __getitem__(self, key):
            calls.append(key)
            if len(calls) < 3:
                raise OSError("connection reset")
            return b"ok"

    waits = []
    store = d.CountingStore(Flaky(), d.Budget(100, 100), sleep=waits.append)
    assert store["k"] == b"ok" and len(calls) == 3 and waits == [2.0, 4.0]
    assert store.budget.retries == 2 and store.budget.attempts == 3 and store.budget.used == 2

    class Down(dict):
        def __getitem__(self, key):
            raise OSError("down")

    down = d.CountingStore(Down(), d.Budget(100, 100), sleep=lambda s: None)
    with pytest.raises(d.DestineTransferError, match="after 3 attempts"):
        down["k"]
    assert down.budget.attempts == 3 and down.budget.used == 0


@pytest.mark.parametrize("status, error, match", [
    (401, "DestineSetupError", "refused the token"),
    (403, "DestineSetupError", "refused the token"),
    (429, "DestineVolumeError", "429"),
])
def test_refused_token_and_quota_are_not_retried(status, error, match):
    exc = Exception("refused")
    exc.status = status
    calls = []

    class Denied(dict):
        def __getitem__(self, key):
            calls.append(key)
            raise exc

    store = d.CountingStore(Denied(), d.Budget(100, 100), sleep=lambda s: None)
    with pytest.raises(getattr(d, error), match=match):
        store["k"]
    assert len(calls) == 1


def test_missing_chunk_is_not_filled_with_nan(store, track, tmp_path):
    key = next(k for k in sorted(store) if k.startswith("t/") and not d._is_metadata(k))
    del store[key]
    with pytest.raises(d.DestineContentError, match="absent from the store"):
        d.get_destine_data(args_for(tmp_path / "m.nc"), track, LOGGER, mapping=store)
    assert not (tmp_path / "m.nc").exists() and not (tmp_path / "m.nc.destine-tmp.nc").exists()


@pytest.mark.parametrize("version, message", [("v3", "updated monthly"), ("v2", "no longer updated")])
def test_time_outside_the_store_is_refused_before_any_chunk(version, message, track, tmp_path):
    mapping = Recorder(build_store(LAT, LON, pd.date_range("2019-01-01", periods=2, freq="1h"), version=version))
    with pytest.raises(d.DestineRequestError, match=message):
        d.get_destine_data(args_for(tmp_path / "t.nc", destine_store=version), track, LOGGER, mapping=mapping)
    assert chunk_reads(mapping.reads) == []


def test_sharded_store_is_refused_before_any_chunk(track, tmp_path):
    """The estimate counts whole chunks, one request each: a sharded store (partial reads) is refused."""
    mapping = Recorder(build_store(LAT, LON, HOURS, chunk_time=12, chunk_lat=16, chunk_lon=16, shards=(24, 1, 32, 32)))
    with pytest.raises(d.DestineContentError, match="sharded"):
        d.get_destine_data(args_for(tmp_path / "s.nc"), track, LOGGER, mapping=mapping)
    assert chunk_reads(mapping.reads) == [] and not (tmp_path / "s.nc").exists()


def test_v2_store_option_and_unknown_store(store, track, tmp_path):
    mapping = build_store(LAT, LON, HOURS, version="v2")
    out = d.get_destine_data(args_for(tmp_path / "v2.nc", destine_store="v2"), track, LOGGER, mapping=mapping)
    prov = json.loads((tmp_path / "v2.nc.provenance.json").read_text())
    assert prov["store_version"] == "v2" and prov["zarr_format"] == 2 and prov["warnings"][1] == d.PRECISION_WARNING
    assert os.path.exists(out)
    assert d.STORES["v3"]["url"].endswith("/era5/era5-pressure-levels-v0.zarr") and d.STORE_URL == d.STORES["v3"]["url"]
    assert d.STORES["v2"]["url"].endswith("/era5/reanalysis-era5-pressure-levels-v0.zarr")
    with pytest.raises(d.DestineRequestError, match="unknown Earth Data Hub store"):
        d.get_destine_data(args_for(tmp_path / "x.nc", destine_store="v9"), track, LOGGER, mapping=store)


def test_at_most_four_reads_at_a_time(track, tmp_path):
    """zarr 3 reads chunks concurrently: with 4 worker threads the counted store still lets at most 4 run."""
    import threading
    import time as clock
    lock, state = threading.Lock(), {"now": 0, "peak": 0}

    class Slow(dict):
        def __getitem__(self, key):
            with lock:
                state["now"] += 1
                state["peak"] = max(state["peak"], state["now"])
            clock.sleep(0.005)
            with lock:
                state["now"] -= 1
            return super().__getitem__(key)

    mapping = Slow(build_store(LAT, LON, HOURS, chunk_lat=16, chunk_lon=16))
    d.get_destine_data(args_for(tmp_path / "par.nc"), track, LOGGER, mapping=mapping)
    prov = json.loads((tmp_path / "par.nc.provenance.json").read_text())
    assert 1 < state["peak"] <= d.WORKERS and prov["max_parallel_reads"] == state["peak"]


# ------------------------------------------------------------------ credential and packages
def test_token_from_environment_and_netrc(tmp_path):
    assert d.read_token({d.TOKEN_ENV: " abc "}) == "abc"
    rc = tmp_path / "netrc"
    rc.write_text("machine data.earthdatahub.destine.eu\n  password from-netrc\n"
                  "machine api.earthdatahub.destine.eu\n  password other\n")  # the Earth Data Hub example
    rc.chmod(0o600)
    assert d.read_token({}, netrc_file=str(rc)) == "from-netrc"
    assert d.read_token({"NETRC": str(rc)}) == "from-netrc"


def test_missing_token_says_how_to_fix(tmp_path):
    with pytest.raises(d.DestineSetupError) as error:
        d.read_token({"NETRC": str(tmp_path / "absent")})
    text = str(error.value)
    assert d.TOKEN_ENV in text and "~/.netrc" in text and d.STORE_HOST in text and d.TOKEN_URL in text


def test_missing_token_stops_before_any_request(track, tmp_path, monkeypatch):
    pytest.importorskip("zarr")
    pytest.importorskip("aiohttp")
    monkeypatch.delenv(d.TOKEN_ENV, raising=False)
    monkeypatch.setenv("NETRC", str(tmp_path / "absent"))
    monkeypatch.setattr(d, "open_remote_mapping", lambda *a, **k: (_ for _ in ()).throw(AssertionError("opened")))
    with pytest.raises(d.DestineSetupError, match="token not found"):
        d.get_destine_data(args_for(tmp_path / "c.nc"), track, LOGGER)
    assert not (tmp_path / "c.nc").exists()


def test_missing_package(monkeypatch):
    real = d.importlib.import_module

    def fake_import(name, *args):
        if name == "zarr":
            raise ImportError("No module named 'zarr'")
        return real(name, *args)

    monkeypatch.setattr(d.importlib, "import_module", fake_import)
    with pytest.raises(d.DestineSetupError, match=r"lorenzcycletoolkit\[destine\]"):
        d.load_packages()


def test_zarr_2_is_refused(monkeypatch):
    fakes = {"zarr": types.SimpleNamespace(__version__="2.18.4"), "fsspec": types.SimpleNamespace(),
             "aiohttp": types.SimpleNamespace()}
    real = d.importlib.import_module
    monkeypatch.setattr(d.importlib, "import_module", lambda name, *x: fakes[name] if name in fakes else real(name, *x))
    with pytest.raises(d.DestineSetupError, match="needs zarr 3"):
        d.load_packages()


def test_basic_auth_and_request_trace_only_in_the_mapper():
    seen = {}

    class Trace:
        def __init__(self):
            self.on_request_end, self.on_request_exception = [], []

    fsspec = types.SimpleNamespace(get_mapper=lambda url, **kw: seen.update(url=url, **kw) or {})
    aiohttp = types.SimpleNamespace(ClientTimeout=lambda total: ("timeout", total), TraceConfig=Trace,
                                    BasicAuth=lambda login, password: ("basic", login, password))
    counter = d.HttpCounter()
    d.open_remote_mapping("TOKEN", counter, packages={"fsspec": fsspec, "aiohttp": aiohttp})
    kwargs = seen["client_kwargs"]
    assert seen["url"] == d.STORE_URL and kwargs["auth"] == ("basic", "edh", "TOKEN")
    assert len(kwargs["trace_configs"][0].on_request_end) == 1
    assert "TOKEN" not in repr(d.CountingStore({}, d.Budget(1, 1)).__dict__)


def test_http_counter():
    import asyncio
    counter = d.HttpCounter()

    class Trace:
        def __init__(self):
            self.on_request_end, self.on_request_exception = [], []

    traces = []
    aiohttp = types.SimpleNamespace(ClientTimeout=lambda total: None, BasicAuth=lambda *a: None,
                                    TraceConfig=lambda: traces.append(Trace()) or traces[-1])
    d.open_remote_mapping("x", counter, packages={"aiohttp": aiohttp,
                                                  "fsspec": types.SimpleNamespace(get_mapper=lambda *a, **k: {})})
    params = types.SimpleNamespace(response=types.SimpleNamespace(status=200))
    asyncio.run(traces[0].on_request_end[0](None, None, params))
    asyncio.run(traces[0].on_request_exception[0](None, None, None))
    assert counter.summary() == {"requests": 2, "status": {"200": 1}, "exceptions": 1}


def test_existing_file_is_reused_by_the_command(tmp_path, monkeypatch):
    """With an existing input file nothing is read (same rule as --cdsapi)."""
    from src.utils import preprocessing
    (tmp_path / "x.nc").write_bytes(b"")
    monkeypatch.setattr(d, "get_destine_data", lambda *x, **k: (_ for _ in ()).throw(AssertionError("read")))
    ns = argparse.Namespace(infile=str(tmp_path / "x.nc"), cdsapi=False, wekeo=False, arco=False, destine=True,
                            ncar=False)
    with pytest.raises(Exception):  # the empty file is not NetCDF: opened, not downloaded
        preprocessing.get_data(ns, LOGGER)


# ------------------------------------------------------------------ command line
LAUNCHER = os.path.join(os.path.dirname(__file__), "run_cli_offline.py")
SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lorenzcycletoolkit.py")


@pytest.mark.parametrize("extra, message", [
    (["--destine", "--cdsapi"], "alternative data sources"),
    (["--destine", "--wekeo"], "alternative data sources"),
    (["--destine", "--ncar"], "alternative data sources"),
    (["--destine", "--arco"], "alternative data sources"),
    (["--max-requests", "10"], "--max-requests needs --destine"),
    (["--destine-store", "v2"], "--destine-store needs --destine"),
    (["--destine", "--destine-store", "v1"], "invalid choice"),
    (["--destine", "--max-requests", "0"], "--max-requests must be greater than 0"),
    (["--delete-downloaded-input"], "--delete-downloaded-input needs --cdsapi, --ncar, --wekeo or --arco (or --destine)"),
])
def test_cli_refusals(tmp_path, extra, message):
    result = subprocess.run([sys.executable, LAUNCHER, SCRIPT, str(tmp_path / "x.nc"), "-t", *extra],
                            capture_output=True, text=True, timeout=300)
    assert result.returncode == 2 and message in result.stderr
    assert not (tmp_path / "x.nc").exists()


def test_cli_accepts_destine_with_its_options():
    import lorenzcycletoolkit
    args = lorenzcycletoolkit.create_arg_parser().parse_args(
        ["x.nc", "-t", "--destine", "--time-resolution", "3", "--max-gib", "5", "--max-requests", "4000",
         "--delete-downloaded-input", "--trackfile", "inputs/track", "--destine-store", "v2"])
    assert args.destine and args.max_gib == 5 and args.max_requests == 4000 and args.delete_downloaded_input
    assert args.destine_store == "v2"
    assert not (args.cdsapi or args.arco)


def test_cli_help_states_the_limitations():
    import lorenzcycletoolkit
    text = " ".join(lorenzcycletoolkit.create_arg_parser().format_help().split())
    assert "--destine" in text and "18 of 240 terms" in text and "19 pressure levels" in text
    assert "--max-requests" in text and "500 000" in text
    assert "--destine-store" in text and "13 bits" in text and "10 bits" in text and "no longer updated" in text

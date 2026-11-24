"""
Offline tests of the CDS ARCO Zarr route (src/utils/arco_zarr.py, --arco).

No network is used. The remote store is a local zarr store (a dict of keys to
bytes) holding a synthetic ERA5 shaped like the product (13 levels, 6-hourly,
0.25-degree grid, chunked by time with whole latitude bands, longitudes in
-180..180 or 0..360), whose values depend only on variable, time, level,
latitude and longitude. The tests that write or read zarr need zarr 2
(``pip install "lorenzcycletoolkit[arco]"``) and are skipped without it; the
others (resolution, credential, byte counting, retries) do not.
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

from src.utils import arco_zarr as a
from src.utils.ncar_gdex import area_from_track
from src.utils.track_alignment import planned_request_times

LOGGER = logging.getLogger("test_arco")
STORE_LEVELS = (50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000)  # ascending, as a store would be
UNITS = {"u": "m s**-1", "v": "m s**-1", "t": "K", "w": "Pa s**-1", "z": "m**2 s**-2"}


def args_for(path, time_resolution=6, max_gib=None):
    return argparse.Namespace(infile=str(path), time_resolution=time_resolution, max_gib=max_gib)


def make_track(start, end, lat=(-30, -28), lon=(-50, -48)):
    index = pd.date_range(start, end, freq="6h")
    return pd.DataFrame({"Lat": np.linspace(*lat, len(index)), "Lon": np.linspace(*lon, len(index))}, index=index)


@pytest.fixture
def track():
    return make_track("2020-01-30 00:00", "2020-01-30 12:00")  # 3 times; area N-13 W-65 S-45 E-33


def field(short, times, levels, lat, lon):
    """Synthetic value of a variable: depends on its name, time, level, latitude and longitude (-180..180)."""
    k = a.SHORT_NAMES.index(short) + 1
    hours = ((pd.DatetimeIndex(times) - pd.Timestamp("2020-01-01")) / pd.Timedelta(hours=1)).values
    lon = ((np.asarray(lon) + 180.0) % 360.0) - 180.0
    return (k * (hours[:, None, None, None] * 1e-2 + np.asarray(levels, float)[None, :, None, None] * 1e-3
                 + np.asarray(lat)[None, None, :, None] * 1e-1 + lon[None, None, None, :] * 1e-3)).astype("float32")


def build_store(lat, lon_signed, times, convention="signed", dims_order=("time", "pressureLevels", "latitude", "longitude"),
                chunk_lat=40, fields="long", chunks=None):
    """
    Synthetic ARCO store as a zarr-2 dict: consolidated metadata, chunks of one time step and whole bands.

    ``fields``: "long" (names of the user guide) or "short" (names of the real store, observed 2026-11-24);
    ``chunks``: sizes by dimension name, to override the default (one time step, all levels, bands of ``chunk_lat``).
    """
    pytest.importorskip("zarr")
    lon_signed = np.asarray(lon_signed, float)
    lon_store = lon_signed % 360.0 if convention == "360" else lon_signed
    order = np.argsort(lon_store)
    lon_store, lon_signed = lon_store[order], lon_signed[order]
    levels = np.asarray(STORE_LEVELS)
    level_dim = dims_order[1] if "pressureLevel" in dims_order[1] else next(d for d in dims_order if "ressure" in d)
    data, encoding = {}, {}
    axes = {"time": 0, level_dim: 1, "latitude": 2, "longitude": 3}
    for short, candidates in a.STORE_VARIABLES.items():
        name = candidates[0] if fields == "short" else candidates[1]
        values = field(short, times, levels, lat, lon_signed)
        perm = [axes[d] for d in dims_order]
        data[name] = (dims_order, values.transpose(perm))
        sizes = {"time": 1, level_dim: len(levels), "latitude": chunk_lat, "longitude": len(lon_store)}
        sizes.update(chunks or {})
        encoding[name] = {"chunks": tuple(sizes[d] for d in dims_order)}
    ds = xr.Dataset(data, coords={"time": pd.DatetimeIndex(times), level_dim: levels,
                                  "latitude": np.asarray(lat, float), "longitude": lon_store})
    store = {}
    ds.to_zarr(store, mode="w", consolidated=True, encoding=encoding)
    return store


class Recorder(dict):
    """Store mapping that remembers the keys read."""

    def __init__(self, *args):
        super().__init__(*args)
        self.reads = []

    def __getitem__(self, key):
        self.reads.append(key)
        return super().__getitem__(key)


def times_6h(start, n):
    return pd.date_range(start, periods=n, freq="6h")


@pytest.fixture
def store(track):
    lat = np.arange(-10.0, -50.0 - 1e-9, -0.25)  # covers the area -13..-45 (descending, as ERA5)
    lon = np.arange(-70.0, -28.0 + 1e-9, 0.25)
    return Recorder(build_store(lat, lon, times_6h("2020-01-29 00:00", 8)))


# ------------------------------------------------------------------ request
def test_plan_as_cds_route(track):
    plan = a.plan_request(track, 6)
    assert plan["area"] == area_from_track(track) == (-13, -65, -45, -33)
    assert list(plan["times"]) == list(planned_request_times(track.index, 6))
    assert plan["latitudes"][0] == -13 and plan["latitudes"][-1] == -45 and np.all(np.diff(plan["latitudes"]) < 0)
    assert plan["longitudes"][0] == -65 and plan["longitudes"][-1] == -33 and np.all(np.diff(plan["longitudes"]) > 0)
    assert list(plan["levels"]) == [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50]


def test_plan_greenwich_is_one_arc_west_to_east():
    plan = a.plan_request(make_track("2020-03-01", "2020-03-02", lon=(-5, 5)), 6)
    assert (plan["area"][1], plan["area"][3]) == (-20, 20)
    assert plan["longitudes"][0] == -20 and plan["longitudes"][-1] == 20 and np.all(np.diff(plan["longitudes"]) > 0)


def test_plan_refuses_180_crossing():
    with pytest.raises(a.ArcoRequestError, match="crosses 180"):
        a.plan_request(make_track("2020-03-01", "2020-03-02", lon=(170, 178)), 6)


@pytest.mark.parametrize("step", [1, 3, 4, 5, 9])
def test_time_resolution_not_multiple_of_6_refused_before_any_read(step, track, tmp_path):
    class Boom(dict):
        def __getitem__(self, key):
            raise AssertionError("the store was read")

    with pytest.raises(a.ArcoRequestError, match="multiple of 6"):
        a.get_arco_data(args_for(tmp_path / "x.nc", time_resolution=step), track, LOGGER, mapping=Boom())
    assert not (tmp_path / "x.nc").exists()


def test_twelve_hours_is_served(tmp_path):
    pytest.importorskip("zarr")
    track = make_track("2020-01-30 00:00", "2020-01-31 00:00")
    lat = np.arange(-10.0, -50.0 - 1e-9, -0.25)
    lon = np.arange(-70.0, -28.0 + 1e-9, 0.25)
    mapping = build_store(lat, lon, times_6h("2020-01-29 00:00", 12))
    out = a.get_arco_data(args_for(tmp_path / "x.nc", 12), track, LOGGER, mapping=mapping)
    with xr.open_dataset(out) as ds:
        assert list(ds.valid_time.values) == [np.datetime64("2020-01-30T00:00"), np.datetime64("2020-01-30T12:00"),
                                              np.datetime64("2020-01-31T00:00")]


# ------------------------------------------------------------------ output contract and values
@pytest.mark.parametrize("convention, lon_signed", [
    ("signed", np.arange(-70.0, -28.0 + 1e-9, 0.25)),
    ("360", np.arange(-70.0, -28.0 + 1e-9, 0.25)),
])
def test_file_follows_the_cds_contract_with_the_values_of_the_store(convention, lon_signed, track, tmp_path):
    lat = np.arange(-10.0, -50.0 - 1e-9, -0.25)
    mapping = build_store(lat, lon_signed, times_6h("2020-01-29 00:00", 8), convention=convention)
    out = a.get_arco_data(args_for(tmp_path / "arco.nc"), track, LOGGER, mapping=mapping)
    plan = a.plan_request(track, 6)
    with xr.open_dataset(out) as ds:
        assert set(ds.data_vars) == set("uvtwz")
        for short in "uvtwz":
            v = ds[short]
            assert v.dims == ("valid_time", "pressure_level", "latitude", "longitude") and v.dtype == np.float32
            assert v.attrs["units"] == UNITS[short]
            expected = field(short, plan["times"], plan["levels"], plan["latitudes"], plan["longitudes"])
            np.testing.assert_array_equal(v.values, expected)  # values unchanged, whatever the store's order
        assert list(ds.pressure_level.values) == list(plan["levels"])  # 13 levels, 1000 first, nothing added
        np.testing.assert_array_equal(ds.latitude.values, plan["latitudes"])  # north to south
        np.testing.assert_array_equal(ds.longitude.values, plan["longitudes"])  # west to east, -180..180 as CDS
        assert ds.valid_time.encoding["units"] == "seconds since 1970-01-01"


def test_greenwich_with_a_store_in_0_360(tmp_path):
    track = make_track("2020-03-01 00:00", "2020-03-01 12:00", lat=(-10, -8), lon=(-5, 5))
    lat = np.arange(8.0, -26.0 - 1e-9, -0.25)
    lon = np.arange(-25.0, 25.0 + 1e-9, 0.25)  # stored as 335..359.75 and 0..25
    mapping = build_store(lat, lon, times_6h("2020-03-01 00:00", 4), convention="360")
    out = a.get_arco_data(args_for(tmp_path / "g.nc"), track, LOGGER, mapping=mapping)
    plan = a.plan_request(track, 6)
    with xr.open_dataset(out) as ds:
        np.testing.assert_array_equal(ds.longitude.values, plan["longitudes"])
        assert ds.longitude.values[0] == -20 and ds.longitude.values[-1] == 20 and 0.0 in ds.longitude.values
        np.testing.assert_array_equal(
            ds["t"].values, field("t", plan["times"], plan["levels"], plan["latitudes"], plan["longitudes"]))


def test_real_store_layout(track, tmp_path):
    """Layout observed in the real store: short names, pressureLevel, latitude ascending, one level per chunk."""
    lat = np.arange(-50.0, -10.0 + 1e-9, 0.25)  # ascending, as the real store (-90..90)
    lon = np.arange(-70.0, -28.0 + 1e-9, 0.25)
    mapping = build_store(lat, lon, times_6h("2020-01-29 00:00", 8), fields="short",
                          dims_order=("time", "pressureLevel", "latitude", "longitude"),
                          chunks={"pressureLevel": 1, "latitude": len(lat), "longitude": 100})
    out = a.get_arco_data(args_for(tmp_path / "real.nc"), track, LOGGER, mapping=mapping)
    plan = a.plan_request(track, 6)
    with xr.open_dataset(out) as ds:
        for short in "uvtwz":
            np.testing.assert_array_equal(
                ds[short].values, field(short, plan["times"], plan["levels"], plan["latitudes"], plan["longitudes"]))
        assert ds.latitude.values[0] == -13 and ds.latitude.values[-1] == -45  # north to south whatever the store
    prov = json.loads((tmp_path / "real.nc.provenance.json").read_text())
    # 3 times x 13 levels x 1 latitude band x the longitude chunks touched (169 values over chunks of 100), 5 fields
    assert prov["estimate_by_variable"]["t"] == 3 * 13 * 1 * 2 * (len(lat) * 100) * 4
    assert prov["bytes_received"] <= prov["estimate_bytes"]


def test_task_graph_does_not_grow_with_the_store(track):
    """One task per time step and level and field, whatever the length of the store (the real one has 126 936 steps)."""
    lat = np.arange(-10.0, -50.0 - 1e-9, -0.25)
    lon = np.arange(-70.0, -28.0 + 1e-9, 0.25)
    mapping = build_store(lat, lon, times_6h("2020-01-01 00:00", 120))
    ds = xr.open_zarr(a.CountingStore(mapping, a.Budget(10 ** 12)), consolidated=True, chunks=None)
    plan = a.plan_request(track, 6)
    result = a.build_output(ds, plan, a.select(ds, plan))
    graph = result["t"].data.__dask_graph__()
    assert result["t"].data.numblocks == (3, 13, 1, 1)
    assert len(graph) < 5 * 3 * 13 + 50  # the selected blocks (plus the indexing and casting layers), not 120 x 13


def test_store_with_another_dimension_order(track, tmp_path):
    lat = np.arange(-10.0, -50.0 - 1e-9, -0.25)
    lon = np.arange(-70.0, -28.0 + 1e-9, 0.25)
    mapping = build_store(lat, lon, times_6h("2020-01-29 00:00", 8),
                          dims_order=("time", "latitude", "longitude", "pressureLevels"))
    out = a.get_arco_data(args_for(tmp_path / "o.nc"), track, LOGGER, mapping=mapping)
    plan = a.plan_request(track, 6)
    with xr.open_dataset(out) as ds:
        assert ds["z"].dims == ("valid_time", "pressure_level", "latitude", "longitude")
        np.testing.assert_array_equal(
            ds["z"].values, field("z", plan["times"], plan["levels"], plan["latitudes"], plan["longitudes"]))


def test_namelist_of_the_cds_contract_is_enough(store, track, tmp_path):
    from src.utils.validation import validate_namelist_file, validate_variable_match
    from src.utils import preprocessing
    out = a.get_arco_data(args_for(tmp_path / "n.nc"), track, LOGGER, mapping=store)
    namelist = validate_namelist_file(preprocessing.ERA5_NAMELIST, LOGGER)
    with xr.open_dataset(out) as ds:
        validate_variable_match(namelist, ds, preprocessing.ERA5_NAMELIST, LOGGER)  # raises if a name is absent
    ns = argparse.Namespace(cdsapi=False, wekeo=False, arco=True, ncar=False)
    assert preprocessing.namelist_for_run(ns) == preprocessing.ERA5_NAMELIST


# ------------------------------------------------------------------ levels warning and provenance
def test_levels_warning_and_provenance(store, track, tmp_path, caplog, monkeypatch):
    monkeypatch.setenv("CDSAPI_KEY", "SECRET-KEY-VALUE")
    with caplog.at_level(logging.INFO, logger="test_arco"):
        out = a.get_arco_data(args_for(tmp_path / "p.nc"), track, LOGGER, mapping=store)
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any("13 levels" in w and "Ca -5 %" in w and "37 levels" in w for w in warnings)
    text = (tmp_path / "p.nc.provenance.json").read_text()
    prov = json.loads(text)
    assert prov["route"] == "arco" and prov["warning"] == a.LEVELS_WARNING and len(prov["levels_hpa"]) == 13
    assert prov["area_nwse"] == [-13, -65, -45, -33] and len(prov["times"]) == 3
    assert prov["lectk_commit"] is None or len(prov["lectk_commit"]) == 40
    assert set(prov["packages"]) == {"zarr", "numcodecs", "fsspec", "aiohttp"}
    assert len(prov["output"]["sha256"]) == 64 and prov["output"]["bytes"] == os.path.getsize(out)
    assert 0 < prov["bytes_received"] <= prov["estimate_bytes"]  # the estimate is an upper bound
    assert prov["bytes_received"] > prov["bytes_metadata"] > 0 and prov["reads"] > 0
    assert "SECRET-KEY-VALUE" not in text and "SECRET-KEY-VALUE" not in caplog.text


def test_estimate_counts_whole_chunks(store, track):
    import zarr  # noqa: F401
    counting = a.CountingStore(store, a.Budget(10 ** 12))
    ds = xr.open_zarr(counting, consolidated=True, chunks=None)
    selection = a.select(ds, a.plan_request(track, 6))
    total, by_variable = a.estimate_bytes(ds, selection)
    # 3 times x 1 chunk of 13 levels x (latitudes 129 -> chunks of 40) x all 169 longitudes, float32
    lat_values = ds["latitude"].values
    positions = selection["positions"]["latitude"]
    n_lat_chunks = len(np.unique(positions // 40))
    assert by_variable["temperature"] == 3 * n_lat_chunks * (13 * 40 * len(ds["longitude"])) * 4
    assert total == 5 * by_variable["temperature"] and len(lat_values) > positions.max()
    assert total > a.output_bytes(a.plan_request(track, 6))  # more than the values kept


# ------------------------------------------------------------------ volume protections
def test_estimate_above_limit_reads_no_chunk(store, track, tmp_path):
    out = tmp_path / "v.nc"
    with pytest.raises(a.ArcoVolumeError, match=r"--max-gib 0\.\d or more") as error:
        a.get_arco_data(args_for(out, max_gib=0.001), track, LOGGER, mapping=store)
    assert "no chunk was read" in str(error.value)
    chunk_reads = [k for k in store.reads if not os.path.basename(k).startswith(".z") and "/" in k
                   and not k.split("/")[0] in ("time", "latitude", "longitude", "pressureLevels")]
    assert chunk_reads == [] and not out.exists()


def test_max_gib_named_in_the_refusal_is_enough(store, track, tmp_path):
    with pytest.raises(a.ArcoVolumeError) as error:
        a.get_arco_data(args_for(tmp_path / "w.nc", max_gib=0.001), track, LOGGER, mapping=store)
    needed = float(str(error.value).split("--max-gib ")[1].split()[0])
    out = a.get_arco_data(args_for(tmp_path / "w.nc", max_gib=needed), track, LOGGER, mapping=store)
    assert os.path.exists(out)


def test_disk_space_checked_before_reading(store, track, tmp_path, monkeypatch):
    monkeypatch.setattr(a.shutil, "disk_usage", lambda d: types.SimpleNamespace(free=1000))
    with pytest.raises(a.ArcoVolumeError, match="of disk"):
        a.get_arco_data(args_for(tmp_path / "d.nc"), track, LOGGER, mapping=store)


def test_real_count_stops_at_the_limit():
    inner = Recorder({"a": b"x" * 600, "b": b"y" * 600, "c": b"z" * 600})
    store = a.CountingStore(inner, a.Budget(1000))
    assert len(store["a"]) == 600
    with pytest.raises(a.ArcoVolumeError, match="passed reading b"):
        store["b"]  # 1200 > 1000: counted, then stopped
    with pytest.raises(a.ArcoVolumeError, match="not requested"):
        store["c"]  # no further request once the limit is reached
    assert inner.reads == ["a", "b"] and store.budget.used == 1200


def test_failed_attempts_are_limited_and_not_counted_as_bytes():
    calls = []

    class Flaky(dict):
        def __getitem__(self, key):
            calls.append(key)
            if len(calls) < 3:
                raise OSError("connection reset")
            return b"ok"

    waits = []
    store = a.CountingStore(Flaky(), a.Budget(100), sleep=waits.append)
    assert store["k"] == b"ok" and len(calls) == 3 and store.budget.retries == 2 and waits == [2.0, 4.0]

    class Down(dict):
        def __getitem__(self, key):
            raise OSError("down")

    down = a.CountingStore(Down(), a.Budget(100), sleep=lambda s: None)
    with pytest.raises(a.ArcoTransferError, match="after 3 attempts"):
        down["k"]
    assert down.budget.retries == 3 and down.budget.used == 0


def test_refused_credential_is_not_retried():
    error = Exception("unauthorized")
    error.status = 401

    class Denied(dict):
        calls = 0

        def __getitem__(self, key):
            Denied.calls += 1
            raise error

    store = a.CountingStore(Denied(), a.Budget(100), sleep=lambda s: None)
    with pytest.raises(a.ArcoSetupError, match="refused the credential"):
        store["k"]
    assert Denied.calls == 1


def test_missing_chunk_is_not_filled_with_nan(store, track, tmp_path):
    del store["temperature/5.0.1.0"]  # time 2020-01-30 06:00, a latitude band the box touches
    with pytest.raises(a.ArcoContentError, match="absent from the store"):
        a.get_arco_data(args_for(tmp_path / "m.nc"), track, LOGGER, mapping=store)
    assert not (tmp_path / "m.nc").exists()


def test_time_outside_the_store_is_refused(track, tmp_path):
    lat = np.arange(-10.0, -50.0 - 1e-9, -0.25)
    lon = np.arange(-70.0, -28.0 + 1e-9, 0.25)
    mapping = build_store(lat, lon, times_6h("2019-01-01 00:00", 4))
    with pytest.raises(a.ArcoContentError, match="not in the store"):
        a.get_arco_data(args_for(tmp_path / "t.nc"), track, LOGGER, mapping=mapping)


# ------------------------------------------------------------------ credential and packages
def test_token_from_environment_and_file(tmp_path):
    assert a.read_token({"CDSAPI_KEY": " abc "}) == "abc"
    rc = tmp_path / "rc"
    rc.write_text("url: https://cds.climate.copernicus.eu/api\nkey: from-file\n")
    assert a.read_token({"CDSAPI_RC": str(rc)}) == "from-file"
    (tmp_path / ".cdsapirc").write_text("key: home-key\n")
    assert a.read_token({}, home=str(tmp_path)) == "home-key"


def test_missing_credential_says_how_to_fix_without_a_value(tmp_path):
    with pytest.raises(a.ArcoSetupError) as error:
        a.read_token({"CDSAPI_RC": str(tmp_path / "absent")})
    text = str(error.value)
    assert "~/.cdsapirc" in text and "CDSAPI_KEY" in text and "CDSAPI_RC" in text and a.CDS_PROFILE_URL in text


def test_missing_credential_stops_before_any_request(track, tmp_path, monkeypatch):
    pytest.importorskip("zarr")
    pytest.importorskip("aiohttp")
    monkeypatch.delenv("CDSAPI_KEY", raising=False)
    monkeypatch.setenv("CDSAPI_RC", str(tmp_path / "absent"))
    with pytest.raises(a.ArcoSetupError, match="credential"):
        a.get_arco_data(args_for(tmp_path / "c.nc"), track, LOGGER)


def test_missing_package(monkeypatch):
    real = a.importlib.import_module

    def fake_import(name, *args):
        if name == "zarr":
            raise ImportError("No module named 'zarr'")
        return real(name, *args)

    monkeypatch.setattr(a.importlib, "import_module", fake_import)
    with pytest.raises(a.ArcoSetupError, match=r"lorenzcycletoolkit\[arco\]"):
        a.load_packages()


def test_zarr_3_is_refused(monkeypatch):
    fakes = {"zarr": types.SimpleNamespace(__version__="3.0.1"), "fsspec": types.SimpleNamespace(),
             "aiohttp": types.SimpleNamespace()}
    real = a.importlib.import_module
    monkeypatch.setattr(a.importlib, "import_module", lambda name, *x: fakes[name] if name in fakes else real(name, *x))
    with pytest.raises(a.ArcoSetupError, match="zarr 2"):
        a.load_packages()


def test_bearer_header_only_in_the_mapper():
    seen = {}
    fsspec = types.SimpleNamespace(get_mapper=lambda url, **kw: seen.update(url=url, **kw) or {})
    aiohttp = types.SimpleNamespace(ClientTimeout=lambda total: ("timeout", total))
    a.open_remote_mapping("TOKEN", packages={"fsspec": fsspec, "aiohttp": aiohttp})
    assert seen["url"] == a.STORE_URL and seen["headers"] == {"Authorization": "Bearer TOKEN"}
    assert "TOKEN" not in repr(a.CountingStore({}, a.Budget(1)).__dict__)  # the counted store holds no header


def test_existing_file_is_reused_by_the_command(tmp_path, monkeypatch):
    """With an existing input file nothing is read (same rule as --cdsapi)."""
    from src.utils import preprocessing
    (tmp_path / "x.nc").write_bytes(b"")
    monkeypatch.setattr(a, "get_arco_data", lambda *x, **k: (_ for _ in ()).throw(AssertionError("read")))
    ns = argparse.Namespace(infile=str(tmp_path / "x.nc"), cdsapi=False, wekeo=False, arco=True, ncar=False)
    with pytest.raises(Exception):  # the empty file is not NetCDF: opened, not downloaded
        preprocessing.get_data(ns, LOGGER)


# ------------------------------------------------------------------ command line
LAUNCHER = os.path.join(os.path.dirname(__file__), "run_cli_offline.py")
SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lorenzcycletoolkit.py")


@pytest.mark.parametrize("extra, message", [
    (["--arco", "--cdsapi"], "alternative data sources"),
    (["--arco", "--wekeo"], "alternative data sources"),
    (["--arco", "--ncar"], "alternative data sources"),
    (["--arco", "--time-resolution", "3"], "not a multiple of 6"),
    (["--arco", "--time-resolution", "9"], "not a multiple of 6"),
    (["--delete-downloaded-input"], "--delete-downloaded-input needs --cdsapi, --ncar, --wekeo or --arco"),
])
def test_cli_refusals(tmp_path, extra, message):
    result = subprocess.run([sys.executable, LAUNCHER, SCRIPT, str(tmp_path / "x.nc"), "-t", *extra],
                            capture_output=True, text=True, timeout=300)
    assert result.returncode == 2 and message in result.stderr
    assert not (tmp_path / "x.nc").exists()


def test_cli_accepts_arco_with_its_options():
    import lorenzcycletoolkit
    args = lorenzcycletoolkit.create_arg_parser().parse_args(
        ["x.nc", "-t", "--arco", "--time-resolution", "6", "--max-gib", "5", "--delete-downloaded-input",
         "--trackfile", "inputs/track"])
    assert args.arco and args.max_gib == 5 and args.delete_downloaded_input and not args.cdsapi

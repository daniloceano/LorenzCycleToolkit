"""
The LEC reads the namelist with which the data were prepared (offline tests).

With ``--cdsapi`` or ``--wekeo`` the data are prepared with the preset
``inputs/namelist_ERA5-cdsapi``. The LEC used to read ``inputs/namelist`` all
the same and failed with ``KeyError: 'time'`` when the two differed. Here the
input file already exists (the routes reuse it, nothing is downloaded) and
``inputs/namelist`` is the ERA5 variant with other names (``time``, ``level``).
With a local input file nothing changes: ``inputs/namelist`` serves both.
"""

import datetime
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import lorenzcycletoolkit
from src.utils import preprocessing, wekeo_hda
from tests import ncss_emulator as emu

REPO_ROOT = Path(__file__).resolve().parent.parent
TRACK = "time;Lat;Lon\n2018-08-09-0000;-23.00;-45.50\n2018-08-09-0600;-22.75;-45.25\n" \
        "2018-08-09-1200;-22.50;-45.00\n2018-08-09-1800;-22.25;-44.75\n"
TIMES = [datetime.datetime(2018, 8, 9, h) for h in (0, 6, 12, 18)]
SHORT = {"t": "T", "z": "Z", "w": "W", "u": "U", "v": "V"}


def era5_file(path, names):
    """ERA5-like file; ``names`` maps the roles to the names of the file's variables and coordinates."""
    lat = np.arange(-16.0, -34.01, -0.5)
    lon = np.arange(-56.0, -34.0 + 1e-9, 0.5)
    dims = (names["time"], names["level"], names["lat"], names["lon"])
    data = {names[k]: (dims, emu.field_values(SHORT[k], TIMES, lat, lon)) for k in SHORT}
    ds = xr.Dataset(data, coords={names["time"]: pd.DatetimeIndex(TIMES), names["level"]: np.asarray(emu.LEVELS),
                                  names["lat"]: lat, names["lon"]: lon})
    ds.to_netcdf(path)


CDS_NAMES = dict(time="valid_time", level="pressure_level", lat="latitude", lon="longitude",
                 t="t", z="z", w="w", u="u", v="v")


@pytest.fixture
def workdir(lectk_workdir):
    """Working directory whose inputs/namelist differs from the ERA5-cdsapi preset."""
    work = lectk_workdir(namelist="namelist_ERA5")
    shutil.copyfile(REPO_ROOT / "inputs" / "namelist_ERA5-cdsapi", work / "inputs" / "namelist_ERA5-cdsapi")
    (work / "inputs" / "track").write_text(TRACK)
    (work / "inputs" / "box_limits").write_text("min_lon;-48\nmax_lon;-42\nmin_lat;-25.5\nmax_lat;-19.5\n")
    return work


def forbid_download(monkeypatch):
    """Any request to a data service fails the test."""
    def boom(*args, **kwargs):
        raise AssertionError("a download was requested")
    monkeypatch.setattr(preprocessing, "get_cdsapi_data", boom)
    monkeypatch.setattr(wekeo_hda, "get_wekeo_data", boom)


@pytest.mark.parametrize("route", ["--cdsapi", "--wekeo"])
def test_lec_uses_the_namelist_of_the_route(route, workdir, monkeypatch):
    forbid_download(monkeypatch)
    era5_file(workdir / "era5.nc", CDS_NAMES)  # reused, as the routes do for an existing file
    monkeypatch.setattr(sys, "argv", ["lorenzcycletoolkit.py", "era5.nc", "-r", "-t", route,
                                      "--trackfile", "inputs/track", "--time-resolution", "6"])
    assert lorenzcycletoolkit.main() == 0
    status = json.loads((workdir / "LEC_Results" / "era5_track" / "run_status.json").read_text())
    assert status["status"] == "success"


@pytest.mark.parametrize("route", ["--cdsapi", "--wekeo"])
def test_missing_route_namelist_stops_before_any_download(route, workdir, monkeypatch):
    forbid_download(monkeypatch)
    (workdir / "inputs" / "namelist_ERA5-cdsapi").unlink()
    monkeypatch.setattr(sys, "argv", ["lorenzcycletoolkit.py", "new_case.nc", "-r", "-t", route,
                                      "--trackfile", "inputs/track", "--time-resolution", "6"])
    with pytest.raises(FileNotFoundError, match="namelist_ERA5-cdsapi"):
        lorenzcycletoolkit.main()
    assert not (workdir / "new_case.nc").exists()


def test_local_file_still_uses_inputs_namelist(workdir, monkeypatch):
    forbid_download(monkeypatch)
    era5_file(workdir / "local.nc", dict(CDS_NAMES, time="time", level="level",
                                         t="T", z="Z", w="W", u="U", v="V"))
    monkeypatch.setattr(sys, "argv", ["lorenzcycletoolkit.py", "local.nc", "-r", "-t"])
    assert lorenzcycletoolkit.main() == 0


def test_namelist_for_run(workdir):
    parser = lorenzcycletoolkit.create_arg_parser()
    parse = lambda *extra: parser.parse_args(["x.nc", "-t", *extra])
    assert preprocessing.namelist_for_run(parse("--cdsapi")) == "inputs/namelist_ERA5-cdsapi"
    assert preprocessing.namelist_for_run(parse("--wekeo")) == "inputs/namelist_ERA5-cdsapi"
    assert preprocessing.namelist_for_run(parse("--ncar")) == "inputs/namelist"
    assert preprocessing.namelist_for_run(parse()) == "inputs/namelist"

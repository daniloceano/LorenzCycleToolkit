"""
Running the toolkit with no command-line arguments.

main() used to replace an empty command line with built-in example
arguments ('test.nc', -t, -r, -p, -v, --cdsapi). Running the program
without arguments therefore created a results directory, read inputs/track
and started an ERA5 download from the CDS API with the user's credentials.
It must instead print the usage and exit with status 2, doing nothing else.

No network access happens in these tests: the download function is
replaced by a stub that records the call and raises, and creating a CDS API
client fails the test.
"""

import hashlib
import shutil
import sys
from pathlib import Path

import cdsapi
import pytest

import lorenzcycletoolkit
import src.utils.preprocessing as preprocessing

REPO_ROOT = Path(__file__).resolve().parents[1]


class _DownloadReached(Exception):
    """Raised by the stub that replaces the CDS API download."""


def _snapshot(root):
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file() and not p.is_symlink() and "samples" not in p.parts
    }


@pytest.fixture
def user_like_workdir(lectk_workdir):
    """A working directory with the files a user keeps: namelists and a track."""
    workdir = lectk_workdir(namelist="namelist_ERA5", track="track_testdata_ERA5")
    shutil.copyfile(REPO_ROOT / "inputs" / "namelist_ERA5-cdsapi", workdir / "inputs" / "namelist_ERA5-cdsapi")
    return workdir


@pytest.fixture
def no_download(monkeypatch):
    calls = []

    def stub(args, track, app_logger):
        calls.append({"infile": args.infile, "cdsapi": args.cdsapi, "trackfile": args.trackfile})
        raise _DownloadReached(f"CDS API download reached with arguments {vars(args)}")

    def no_client(*args, **kwargs):
        raise AssertionError("a CDS API client must not be created")

    monkeypatch.setattr(preprocessing, "get_cdsapi_data", stub)
    monkeypatch.setattr(cdsapi, "Client", no_client)
    return calls


def test_no_arguments_prints_usage_and_does_nothing(user_like_workdir, no_download, monkeypatch, capsys):
    before = _snapshot(user_like_workdir)
    monkeypatch.setattr(sys, "argv", ["lorenzcycletoolkit.py"])

    with pytest.raises(SystemExit) as exit_info:
        lorenzcycletoolkit.main()

    assert exit_info.value.code == 2
    assert "usage:" in capsys.readouterr().err
    assert no_download == []
    assert not (user_like_workdir / "LEC_Results").exists()
    assert _snapshot(user_like_workdir) == before


def test_explicit_cdsapi_flag_is_still_accepted():
    args = lorenzcycletoolkit.create_arg_parser().parse_args(
        ["samples/new_case.nc", "-t", "-r", "--cdsapi", "--trackfile", "inputs/track"]
    )
    assert args.cdsapi is True
    assert args.track is True
    assert args.infile == "samples/new_case.nc"


def test_explicit_cdsapi_run_reaches_the_download_step(user_like_workdir, no_download, monkeypatch):
    """An explicit --cdsapi run still goes to the download (stubbed here)."""
    monkeypatch.setattr(sys, "argv", ["lorenzcycletoolkit.py", "new_case.nc", "-t", "-r", "--cdsapi"])
    with pytest.raises(_DownloadReached):
        lorenzcycletoolkit.main()
    assert no_download == [{"infile": "new_case.nc", "cdsapi": True, "trackfile": "inputs/track"}]

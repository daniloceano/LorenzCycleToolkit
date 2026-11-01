"""
D-Bench runner (benchmarks/d_bench/dbench.py), offline.

They check the pre-registered order, the private reading of the CDS key
(last line, one token, never written), the exclusive run directory with
the same namelist for both providers, the credential given only to the CDS
run, the outcome of a run (success, censored, failed), the stop rules of
the campaign (CDS stops at its first failure or when the account cannot be
checked; NCAR goes on), and the pre-registered comparison of two files.
The toolkit process and the CDS listing are replaced by fakes: nothing is
downloaded and no connection is opened.
"""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("dbench", REPO / "benchmarks" / "d_bench" / "dbench.py")
dbench = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dbench)

DUMMY = "00000000-1111-2222-3333-444444444444"
TRACK = "time;Lat;Lon\n1980-04-09-0000;-44.0;-68.0\n1980-04-09-0300;-45.0;-66.0\n1980-04-09-0600;-45.5;-65.0\n"


@pytest.fixture
def tracks(tmp_path, monkeypatch):
    directory = tmp_path / "tracks"
    directory.mkdir()
    for case in dbench.CASES:
        (directory / f"track_{case}.txt").write_text(TRACK)
    # the test track stands for every case: register its hash instead of the real ones
    monkeypatch.setattr(dbench, "TRACK_SHA256", {c: dbench.sha256(directory / f"track_{c}.txt") for c in dbench.CASES})
    return directory


@pytest.fixture
def key_file(tmp_path):
    path = tmp_path / "keys"
    path.write_text(f"old label 99999999-1111-2222-3333-444444444444\nnew label {DUMMY}\n")
    os.chmod(path, 0o600)
    return path


def _dataset(times, lon_offset=0.0, perturb=None):
    lat = np.arange(-30.0, -32.5, -0.25)
    lon = np.arange(-50.0, -47.5, 0.25) % 360.0 + lon_offset
    shape = (len(times), len(dbench.LEVELS_HPA), len(lat), len(lon))
    rng = np.random.default_rng(1)
    data = {}
    for var, units in (("t", "K"), ("z", "m**2 s**-2"), ("w", "Pa s**-1"), ("u", "m s**-1"), ("v", "m s**-1")):
        values = rng.normal(size=shape).astype("float32")
        if perturb and var in perturb:
            values = values + np.float32(perturb[var])
        data[var] = (("valid_time", "pressure_level", "latitude", "longitude"), values, {"units": units})
    return xr.Dataset(data, coords={"valid_time": times, "pressure_level": np.array(dbench.LEVELS_HPA, float),
                                    "latitude": lat, "longitude": lon})


class FakeProcess:
    """Stands for the toolkit run: writes the data file and run_status.json as the real CLI would."""

    calls = []

    def __init__(self, command, cwd, env, stdout, stderr, start_new_session):
        FakeProcess.calls.append({"command": command, "cwd": Path(cwd), "env": env})
        self.cwd = Path(cwd)
        self.returncode = FakeProcess.exit_code
        outname = command[max(i for i, c in enumerate(command) if c == "-o") + 1]  # not /usr/bin/time's -o
        infile = command[command.index("-t") - 1]
        _dataset(pd.date_range("1980-04-09", periods=3, freq="3h")).to_netcdf(self.cwd / infile)
        results = self.cwd / "LEC_Results" / f"{outname}_track"
        results.mkdir(parents=True)
        (results / "run_status.json").write_text(json.dumps(FakeProcess.status))

    def wait(self, timeout=None):
        return self.returncode


@pytest.fixture
def fake_run(monkeypatch):
    FakeProcess.calls = []
    FakeProcess.exit_code = 0
    FakeProcess.status = {"status": "success", "exit_code": 0, "error": None,
                          "timing": {"provider": "x", "input_data": "downloaded"}}
    monkeypatch.setattr(dbench, "_popen", FakeProcess)
    monkeypatch.setattr(dbench, "git_sha", lambda: "0" * 40)
    return FakeProcess


# ------------------------------------------------------------------- plan
def test_order_is_preregistered_and_complete():
    order = dbench.plan_order()
    assert order == dbench.plan_order()  # fixed seed
    assert order[0][0] == dbench.PILOT and sorted(c for c, _ in order) == sorted(dbench.CASES)
    assert all(sorted(p) == ["cds", "ncar"] for _, p in order)


def test_track_differing_from_the_registered_file_is_refused(tmp_path, tracks):
    (tracks / f"track_{dbench.PILOT}.txt").write_text(TRACK + "1980-04-09-0900;-46.0;-64.0\n")
    with pytest.raises(SystemExit, match="differs from the registered file"):
        dbench.track_file(tracks, dbench.PILOT)


def test_registered_hashes_are_the_ten_cases():
    assert sorted(dbench.TRACK_SHA256) == sorted(dbench.CASES)


def test_estimated_bytes_of_the_pilot_like_ncar_plan(tracks):
    track = dbench.read_track(dbench.track_file(tracks, dbench.PILOT))
    assert dbench.estimated_bytes(track) > 0


# -------------------------------------------------------------- credential
def test_key_is_the_last_line_token(key_file):
    assert dbench.cds_credential(key_file) == DUMMY


def test_key_file_readable_by_others_is_refused(key_file):
    os.chmod(key_file, 0o644)
    with pytest.raises(SystemExit, match="readable by other users"):
        dbench.cds_credential(key_file)


def test_key_line_without_exactly_one_token_is_refused(tmp_path):
    path = tmp_path / "keys"
    path.write_text("label only\n")
    os.chmod(path, 0o600)
    with pytest.raises(SystemExit, match="exactly one key"):
        dbench.cds_credential(path)


# --------------------------------------------------------------------- run
def test_ncar_run_record_and_isolated_directory(tmp_path, tracks, fake_run):
    record = dbench.run(dbench.PILOT, "ncar", 1, tmp_path / "root", tracks)
    call = fake_run.calls[0]
    assert "--ncar" in call["command"] and "--cdsapi" not in call["command"]
    assert not any(k.startswith("CDSAPI_") for k in call["env"])
    inputs = call["cwd"] / "inputs"
    assert (inputs / "namelist").read_bytes() == dbench.NAMELIST.read_bytes()
    assert (inputs / "namelist_ERA5-cdsapi").read_bytes() == dbench.NAMELIST.read_bytes()
    assert record["outcome"] == "success" and record["cold"] is True
    assert record["input_file"]["bytes"] > 0 and len(record["input_file"]["sha256"]) == 64
    saved = tmp_path / "root" / "records" / f"{dbench.PILOT}_ncar_r1.json"
    assert json.loads(saved.read_text())["outcome"] == "success"
    assert oct(os.stat(saved).st_mode & 0o777) == "0o600"


def test_run_directory_is_exclusive(tmp_path, tracks, fake_run):
    dbench.run(dbench.PILOT, "ncar", 1, tmp_path / "root", tracks)
    with pytest.raises(FileExistsError):
        dbench.run(dbench.PILOT, "ncar", 1, tmp_path / "root", tracks)
    assert len(fake_run.calls) == 1


def test_cds_run_gets_credential_only_in_its_environment(tmp_path, tracks, key_file, fake_run):
    record = dbench.run(dbench.PILOT, "cds", 1, tmp_path / "root", tracks, key_file=key_file)
    call = fake_run.calls[0]
    assert call["env"]["CDSAPI_KEY"] == DUMMY and call["env"]["CDSAPI_URL"] == dbench.CDS_URL
    command = call["command"]
    for flag in ("--cdsapi", "--cds-max-seconds", "--cds-ledger", "--cds-retries", "--cds-timeout", "--cds-case-id"):
        assert flag in command
    assert DUMMY not in " ".join(command)
    root = tmp_path / "root"
    for path in root.rglob("*"):
        if path.is_file() and path.suffix in (".json", ".jsonl", ".txt"):
            assert DUMMY not in path.read_text(errors="replace"), path
    assert record["outcome"] == "success"


def test_cds_run_needs_key_file(tmp_path, tracks, fake_run):
    with pytest.raises(SystemExit, match="--key-file"):
        dbench.run(dbench.PILOT, "cds", 1, tmp_path / "root", tracks)


@pytest.mark.parametrize("exit_code, status, error, expected", [
    (0, "success", None, "success"),
    (1, "failed", "CdsDeadlineExceeded", "censored"),
    (1, "failed", "NcarRetrievalError", "failed"),
    (0, "partial", None, "failed"),
])
def test_outcome(tmp_path, tracks, fake_run, exit_code, status, error, expected):
    fake_run.exit_code = exit_code
    fake_run.status = {"status": status, "exit_code": exit_code, "error": error and {"type": error}, "timing": {}}
    assert dbench.run(dbench.PILOT, "ncar", 1, tmp_path / "root", tracks)["outcome"] == expected


def test_not_enough_disk_is_refused(tmp_path, tracks, fake_run, monkeypatch):
    monkeypatch.setattr(dbench.shutil, "disk_usage", lambda path: type("U", (), {"free": 1})())
    with pytest.raises(SystemExit, match="not enough free space"):
        dbench.run(dbench.PILOT, "ncar", 1, tmp_path / "root", tracks)
    assert fake_run.calls == []


# ---------------------------------------------------------------- campaign
def _fake_campaign(monkeypatch, outcomes, pending=0):
    runs = []

    def fake_run(case, provider, rep, root, tracks, key_file=None, order=None):
        runs.append((case, provider))
        Path(dbench.run_directory(root, case, provider, rep)).mkdir(parents=True)
        return {"outcome": outcomes.get((case, provider), "success")}

    def fake_pending(key_file):
        if isinstance(pending, Exception):
            raise pending
        return pending

    monkeypatch.setattr(dbench, "run", fake_run)
    monkeypatch.setattr(dbench, "cds_pending", fake_pending)
    monkeypatch.setattr(dbench, "compare", lambda *a, **k: {"equivalent": True})
    monkeypatch.setattr(dbench, "pilot_gate", lambda *a, **k: [])
    monkeypatch.setattr(dbench, "campaign_elapsed", lambda root: 0.0)
    return runs


def test_campaign_runs_every_pair_in_order(tmp_path, tracks, monkeypatch):
    runs = _fake_campaign(monkeypatch, {})
    result = dbench.campaign(tmp_path, tracks, "keys")
    assert runs == [(c, p) for c, providers in dbench.plan_order() for p in providers]
    assert result["cds_stopped"] is None


def test_campaign_stops_cds_after_first_failure_and_keeps_ncar(tmp_path, tracks, monkeypatch):
    runs = _fake_campaign(monkeypatch, {(dbench.PILOT, "cds"): "censored"})
    result = dbench.campaign(tmp_path, tracks, "keys")
    assert [r for r in runs if r[1] == "cds"] == [(dbench.PILOT, "cds")]
    assert len([r for r in runs if r[1] == "ncar"]) == 10
    assert "censored" in result["cds_stopped"]


@pytest.mark.parametrize("pending, message", [(1, "still pending"), (RuntimeError("down"), "could not be listed")])
def test_campaign_submits_no_cds_when_account_not_clear(tmp_path, tracks, monkeypatch, pending, message):
    runs = _fake_campaign(monkeypatch, {}, pending=pending)
    result = dbench.campaign(tmp_path, tracks, "keys")
    assert not [r for r in runs if r[1] == "cds"]
    assert message in result["cds_stopped"]


def test_campaign_stops_ncar_after_two_consecutive_failures(tmp_path, tracks, monkeypatch):
    order = [c for c, _ in dbench.plan_order()]
    runs = _fake_campaign(monkeypatch, {(order[0], "ncar"): "failed", (order[1], "ncar"): "failed"})
    dbench.campaign(tmp_path, tracks, "keys")
    assert [r for r in runs if r[1] == "ncar"] == [(order[0], "ncar"), (order[1], "ncar")]


def test_campaign_never_repeats_a_done_run(tmp_path, tracks, monkeypatch):
    runs = _fake_campaign(monkeypatch, {})
    dbench.run_directory(tmp_path, dbench.PILOT, "cds", 1).mkdir(parents=True)
    dbench.campaign(tmp_path, tracks, "keys", cases=[dbench.PILOT])
    assert runs == [(dbench.PILOT, "ncar")]


# ----------------------------------------------------------------- compare
def _write_pair(root, cds, ncar, case=dbench.PILOT):
    for provider, ds in (("cds", cds), ("ncar", ncar)):
        directory = dbench.run_directory(root, case, provider, 1)
        directory.mkdir(parents=True)
        ds.to_netcdf(directory / f"era5_{case}_{provider}.nc")


TIMES = pd.date_range("1980-04-09", periods=3, freq="3h")


def test_identical_files_are_equivalent_except_missing_lec(tmp_path, tracks):
    # CDS-like longitudes in -180..180, NCAR-like in 0..360
    cds = _dataset(TIMES).assign_coords(longitude=lambda d: ((d.longitude + 180) % 360) - 180)
    _write_pair(tmp_path, cds, _dataset(TIMES))
    report = dbench.compare(dbench.PILOT, 1, tmp_path, tracks)
    assert report["problems"] == ["a provider has no LEC results; scientific equivalence not assessed"]
    assert all(v["max_abs_diff"] == 0 and v["fraction_identical"] == 1.0 for v in report["variables"].values())


def test_difference_beyond_tolerance_is_not_equivalent(tmp_path, tracks):
    _write_pair(tmp_path, _dataset(TIMES), _dataset(TIMES, perturb={"t": 0.05}))
    report = dbench.compare(dbench.PILOT, 1, tmp_path, tracks)
    assert report["equivalent"] is False
    assert any(p.startswith("t differs") for p in report["problems"])


def test_different_times_and_missing_track_times_are_reported(tmp_path, tracks):
    _write_pair(tmp_path, _dataset(TIMES), _dataset(TIMES[:2]))
    problems = dbench.compare(dbench.PILOT, 1, tmp_path, tracks)["problems"]
    assert "different time sets" in problems and "ncar misses track times" in problems


def test_lec_tables_are_compared(tmp_path, tracks):
    _write_pair(tmp_path, _dataset(TIMES), _dataset(TIMES))
    for provider, value in (("cds", 1.0), ("ncar", 1.0 + 1e-6)):
        results = (dbench.run_directory(tmp_path, dbench.PILOT, provider, 1) / "LEC_Results"
                   / f"{dbench.PILOT}_{provider}_r1_track")
        results.mkdir(parents=True)
        pd.DataFrame({"Ca": [value, 2.0], "Az": [10.0, 20.0]}, index=["a", "b"]).to_csv(
            results / f"{dbench.PILOT}_{provider}_r1_track_results.csv")
    report = dbench.compare(dbench.PILOT, 1, tmp_path, tracks)
    assert report["equivalent"] is True and report["lec_columns"]["Ca"]["max_abs_diff"] == pytest.approx(1e-6)


def test_missing_provider_file_is_no_pair(tmp_path, tracks):
    report = dbench.compare(dbench.PILOT, 1, tmp_path, tracks)
    assert report["equivalent"] is False and "no pair" in report["problems"][0]


# ----------------------------------------------------------------- summary
def test_median_ci_order_statistics():
    median, low, high, coverage = dbench.median_ci(list(range(1, 11)))
    assert (median, low, high) == (5.5, 2, 9) and coverage == pytest.approx(1 - 2 * 11 / 1024)
    assert dbench.median_ci([1, 2, 3, 4, 5])[1:] == (None, None, None)  # n=5: no 95% CI


def _records(root, ratios, outcome="success", equivalent=True):
    (root / "records").mkdir(parents=True, exist_ok=True)
    (root / "compare").mkdir(exist_ok=True)
    for case, ratio in zip(dbench.CASES, ratios):
        for provider, seconds in (("cds", 100.0), ("ncar", 100.0 * ratio)):
            record = {"outcome": outcome, "wall_seconds_external": seconds + 5, "input_file": {"bytes": 1},
                      "run_status": {"timing": {"total_seconds": seconds + 2, "phases": {
                          "acquisition": {"seconds": seconds}, "preparation": {"seconds": 1}, "lec": {"seconds": 1}}}}}
            (root / "records" / f"{case}_{provider}_r1.json").write_text(json.dumps(record))
        (root / "compare" / f"{case}_r1.json").write_text(json.dumps({"equivalent": equivalent}))


@pytest.mark.parametrize("ratios, verdict", [
    ([0.5] * 10, "advantage of NCAR"),
    ([2.0] * 10, "advantage of CDS"),
    ([0.7, 1.3] * 5, "no superiority demonstrated"),
    ([0.5] * 5, "insufficient valid pairs"),
])
def test_summary_verdicts(tmp_path, ratios, verdict):
    _records(tmp_path, ratios)
    result = dbench.summary(tmp_path)
    assert verdict in result["verdict"] and result["n_valid_pairs"] == len(ratios)


def test_summary_excludes_censored_and_non_equivalent_pairs(tmp_path):
    _records(tmp_path, [0.5] * 10, outcome="censored")
    assert dbench.summary(tmp_path)["n_valid_pairs"] == 0
    _records(tmp_path, [0.5] * 10, equivalent=False)
    assert dbench.summary(tmp_path)["n_valid_pairs"] == 0


# -------------------------------------------------------------- pilot gate
def _pilot_record(root, provider, outcome="success", rss=1000, acquisition=10.0):
    (root / "records").mkdir(parents=True, exist_ok=True)
    record = {"outcome": outcome, "peak_rss_kib": rss, "run_status": {"timing": {"phases": {
        "acquisition": {"seconds": acquisition}, "preparation": {"seconds": 1.0}, "lec": {"seconds": 2.0}}}}}
    (root / "records" / f"{dbench.PILOT}_{provider}_r1.json").write_text(json.dumps(record))


def test_pilot_gate_passes_and_fails(tmp_path):
    for provider in dbench.PROVIDERS:
        _pilot_record(tmp_path, provider)
    assert dbench.pilot_gate(tmp_path, 1, {"equivalent": True}) == []
    assert any("not equivalent" in r for r in dbench.pilot_gate(tmp_path, 1, {"equivalent": False, "problems": ["x"]}))
    _pilot_record(tmp_path, "cds", outcome="censored", acquisition=None)
    reasons = dbench.pilot_gate(tmp_path, 1, {"equivalent": True})
    assert "cds: censored" in reasons and "cds: a phase was not measured" in reasons
    _pilot_record(tmp_path, "cds", rss=10 ** 12)
    assert any("peak RSS" in r for r in dbench.pilot_gate(tmp_path, 1, {"equivalent": True}))


def test_campaign_stops_after_a_failed_pilot_gate(tmp_path, tracks, monkeypatch):
    runs = _fake_campaign(monkeypatch, {})
    monkeypatch.setattr(dbench, "pilot_gate", lambda *a, **k: ["cds: censored"])
    result = dbench.campaign(tmp_path, tracks, "keys")
    assert result["stopped"] == "pilot gate not passed" and {c for c, _ in runs} == {dbench.PILOT}


# --------------------------------------------------------- CSV tolerance
@pytest.mark.parametrize("cds, ncar, ok", [
    ([1e5, 2e5], [1e5 + 5, 2e5], True),             # 5 <= 1e-4 x 2e5 = 20
    ([1e5, 2e5], [1e5 + 50, 2e5], False),
    ([-3.0, 1.0], [-3.0 + 2e-4, 1.0], True),         # scale uses |values|: 1e-4 x 3
    ([-3.0, 1.0], [-3.0 + 4e-4, 1.0], False),
    ([0.0, 0.0], [0.0, 5e-10], True),                # near zero: absolute floor 1e-9
    ([0.0, 0.0], [0.0, 5e-9], False),
    ([1e-12, 0.0], [0.0, 1e-12], True),
    ([1.0, float("nan")], [1.0, float("nan")], True),
    ([1.0, float("nan")], [1.0, 1.0], False),        # missing at different times
    ([1.0, float("inf")], [1.0, float("inf")], False),
])
def test_csv_column_check(cds, ncar, ok):
    assert dbench.csv_column_check(cds, ncar)["ok"] is ok


def test_csv_scale_is_symmetric():
    one = dbench.csv_column_check([1.0, 0.0], [10.0, 0.0])
    other = dbench.csv_column_check([10.0, 0.0], [1.0, 0.0])
    assert one["scale"] == other["scale"] == 10.0 and one["tolerance"] == other["tolerance"]


# --------------------------------------------------------- campaign limits
def test_campaign_ceiling_counts_from_the_first_start(tmp_path, monkeypatch):
    clock = {"now": 1000.0}
    monkeypatch.setattr(dbench.time, "time", lambda: clock["now"])
    assert dbench.campaign_elapsed(tmp_path) == 0.0
    clock["now"] += 3600
    assert dbench.campaign_elapsed(tmp_path) == 3600.0  # a resumed campaign keeps its start


def test_no_run_starts_without_room_for_its_deadline(tmp_path, tracks, monkeypatch):
    runs = _fake_campaign(monkeypatch, {})
    room = dbench.CAMPAIGN_WALL_S - dbench.CASE_DEADLINE_S - dbench.KILL_MARGIN_S
    monkeypatch.setattr(dbench, "campaign_elapsed", lambda root: room + 1)
    result = dbench.campaign(tmp_path, tracks, "keys")
    assert runs == [] and result["partial"] is True and "ceiling" in result["stopped"]


def test_r2_only_after_ten_successful_r1_pairs(tmp_path, tracks, monkeypatch):
    runs = _fake_campaign(monkeypatch, {})
    result = dbench.campaign(tmp_path, tracks, "keys", rep=2)
    assert runs == [] and "r2 refused" in result["stopped"]
    (tmp_path / "records").mkdir(exist_ok=True)
    for case in dbench.CASES:
        for provider in dbench.PROVIDERS:
            (tmp_path / "records" / f"{case}_{provider}_r1.json").write_text(json.dumps({"outcome": "success"}))
    assert dbench.r1_complete(tmp_path)
    dbench.campaign(tmp_path, tracks, "keys", rep=2)
    assert len(runs) == 20


def test_campaign_with_cds_stopped_is_partial(tmp_path, tracks, monkeypatch):
    _fake_campaign(monkeypatch, {(dbench.PILOT, "cds"): "failed"})
    assert dbench.campaign(tmp_path, tracks, "keys")["partial"] is True

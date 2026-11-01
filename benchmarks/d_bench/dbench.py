"""
D-Bench runner: controlled CDS x NCAR acquisition benchmark of LorenzCycleToolkit.

Protocol: docs/evolution/evidence/records/2026-10-29-dbench-preregistration.md
(criteria, tolerances, order, budget and stop rules fixed before any result).

Every run is one cold acquisition of one track by one provider, followed by
the LEC computation (``--track -r``), in a new private working directory
``ROOT/runs/<case>/<provider>/r<rep>/`` that must not exist yet: no file of
an earlier run can be taken for a new download. The working directory holds
its own ``inputs/namelist`` and ``inputs/namelist_ERA5-cdsapi`` (both the
versioned CDS preset), so that the namelist read for the data, by --ncar for
the output names and by the LEC computation is the same for both providers.

Nothing secret is written: the CDS credential is read from a private file,
given to the CDS run through its environment only (CDSAPI_URL/CDSAPI_KEY),
and never printed, logged or recorded. Records, logs, data and the CDS
ledger stay under ROOT, outside Git.

Commands (run with the Python environment of the toolkit, from anywhere):

    dbench.py plan
    dbench.py cds-pending --key-file F                   # read-only job listing
    dbench.py run --case ID --provider {ncar,cds} --rep N --root ROOT --tracks DIR [--key-file F]
    dbench.py compare --case ID --rep N --root ROOT --tracks DIR
    dbench.py campaign --root ROOT --tracks DIR --key-file F [--cases ...]
"""

import argparse
import datetime
import hashlib
import json
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

CASES = ("19790166", "19790205", "19790535", "19790585", "19790644",
         "19791008", "19791150", "19800256", "19800296", "19800318")
PILOT = "19800296"
# SHA-256 of the ten supplied track files (docs/evolution/fronts/front-D-bench-cds-ncar.md)
TRACK_SHA256 = {
    "19790166": "476d0197842408c36ae6e785b913e2ae3e22a74bb25b6910409968eaf0bef992",
    "19790205": "9da9d8092731606c58c31709469b9bb6e1ae6e3ee2750ee000cd3032a5f4f0c9",
    "19790535": "7649033bf5a22aaa4d1c57833ec69e72c603cc2b5d19d43b941862176f59a2fc",
    "19790585": "d6b6991da9c85fb406422f9c5f8b91c892fe1027641d7c53e0b26f0a1332a615",
    "19790644": "f2a1c53ec29bd66290c1e443c60bcaec45da51ceda69602801716f3a0a0861f9",
    "19791008": "22b8055e96ce7491d86771a45532a57e9dbeb1dce9c485ae6fa566d7401376ff",
    "19791150": "6134b532561fb91df35b75ede098f399340927ebcae0c8b5678ff4ac854f3f92",
    "19800256": "ca838ff9ba98faf5c61d8cc31ec6509fc25dd1ad4fd9bb819e61208322146c36",
    "19800296": "d3ae66afff52ada88b54f9861b702451c291df1366fbbbf817f132fcd342b857",
    "19800318": "616dec4389968e3ef70751be2b3158b81761ff22528a365048aec3770c2cc3b5",
}
SEED = 20261029
PROVIDERS = ("ncar", "cds")
CDS_URL = "https://cds.climate.copernicus.eu/api"
NAMELIST = REPO / "inputs" / "namelist_ERA5-cdsapi"
TIME_RESOLUTION = 3

# Pre-registered limits (see the pre-registration record)
CASE_DEADLINE_S = 12 * 3600          # per case and provider: acquisition deadline (censoring)
KILL_MARGIN_S = 30 * 60               # external kill only if the run overruns the deadline by this much
CDS_RETRIES = 3
CDS_TIMEOUT_S = 120
CDS_DELETE_ON_STOP = True             # best-effort delete of our own job after a deadline (subject to approval)
CAMPAIGN_WALL_S = 96 * 3600
MIN_FREE_FACTOR = 3                   # free space >= 3 x estimated file + 5 GiB before each run
MIN_FREE_EXTRA = 5 * 1024 ** 3
NCAR_MAX_CONSECUTIVE_FAILURES = 2
PILOT_MAX_RSS_FRACTION = 0.25         # pilot gate: peak RSS below 25% of the RAM
TIME_V = "/usr/bin/time"              # peak RSS and CPU of the run (null when absent)
_popen = subprocess.Popen             # replaced by a fake in the tests

# Pre-registered integrity tolerances (absolute, per variable, in the file units)
TOLERANCES = {"t": 0.005, "z": 2.0, "w": 0.001, "u": 0.01, "v": 0.01}
LEVELS_HPA = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450,
              500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
# LEC CSV, per column: |delta| <= max(1e-4 x scale, 1e-9), scale = max(max|CDS|, max|NCAR|) >= 0;
# the absolute floor (column units) is the explicit rule for columns at or near zero
CSV_RELATIVE_TOLERANCE = 1e-4
CSV_ABSOLUTE_FLOOR = 1e-9

# Pre-registered decision rule on the per-case ratio NCAR/CDS of the acquisition time (rep 1)
MIN_PAIRS = 6                         # below this no 95% distribution-free CI of the median exists
NCAR_ADVANTAGE_MEDIAN = 0.80          # NCAR: median ratio <= 0.80 and CI upper bound < 1
CDS_ADVANTAGE_MEDIAN = 1.25           # CDS: median ratio >= 1.25 and CI lower bound > 1


# ------------------------------------------------------------------ helpers
def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_sha():
    try:
        return subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def private_dir(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    os.chmod(path, 0o700)
    return path


def write_json(path, content):
    tmp = Path(f"{path}.tmp")
    with open(os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600), "w", encoding="utf-8") as handle:
        json.dump(content, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(tmp, path)


def append_jsonl(path, entry):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
    try:
        os.write(fd, (json.dumps(entry, sort_keys=True) + "\n").encode("utf-8"))
        os.fsync(fd)
    finally:
        os.close(fd)


def track_file(tracks, case):
    path = Path(tracks) / f"track_{case}.txt"
    if not path.is_file():
        raise SystemExit(f"track of case {case} not found in the tracks directory")
    if sha256(path) != TRACK_SHA256[case]:
        raise SystemExit(f"track of case {case} differs from the registered file (SHA-256)")
    return path


def read_track(path):
    import pandas as pd

    return pd.read_csv(path, sep=";", index_col="time", parse_dates=["time"], date_format="%Y-%m-%d-%H%M")


def estimated_bytes(track):
    """float32 content of the five fields on the 37 levels (the NCAR plan of the same request)."""
    from src.utils.ncar_gdex import area_from_track, daily_hours, plan_retrieval

    return plan_retrieval(area_from_track(track), daily_hours(track.index, TIME_RESOLUTION)).output_bytes


def cds_credential(key_file):
    """The token of the last non-empty line of the private key file (never printed)."""
    mode = os.stat(key_file).st_mode
    if mode & 0o077:
        raise SystemExit("the CDS key file is readable by other users; refusing to use it")
    lines = [line for line in Path(key_file).read_text(encoding="utf-8").splitlines() if line.strip()]
    tokens = re.findall(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", lines[-1]) if lines else []
    if len(tokens) != 1:
        raise SystemExit("the last line of the CDS key file does not hold exactly one key")
    return tokens[0]


def plan_order():
    """Pre-registered order: the pilot first, then the other nine cases shuffled; provider order per case."""
    rng = random.Random(SEED)
    others = [c for c in CASES if c != PILOT]
    rng.shuffle(others)
    return [(case, tuple(rng.sample(PROVIDERS, 2))) for case in [PILOT, *others]]


# ------------------------------------------------------------- CDS pending
def cds_pending(key_file):
    """
    Number of the account's requests still accepted or running (read-only).

    Returns:
        int, or raises when the listing is not possible (then nothing may be submitted).
    """
    from ecmwf.datastores import Client

    client = Client(url=CDS_URL, key=cds_credential(key_file), progress=False, maximum_tries=2, timeout=60)
    jobs = client.get_jobs(status=["accepted", "running"], limit=100)
    return len(jobs.request_ids)


# --------------------------------------------------------------------- run
def parse_time_v(path):
    """Peak RSS (KiB) and CPU seconds from the report of /usr/bin/time -v."""
    out = {"peak_rss_kib": None, "user_cpu_s": None, "system_cpu_s": None}
    try:
        text = Path(path).read_text(errors="replace")
    except OSError:
        return out
    for key, pattern in (("peak_rss_kib", r"Maximum resident set size \(kbytes\): (\d+)"),
                         ("user_cpu_s", r"User time \(seconds\): ([\d.]+)"),
                         ("system_cpu_s", r"System time \(seconds\): ([\d.]+)")):
        match = re.search(pattern, text)
        if match:
            out[key] = float(match.group(1)) if "." in match.group(1) else int(match.group(1))
    return out


def run_directory(root, case, provider, rep):
    return Path(root) / "runs" / case / provider / f"r{rep}"


def run(case, provider, rep, root, tracks, key_file=None, order=None, deadline=CASE_DEADLINE_S):
    """One cold acquisition + LEC computation; returns its record (also written under ROOT/records)."""
    root = private_dir(root)
    track_path = track_file(tracks, case)
    track = read_track(track_path)
    estimate = estimated_bytes(track)
    free = shutil.disk_usage(root).free
    if free < MIN_FREE_FACTOR * estimate + MIN_FREE_EXTRA:
        raise SystemExit(f"not enough free space under ROOT for case {case}: {free / 1024 ** 3:.1f} GiB free")

    workdir = run_directory(root, case, provider, rep)
    workdir.parent.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(mode=0o700)  # exclusive destination: fails if the run directory exists
    inputs = workdir / "inputs"
    inputs.mkdir()
    for name in ("namelist", "namelist_ERA5-cdsapi"):
        shutil.copyfile(NAMELIST, inputs / name)
    shutil.copyfile(track_path, inputs / "track")
    tmpdir = private_dir(workdir / "tmp")

    infile = f"era5_{case}_{provider}.nc"
    outname = f"{case}_{provider}_r{rep}"
    command = [sys.executable, str(REPO / "lorenzcycletoolkit.py"), infile, "-t", "-r",
               "--trackfile", "inputs/track", "--time-resolution", str(TIME_RESOLUTION), "-o", outname]
    env = {k: v for k, v in os.environ.items() if not k.startswith("CDSAPI_")}
    env.update(MPLBACKEND="Agg", TMPDIR=str(tmpdir))
    ledger = private_dir(root / "private") / "cds_ledger.jsonl"
    if provider == "cds":
        if key_file is None:
            raise SystemExit("a CDS run needs --key-file")
        command += ["--cdsapi", "--cds-max-seconds", str(deadline), "--cds-ledger", str(ledger),
                    "--cds-retries", str(CDS_RETRIES), "--cds-timeout", str(CDS_TIMEOUT_S),
                    "--cds-case-id", f"{case}-cds-r{rep}"]
        if CDS_DELETE_ON_STOP:
            command.append("--cds-delete-on-stop")
        env.update(CDSAPI_URL=CDS_URL, CDSAPI_KEY=cds_credential(key_file))
    else:
        command += ["--ncar"]
    time_report = workdir / "time_v.txt"
    if os.path.exists(TIME_V):
        command = [TIME_V, "-v", "-o", str(time_report), *command]

    record = {
        "case": case, "provider": provider, "rep": rep, "order": order, "cold": True,
        "lectk_commit": git_sha(), "python": sys.version.split()[0], "command": [
            c if not c.startswith(str(root)) else c.replace(str(root), "ROOT") for c in command],
        "track_sha256": sha256(track_path), "track_positions": int(len(track)),
        "namelist_sha256": sha256(NAMELIST), "estimated_float32_bytes": int(estimate),
        "free_bytes_before": int(free), "load_avg_before": os.getloadavg(), "started_utc": utc(),
    }
    started = time.perf_counter()
    timed_out = False
    with open(workdir / "stdout.txt", "w") as out, open(workdir / "stderr.txt", "w") as err:
        # own process group, so that the toolkit itself is stopped, not only /usr/bin/time
        process = _popen(command, cwd=workdir, env=env, stdout=out, stderr=err, start_new_session=True)
        try:
            process.wait(timeout=deadline + KILL_MARGIN_S)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=60)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
    record["wall_seconds_external"] = round(time.perf_counter() - started, 3)
    record.update(finished_utc=utc(), exit_code=process.returncode, killed_after_deadline=timed_out,
                  load_avg_after=os.getloadavg(), **parse_time_v(time_report))

    status_path = workdir / "LEC_Results" / f"{outname}_track" / "run_status.json"
    status = json.loads(status_path.read_text()) if status_path.is_file() else None
    record["run_status"] = None if status is None else {
        "status": status.get("status"), "exit_code": status.get("exit_code"), "error": status.get("error"),
        "timing": status.get("timing")}
    data = workdir / infile
    record["input_file"] = {"bytes": data.stat().st_size, "sha256": sha256(data)} if data.is_file() else None
    provenance = workdir / f"{infile}.provenance.json"
    if provenance.is_file():
        totals = json.loads(provenance.read_text()).get("totals", {})
        record["ncar_totals"] = totals
    record["outcome"] = outcome(record)
    shutil.rmtree(tmpdir, ignore_errors=True)

    records = private_dir(root / "records")
    write_json(records / f"{case}_{provider}_r{rep}.json", record)
    append_jsonl(root / "records" / "campaign.jsonl", {k: record[k] for k in (
        "case", "provider", "rep", "order", "started_utc", "finished_utc", "exit_code", "outcome",
        "wall_seconds_external")})
    return record


def outcome(record):
    """success, censored (deadline) or failed; never a success without the complete scientific result."""
    status = record.get("run_status") or {}
    error = (status.get("error") or {}).get("type")
    if record.get("killed_after_deadline") or error == "CdsDeadlineExceeded":
        return "censored"
    if record.get("exit_code") == 0 and status.get("status") == "success":
        return "success"
    return "failed"


# ----------------------------------------------------------------- compare
def compare(case, rep, root, tracks):
    """Pre-registered integrity and equivalence checks of the CDS and NCAR files and LEC results."""
    import numpy as np
    import pandas as pd
    import xarray as xr

    root = Path(root)
    files = {p: run_directory(root, case, p, rep) / f"era5_{case}_{p}.nc" for p in PROVIDERS}
    report = {"case": case, "rep": rep, "checked_utc": utc(), "tolerances": TOLERANCES,
              "csv_relative_tolerance": CSV_RELATIVE_TOLERANCE, "csv_absolute_floor": CSV_ABSOLUTE_FLOOR,
              "problems": []}
    problems = report["problems"]
    if not all(f.is_file() for f in files.values()):
        problems.append("a provider has no data file; no pair")
        report["equivalent"] = False
        return _save_compare(root, report)

    track_times = set(pd.DatetimeIndex(read_track(track_file(tracks, case)).index))
    ds = {p: xr.open_dataset(f) for p, f in files.items()}
    try:
        times = {p: pd.DatetimeIndex(d["valid_time"].values) for p, d in ds.items()}
        report["n_times"] = {p: len(t) for p, t in times.items()}
        if set(times["cds"]) != set(times["ncar"]):
            problems.append("different time sets")
        for p, t in times.items():
            if not track_times <= set(t):
                problems.append(f"{p} misses track times")
        for p, d in ds.items():
            if [float(x) for x in d["pressure_level"].values] != [float(x) for x in LEVELS_HPA] and \
                    sorted(float(x) for x in d["pressure_level"].values) != [float(x) for x in LEVELS_HPA]:
                problems.append(f"{p} levels are not the 37 levels")
        lat = {p: np.sort(d["latitude"].values) for p, d in ds.items()}
        lon = {p: np.sort(d["longitude"].values % 360.0) for p, d in ds.items()}
        if lat["cds"].shape != lat["ncar"].shape or not np.allclose(lat["cds"], lat["ncar"], atol=1e-6):
            problems.append("different latitudes")
        if lon["cds"].shape != lon["ncar"].shape or not np.allclose(lon["cds"], lon["ncar"], atol=1e-6):
            problems.append("different longitudes")
        report["grid"] = {p: [len(lat[p]), len(lon[p])] for p in PROVIDERS}

        common = sorted(set(times["cds"]) & set(times["ncar"]))
        report["variables"] = {}
        for var, tol in TOLERANCES.items():
            if any(var not in d for d in ds.values()):
                problems.append(f"variable {var} missing")
                continue
            units = {p: d[var].attrs.get("units") for p, d in ds.items()}
            if units["cds"] != units["ncar"]:
                problems.append(f"{var} units differ: {units}")
            arrays = {}
            for p, d in ds.items():
                field = d[var].assign_coords(longitude=d["longitude"] % 360.0)
                field = field.sel(valid_time=common).sortby(["valid_time", "pressure_level", "latitude", "longitude"])
                arrays[p] = field.transpose("valid_time", "pressure_level", "latitude", "longitude").values
            entry = {"units": units, "dtype": {p: str(a.dtype) for p, a in arrays.items()}}
            if arrays["cds"].shape != arrays["ncar"].shape:
                problems.append(f"{var} shapes differ")
                report["variables"][var] = entry
                continue
            finite = {p: bool(np.isfinite(a).all()) for p, a in arrays.items()}
            diff = np.abs(arrays["cds"].astype("float64") - arrays["ncar"].astype("float64"))
            entry.update(finite=finite, max_abs_diff=float(np.nanmax(diff)), mean_abs_diff=float(np.nanmean(diff)),
                         fraction_identical=float((diff == 0).mean()), tolerance=tol)
            if not all(finite.values()):
                problems.append(f"{var} has non-finite values")
            if entry["max_abs_diff"] > tol:
                problems.append(f"{var} differs by {entry['max_abs_diff']:.3g} > {tol}")
            report["variables"][var] = entry
    finally:
        for d in ds.values():
            d.close()

    csvs = {p: run_directory(root, case, p, rep) / "LEC_Results" / f"{case}_{p}_r{rep}_track"
            / f"{case}_{p}_r{rep}_track_results.csv" for p in PROVIDERS}
    if all(c.is_file() for c in csvs.values()):
        a, b = (pd.read_csv(csvs[p], index_col=0) for p in PROVIDERS)
        if list(a.columns) != list(b.columns) or len(a) != len(b) or list(a.index) != list(b.index):
            problems.append("LEC result tables differ in shape, columns or times")
        else:
            report["lec_columns"] = {}
            for column in a.columns:
                entry = csv_column_check(a[column].values, b[column].values)
                report["lec_columns"][column] = entry
                if not entry["ok"]:
                    problems.append(f"LEC column {column}: {entry['why']}")
    else:
        problems.append("a provider has no LEC results; scientific equivalence not assessed")
    report["equivalent"] = not problems
    return _save_compare(root, report)


def csv_column_check(cds, ncar):
    """Pre-registered check of one LEC column (see CSV_RELATIVE_TOLERANCE and CSV_ABSOLUTE_FLOOR)."""
    import numpy as np

    cds, ncar = np.asarray(cds, dtype="float64"), np.asarray(ncar, dtype="float64")
    nan_cds, nan_ncar = np.isnan(cds), np.isnan(ncar)
    entry = {"ok": True, "why": None, "max_abs_diff": None, "scale": None, "tolerance": None,
             "nan_count": int(nan_cds.sum())}
    if not np.array_equal(nan_cds, nan_ncar):
        entry.update(ok=False, why="missing values at different times")
        return entry
    if np.isinf(cds).any() or np.isinf(ncar).any():
        entry.update(ok=False, why="infinite values")
        return entry
    valid = ~nan_cds
    if not valid.any():
        return entry  # all missing in both, at the same times
    scale = float(max(np.abs(cds[valid]).max(), np.abs(ncar[valid]).max()))
    tolerance = max(CSV_RELATIVE_TOLERANCE * scale, CSV_ABSOLUTE_FLOOR)
    delta = float(np.abs(cds[valid] - ncar[valid]).max())
    entry.update(max_abs_diff=delta, scale=scale, tolerance=tolerance)
    if delta > tolerance:
        entry.update(ok=False, why=f"differs by {delta:.3g} > {tolerance:.3g} (scale {scale:.3g})")
    return entry


def _save_compare(root, report):
    out = private_dir(Path(root) / "compare")
    write_json(out / f"{report['case']}_r{report['rep']}.json", report)
    return report


# ----------------------------------------------------------------- summary
def median_ci(values, level=0.95):
    """
    Median and distribution-free CI of the median from order statistics
    (binomial, as for the sign test); (median, low, high, coverage), with
    low/high None when the sample is too small for the level.
    """
    from math import comb

    x = sorted(values)
    n = len(x)
    if n == 0:
        return None, None, None, None
    median = (x[(n - 1) // 2] + x[n // 2]) / 2
    cdf = [sum(comb(n, i) for i in range(k + 1)) / 2 ** n for k in range(n + 1)]
    best = None
    for k in range(1, n // 2 + 1):  # CI = (x_(k), x_(n-k+1)), coverage 1 - 2 P(B <= k-1)
        coverage = 1 - 2 * cdf[k - 1]
        if coverage >= level:
            best = (x[k - 1], x[n - k], coverage)
    if best is None:
        return median, None, None, None
    return median, *best


def summary(root, rep=1):
    """Pre-registered comparison of the records (written to ROOT/summary/summary_r<rep>.json)."""
    root = Path(root)
    rows, ratios = [], []
    for case in CASES:
        row = {"case": case}
        for provider in PROVIDERS:
            path = root / "records" / f"{case}_{provider}_r{rep}.json"
            record = json.loads(path.read_text()) if path.is_file() else None
            phases = (((record or {}).get("run_status") or {}).get("timing") or {}).get("phases") or {}
            row[provider] = None if record is None else {
                "outcome": record["outcome"],
                "acquisition_s": (phases.get("acquisition") or {}).get("seconds"),
                "preparation_s": (phases.get("preparation") or {}).get("seconds"),
                "lec_s": (phases.get("lec") or {}).get("seconds"),
                "total_s": (((record.get("run_status") or {}).get("timing")) or {}).get("total_seconds"),
                "external_s": record.get("wall_seconds_external"),
                "bytes": (record.get("input_file") or {}).get("bytes"),
                "peak_rss_kib": record.get("peak_rss_kib"),
            }
        compare_path = root / "compare" / f"{case}_r{rep}.json"
        row["equivalent"] = json.loads(compare_path.read_text())["equivalent"] if compare_path.is_file() else None
        a, b = row["ncar"], row["cds"]
        valid = bool(a and b and a["outcome"] == b["outcome"] == "success" and row["equivalent"]
                     and a["acquisition_s"] and b["acquisition_s"])
        row["valid_pair"] = valid
        if valid:
            row["ratio_ncar_cds"] = a["acquisition_s"] / b["acquisition_s"]
            ratios.append(row["ratio_ncar_cds"])
        rows.append(row)

    median, low, high, coverage = median_ci(ratios)
    if len(ratios) < MIN_PAIRS or low is None:
        verdict = f"insufficient valid pairs ({len(ratios)} < {MIN_PAIRS}): no comparison of speed"
    elif median <= NCAR_ADVANTAGE_MEDIAN and high < 1:
        verdict = "practical advantage of NCAR under the observed conditions"
    elif median >= CDS_ADVANTAGE_MEDIAN and low > 1:
        verdict = "practical advantage of CDS under the observed conditions"
    else:
        verdict = "no superiority demonstrated (not a proof of equal speed)"
    result = {"rep": rep, "created_utc": utc(), "rows": rows, "n_valid_pairs": len(ratios),
              "complete": len(ratios) == len(CASES),
              "median_ratio_ncar_cds": median, "ci": [low, high], "ci_coverage": coverage, "verdict": verdict,
              "ci_note": ("distribution-free interval of the median from binomial order statistics; its exact "
                          "coverage is reported (97.9% for 10 pairs) and is at least 95%; it describes the median "
                          "ratio of these tracks under the observed service conditions and is not a significance "
                          "test; with few pairs it is wide"),
              "rule": {"min_pairs": MIN_PAIRS, "ncar_median_max": NCAR_ADVANTAGE_MEDIAN,
                       "cds_median_min": CDS_ADVANTAGE_MEDIAN}}
    out = private_dir(root / "summary")
    write_json(out / f"summary_r{rep}.json", result)
    return result


# ---------------------------------------------------------------- campaign
def pilot_gate(root, rep=1, compare_report=None):
    """Reasons why the pilot does not open the campaign (empty when it does)."""
    root = Path(root)
    ram_kib = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1024
    reasons = []
    for provider in PROVIDERS:
        path = root / "records" / f"{PILOT}_{provider}_r{rep}.json"
        if not path.is_file():
            reasons.append(f"{provider}: no pilot run")
            continue
        record = json.loads(path.read_text())
        if record["outcome"] != "success":
            reasons.append(f"{provider}: {record['outcome']}")
        phases = (((record.get("run_status") or {}).get("timing")) or {}).get("phases") or {}
        if any((phases.get(name) or {}).get("seconds") is None for name in ("acquisition", "preparation", "lec")):
            reasons.append(f"{provider}: a phase was not measured")
        rss = record.get("peak_rss_kib")
        if rss is not None and rss > PILOT_MAX_RSS_FRACTION * ram_kib:
            reasons.append(f"{provider}: peak RSS above {PILOT_MAX_RSS_FRACTION:.0%} of the RAM")
    if not (compare_report or {}).get("equivalent"):
        reasons.append("pilot pair not equivalent: " + "; ".join((compare_report or {}).get("problems", ["not compared"])))
    return reasons


def campaign_elapsed(root):
    """
    Seconds since the first start of the campaign under ROOT (persisted, so a
    resumed campaign keeps the same 96 h ceiling; UTC wall clock across restarts).
    """
    path = private_dir(Path(root) / "records") / "campaign_start.json"
    if not path.is_file():
        write_json(path, {"started_utc": utc(), "epoch": time.time()})
    return time.time() - json.loads(path.read_text())["epoch"]


def r1_complete(root):
    """True when the ten r1 pairs ended in success for both providers."""
    for case in CASES:
        for provider in PROVIDERS:
            path = Path(root) / "records" / f"{case}_{provider}_r1.json"
            if not path.is_file() or json.loads(path.read_text())["outcome"] != "success":
                return False
    return True


def campaign(root, tracks, key_file, cases=None, rep=1):
    """
    Pre-registered order, one cold pair per case. CDS stops for the rest of
    the campaign at its first failure, censoring or unverifiable account
    state (and the ledger blocks it anyway); NCAR stops after two consecutive
    failures. A run starts only if its deadline (and the kill margin) still
    fits in the 96 h campaign ceiling counted from the first start; limits
    are never raised. r2 starts only after the ten r1 pairs succeeded. The
    result says when the campaign is partial.
    """
    if rep > 1 and not r1_complete(root):
        return {"stopped": f"r{rep} refused: the ten r1 pairs are not all successful", "partial": True}
    cds_stopped, ncar_failures = None, 0
    for index, (case, providers) in enumerate(plan_order()):
        if cases and case not in cases:
            continue
        for provider in providers:
            if campaign_elapsed(root) + CASE_DEADLINE_S + KILL_MARGIN_S > CAMPAIGN_WALL_S:
                return {"stopped": "campaign ceiling: no room for another run within 96 h", "partial": True,
                        "cds_stopped": cds_stopped}
            if run_directory(root, case, provider, rep).exists():
                continue  # done earlier (resumed campaign); never repeated as a new cold run
            if provider == "cds":
                if cds_stopped:
                    continue
                try:
                    pending = cds_pending(key_file)
                except Exception as exc:  # cannot verify: submit nothing
                    cds_stopped = f"pending requests could not be listed ({type(exc).__name__})"
                    continue
                if pending:
                    cds_stopped = f"{pending} CDS request(s) still pending on the account"
                    continue
            if provider == "ncar" and ncar_failures >= NCAR_MAX_CONSECUTIVE_FAILURES:
                continue
            record = run(case, provider, rep, root, tracks, key_file=key_file, order=index)
            if provider == "cds" and record["outcome"] != "success":
                cds_stopped = f"case {case}: {record['outcome']}"
            if provider == "ncar":
                ncar_failures = 0 if record["outcome"] == "success" else ncar_failures + 1
        report = compare(case, rep, root, tracks)
        if case == PILOT and rep == 1:
            reasons = pilot_gate(root, rep, report)
            if reasons:
                return {"stopped": "pilot gate not passed", "reasons": reasons, "cds_stopped": cds_stopped,
                        "partial": True}
    complete = all((Path(root) / "records" / f"{c}_{p}_r{rep}.json").is_file() for c in CASES for p in PROVIDERS)
    return {"stopped": None, "cds_stopped": cds_stopped, "ncar_consecutive_failures": ncar_failures,
            "partial": bool(cds_stopped) or not complete}


# ----------------------------------------------------------------- command
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("plan")
    pending = sub.add_parser("cds-pending")
    pending.add_argument("--key-file", required=True)
    total = sub.add_parser("summary")
    total.add_argument("--root", required=True)
    total.add_argument("--rep", type=int, default=1)
    for name in ("run", "compare", "campaign"):
        p = sub.add_parser(name)
        p.add_argument("--root", required=True)
        p.add_argument("--tracks", required=True)
        p.add_argument("--rep", type=int, default=1)
        if name in ("run", "compare"):
            p.add_argument("--case", required=True, choices=CASES)
        if name == "run":
            p.add_argument("--provider", required=True, choices=PROVIDERS)
        if name in ("run", "campaign"):
            p.add_argument("--key-file")
        if name == "campaign":
            p.add_argument("--cases", nargs="*", choices=CASES)
    args = parser.parse_args(argv)

    if args.command == "plan":
        for index, (case, providers) in enumerate(plan_order()):
            print(f"{index:2d} {case} {' -> '.join(providers)}")
        return 0
    if args.command == "summary":
        result = summary(args.root, args.rep)
        print(f"{result['n_valid_pairs']} valid pairs; median NCAR/CDS {result['median_ratio_ncar_cds']}; "
              f"CI {result['ci']}; {result['verdict']}")
        return 0
    if args.command == "cds-pending":
        print(f"pending CDS requests (accepted or running): {cds_pending(args.key_file)}")
        return 0
    if args.command == "run":
        record = run(args.case, args.provider, args.rep, args.root, args.tracks, key_file=args.key_file)
        print(f"{args.case} {args.provider} r{args.rep}: {record['outcome']} "
              f"({record['wall_seconds_external']:.1f} s)")
        return 0 if record["outcome"] == "success" else 1
    if args.command == "compare":
        report = compare(args.case, args.rep, args.root, args.tracks)
        print(f"{args.case} r{args.rep}: equivalent={report['equivalent']} problems={report['problems']}")
        return 0 if report["equivalent"] else 1
    print(json.dumps(campaign(args.root, args.tracks, args.key_file, cases=args.cases, rep=args.rep)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

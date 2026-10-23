"""
Checks on the files written by a LorenzCycleToolkit run, used by the
integration tests.

Every check raises ``AssertionError`` with an explicit message when an
expected file is missing, a time step or a term is absent, a value is not
finite, or a value differs from the reference outputs beyond the tolerance
registered for the audit baseline (see ``tests/data/reference/PROVENANCE.md``).
The files are parsed with the ``csv`` module and ``float`` so that numbers
are compared exactly as written, without another parser in between.
"""

import csv
import json
import math
from datetime import datetime
from pathlib import Path

import pandas as pd
import xarray as xr

REFERENCE_DIR = Path(__file__).resolve().parent / "data" / "reference"

MAIN_TERMS_WITH_RESIDUALS = [
    "Az", "Ae", "Kz", "Ke", "Cz", "Ca", "Ck", "Ce", "C_overturning",
    "BAz", "BAe", "BKz", "BKe", "BΦZ", "BΦE", "M", "Gz", "Ge",
    "∂Az/∂t (finite diff.)", "∂Ae/∂t (finite diff.)",
    "∂Kz/∂t (finite diff.)", "∂Ke/∂t (finite diff.)",
    "RGz", "RKz", "RGe", "RKe",
]

VERTICAL_TERMS = [
    "Az", "Ae", "Kz", "Ke", "Ge", "Gz",
    "Cz", "Cz_1", "Cz_2", "Ca", "Ca_1", "Ca_2",
    "Ce", "Ce_1", "Ce_2", "C_overturning", "M",
    "Ck", "Ck_1", "Ck_2", "Ck_3", "Ck_4", "Ck_5",
]

# Profiles that do not depend on time; in the fixed framework they are
# written as two rows (level values, then the profile).
TIME_INVARIANT_PROFILES = ("Cz_1", "Ce_1")

TRACKFILE_COLUMNS = [
    "time", "Lat", "Lon", "length", "width",
    "min_lon", "max_lon", "min_lat", "max_lat",
    "min_max_zeta_850_lat", "min_max_zeta_850_lon", "min_max_zeta_850",
    "min_hgt_850_lat", "min_hgt_850_lon", "min_hgt_850",
    "max_wind_850_lat", "max_wind_850_lon", "max_wind_850",
]

# Tolerances registered for the audit baseline: (atol, rtol).
ENERGY_TERMS = ("Az", "Ae", "Kz", "Ke")
DEGREE_COLUMNS = (
    "Lat", "Lon", "length", "width", "min_lon", "max_lon", "min_lat", "max_lat",
    "min_max_zeta_850_lat", "min_max_zeta_850_lon",
    "min_hgt_850_lat", "min_hgt_850_lon", "max_wind_850_lat", "max_wind_850_lon",
)
RTOL = 1e-10


def read_rows(path, delimiter=","):
    """Read a CSV file into a list of rows of strings."""
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle, delimiter=delimiter))


def _as_float(text):
    try:
        return float(text)
    except ValueError:
        return None


def _assert_finite_row(row, where):
    for column, cell in enumerate(row):
        value = _as_float(cell)
        assert value is not None, f"{where}: column {column} is not numeric: {cell!r}"
        assert math.isfinite(value), f"{where}: column {column} is not finite: {cell!r}"


def input_times(sample_path, time_name):
    """Time stamps of a NetCDF input as 'YYYY-MM-DD HH:MM:SS' strings."""
    with xr.open_dataset(sample_path) as dataset:
        stamps = pd.to_datetime(dataset[time_name].values)
    return [stamp.strftime("%Y-%m-%d %H:%M:%S") for stamp in stamps]


def track_times(track_path):
    """Time stamps of a track file ('YYYY-MM-DD-HHMM') as 'YYYY-MM-DD HH:MM:SS'."""
    rows = read_rows(track_path, delimiter=";")
    column = rows[0].index("time")
    return [
        datetime.strptime(row[column], "%Y-%m-%d-%H%M").strftime("%Y-%m-%d %H:%M:%S")
        for row in rows[1:] if row
    ]


def check_main_results(results_dir, stem, expected_times, terms=MAIN_TERMS_WITH_RESIDUALS):
    """The main CSV exists, has every term, every expected time and finite values."""
    path = Path(results_dir) / f"{stem}_results.csv"
    assert path.is_file(), f"main results file is missing: {path}"
    rows = read_rows(path)
    header, body = rows[0], rows[1:]
    assert header[0] == "", f"{path.name}: unexpected index header {header[0]!r}"
    assert header[1:] == list(terms), (
        f"{path.name}: terms differ from the expected list; "
        f"missing {sorted(set(terms) - set(header[1:]))}, "
        f"unexpected {sorted(set(header[1:]) - set(terms))}"
    )
    times = [row[0] for row in body]
    assert len(times) == len(expected_times), (
        f"{path.name}: {len(times)} time steps, expected {len(expected_times)}"
    )
    assert times == list(expected_times), f"{path.name}: time axis differs from the input"
    for row in body:
        assert len(row) == len(header), f"{path.name}: row {row[0]} has {len(row)} fields"
        _assert_finite_row(row[1:], f"{path.name} at {row[0]}")
    return path


def check_vertical_levels(results_dir, level_name, time_name, expected_times, method):
    """The 23 per-level CSVs exist, carry the levels and finite values for every time."""
    vertical_dir = Path(results_dir) / "results_vertical_levels"
    assert vertical_dir.is_dir(), f"per-level directory is missing: {vertical_dir}"
    paths = []
    for term in VERTICAL_TERMS:
        path = vertical_dir / f"{term}_{level_name}.csv"
        assert path.is_file(), f"per-level file is missing: {path.name}"
        rows = read_rows(path)
        header, body = rows[0], rows[1:]
        assert header[0] == time_name, f"{path.name}: first header cell {header[0]!r}"
        levels = [_as_float(cell) for cell in header[1:]]
        assert levels and all(v is not None for v in levels), f"{path.name}: bad level header"
        assert levels == sorted(set(levels)), f"{path.name}: levels not strictly increasing"
        if method == "fixed" and term in TIME_INVARIANT_PROFILES:
            assert [row[0] for row in body] == [level_name, term], (
                f"{path.name}: unexpected layout of a time-invariant profile"
            )
            assert [_as_float(c) for c in body[0][1:]] == levels, f"{path.name}: level row"
            _assert_finite_row(body[1][1:], path.name)
        else:
            times = [row[0] for row in body]
            assert times == list(expected_times), (
                f"{path.name}: {len(times)} rows, expected the {len(expected_times)} time steps"
            )
            for row in body:
                assert len(row) == len(header), f"{path.name}: row {row[0]} has {len(row)} fields"
                _assert_finite_row(row[1:], f"{path.name} at {row[0]}")
        paths.append(path)
    return paths


def check_trackfile(results_dir, stem, expected_times):
    """The processed track file exists, has every column, time step and finite values."""
    path = Path(results_dir) / f"{stem}_trackfile"
    assert path.is_file(), f"trackfile is missing: {path}"
    rows = read_rows(path, delimiter=";")
    header, body = rows[0], rows[1:]
    assert header == TRACKFILE_COLUMNS, f"{path.name}: columns differ: {header}"
    times = [
        datetime.strptime(row[0], "%Y-%m-%d-%H%M").strftime("%Y-%m-%d %H:%M:%S")
        for row in body
    ]
    assert times == list(expected_times), f"{path.name}: time steps differ from the track"
    for row in body:
        assert len(row) == len(header), f"{path.name}: row {row[0]} has {len(row)} fields"
        _assert_finite_row(row[1:], f"{path.name} at {row[0]}")
    return path


def check_log(results_dir, infile_stem):
    """
    The run log exists, is not empty and holds the messages of this run only.

    The integration tests log "Starting LEC analysis" and the parsed arguments
    once per run; a second start, or arguments naming another input file,
    mean that messages of another run were written to this log.
    """
    path = Path(results_dir) / f"log.{infile_stem}"
    assert path.is_file(), f"log file is missing: {path}"
    assert path.stat().st_size > 0, f"log file is empty: {path}"
    lines = path.read_text(encoding="utf-8").splitlines()
    starts = sum("Starting LEC analysis" in line for line in lines)
    assert starts == 1, f"{path.name}: {starts} runs logged, expected 1"
    infiles = [line.split("infile='", 1)[1].split("'", 1)[0] for line in lines if "infile='" in line]
    own = f"samples/{infile_stem}.nc"
    assert infiles == [own], f"{path.name}: arguments of other runs logged: {infiles}"
    return path


# Products of run_status.json (pre-registered inventory, Onda 1B; Onda 3 adds
# lec_metadata.json and domain_selection.csv): required scientific outputs and
# optional figures/diagnostics requested with --plots.
SCIENTIFIC_PRODUCTS = {
    "fixed": ("results_csv", "vertical_level_csvs", "lec_metadata", "domain_selection"),
    "track": ("results_csv", "vertical_level_csvs", "trackfile", "lec_metadata", "domain_selection"),
}
FIGURE_PRODUCTS = {
    "fixed": ("timeseries", "box_limits_map", "boxplots", "hovmoller", "lec_cycle"),
    "track": ("domain_maps", "track_map", "zeta_hgt_timeseries", "timeseries", "hovmoller",
              "boxplots", "periods", "lec_cycle", "lec_periods", "lps"),
}


def products(method, figures="success", **overrides):
    """Expected state of every product: scientific outputs succeed, figures as given."""
    states = {name: "success" for name in SCIENTIFIC_PRODUCTS[method]}
    states.update({name: figures for name in FIGURE_PRODUCTS[method]})
    states.update(overrides)
    return states


def load_run_status(results_dir):
    """Return the parsed run_status.json of a run."""
    path = Path(results_dir) / "run_status.json"
    assert path.is_file(), f"run_status.json is missing in {results_dir}"
    return json.loads(path.read_text(encoding="utf-8"))


def check_run_status(results_dir, status, exit_code, products):
    """
    run_status.json agrees with the expected run and with the files on disk.

    ``products`` maps every product name to its expected state. A product in
    "success" must list expected files, all of them written by this run and
    present; a "failed" one must give a reason and must not list as written
    any file it is missing.
    """
    data = load_run_status(results_dir)
    assert data["schema_version"] == 1, f"schema_version {data['schema_version']}"
    assert data["status"] == status, f"run status {data['status']!r}, expected {status!r}"
    assert data["exit_code"] == exit_code, f"exit_code {data['exit_code']!r}, expected {exit_code!r}"
    by_name = {p["name"]: p for p in data["products"]}
    assert set(by_name) == set(products), f"products {sorted(by_name)}, expected {sorted(products)}"
    for name, expected in products.items():
        product = by_name[name]
        assert product["status"] == expected, (
            f"{name}: status {product['status']!r}, expected {expected!r} (reason: {product['reason']})"
        )
        if expected == "success":
            assert product["expected"], f"{name}: no expected files declared"
            for rel in product["expected"]:
                assert rel in product["written"], f"{name}: {rel} declared but not written by this run"
                assert (Path(results_dir) / rel).is_file(), f"{name}: {rel} declared and missing"
        if expected == "failed":
            assert product["reason"], f"{name}: failed without a reason"
            assert not set(product["missing"]) & set(product["written"]), f"{name}: missing files listed as written"
    for state, count in data["summary"].items():
        assert count == sum(p["status"] == state for p in data["products"]), f"summary count of {state!r}"
    return data


def check_final_messages(results_dir, infile_stem, complete):
    """
    The closing messages of the log tell the truth: the success messages only
    when every requested product was generated, the PARTIAL message otherwise.
    """
    text = (Path(results_dir) / f"log.{infile_stem}").read_text(encoding="utf-8")
    if complete:
        assert "🎉 Analysis complete!" in text, "success message missing"
        assert "PARTIAL" not in text, "PARTIAL message in a complete run"
    else:
        assert "Analysis complete with PARTIAL output" in text, "PARTIAL message missing"
        assert "🎉 Analysis complete!" not in text, "success message in a partial run"
        assert "All plots generated successfully!" not in text, "plots declared successful in a partial run"


def _tolerance(kind, column_name):
    if kind == "main":
        if column_name in ENERGY_TERMS:
            return 1e-6, RTOL
        if column_name == "M":
            return 1e-10, RTOL
        return 1e-9, RTOL
    if kind == "trackfile":
        if column_name in DEGREE_COLUMNS:
            return 1e-9, 0.0
        return 1e-12, RTOL
    return 1e-12, RTOL  # per-level profiles, including M (s-1)


def compare_with_reference(new_path, reference_path, kind):
    """
    Compare a CSV with its reference cell by cell.

    Text cells (times, labels, level headers) must be identical; numeric cells
    must satisfy |a - b| <= atol + rtol * max(|a|, |b|) with the registered
    tolerance of their column; NaN must sit in the same cells.
    """
    delimiter = ";" if kind == "trackfile" else ","
    new_rows = read_rows(new_path, delimiter)
    ref_rows = read_rows(reference_path, delimiter)
    name = Path(new_path).name
    assert new_rows[0] == ref_rows[0], f"{name}: header differs from the reference"
    assert len(new_rows) == len(ref_rows), (
        f"{name}: {len(new_rows) - 1} data rows, reference has {len(ref_rows) - 1}"
    )
    header = ref_rows[0]
    failures = []
    for new_row, ref_row in zip(new_rows[1:], ref_rows[1:]):
        assert len(new_row) == len(ref_row), f"{name}: row {ref_row[0]} has a different width"
        assert new_row[0] == ref_row[0], f"{name}: row label {new_row[0]!r} != {ref_row[0]!r}"
        for column, (a_text, b_text) in enumerate(zip(new_row[1:], ref_row[1:]), start=1):
            a, b = _as_float(a_text), _as_float(b_text)
            if a is None or b is None:
                if a_text != b_text:
                    failures.append(f"{ref_row[0]} / {header[column]}: {a_text!r} != {b_text!r}")
                continue
            if math.isnan(a) or math.isnan(b):
                if not (math.isnan(a) and math.isnan(b)):
                    failures.append(f"{ref_row[0]} / {header[column]}: NaN position differs")
                continue
            atol, rtol = _tolerance(kind, header[column])
            if abs(a - b) > atol + rtol * max(abs(a), abs(b)):
                failures.append(
                    f"{ref_row[0]} / {header[column]}: {a!r} vs reference {b!r} "
                    f"(|diff| = {abs(a - b):.3e}, atol = {atol:g}, rtol = {rtol:g})"
                )
    assert not failures, f"{name}: {len(failures)} values differ from the reference:\n" + "\n".join(failures[:20])


def assert_matches_reference(results_dir, reference_case_dir, stem, level_name):
    """Compare the main CSV, the trackfile (if any) and every per-level CSV with the reference."""
    results_dir, reference_case_dir = Path(results_dir), Path(reference_case_dir)
    compare_with_reference(
        results_dir / f"{stem}_results.csv", reference_case_dir / f"{stem}_results.csv", "main"
    )
    reference_track = reference_case_dir / f"{stem}_trackfile"
    if reference_track.is_file():
        compare_with_reference(results_dir / f"{stem}_trackfile", reference_track, "trackfile")
    for term in VERTICAL_TERMS:
        file_name = f"{term}_{level_name}.csv"
        compare_with_reference(
            results_dir / "results_vertical_levels" / file_name,
            reference_case_dir / "results_vertical_levels" / file_name,
            f"level:{term}",
        )

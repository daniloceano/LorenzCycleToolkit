"""
Scientific metadata of one run: ``lec_metadata.json`` and ``domain_selection.csv``.

Both files are additional outputs written next to the results; no existing
file, column or header changes. They describe what the computation actually
used, read from the objects of the run (nothing is recomputed and no
selection is changed here):

``domain_selection.csv``
    One row per time step (``track``/``choose``) or one row with time
    ``all`` (``fixed``): the requested box limits (in the longitude
    convention of the dataset), the coordinates actually selected (first and
    last longitude and latitude of the extracted fields), the difference
    between them, the number of points, the extent of the input data and
    whether each selected edge is the edge of the input data.

``lec_metadata.json``
    Provenance (version, command, SHA-256 of the input file, namelist), the
    time axis, the retained pressure layer, the units of every column of the
    main results file and, for every per-level file, its units, its relation
    to the column value and its layout.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import contextlib
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from ..analysis.energy_contents import EnergyContents
from ..analysis.mass_continuity import MassContinuity
from .output_spec import TENDENCY_UNITS, VERTICAL_PROFILE_SPEC
from .run_status import _version
from .term_output import TermOutputMixin, header_levels

SCHEMA_VERSION = 1
METADATA_FILE = "lec_metadata.json"
DOMAIN_SELECTION_FILE = "domain_selection.csv"

# Attributes set on the dataset by select_area.slice_domain before the domain
# is cut: the extent of the input data in the dataset's own convention.
INPUT_EXTENT_ATTRS = {
    "west": "lectk_input_lon_min",
    "east": "lectk_input_lon_max",
    "south": "lectk_input_lat_min",
    "north": "lectk_input_lat_max",
}

# Attributes set on the dataset by preprocessing.process_data at the cut of the
# levels above 10 hPa: the grid before the cut and the levels it discarded (DEC-09).
TOP_CUT_ATTRS = {
    "input_grid": "lectk_input_grid_before_top_cut_Pa",
    "discarded": "lectk_levels_discarded_above_10hPa_Pa",
}
TOP_CUT_RULE = (
    "levels with a pressure lower than 1000 Pa (10 hPa) are discarded when the data are "
    "read; 1000 Pa itself is kept (DEC-09)"
)

# How the messages name each framework (as in vertical_harmonization.MODES).
FRAMEWORK_LABELS = {
    "fixed": "fixed framework (-f)",
    "track": "track mode (-t)",
    "choose": "--choose mode",
}

COLUMN_STATEMENT = (
    "Column values are integrals over the retained pressure layer "
    "(p_top_Pa to p_bottom_Pa), not over the entire atmospheric column."
)

ENERGY_COLUMNS = ("Az", "Ae", "Kz", "Ke")
MASS_COLUMNS = ("M",)
FLUX_COLUMNS = (
    "Cz", "Ca", "Ck", "Ce", "C_overturning",
    "BAz", "BAe", "BKz", "BKe", "BΦZ", "BΦE",
    "Gz", "Ge", "Dz", "De", "RGz", "RKz", "RGe", "RKe",
)
NOT_VERIFIED_COLUMNS = ("Dz", "De")

DOMAIN_COLUMNS = [
    "time",
    "requested_west", "requested_east", "requested_south", "requested_north",
    "selected_west", "selected_east", "selected_south", "selected_north",
    "offset_west", "offset_east", "offset_south", "offset_north",
    "n_lon", "n_lat",
    "input_west", "input_east", "input_south", "input_north",
    "west_at_input_edge", "east_at_input_edge", "south_at_input_edge", "north_at_input_edge",
    "requested_outside_input",
    "fields_share_selection",
]

# Fields extracted by BoxData whose horizontal selection is recorded.
_BOX_FIELDS = ("tair", "u", "v", "omega", "geopt")


def _units(text):
    return text.replace("^", "**")


def column_units(column):
    """
    Units of a column of the main results file, from the term classes.

    Args:
        column (str): Column name (e.g. ``"Az"``, ``"∂Kz/∂t (finite diff.)"``).

    Returns:
        str or None: Pint-style units (``J/m**2``, ``W/m**2``,
        ``kg/m**2/s``), or None for an unknown column.
    """
    if column in ENERGY_COLUMNS:
        return _units(EnergyContents.OUTPUT_UNITS)
    if column in MASS_COLUMNS:
        return _units(MassContinuity.OUTPUT_UNITS)
    if column in TENDENCY_UNITS:
        return TENDENCY_UNITS[column]
    if column in FLUX_COLUMNS:
        return _units(TermOutputMixin.OUTPUT_UNITS)
    return None


def sha256_file(path, chunk_size=1 << 20):
    """Return the SHA-256 of a file, read in chunks (constant memory)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record_input_extent(dataset, lon_name, lat_name):
    """
    Attach the extent of the input data to the dataset attributes.

    Called before the domain is cut, so that the edges of the selected box can
    later be compared with the edges of the data actually available.

    Args:
        dataset (xr.Dataset): The processed input data (before the cut).
        lon_name (str): Name of the longitude coordinate.
        lat_name (str): Name of the latitude coordinate.

    Returns:
        dict: {"west", "east", "south", "north"} in the dataset convention.
    """
    lon = np.asarray(dataset[lon_name].values, dtype=float)
    lat = np.asarray(dataset[lat_name].values, dtype=float)
    extent = {
        "west": float(np.nanmin(lon)), "east": float(np.nanmax(lon)),
        "south": float(np.nanmin(lat)), "north": float(np.nanmax(lat)),
    }
    for edge, attr in INPUT_EXTENT_ATTRS.items():
        dataset.attrs[attr] = extent[edge]
    return extent


def input_extent(dataset):
    """Return the extent recorded by record_input_extent, or None for each missing edge."""
    return {edge: dataset.attrs.get(attr) for edge, attr in INPUT_EXTENT_ATTRS.items()}


def _pa(levels):
    return ", ".join(f"{p:.0f} Pa" for p in levels)


def record_top_cut(dataset, input_grid, level_name, app_logger=None):
    """
    Attach to the dataset, and report, what the cut of the levels above 10 hPa did.

    The cut is made by ``preprocessing.process_data``; the levels kept are read
    here from the dataset after the cut, so the report is exactly what the cut
    did. A warning is given only when some level was discarded.

    Args:
        dataset (xr.Dataset): The data after the cut.
        input_grid (list): Pressure levels (Pa) of the input, before the cut.
        level_name (str): Name of the vertical coordinate.
        app_logger (logging.Logger, optional): Logger of the run.

    Returns:
        list: The levels discarded (Pa).
    """
    kept = [float(p) for p in dataset[level_name].values]
    discarded = [p for p in input_grid if p not in kept]
    dataset.attrs[TOP_CUT_ATTRS["input_grid"]] = list(input_grid)
    dataset.attrs[TOP_CUT_ATTRS["discarded"]] = discarded
    if discarded and app_logger is not None:
        app_logger.warning(
            f"⚠️ Vertical levels above 10 hPa discarded: {len(discarded)} of {len(input_grid)} "
            f"input levels have a pressure lower than 1000 Pa (10 hPa) and are not used; 1000 Pa "
            f"itself is kept (DEC-09). Discarded: {_pa(discarded)}. Input grid: {len(input_grid)} "
            f"levels, {_pa(input_grid)}; grid after preprocessing: {len(kept)} levels, {_pa(kept)}."
        )
    return discarded


def top_cut(dataset):
    """Return the record of record_top_cut ({input_grid_Pa, discarded_Pa}), or None if absent."""
    if TOP_CUT_ATTRS["input_grid"] not in dataset.attrs:
        return None
    return {
        "input_grid_Pa": [float(p) for p in dataset.attrs[TOP_CUT_ATTRS["input_grid"]]],
        "discarded_Pa": [float(p) for p in dataset.attrs[TOP_CUT_ATTRS["discarded"]]],
    }


class RunMetadata:
    """
    Collects, during the run, what the computation actually used, and writes
    ``domain_selection.csv`` and ``lec_metadata.json``.

    Args:
        args (argparse.Namespace): Command-line arguments.
        method (str): ``fixed``, ``track`` or ``choose``.
        variable_list_df (pd.DataFrame): Parsed namelist.
        extent (dict): Extent of the input data (see ``input_extent``).
        top_cut (dict, optional): What the cut at 10 hPa did (see ``top_cut``).
    """

    def __init__(self, args, method, variable_list_df, extent, top_cut=None):
        self.args = args
        self.method = method
        self.variable_list_df = variable_list_df
        self.extent = extent
        self.top_cut = top_cut
        self.domain_rows = []
        self.vertical_rows = []
        self.input_levels = None
        self.harmonization = None

    def record_box(self, box_obj, requested, time_label, input_levels, app_logger=None):
        """
        Record the horizontal selection and the retained levels of one box.

        A box requested beyond the input data is cut at the data edge, as
        before; it is reported with a warning that gives the limits requested
        and the limits of the fields actually used (DEC-11).

        Args:
            box_obj (BoxData): The box of this time step (or of the fixed run).
            requested (dict): Requested limits {"west", "east", "south", "north"}
                as passed to BoxData.
            time_label (str): ``"all"`` (fixed) or the time of the step.
            input_levels (sequence): Pressure levels (Pa) of the data given to BoxData.
            app_logger (logging.Logger, optional): Logger of the run.
        """
        lon_name, lat_name = box_obj.LonIndexer, box_obj.LatIndexer
        coords = []
        for name in _BOX_FIELDS:
            field = getattr(box_obj, name, None)
            if field is None:
                continue
            coords.append((
                tuple(np.asarray(field[lon_name].values, dtype=float)),
                tuple(np.asarray(field[lat_name].values, dtype=float)),
            ))
        lon, lat = (np.asarray(c, dtype=float) for c in coords[0])
        selected = {
            "west": float(lon.min()), "east": float(lon.max()),
            "south": float(lat.min()), "north": float(lat.max()),
        }
        row = {"time": time_label}
        for edge in ("west", "east", "south", "north"):
            row[f"requested_{edge}"] = float(requested[edge])
        for edge in ("west", "east", "south", "north"):
            row[f"selected_{edge}"] = selected[edge]
        for edge in ("west", "east", "south", "north"):
            row[f"offset_{edge}"] = selected[edge] - float(requested[edge])
        row["n_lon"], row["n_lat"] = int(lon.size), int(lat.size)
        for edge in ("west", "east", "south", "north"):
            row[f"input_{edge}"] = self.extent.get(edge)
        for edge in ("west", "east", "south", "north"):
            value = self.extent.get(edge)
            row[f"{edge}_at_input_edge"] = None if value is None else bool(selected[edge] == value)
        outside = []
        if self.extent.get("west") is not None and requested["west"] < self.extent["west"]:
            outside.append("west")
        if self.extent.get("east") is not None and requested["east"] > self.extent["east"]:
            outside.append("east")
        if self.extent.get("south") is not None and requested["south"] < self.extent["south"]:
            outside.append("south")
        if self.extent.get("north") is not None and requested["north"] > self.extent["north"]:
            outside.append("north")
        row["requested_outside_input"] = "+".join(outside) if outside else "none"
        row["fields_share_selection"] = all(c == coords[0] for c in coords)
        self.domain_rows.append(row)
        if outside and app_logger is not None:
            when = "" if time_label == "all" else f", {time_label}"
            e = {edge: ("n/a" if v is None else f"{v:g}") for edge, v in self.extent.items()}
            app_logger.warning(
                f"⚠️ Domain clipping ({FRAMEWORK_LABELS.get(self.method, self.method)}{when}): the "
                f"requested box reaches beyond the input data on the {', '.join(outside)} edge(s) and "
                "is cut at the data edge (DEC-11); the computation goes on with the box actually "
                f"used. Requested: lon {requested['west']:g} to {requested['east']:g}, lat "
                f"{requested['south']:g} to {requested['north']:g}; used: lon {selected['west']:g} to "
                f"{selected['east']:g}, lat {selected['south']:g} to {selected['north']:g} "
                f"({row['n_lon']} x {row['n_lat']} grid points); input data: lon {e['west']} to "
                f"{e['east']}, lat {e['south']} to {e['north']}."
            )

        # Dropped = input levels not used, also when the levels were removed
        # before BoxData (common vertical volume of the moving framework).
        self.input_levels = [float(p) for p in input_levels]
        used = [float(p) for p in box_obj.valid_levels]
        self.vertical_rows.append({
            "time": time_label,
            "levels_used_Pa": used,
            "dropped_levels_Pa": [p for p in self.input_levels if p not in used],
        })

    def record_harmonization(self, harmonization):
        """
        Record the levels common to all time steps of the moving framework.

        Args:
            harmonization (VerticalHarmonization): Result of ``harmonize_levels``.
        """
        self.harmonization = harmonization.as_metadata()

    # ------------------------------------------------------------------ files
    def write_domain_selection(self, results_directory):
        """Write ``domain_selection.csv`` and return its path."""
        path = Path(results_directory) / DOMAIN_SELECTION_FILE
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=DOMAIN_COLUMNS, lineterminator="\n")
            writer.writeheader()
            for row in self.domain_rows:
                writer.writerow({k: _csv_value(v) for k, v in row.items()})
        return path

    def _vertical(self, level_files):
        used = [r["levels_used_Pa"] for r in self.vertical_rows]
        constant = all(u == used[0] for u in used)
        retained = sorted(set.intersection(*(set(u) for u in used))) if used else []
        header = header_levels(level_files[0]) if level_files else None
        vertical = {
            "coordinate": self.variable_list_df.loc["Vertical Level"]["Variable"],
            "units": "Pa",
            "input_levels_Pa": self.input_levels,
            "levels_used_Pa": used[0] if constant else None,
            "levels_used_at_every_time_Pa": retained,
            "levels_used_vary_with_time": not constant,
            "dropped_levels_Pa": self.vertical_rows[0]["dropped_levels_Pa"] if constant else None,
            "p_top_Pa": min(used[0]) if constant else None,
            "p_bottom_Pa": max(used[0]) if constant else None,
            "per_level_file_header_levels_Pa": header,
            "column_statement": COLUMN_STATEMENT,
        }
        if self.method != "fixed":
            vertical["by_time"] = self.vertical_rows
        if self.harmonization is not None:
            vertical["harmonization"] = self.harmonization
        if self.top_cut is not None:
            vertical["top_cut"] = {
                "rule": TOP_CUT_RULE,
                "input_grid_before_cut_Pa": self.top_cut["input_grid_Pa"],
                "discarded_levels_Pa": self.top_cut["discarded_Pa"],
            }
        return vertical

    def write_metadata(self, results_directory, results_file, level_files, columns, times,
                       trackfile=None, run_status_file="run_status.json"):
        """
        Write ``lec_metadata.json`` atomically and return its path.

        Args:
            results_directory (str): Directory of the results of this run.
            results_file (str): Main results CSV.
            level_files (list): Per-level CSV files, as written.
            columns (sequence): Columns of the main results file.
            times (sequence): Times of the main results file.
            trackfile (str, optional): System track file (moving framework).
            run_status_file (str): Name of the status file.
        """
        results_directory = Path(results_directory)
        times = np.asarray(times, dtype="datetime64[ns]")
        steps = sorted({int(s) for s in np.diff(times).astype("timedelta64[s]").astype(int)}) if times.size > 1 else []
        names = self.variable_list_df
        namelist = {
            str(row): {
                "variable": str(names.loc[row]["Variable"]),
                "units": None if "Units" not in names.columns else _plain(names.loc[row]["Units"]),
            }
            for row in names.index
        }
        infile = self.args.infile
        metadata = {
            "schema_version": SCHEMA_VERSION,
            "lectk_version": _version(),
            "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "command": {
                "method": self.method,
                "residuals": bool(getattr(self.args, "residuals", False)),
                "outname": getattr(self.args, "outname", None),
                "box_limits": getattr(self.args, "box_limits", None) if self.method == "fixed" else None,
                "trackfile": getattr(self.args, "trackfile", None) if self.method == "track" else None,
                "zeta": bool(getattr(self.args, "zeta", False)),
                "mpas": bool(getattr(self.args, "mpas", False)),
            },
            "input": {
                "file": os.path.basename(infile),
                "sha256": sha256_file(infile) if os.path.isfile(infile) else None,
                "size_bytes": os.path.getsize(infile) if os.path.isfile(infile) else None,
                "namelist": namelist,
            },
            "time": {
                "n": int(times.size),
                "first": _time(times[0]) if times.size else None,
                "last": _time(times[-1]) if times.size else None,
                "distinct_steps_s": steps,
            },
            "horizontal": {
                "domain_selection_file": DOMAIN_SELECTION_FILE,
                "input_extent": self.extent,
                "longitude_convention": "requested limits are expressed in the longitude convention of the dataset",
                "boxes_cut_at_input_edge": [
                    {"time": r["time"], "edges": r["requested_outside_input"].split("+")}
                    for r in self.domain_rows if r["requested_outside_input"] != "none"
                ],
            },
            "vertical": self._vertical(level_files),
            "main_results": {
                "file": _relative(results_file, results_directory),
                "time_column_header": "",
                "time_column_note": "the first column holds the time; its header cell is empty",
                "columns": [
                    {"name": c, "units": column_units(c),
                     **({"verification": "NOT VERIFIED"} if c in NOT_VERIFIED_COLUMNS else {})}
                    for c in columns
                ],
            },
            "vertical_level_files": [
                _level_file_entry(p, results_directory, names.loc["Vertical Level"]["Variable"]) for p in level_files
            ],
            "trackfile": _relative(trackfile, results_directory) if trackfile else None,
            "run_status_file": run_status_file,
        }
        path = results_directory / METADATA_FILE
        tmp = path.with_name(f".{METADATA_FILE}.{os.getpid()}.tmp")
        try:
            with open(tmp, "w", encoding="utf-8") as handle:
                json.dump(metadata, handle, indent=2, ensure_ascii=False)
                handle.write("\n")
            os.replace(tmp, path)
        except BaseException:
            with contextlib.suppress(OSError):
                tmp.unlink()
            raise
        return path


def _csv_value(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    return value


def _plain(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    return str(value)


def _time(value):
    return str(np.datetime_as_string(value, unit="s")).replace("T", " ")


def _relative(path, root):
    try:
        return Path(path).resolve().relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return str(path)


def _layout(path):
    """
    Layout of a per-level file, read from the file itself.

    ``time_by_level``: one row per time, one column per level.
    ``time_invariant_profile``: a row ``level`` with the levels and a row
    with the term name and its values (Cz_1/Ce_1 in the fixed framework,
    whose values do not depend on time).
    """
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    if len(rows) >= 2 and rows[1] and rows[1][0] == "level":
        return "time_invariant_profile"
    return "time_by_level"


def _level_file_entry(path, root, coordinate):
    suffix = f"_{coordinate}.csv"
    name = Path(path).name
    term = name[: -len(suffix)] if name.endswith(suffix) else Path(path).stem
    spec = VERTICAL_PROFILE_SPEC.get(term, {})
    entry = {
        "file": _relative(path, root),
        "term": term,
        "layout": _layout(path),
        "units": spec.get("units"),
        "role": spec.get("role"),
        "column_term": spec.get("column_term"),
        "column_units": spec.get("column_units"),
        "column_factor": spec.get("column_factor"),
    }
    if "part_of" in spec:
        entry["part_of"] = spec["part_of"]
    return entry

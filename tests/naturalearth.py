"""
Natural Earth data for the tests that draw maps (offline, verified).

Cartopy draws coastlines, borders, land, ocean, states and cities from
Natural Earth shapefiles and downloads a missing file the first time it is
needed. The tests never download: they require ``CARTOPY_DATA_DIR`` (the
Cartopy variable for pre-existing data) to point to a copy of the files
listed, with their sha256, in ``tests/data/naturalearth_manifest.tsv``.
If the variable is not set, a file is missing or a file differs, the tests
that need maps stop with an explicit error instead of being skipped.

To prepare the copy, use an existing verified copy of these files, or let
Cartopy download them once outside the tests and check them against the
manifest; keep the layout ``shapefiles/natural_earth/{physical,cultural}/``.
"""

import csv
import hashlib
import os
from pathlib import Path


class NaturalEarthCacheError(RuntimeError):
    """The Natural Earth data needed by the map tests are not available."""


def read_manifest(manifest):
    """Return the rows of the manifest (comment lines start with '#')."""
    with open(manifest, newline="", encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_cache(root, manifest):
    """
    Check every file of the manifest under ``root``.

    Returns:
        list: The manifest rows (all verified).

    Raises:
        NaturalEarthCacheError: Listing the missing and differing files.
    """
    root = Path(root)
    rows = read_manifest(manifest)
    problems = []
    for row in rows:
        path = root / row["relpath"]
        if not path.is_file():
            problems.append(f"missing {row['relpath']}")
        elif path.stat().st_size != int(row["bytes"]) or _sha256(path) != row["sha256"]:
            problems.append(f"differs from the manifest: {row['relpath']}")
    if problems:
        shown = "; ".join(problems[:8]) + ("; ..." if len(problems) > 8 else "")
        raise NaturalEarthCacheError(
            f"The Natural Earth data in {root} do not match {Path(manifest).name} "
            f"({len(problems)} of {len(rows)} files): {shown}. "
            "The map tests need a verified copy and never download it (see tests/naturalearth.py)."
        )
    return rows


def cache_from_environment(manifest):
    """Verify the directory named by CARTOPY_DATA_DIR and return it."""
    root = os.environ.get("CARTOPY_DATA_DIR")
    if not root:
        raise NaturalEarthCacheError(
            "CARTOPY_DATA_DIR is not set. The tests that draw maps need the Natural Earth "
            f"files listed in {Path(manifest).name} and never download them (see tests/naturalearth.py)."
        )
    verify_cache(root, manifest)
    return Path(root)

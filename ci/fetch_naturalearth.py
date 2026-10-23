"""
Download, before the test run, the Natural Earth files needed by the map tests
and check them against tests/data/naturalearth_manifest.tsv.

The tests never download: they require CARTOPY_DATA_DIR to point to a copy of
the files of the manifest (see tests/naturalearth.py). This script prepares
that copy for the CI. The files come from the Natural Earth bucket used by
Cartopy's own downloader (cartopy.io.shapereader.NEShpDownloader), in
Cartopy's layout, under the target directory. Every file is then checked
against the sha256 and size of the manifest; a missing or differing file
stops the script with an error, and the tests are not run.

Usage: python ci/fetch_naturalearth.py <target directory>
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import cartopy  # noqa: E402
from cartopy.io import shapereader  # noqa: E402

from tests.naturalearth import read_manifest, verify_cache  # noqa: E402

MANIFEST = ROOT / "tests" / "data" / "naturalearth_manifest.tsv"


def main(target):
    target = Path(target).resolve()
    target.mkdir(parents=True, exist_ok=True)
    cartopy.config["data_dir"] = str(target)
    cartopy.config["pre_existing_data_dir"] = ""
    downloader = shapereader.NEShpDownloader.default_downloader()
    sets = sorted({(row["scale"], row["category"], row["name"]) for row in read_manifest(MANIFEST)})
    for scale, category, name in sets:
        url = downloader.url_template.format(resolution=scale, category=category, name=name)
        print(f"Natural Earth {scale} {category} {name}: {url}", flush=True)
        shapereader.natural_earth(resolution=scale, category=category, name=name)
    rows = verify_cache(target, MANIFEST)
    print(f"{len(rows)} files verified against {MANIFEST.name} (sha256 and size) in {target}")


if __name__ == "__main__":
    main(sys.argv[1])

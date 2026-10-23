"""
Offline maps (--offline-maps): draw maps from local Natural Earth data only.

Cartopy downloads a Natural Earth shapefile the first time a map needs it.
With --offline-maps, the Natural Earth downloaders of Cartopy are replaced,
for the duration of the run only, by ones that return the local file
(``CARTOPY_DATA_DIR`` or Cartopy's data directory) and raise
NaturalEarthOfflineError when it is missing, before any network access
(no name resolution, no connection, no SSL context). The map is then
reported as failed and the run as partial; the computation is unaffected.
Without the option, Cartopy's behaviour is unchanged.
"""

import contextlib

NE_KEY = ("shapefiles", "natural_earth")


class NaturalEarthOfflineError(FileNotFoundError):
    """A Natural Earth file needed by a map is missing and downloads are off."""


class _OfflineDownloader:
    """Local lookup of a Cartopy Downloader, without its download step."""

    def __init__(self, downloader):
        self.downloader = downloader

    def path(self, format_dict):
        pre = self.downloader.pre_downloaded_path(format_dict)
        target = self.downloader.target_path(format_dict)
        if pre is not None and pre.exists():
            return pre
        if target.exists():
            return target
        raise NaturalEarthOfflineError(
            f"Natural Earth file {target.name} is not available locally "
            f"(looked in {pre.parent if pre is not None else '-'} and {target.parent}) and "
            "--offline-maps forbids downloading it; set CARTOPY_DATA_DIR to a directory "
            "with the Natural Earth shapefiles"
        )


@contextlib.contextmanager
def offline_maps(enabled):
    """Within the block, Natural Earth data are read locally or not at all."""
    if not enabled:
        yield
        return
    import cartopy
    import cartopy.io.shapereader  # noqa: F401  (registers the Natural Earth downloader)

    downloaders = cartopy.config["downloaders"]
    saved = {key: value for key, value in downloaders.items() if tuple(key[:2]) == NE_KEY}
    for key, value in saved.items():
        downloaders[key] = _OfflineDownloader(value)
    try:
        yield
    finally:
        for key, value in saved.items():
            downloaders[key] = value

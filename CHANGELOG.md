# Changelog

All notable changes to this project will be documented in this file.

## [1.1.5] - 2026-01-21

### Improved

- **CDS API Download Logging**: Significantly enhanced logging for ERA5 data downloads:
  - Shows list of expected files before download starts
  - Displays progress with checkmarks (✓) for successful downloads
  - Shows time ranges being downloaded for each day (start time to end time)
  - Provides clear success summary with file size and time coverage
  - More informative error messages if downloads fail

- **Smart Time Range Detection**: Fixed issue where downloads always requested full days (00Z-23Z):
  - Now detects actual start and end times from track file
  - For first day: downloads only from track start time onwards
  - For last day: downloads only until track end time
  - For middle days: downloads complete day as before
  - Respects the `--time-resolution` parameter for all calculations
  - Example: Track from 18Z to 23Z on same day will only download those hours

### Technical Details

- Track files can now start at any hour (e.g., 18Z) and end at any hour (e.g., 15Z)
- Time ranges are automatically rounded to the nearest `time_resolution` interval
- Reduces unnecessary data downloads and API costs
- Better handles single-day tracks with partial day coverage

## [1.1.4] - 2026-01-21

### Changed

- **Dependencies Update**: Updated CDS API and multiurl dependencies to use flexible version constraints:
  - `cdsapi`: Changed from `==0.7.0` to `>=0.7.0` (now compatible with version 0.7.7)
  - `multiurl`: Changed from `==0.3.1` to `>=0.3.1` (now compatible with version 0.3.7)
  - This resolves dependency conflicts when users upgrade cdsapi to the latest version
  - The package now automatically adapts to newer versions of these dependencies

### Notes

- Users can now safely upgrade cdsapi without encountering dependency conflicts
- The deprecated CDS API endpoint warning will persist until Copernicus fully migrates to the new API

## [1.1.3] - 2026-01-21

### Fixed

- **CDS API Download Issues**: Modified `get_cdsapi_data` function to download ERA5 data day-by-day instead of all at once. This resolves issues with large file downloads that were causing failures. The function now:
  - Downloads data separately for each day
  - Stores temporary daily files in a system temporary directory
  - Concatenates all daily files into the final output file
  - Automatically cleans up temporary files after concatenation
  - Provides better progress logging for multi-day downloads
- **Test Suite Warnings**: Fixed FutureWarnings in test suite by updating pandas `freq` parameter from `'H'` to `'h'`. Added `pytest.ini` to filter deprecation warnings from external libraries (cdsapi, pkg_resources, metpy).

### Changed

- **CDS API Workflow**: The CDS API download process is now more robust and reliable for large datasets spanning multiple days. Users will see progress updates for each day being downloaded.

## [1.1.2] - 2026-01-21

### Fixed / Improved

- **Environment & Installation**: Added `environment.yml` and clarified recommended Conda-based installation method for developers. Ensured `pip install -e .` is documented and included in developer install steps.
- **Package installation**: Confirmed `setup.py` includes the root module (`py_modules=['lorenzcycletoolkit']`) and bumped package version to `1.1.2` to reflect this release.
- **Documentation**: Rewrote `installation.rst` and updated `README.md` to keep the README concise while directing full install and usage instructions to the documentation. Fixed RST formatting issues and rebuilt docs to resolve related warnings.
- **Docs version**: Updated Sphinx `release` to `1.1.2` in `docs/source/conf.py`.
- 
## [1.1.1] - 2026-01-20

### Added

- **--time-resolution Flag**: Added `--time-resolution` command-line flag to control the temporal resolution of ERA5 data downloads via CDS API. Default is 3 hours to avoid exceeding API cost limits.
- **environment.yml**: Added `environment.yml` file to simplify conda environment creation for users.

### Fixed

- **CDS API Cost Limits**: Fixed issue where automatic 1-hour temporal resolution detection often exceeded CDS API cost limits. Now defaults to 3-hour resolution with user-configurable options (3, 6, 12, or 24 hours).
- **Test Suite**: Updated `test_cdsapi.py` to include `time_resolution` parameter in test fixtures and adjusted time step calculation test expectations.
- **Dependencies**: Fixed `docutils` version conflict between `sphinx` and `sphinx_rtd_theme` by downgrading to `docutils==0.20.1`.
- **Package Installation**: Fixed `setup.py` to properly include `lorenzcycletoolkit.py` module by adding `py_modules=['lorenzcycletoolkit']` parameter. This resolves the `ModuleNotFoundError` when installing the package.

### Changed

- **CDS API Temporal Resolution**: Changed from automatically detecting temporal resolution from track file to using a configurable `--time-resolution` parameter with a safe default of 3 hours.
- **Documentation Theme**: Updated documentation theme from Alabaster to ReadTheDocs (`sphinx_rtd_theme`) for better readability and navigation.
- **Repository URL**: Updated installation documentation to reflect new repository URL: `https://github.com/daniloceano/LorenzCycleToolkit`.
- **Documentation Version**: Updated documentation version to 1.1.1 in `conf.py`.

### Documentation

- Completely rewrote `installation.rst` with clearer instructions for:
  - Installing from PyPI for regular users
  - Installing from source with Conda (recommended for developers)
  - Installing from source with venv
  - Added troubleshooting section for common installation issues
- Updated `usage.rst` with information about the new `--time-resolution` flag and cost limit guidelines.
- Updated `cdsapi_example.rst` with detailed examples and guidelines for choosing appropriate temporal resolutions.
- Updated `CDSAPI_USAGE.md` with troubleshooting tips for cost limit errors and temporal resolution examples.
- Updated `README.md` with comprehensive installation instructions and examples for both users and developers.
- Updated `installation.rst` with corrected repository URLs.
- Added `sphinx_rtd_theme` to `requirements.txt`.


## [1.1.0] - 2025-11-18

### Updated

- **get_cdsapi_data Function**: Modernized the `get_cdsapi_data` function to align with the latest CDS API syntax (v0.7.x). This includes improved error handling, logging, and support for new request parameters such as `retry_max`.

### Added

- **Test Suite for get_cdsapi_data**: Created a comprehensive test suite (`tests/test_cdsapi.py`) to validate the `get_cdsapi_data` function. The suite includes 13 test cases covering parameter validation, error handling, and edge cases.



## [1.0.9] - 2024-10-21

### Fixed

- **LEC Fixed Framework**: Fixed a bug where the vertical levels output file was not saving the correct pressure levels as columns. Now, the output correctly reflects the pressure levels.
  
- **LEC Moving Framework**: Resolved a bug that prevented the framework from correctly handling maximum vorticity for the Northern Hemisphere track files. The framework now properly processes maximum vorticity values for the Northern Hemisphere.

- **Plot Periods**: Enhanced the functionality to check if the track is for the Northern Hemisphere. If so, the vorticity values are multiplied by `-1` to ensure correct behavior in the CycloPhaser program.

### Added

- **Tools.py**: Added a preprocessing stage to remove the new dimensions (`expver` and `number`) introduced in recent ERA5 datasets, ensuring compatibility with previous workflows.

## [1.0.8] - 2024-09-05

### Added
- Added badges for JOSS, Zenodo, documentation status, and license (GPL v3) to the README.

## [1.0.7] - 2024-09-05

### Fixed
- Fixed the broken link to the **Contributing** section in the `README.md` file.

## [1.0.6] - 2024-09-05

### Added
- Expanded the **Contributing** section to include detailed guidelines on branch naming conventions, code style (autopep8, isort, flake8), testing, and PR submission.
- Added a suggestion for running tests (`pytest`) and code formatting/linting commands (`autopep8`, `isort`, `flake8`) before submitting pull requests.

### Changed
- Improved the **Contributing** section in the documentation to provide better clarity on contribution processes and coding standards.
- Updated the **Contributing** section to include information about the Continuous Integration (CI) pipeline and the required coding practices.
- Removed "flags" from `index.rst` as the content was moved to the **Usage** section in version 1.0.5 but not reflected in the index.

### Fixed
- Resolved a merge conflict in the `README.md` file caused during cherry-picking of changes from the `joss-submission` branch to the `main` branch.


## [1.0.5] - 2024-09-04

### Added
- Updated README to clarify that installation instructions are included in the documentation.
- Revised README and documentation to reflect the GNU license, correcting the previously stated MIT license.
- Moved the "Flags" section to the "Usage" section to improve tutorial clarity and make the process more linear for users.
- Added a note in the documentation explaining that results generated using the LorenzCycleToolkit will be referenced in a future journal publication to aid users in interpreting outputs.

### Fixed
- Corrected links to the Contributing Guide and License in the README and documentation.

## [1.0.4] - 2024-08-28

### Added
- Added a `support` section in the documentation to provide guidelines for users seeking help or reporting issues. Users are encouraged to use the GitHub repository for reporting, with the option to contact via email if needed.

### Updated
- Updated the "Moving (Semi-Lagrangian) Framework Example" section in the documentation to include information about using the default track file for testing and the requirement to overwrite the `inputs/track` file with custom track data.
- Clarified that the preferred method for running the LorenzCycleToolkit is through command line arguments executed from the top-level directory of the project.

### Fixed
- Improved clarity in the documentation by specifying the appropriate directories and command usage for various examples and tutorials.

## [1.0.3] - 2024-07-27
### Bug Fixes
- Deploying documentation

## [1.0.2] - 2024-07-27
### Added
- Improved CI/CD pipeline by ensuring up-to-date gh-pages branch before deployment.
- Updated installation instructions to include pip installation from PyPI.

## [1.0.1] - 2024-07-27
### Added
- Added `CHANGELOG.md` file to document changes.
- Published project to PyPI.
- Integrated CircleCI for continuous integration and deployment.

## [1.0.0] - 2024-07-16
### Initial Release
- Initial release of LorenzCycleToolkit.
- Implemented core functionalities for calculating the Lorenz Energy Cycle (LEC) in specific atmospheric regions.
- Added support for both Eulerian and Semi-Lagrangian frameworks.
- Provided options for fixed and moving computational domains.
- Generated detailed results in CSV format for energy, conversion, and generation terms.
- Included automatic plot generation for spatial and temporal analysis.
- Utilized Xarray, Dask, Pandas, and MetPy for data handling and computations.
- Comprehensive documentation and usage examples.

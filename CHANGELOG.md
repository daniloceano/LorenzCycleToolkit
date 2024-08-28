# Changelog

All notable changes to this project will be documented in this file.

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

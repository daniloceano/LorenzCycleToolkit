# LorenzCycleToolkit

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13765959.svg)](https://doi.org/10.5281/zenodo.13765959)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.07139/status.svg)](https://doi.org/10.21105/joss.07139)
![Documentation Status](https://img.shields.io/badge/docs-passing-brightgreen.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## Overview

The LorenzCycleToolkit is a tool for calculating the Lorenz Energy Cycle (LEC) in specific atmospheric regions. Introduced by Edward Lorenz in 1965, the LEC estimates atmospheric energy, including zonal and eddy components of Kinetic and Available Potential Energy, and the conversions between these forms.

### Importance

- **Climate Studies:** Improves understanding of energy exchanges, aiding climate predictions.
- **Weather Prediction:** Comparison tool for weather models and forecasts.
- **Research and Education:** Provides a structured approach for analyzing atmospheric energy flows.
- **Environmental Diagnostics:** Serves as a diagnostic tool for atmospheric dynamics.

### Applications

- **Synoptic-Scale Phenomena:** Study extratropical cyclones and convergence zones.
- **Mesoscale Phenomena:** Analyze energy transformations in tropical cyclones.
- **Model Comparison:** Compare model dynamics with reanalysis data.

### Features

- **Flexible Frameworks:** Analyze fixed regions using Eulerian or Semi-Lagrangian frameworks.
- **User-Friendly:** Simple command-line interface without needing to program scripts.
- **Customizable Inputs:** Define custom domains and variable configurations.
- **Visualization Options:** Generate detailed plots to visualize the energy cycle components.

## Installation

For detailed installation instructions and troubleshooting, please refer to the [Installation Guide](https://daniloceano.github.io/LorenzCycleToolkit/installation).

## Documentation

For detailed documentation and usage instructions, visit the [LorenzCycleToolkit Documentation](https://daniloceano.github.io/LorenzCycleToolkit/).

## Contributing
Contributions are welcome! Please see the [contributing guide](https://daniloceano.github.io/LorenzCycleToolkit/contributing) for more details.

### Continuous Integration and Deployment

This project uses GitHub Actions and CircleCI for Continuous Integration (CI) and Continuous Deployment (CD). 

- **CI Pipeline with GitHub Actions**: The CI pipeline is triggered on each push and pull request to the repository. It runs the following steps:
  - Sets up the Python environment.
  - Installs dependencies.
  - Builds and installs the package wheel, and checks that the installed package is the checked-out code.
  - Lints code using flake8 (read only: no tool rewrites the code before the tests).
  - Provides the Natural Earth data used by the map tests, verified against their sha256.
  - Runs tests using pytest.

- **CI Pipeline with CircleCI**: CircleCI builds and tests the package on each push. The pipeline runs the following steps:
  - Builds the package.
  - Checks if the wheel file exists.
  - Installs dependencies.
  - Installs the package and checks that it is the checked-out code.
  - Provides the verified Natural Earth data.
  - Runs tests using pytest.

- **Publishing**: Nothing is published automatically; pushes and pull requests only build and test. Publishing is a separate, manual step:
  - TestPyPI (from the `develop` branch) and PyPI (from the `main` branch): a CircleCI pipeline triggered explicitly with the parameter `publish_testpypi` or `publish_pypi`.
  - Documentation on GitHub Pages: a GitHub Actions workflow run explicitly ("Run workflow").

You can view the CI/CD configuration in the following files:
- [python-app.yml](.github/workflows/python-app.yml)
- [deploy.yml](.github/workflows/deploy.yml) and [gh-pages.yml](.github/workflows/gh-pages.yml)
- [config.yml](.circleci/config.yml)


This project is licensed under the GNU License Version 3. See the [LICENSE](LICENSE) file for details.
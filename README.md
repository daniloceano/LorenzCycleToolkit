# LorenzCycleToolkit

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

For detailed installation instructions, please refer to the [Installation Guide](https://daniloceano.github.io/LorenzCycleToolkit/installation) in the official documentation.

## Documentation

For detailed documentation and usage instructions, visit the [LorenzCycleToolkit Documentation](https://daniloceano.github.io/LorenzCycleToolkit/).

## Contributing
Contributions are welcome! Please see the [contributing guide](CONTRIBUTING.md) for more details.

### Continuous Integration and Deployment

This project uses GitHub Actions and CircleCI for Continuous Integration (CI) and Continuous Deployment (CD). 

- **CI Pipeline with GitHub Actions**: The CI pipeline is triggered on each push and pull request to the repository. It runs the following steps:
  - Sets up the Python environment.
  - Installs dependencies.
  - Formats code using autopep8.
  - Sorts imports using isort.
  - Lints code using flake8.
  - Runs tests using pytest.

- **CI Pipeline with CircleCI**: CircleCI is used to build, test, and publish the package. The pipeline runs the following steps:
  - Builds the package.
  - Checks if the wheel file exists.
  - Installs dependencies.
  - Installs the package.
  - Runs tests using pytest.
  - Publishes the package to TestPyPI for the `develop` branch.
  - Publishes the package to PyPI for the `main` branch.

- **CD Pipeline**: The CD pipeline deploys the documentation to GitHub Pages whenever changes are pushed to the `main` branch.

You can view the CI/CD configuration in the following files:
- [python-app.yml](.github/workflows/python-app.yml)
- [deploy-docs.yml](.github/workflows/deploy-docs.yml)
- [config.yml](.circleci/config.yml)


This project is licensed under the GNU License Version 3. See the [LICENSE](LICENSE) file for details.
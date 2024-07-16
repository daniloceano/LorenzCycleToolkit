# LorenzCycleToolkit

## Overview
The LorenzCycleToolkit is a tool for calculating the Lorenz Energy Cycle (LEC) in specific atmospheric regions. The LEC, introduced by Edward Lorenz in 1965, is an analytical framework used to estimate atmospheric energy, including zonal and eddy components of Kinetic and Available Potential Energy, and conversions between these forms.

## Features
- Fixed (Eulerian) Framework
- Moving (Semi-Lagrangian) Framework
- Interactive Domain Selection
- Comprehensive energy calculations and visualizations

## Installation
To install the required dependencies, use the provided `requirements.txt` file:
```sh
pip install -r requirements.txt
```

## Basic Usage
Run the toolkit with a sample NetCDF file:
```sh
python lorenzcycletoolkit.py samples/testdata_ERA5.nc -r -f
```

For more detailed usage instructions, including configuration options and examples, please refer to the documentation.

## Documentation
Comprehensive documentation is available in the `docs/` directory, covering:
- Installation and setup
- Configuration files
- Detailed examples and tutorials

## Contributing
Contributions are welcome! Please see the [contributing guide](CONTRIBUTING.md) for more details.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

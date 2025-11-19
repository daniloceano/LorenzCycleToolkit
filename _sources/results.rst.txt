Results
=======

The `LorenzCycleToolkit` generates a directory structure to store the results of the Lorenz Energy Cycle (LEC) calculations. This structure includes CSV files containing the calculated terms, a log file, and optional figures if the plotting flag is enabled. Here's an example of the results directory structure:

.. code-block:: bash

    LEC_Results/
    ├── testdata_ERA5_fixed
    │   ├── Figures
    │   │   ├── LEC
    │   │   ├── boxplots
    │   │   ├── hovmollers
    │   │   ├── timeseries
    │   ├── log.testdata_ERA5
    │   ├── results_vertical_levels
    │   │   ├── Ae_level.csv
    │   │   ├── Az_level.csv
    │   │   ├── Ca_1_level.csv
    │   │   ├── Ck_1_level.csv
    │   │   ├── Cz_1_level.csv
    │   │   ├── Ge_level.csv
    │   │   ├── Gz_level.csv
    │   │   ├── Ke_level.csv
    │   │   ├── Kz_level.csv
    │   ├── testdata_ERA5_fixed_results.csv

Each subdirectory within the `LEC_Results` directory corresponds to a different dataset or analysis run. The structure of each subdirectory includes:

- **Figures/**: Contains figures generated when the `-p` flag is used.
  - **LEC/**: Contains LEC diagrams for each time step.
  - **boxplots/**: Contains boxplots of the LEC terms for each time step.
  - **hovmollers/**: Contains Hovmöller diagrams of the LEC terms.
  - **timeseries/**: Contains time series plots of the LEC terms.
  
- **log.*: This log file contains detailed information about the processing and calculations performed during the analysis.

- **results_vertical_levels/**: Contains CSV files for each LEC term at various vertical levels. Each conversion term is split into multiple sub-terms. For example:
  - `Ae_level.csv`: Contains the values of the eddy component of Available Potential Energy (Ae) at different levels.
  - `Az_level.csv`: Contains the values of the zonal component of Available Potential Energy (Az) at different levels.
  - `Ca_1_level.csv`, `Ca_2_level.csv`: Contains the values of the conversion term Ca at different levels, split into multiple sub-terms.
  - `Ck_1_level.csv`, `Ck_2_level.csv`, `Ck_3_level.csv`, `Ck_4_level.csv`, `Ck_5_level.csv`: Contains the values of the conversion term Ck at different levels, split into multiple sub-terms.
  - `Cz_1_level.csv`, `Cz_2_level.csv`: Contains the values of the conversion term Cz at different levels, split into multiple sub-terms.
  - `Ge_level.csv`: Contains the values of the generation term Ge at different levels.
  - `Gz_level.csv`: Contains the values of the generation term Gz at different levels.
  - `Ke_level.csv`: Contains the values of the eddy component of Kinetic Energy (Ke) at different levels.
  - `Kz_level.csv`: Contains the values of the zonal component of Kinetic Energy (Kz) at different levels.

- **testdata_*_fixed_results.csv**: This CSV file contains the overall results of the LEC calculations for the given dataset.

### Log File

The log file provides detailed information about the processing and calculations performed during the analysis. It includes information about the input data, the configuration used for the analysis, and any errors or warnings encountered during the process.
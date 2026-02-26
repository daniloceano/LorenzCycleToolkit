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

- **Figures/**: Contains figures generated when the ``-p`` flag is used.
  
  - **LEC/**: Contains LEC diagrams for each time step.
  - **boxplots/**: Contains boxplots of the LEC terms for each time step.
  - **hovmollers/**: Contains Hovmöller diagrams of the LEC terms.
  - **timeseries/**: Contains time series plots of the LEC terms.
  
- **log.*** (log file): Contains detailed information about the processing and calculations performed during the analysis.

- **results_vertical_levels/**: Contains CSV files for each LEC term at various vertical levels. Each conversion term may be split into multiple sub-terms for detailed analysis. For example:
  
  - ``Ae_level.csv``: Contains the values of the eddy component of Available Potential Energy (Ae) at different levels.
  - ``Az_level.csv``: Contains the values of the zonal component of Available Potential Energy (Az) at different levels.
  - ``Ca_1_level.csv``, ``Ca_2_level.csv``: Contains the values of the conversion term Ca at different levels, split into multiple sub-terms for detailed horizontal and vertical advection analysis.
  - ``Ce_level.csv``: Contains the values of the conversion term Ce at different levels.
  - ``Ck_1_level.csv``, ``Ck_2_level.csv``, ``Ck_3_level.csv``, ``Ck_4_level.csv``, ``Ck_5_level.csv``: Contains the values of the conversion term Ck at different levels, split into five sub-terms representing different physical processes (horizontal shear, meridional momentum transport, Coriolis effect, vertical shear).
  - ``Cz_1_level.csv``, ``Cz_2_level.csv``: Contains the values of the conversion term Cz at different levels, split into horizontal and vertical components.
  - ``Ge_level.csv``: Contains the values of the generation term Ge (eddy) at different levels.
  - ``Gz_level.csv``: Contains the values of the generation term Gz (zonal) at different levels.
  - ``Ke_level.csv``: Contains the values of the eddy component of Kinetic Energy (Ke) at different levels.
  - ``Kz_level.csv``: Contains the values of the zonal component of Kinetic Energy (Kz) at different levels.
  
  **Note**: The split of conversion terms into sub-components allows for detailed analysis of the physical processes contributing to energy transformations at each vertical level.

Main Results File
~~~~~~~~~~~~~~~~~

The main CSV file (e.g., ``testdata_ERA5_fixed_results.csv``) contains the integrated LEC terms over the entire atmospheric column for each timestep. This file provides a time series of the energy cycle and is the primary output for analysis.

**File Columns**:

- **Time**: Timestamp of each analysis step
- **Az, Ae, Kz, Ke**: Energy reservoirs (J/m²), integrated over all pressure levels:
  
  - ``Az``: Zonal Available Potential Energy
  - ``Ae``: Eddy Available Potential Energy
  - ``Kz``: Zonal Kinetic Energy
  - ``Ke``: Eddy Kinetic Energy

- **Ca, Ce, Cz, Ck**: Conversion terms (W/m²), representing energy transformations:
  
  - ``Ca``: Conversion between Az and Ae
  - ``Ce``: Conversion from Ae to Ke
  - ``Cz``: Conversion from Az to Kz
  - ``Ck``: Conversion from Ke to Kz

- **BAz, BAe, BKz, BKe**: Boundary terms (W/m²), energy transport across domain boundaries (only for regional analyses)

- **Generation and Dissipation terms** (W/m²):
  
  When using ``-r/--residuals`` flag:
  
  - ``RGz, RGe``: Residual generation terms (include generation + numerical errors)
  - ``RKz, RKe``: Residual dissipation terms (include boundary flux + dissipation + numerical errors)
  
  When **not** using ``-r`` flag:
  
  - ``Gz, Ge``: Generation terms (computed from diabatic heating)
  - ``Dz, De``: Dissipation terms (computed from friction)
  - ``BΦz, BΦe``: Geopotential flux terms

**Usage Example**:

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   
   # Load results
   df = pd.read_csv('LEC_Results/testdata_ERA5_fixed/testdata_ERA5_fixed_results.csv', 
                    parse_dates=['Time'], index_col='Time')
   
   # Plot energy reservoirs
   fig, ax = plt.subplots(figsize=(12, 6))
   df[['Az', 'Ae', 'Kz', 'Ke']].plot(ax=ax)
   ax.set_ylabel('Energy (J/m²)')
   ax.set_title('Lorenz Energy Cycle - Energy Reservoirs')
   plt.show()

Log File
~~~~~~~~

The log file (``log.<filename>``) provides detailed information about the processing and calculations performed during the analysis. It includes:

- Input data information and validation checks
- Configuration used for the analysis
- Processing steps and timing information
- Any errors or warnings encountered during execution
- Download progress (when using ``--cdsapi``)

Use the ``-v/--verbosity`` flag for more detailed debug-level logging.
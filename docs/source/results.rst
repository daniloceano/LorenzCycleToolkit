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
    │   │   ├── C_overturning_level.csv
    │   │   ├── Cz_1_level.csv
    │   │   ├── Ge_level.csv
    │   │   ├── Gz_level.csv
    │   │   ├── Ke_level.csv
    │   │   ├── Kz_level.csv
    │   │   ├── M_level.csv
    │   ├── run_status.json
    │   ├── testdata_ERA5_fixed_results.csv

Each subdirectory within the `LEC_Results` directory corresponds to a different dataset or analysis run. The structure of each subdirectory includes:

- **Figures/**: Contains figures generated when the ``-p`` flag is used.
  
  - **LEC/**: Contains LEC diagrams for each time step.
  - **boxplots/**: Contains boxplots of the LEC terms for each time step.
  - **hovmollers/**: Contains Hovmöller diagrams of the LEC terms.
  - **timeseries/**: Contains time series plots of the LEC terms.
  
- **log.*** (log file): Contains detailed information about the processing and calculations performed during the analysis.

- **run_status.json**: Status of the run: whether the computation completed and whether each requested product was generated (see `Run Status File`_).

- **results_vertical_levels/**: Contains CSV files for each LEC term at various vertical levels. Each conversion term may be split into multiple sub-terms for detailed analysis. For example:
  
  - ``Ae_level.csv``: Contains the values of the eddy component of Available Potential Energy (Ae) at different levels.
  - ``Az_level.csv``: Contains the values of the zonal component of Available Potential Energy (Az) at different levels.
  - ``Ca_1_level.csv``, ``Ca_2_level.csv``: Contains the values of the conversion term Ca at different levels, split into multiple sub-terms for detailed horizontal and vertical advection analysis.
  - ``Ce_level.csv``: Contains the values of the conversion term Ce at different levels.
  - ``Ck_1_level.csv``, ``Ck_2_level.csv``, ``Ck_3_level.csv``, ``Ck_4_level.csv``, ``Ck_5_level.csv``: Contains the values of the conversion term Ck at different levels, split into five sub-terms representing different physical processes (horizontal shear, meridional momentum transport, Coriolis effect, vertical shear).
  - ``Cz_1_level.csv``, ``Cz_2_level.csv``: Contains the values of the conversion term Cz at different levels, split into horizontal and vertical components.
  - ``C_overturning_level.csv``: Contains the domain-mean overturning pressure work at each level. This is a diagnostic and does not enter the budget.
  - ``Ge_level.csv``: Contains the values of the generation term Ge (eddy) at different levels.
  - ``Gz_level.csv``: Contains the values of the generation term Gz (zonal) at different levels.
  - ``Ke_level.csv``: Contains the values of the eddy component of Kinetic Energy (Ke) at different levels.
  - ``Kz_level.csv``: Contains the values of the zonal component of Kinetic Energy (Kz) at different levels.
  - ``M_level.csv``: Contains the mass-continuity residual M(p), in s⁻¹, which sets the numerical noise floor of the budget.
  
  **Note**: The split of conversion terms into sub-components allows for detailed analysis of the physical processes contributing to energy transformations at each vertical level.

  **Note**: Hovmöller diagrams and boxplots are produced for the energy, conversion and generation/dissipation groups only, so ``C_overturning_level.csv`` and ``M_level.csv`` are archived for analysis but are not plotted at vertical levels. Their column integrals do appear in the time series figures.

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

- **BΦZ, BΦE**: Boundary pressure work terms (W/m²), always exported

- **C_overturning**: Domain-mean overturning pressure work (W/m²). A diagnostic; it does **not** enter ``RGz`` or ``RKz``.

- **M**: Mass-continuity residual (kg·m⁻²·s⁻¹). Sets the numerical noise floor of the budget; it does **not** enter the residuals.

- **Generation and Dissipation terms** (W/m²):
  
  When using ``-r/--residuals`` flag:
  
  - ``RGz, RGe``: Residual generation terms (include generation + numerical errors)
  - ``RKz, RKe``: Residual dissipation terms (include boundary flux + dissipation + numerical errors)
  
  When **not** using ``-r`` flag:
  
  - ``Dz, De``: Dissipation terms (computed from friction)

  ``Gz, Ge`` (generation, computed from diabatic heating) are exported in both
  modes.

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

Run Status File
~~~~~~~~~~~~~~~

``run_status.json`` is written by every command-line run, first with the state ``running`` before the data are read and then with the final state. It is an additional file: the CSV files are not changed. The file is replaced in one step (temporary file and rename), so a run that is killed leaves ``running``, never a complete-looking status. Schema version 1:

- ``schema_version`` (``1``), ``lectk_version``, ``command`` (arguments, method, input file), ``started_utc``, ``finished_utc``, ``log_file``.
- ``status``: ``success`` (computation complete, every requested product generated), ``partial`` (computation complete, at least one requested figure failed), ``failed`` (the computation or the run failed) or ``running``.
- ``exit_code``: the exit status of the command line (0, 3 or 1; ``null`` while running or after an interruption). ``allow_partial_plots`` and ``plots_requested`` record the options.
- ``scientific``: ``status`` of the computation and the ``error`` (type and message) when it failed. ``interrupted`` and ``error`` describe a run stopped by an exception.
- ``summary``: number of products in each state.
- ``products``: one entry per product with ``name``, ``kind`` (``scientific``, ``figure`` or ``diagnostic``), ``required``, ``status``, ``expected`` (files the product must write), ``written`` (files written by this run), ``missing``, ``reason``, ``error``, ``failures`` (one per failed call, e.g. per time step), ``depends_on`` and ``consequences``. Paths are relative to the results directory.

Product states: ``success`` (requested and every expected file written by this run), ``failed`` (requested and not, or not completely, generated), ``unavailable`` (not applicable to the data for a documented reason; not used in this version), ``not_requested`` (optional product without ``-p``). ``pending`` appears only while the run is in progress.

Required products: ``results_csv``, ``vertical_level_csvs`` and, in the moving framework, ``trackfile``. Optional products with ``-p`` — fixed framework: ``timeseries``, ``box_limits_map``, ``boxplots``, ``hovmoller``, ``lec_cycle``; moving framework: ``domain_maps`` (one map per time step), ``track_map``, ``zeta_hgt_timeseries``, ``timeseries``, ``hovmoller``, ``boxplots``, ``periods`` (cyclophaser life-cycle periods and ``periods.csv``), ``lec_cycle``, ``lec_periods`` and ``lps`` (these two need ``periods``).

Example (moving framework, a series too short for the cyclophaser filter):

.. code-block:: json

   {
     "schema_version": 1,
     "status": "partial",
     "exit_code": 3,
     "scientific": {"status": "success", "error": null},
     "products": [
       {"name": "results_csv", "kind": "scientific", "required": true, "status": "success",
        "written": ["testdata_ERA5_track_results.csv"]},
       {"name": "periods", "kind": "diagnostic", "required": false, "status": "failed",
        "reason": "ValueError: polyorder must be less than window_length.",
        "consequences": ["lec_periods", "lps"]},
       {"name": "lps", "kind": "figure", "required": false, "status": "failed",
        "depends_on": ["periods"],
        "reason": "not generated: requires 'periods', which failed (ValueError: polyorder must be less than window_length.)"}
     ]
   }

(Abridged: the file lists every product and every field.)

Use the ``-v/--verbosity`` flag for more detailed debug-level logging.
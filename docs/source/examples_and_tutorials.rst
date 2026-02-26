Examples and Tutorials
======================

This page provides step-by-step tutorials progressing from basic to advanced usage. Each tutorial includes objectives, estimated time, and detailed instructions.

Tutorial 1: Basic Fixed Framework Analysis
-------------------------------------------

**Objective**: Calculate the Lorenz Energy Cycle for a fixed region using sample data

**Estimated time**: 10 minutes

**Prerequisites**: LorenzCycleToolkit installed and activated

**Step 1: Prepare the namelist file**

The ``namelist`` file specifies the variable names and units in your NetCDF file. For this example using NCEP-R2 data, copy the appropriate preset:

.. code-block:: bash

   cp inputs/namelist_NCEP-R2 inputs/namelist

The namelist should look like this:

.. code-block:: text

   ;standard_name;Variable;Units
   Air Temperature;air_temperature;TMP_2_ISBL;K
   Geopotential Height;geopotential_height;HGT_2_ISBL;m
   Omega Velocity;omega;V_VEL_2_ISBL;Pa/s
   Eastward Wind Component;eastward_wind;U_GRD_2_ISBL;m/s
   Northward Wind Component;northward_wind;V_GRD_2_ISBL;m/s
   Longitude;;lon_2
   Latitude;;lat_2
   Time;;initial_time0_hours
   Vertical Level;;lv_ISBL3

**Step 2: Prepare the box_limits file**

Create or edit ``inputs/box_limits`` to define your spatial domain. For this example, we'll analyze a region over the South Atlantic:

.. code-block:: bash

   cat > inputs/box_limits << EOF
   min_lon;-60
   max_lon;-30
   min_lat;-42.5
   max_lat;-17.5
   EOF

**Step 3: Run the analysis**

Execute the LorenzCycleToolkit with the fixed framework. The ``-r`` flag computes generation and dissipation as residuals, and ``-p`` generates plots:

.. code-block:: bash

   python lorenzcycletoolkit.py samples/testdata_NCEP-R2.nc -r -f -p

**Step 4: Check the results**

Results are saved in ``LEC_Results/testdata_NCEP-R2_fixed/``:

.. code-block:: text

   LEC_Results/testdata_NCEP-R2_fixed/
   ├── testdata_NCEP-R2_fixed_results.csv    # Main time series results
   ├── log.testdata_NCEP-R2                   # Processing log
   ├── results_vertical_levels/               # Level-by-level results
   │   ├── Az_level.csv
   │   ├── Ae_level.csv
   │   └── ...
   └── Figures/                               # Visualizations
       ├── LEC/                               # LEC diagrams
       ├── boxplots/                          # Statistical distributions
       ├── hovmollers/                        # Vertical-time diagrams
       └── timeseries/                        # Time series plots

**Step 5: Analyze the output**

Open the main results file to examine the energy cycle:

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   
   # Load results
   df = pd.read_csv('LEC_Results/testdata_NCEP-R2_fixed/testdata_NCEP-R2_fixed_results.csv',
                    parse_dates=['Time'], index_col='Time')
   
   # Display first few rows
   print(df.head())
   
   # Quick visualization
   fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
   
   # Energy reservoirs
   df[['Az', 'Ae', 'Kz', 'Ke']].plot(ax=ax1)
   ax1.set_ylabel('Energy (J/m²)')
   ax1.set_title('Energy Reservoirs')
   ax1.legend(loc='best')
   ax1.grid(True)
   
   # Conversion terms
   df[['Ca', 'Ce', 'Cz', 'Ck']].plot(ax=ax2)
   ax2.set_ylabel('Conversion (W/m²)')
   ax2.set_title('Conversion Terms')
   ax2.legend(loc='best')
   ax2.grid(True)
   
   plt.tight_layout()
   plt.show()

**Understanding the results**:

- **Positive Az, Ae**: Indicates potential energy stored in temperature gradients
- **Positive Kz, Ke**: Indicates kinetic energy of atmospheric motions
- **Positive Ca**: Energy conversion from zonal to eddy APE (baroclinic instability)
- **Positive Ce**: Energy conversion from eddy APE to eddy KE (cyclone intensification)
- **Positive Ck**: Energy conversion from eddy to zonal KE (barotropic conversion)
- **Positive Cz**: Energy conversion from zonal APE to zonal KE

Tutorial 2: Moving Framework with Track File
---------------------------------------------

**Objective**: Analyze a moving atmospheric system (e.g., cyclone) using a track file

**Estimated time**: 15 minutes

**Prerequisites**: Tutorial 1 completed

**Step 1: Prepare the track file**

Create a track file defining your system's center over time. The file must use semicolons as delimiters:

.. code-block:: bash

   cat > inputs/track << EOF
   time;Lat;Lon
   2005-08-08 00:00:00;-22.5;-45.0
   2005-08-08 06:00:00;-22.6;-44.9
   2005-08-08 12:00:00;-22.7;-44.8
   2005-08-08 18:00:00;-22.8;-44.7
   2005-08-09 00:00:00;-22.9;-44.6
   EOF

For this tutorial, you can use the provided test track:

.. code-block:: bash

   python lorenzcycletoolkit.py samples/testdata_NCEP-R2.nc -r -t -p --trackfile inputs/track_testdata_NCEP-R2

**Step 2: Configure namelist (same as Tutorial 1)**

.. code-block:: bash

   cp inputs/namelist_NCEP-R2 inputs/namelist

**Step 3: Run the analysis**

.. code-block:: bash

   python lorenzcycletoolkit.py samples/testdata_NCEP-R2.nc -r -t -p --trackfile inputs/track

The moving framework automatically follows the system center, using a domain that moves with the tracked feature.

**Step 4: Examine track-specific output**

The moving framework produces additional visualizations in ``Figures/``:

- **Track maps**: Show the system trajectory with domain boxes
- **Time-evolving LEC diagrams**: Show how energy cycle changes as system moves

**Differences from fixed framework**:

- Domain follows the system (Semi-Lagrangian approach)
- Better for studying individual cyclones or moving systems
- Reduces influence of features entering/leaving domain
- Boundary terms represent actual system interaction with environment

Tutorial 3: Downloading and Analyzing ERA5 Data
------------------------------------------------

**Objective**: Download ERA5 data using CDS API and perform LEC analysis

**Estimated time**: 30-60 minutes (depending on download speed)

**Prerequisites**: CDS API account and credentials configured

**Step 1: Setup CDS API credentials**

Create or edit ``~/.cdsapirc``:

.. code-block:: bash

   cat > ~/.cdsapirc << EOF
   url: https://cds.climate.copernicus.eu/api
   key: YOUR-UID:YOUR-API-KEY
   EOF

Replace ``YOUR-UID:YOUR-API-KEY`` with your credentials from https://cds.climate.copernicus.eu/user

**Step 2: Prepare track file**

Create a track file for the period you want to analyze:

.. code-block:: bash

   cat > inputs/track_era5_example << EOF
   time;Lat;Lon
   2020-06-15 00:00:00;-30.0;-50.0
   2020-06-15 06:00:00;-30.5;-49.5
   2020-06-15 12:00:00;-31.0;-49.0
   2020-06-15 18:00:00;-31.5;-48.5
   2020-06-16 00:00:00;-32.0;-48.0
   EOF

**Step 3: Run with automatic download**

The ``--cdsapi`` flag triggers automatic ERA5 download. The namelist is configured automatically:

.. code-block:: bash

   python lorenzcycletoolkit.py era5_cyclone_2020.nc -t -r -p --cdsapi --trackfile inputs/track_era5_example

**What happens**:

1. Reads track file to determine time range and spatial domain
2. Adds 15° buffer around track boundaries
3. Downloads ERA5 data day-by-day from CDS
4. Saves combined data to ``era5_cyclone_2020.nc``
5. Performs LEC analysis automatically

**Step 4: Re-run without re-downloading**

If the NetCDF file already exists, the toolkit skips the download:

.. code-block:: bash

   # This will use existing era5_cyclone_2020.nc
   python lorenzcycletoolkit.py era5_cyclone_2020.nc -t -r -p --cdsapi --trackfile inputs/track_era5_example

**Controlling temporal resolution**:

.. code-block:: bash

   # Download with 6-hour resolution (faster, less data)
   python lorenzcycletoolkit.py era5_6h.nc -t -r -p --cdsapi --time-resolution 6 --trackfile inputs/track_era5_example

**For more details**, see :doc:`cdsapi_example`.

Tutorial 4: Interactive Domain Selection
-----------------------------------------

**Objective**: Manually select analysis domain at each timestep

**Estimated time**: Variable (depends on number of timesteps)

**When to use**: When you want to manually adjust the domain based on visual inspection of meteorological fields at each time

**Step 1: Prepare namelist**

.. code-block:: bash

   cp inputs/namelist_NCEP-R2 inputs/namelist

**Step 2: Run with interactive mode**

.. code-block:: bash

   python lorenzcycletoolkit.py samples/testdata_NCEP-R2.nc -r -c

**Step 3: Select domain interactively**

For each timestep, a map appears showing meteorological fields. Click on the map to define:

1. First click: Northwestern corner
2. Second click: Southeastern corner

.. image:: https://user-images.githubusercontent.com/56005607/214922008-5b7c094f-c160-4415-a528-07cc58730827.png
   :width: 350

The toolkit then computes LEC for that domain and proceeds to the next timestep.

**Use cases**:

- Analyzing systems with irregular shapes
- Studies requiring domain adjustment based on synoptic patterns
- Exploratory analysis to understand system evolution

Advanced Example: Batch Processing Multiple Cases
--------------------------------------------------

**Objective**: Process multiple cyclone cases efficiently

**Estimated time**: Variable

**Scenario**: You have multiple cyclone track files and want to process them all

**Step 1: Organize your data**

.. code-block:: text

   inputs/
   ├── tracks/
   │   ├── cyclone_20200115.csv
   │   ├── cyclone_20200223.csv
   │   └── cyclone_20200315.csv
   └── data/
       ├── cyclone_20200115_ERA5.nc
       ├── cyclone_20200223_ERA5.nc
       └── cyclone_20200315_ERA5.nc

**Step 2: Create a batch processing script**

.. code-block:: bash

   #!/bin/bash
   # batch_process.sh
   
   # Prepare namelist
   cp inputs/namelist_ERA5-cdsapi inputs/namelist
   
   # Process each case
   for case in cyclone_20200115 cyclone_20200223 cyclone_20200315; do
       echo "Processing $case..."
       python lorenzcycletoolkit.py \
           inputs/data/${case}_ERA5.nc \
           -t -r -p \
           --trackfile inputs/tracks/${case}.csv \
           -o ${case}_analysis
       echo "Completed $case"
   done
   
   echo "All cases processed!"

**Step 3: Run the batch script**

.. code-block:: bash

   chmod +x batch_process.sh
   ./batch_process.sh

**Step 4: Compare results**

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   
   cases = ['cyclone_20200115', 'cyclone_20200223', 'cyclone_20200315']
   
   fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
   
   for i, case in enumerate(cases):
       df = pd.read_csv(f'LEC_Results/{case}_analysis_track/{case}_analysis_track_results.csv',
                        parse_dates=['Time'], index_col='Time')
       
       axes[i].plot(df.index, df['Ke'], label='Ke', color='red')
       axes[i].plot(df.index, df['Ae'], label='Ae', color='blue')
       axes[i].set_ylabel('Energy (J/m²)')
       axes[i].set_title(f'{case}')
       axes[i].legend()
       axes[i].grid(True)
   
   plt.tight_layout()
   plt.savefig('comparison_all_cases.png')
   plt.show()

Next Steps
----------

After completing these tutorials:

- Explore the :doc:`results` section to understand all output files in detail
- Review :doc:`math` for detailed mathematical formulations of each term
- Check :doc:`contributing` if you want to contribute to the project
- Visit :doc:`support` for help with specific issues

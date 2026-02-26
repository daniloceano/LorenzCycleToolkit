Usage
=====

Before You Start
----------------

Before using the LorenzCycleToolkit, you need to configure the required input files. The specific files depend on which framework you choose:

**Required for all frameworks:**

- **Namelist file** (``inputs/namelist``): Specifies the variable names and units in your NetCDF file. See the :doc:`configuration` page for details on creating this file from the provided presets.

**Required for specific frameworks:**

- **Fixed framework** (``-f`` flag): Requires ``inputs/box_limits`` file defining the spatial domain
- **Moving framework** (``-t`` flag): Requires ``inputs/track`` file (or custom path via ``--trackfile``) defining the system's center over time
- **Interactive framework** (``-c`` flag): No additional file required; domain is selected interactively

For detailed information on creating these configuration files, see the :doc:`configuration` page.

Input Data Requirements
-----------------------

The LorenzCycleToolkit requires NetCDF files with specific characteristics:

**Required Format:**

- **Grid type**: Regular latitude-longitude grid (not supported: irregular grids, rotated poles)
- **Vertical levels**: Isobaric (pressure) levels in Pa or hPa
- **Required variables**: 
  
  - Air temperature
  - Geopotential or geopotential height
  - Omega (vertical velocity)
  - Eastward wind component (u)
  - Northward wind component (v)

**Data Preparation Best Practices:**

To avoid memory issues and improve processing speed, it is **strongly recommended** to pre-process your data:

1. **Temporal subsetting**: Extract only the time period needed for your analysis
2. **Spatial subsetting**: Crop to a domain slightly larger than your analysis region
3. **Variable selection**: Keep only the required variables listed above
4. **Vertical subsetting**: Include only the pressure levels needed (typically 1000-100 hPa for tropospheric studies)

**Example: Preparing data with Python (xarray)**

.. code-block:: python

   import xarray as xr
   
   # Open original dataset
   ds = xr.open_dataset('original_data.nc')
   
   # Subset in time (example: January 2020)
   ds_subset = ds.sel(time=slice('2020-01-01', '2020-01-31'))
   
   # Subset in space (example: South Atlantic region)
   ds_subset = ds_subset.sel(
       latitude=slice(-60, 0),
       longitude=slice(-80, 20)
   )
   
   # Select only required variables
   required_vars = ['t', 'z', 'u', 'v', 'w']  # Adjust names to your dataset
   ds_subset = ds_subset[required_vars]
   
   # Select pressure levels (example: 1000 to 100 hPa)
   ds_subset = ds_subset.sel(level=slice(1000, 100))
   
   # Save processed dataset
   ds_subset.to_netcdf('processed_data.nc')

**Example: Preparing data with CDO (Climate Data Operators)**

.. code-block:: bash

   # Select time period, spatial domain, and variables
   cdo -select,name=t,z,u,v,w \
       -sellonlatbox,-80,20,-60,0 \
       -seldate,2020-01-01,2020-01-31 \
       -sellevel,1000,925,850,700,500,300,250,200,150,100 \
       original_data.nc processed_data.nc

**Note**: Adjust variable names, coordinates, and domains according to your specific dataset.

System Requirements and Performance
------------------------------------

**Memory Requirements**

The toolkit loads entire NetCDF files into memory using xarray with dask. Typical requirements:

- **Small domain** (10° × 10°, 7 days, 3-hour resolution): ~500 MB RAM
- **Medium domain** (30° × 30°, 30 days, 6-hour resolution): ~2-4 GB RAM
- **Large domain** (60° × 60°, 60 days, 3-hour resolution): ~8-16 GB RAM

**Performance Tips**

1. **Pre-process your data** to include only the necessary components:
   
   - Time period of interest
   - Spatial domain (with small buffer around analysis region)
   - Pressure levels relevant to your study (e.g., 1000-100 hPa for tropospheric analyses)
   - Required variables only (u, v, ω, T, geopotential)

2. **Choose appropriate temporal resolution**:
   
   - Synoptic-scale analysis: 6-hour intervals
   - Mesoscale systems: 3-hour intervals
   - High-resolution studies: 1-hour intervals (requires significantly more memory)

3. **For very large datasets**, consider:
   
   - Splitting the analysis into smaller time periods
   - Using ERA5 on a coarser grid (0.5° instead of 0.25°)
   - Processing days or weeks individually

**Computational Time**

Typical processing times (on a modern desktop with 16 GB RAM):

- Small domain, 7 days: ~2-5 minutes
- Medium domain, 30 days: ~10-20 minutes  
- Large domain, 60 days: ~30-60 minutes

Time increases with:

- Domain size
- Number of timesteps
- Number of vertical levels
- Plot generation (``-p`` flag)
- Interactive framework (``-c`` flag, requires user input at each timestep)

Fixed Framework
---------------

The fixed framework analyzes a stationary spatial domain over time.

**Prerequisites**

- Prepared NetCDF file (see "Input Data Requirements" above)
- ``inputs/namelist`` file (see :doc:`configuration`)
- ``inputs/box_limits`` file defining the domain (see :doc:`configuration`)

**Command**
::

   python lorenzcycletoolkit.py path/to/infile.nc -r -f

Moving Framework
----------------

The moving framework follows a moving system (e.g., a cyclone) over time.

**Prerequisites**

- Prepared NetCDF file (see "Input Data Requirements" above)
- ``inputs/namelist`` file (see :doc:`configuration`)
- ``inputs/track`` file defining the system's center over time (see :doc:`configuration`)

**Command**
::

   python lorenzcycletoolkit.py path/to/infile.nc -r -t

Interactive Domain Selection
----------------------------

The interactive (choose) framework allows you to manually select the analysis domain at each timestep.

**Prerequisites**

- Prepared NetCDF file (see "Input Data Requirements" above)
- ``inputs/namelist`` file (see :doc:`configuration`)

**Command**
::

   python lorenzcycletoolkit.py path/to/infile.nc -r -c

Advanced Options
================

Custom Output Names
-------------------

By default, results are saved in ``LEC_Results/<inputfile>_<framework>/``. You can customize the output directory name using the ``-o/--outname`` flag:

**Example**:

.. code-block:: bash

   python lorenzcycletoolkit.py data.nc -r -f -o my_cyclone_analysis

This creates the directory: ``LEC_Results/my_cyclone_analysis_fixed/``

**Use cases**:

- Organizing multiple analyses of the same dataset with different configurations
- Creating descriptive names for specific case studies
- Avoiding overwriting previous results

Using Pre-computed Vorticity
-----------------------------

**Flag**: ``-z`` or ``--zeta``

By default, the toolkit computes relative vorticity at 850 hPa for system tracking and visualization. If your track file already contains accurate vorticity values, use the ``-z`` flag to use those values instead of computing them:

.. code-block:: bash

   python lorenzcycletoolkit.py data.nc -r -t -z

**When to use this flag**:

- Your track file includes a vorticity column with pre-computed values
- Working with high-resolution data where pre-computed vorticity is more accurate
- Analyzing systems where a different level or method for vorticity computation is preferred
- Computational efficiency (skips vorticity calculation)

**Track file format with vorticity**:

.. code-block:: text

   time;Lat;Lon;vorticity
   2020-01-01 00:00:00;-25.0;-45.0;2.5e-5
   2020-01-01 06:00:00;-25.5;-44.5;3.2e-5
   2020-01-01 12:00:00;-26.0;-44.0;4.1e-5

Working with MPAS-A Data
-------------------------

**Flag**: ``-m`` or ``--mpas``

The Model for Prediction Across Scales - Atmosphere (MPAS-A) uses an unstructured mesh. When using MPAS-A data that has been processed with MPAS-BR routines (converted to pressure levels on a regular grid), use the ``-m`` flag:

.. code-block:: bash

   python lorenzcycletoolkit.py mpas_data.nc -r -f -m

**MPAS-A Specific Requirements**:

1. **Data preprocessing**: MPAS-A data must be processed with MPAS-BR routines to:
   
   - Convert from native unstructured mesh to regular lat-lon grid
   - Interpolate from model levels to isobaric (pressure) levels
   - Extract required variables

2. **Namelist configuration**: Use the MPAS-A specific namelist:
   
   .. code-block:: bash
   
      cp inputs/namelist_MPAS-A inputs/namelist

3. **Dimension handling**: The ``-m`` flag automatically handles the ``standard_height`` dimension present in MPAS-BR processed data

**Example workflow**:

.. code-block:: bash

   # Step 1: Prepare namelist
   cp inputs/namelist_MPAS-A inputs/namelist
   
   # Step 2: Prepare box_limits
   cat > inputs/box_limits << EOF
   min_lon;-55
   max_lon;-35
   min_lat;-35
   max_lat;-20
   EOF
   
   # Step 3: Run analysis
   python lorenzcycletoolkit.py mpas_processed_data.nc -r -f -m -p

**Note**: The MPAS-A flag is only needed for data processed with MPAS-BR routines. Regular MPAS output requires additional preprocessing steps not covered by the toolkit.

Automatic ERA5 Data Download
============================

The LorenzCycleToolkit supports automatic downloading of ERA5 reanalysis data from the Copernicus Climate Data Store (CDS) using the ``--cdsapi`` flag. This feature is particularly useful when you have a track file but haven't yet downloaded the corresponding atmospheric data.

**Prerequisites**

1. **CDS API Account**: Create a free account at the `CDS registration page <https://cds.climate.copernicus.eu/user/register>`_.

2. **API Key Configuration**: After registration, obtain your API key from your `user profile <https://cds.climate.copernicus.eu/user>`_ and create a ``.cdsapirc`` file in your home directory:

   .. code-block:: text

      url: https://cds.climate.copernicus.eu/api
      key: YOUR-UID:YOUR-API-KEY

   Replace ``YOUR-UID`` and ``YOUR-API-KEY`` with your actual credentials.

3. **Namelist Configuration**: The ``--cdsapi`` flag now **automatically** uses the ERA5-compatible namelist (``inputs/namelist_ERA5-cdsapi``). You no longer need to manually copy or configure the namelist file.

**How It Works**

When using the ``--cdsapi`` flag:

1. The program reads your track file to determine:
   - Temporal range (exact start and end times, not just dates)
   - Spatial domain (automatically adds a 15° buffer around track boundaries)

2. Downloads ERA5 pressure level data with smart time detection:
   - **First day**: Downloads only from track start time onwards
   - **Last day**: Downloads only until track end time
   - **Middle days**: Downloads complete days (all hours)
   - **Temporal resolution**: Configurable via ``--time-resolution`` flag (default: 3 hours)
   - **Variables**: U and V wind components, temperature, vertical velocity (omega), geopotential
   - **Pressure levels**: All standard pressure levels (1000 hPa to 1 hPa)

3. Saves the data to the filename specified as the first argument

This smart time detection minimizes unnecessary downloads and reduces API costs.

**Usage Example**

For a track file covering January 1-3, 2020:

.. code-block:: bash

   python lorenzcycletoolkit.py my_output_file.nc -t -r -p --cdsapi --trackfile inputs/track_20200101.csv

**Important Notes**

- The filename provided (``my_output_file.nc``) doesn't need to exist beforehand—it will be created by the download process
- **Temporal resolution**: Default is 3 hours. Using 1-hour intervals may exceed CDS API cost limits for large domains/long periods
- Recommended time resolutions: 3, 6, or 12 hours
- Download times depend on data volume and CDS server load (typically several minutes)
- A 15° buffer is automatically applied to ensure adequate spatial coverage
- If the file already exists, the download is skipped and the existing file is used

**Troubleshooting**

- **Authentication errors**: Verify your ``.cdsapirc`` file is correctly formatted and contains valid credentials
- **Slow downloads**: CDS API can be slow during peak hours; use the ``-v`` flag for detailed progress logging
- **File size**: ERA5 files can be large (>1 GB); ensure adequate disk space

Flags
=====

- `-r`, `--residuals`: Compute the Dissipation and Generation terms as residuals.
- `-f`, `--fixed`: Compute the energetics for a fixed domain specified by the 'box_limits' file.
- `-t`, `--track`: Define the domain using a track file.
- `-c`, `--choose`: Interactively select the domain for each time step.
- `-o`, `--outname`: Specify an output name for the results.
- `-z`, `--zeta`: Use the vorticity from the track file instead of computing it at 850 hPa.
- `-m`, `--mpas`: Specify this flag if working with MPAS-A data processed with MPAS-BR routines.
- `-p`, `--plots`: Generate plots.
- `-v`, `--verbosity`: Enable debug logging for detailed output.
- `--cdsapi`: Automatically download ERA5 data from CDS API based on track file (experimental).
- `--time-resolution`: Temporal resolution in hours for CDS API downloads (default: 3). Recommended values: 3, 6, or 12 hours to avoid exceeding API cost limits.
- `--trackfile`: Specify a custom track file path (default: ``inputs/track``).
- `--box_limits`: Specify a custom box limits file path (default: ``inputs/box_limits``).

Ensure that the provided NetCDF file and the `namelist` configuration align with the selected flags.

Common Issues and Questions
===========================

If you encounter problems:

- **Track file timesteps don't match data**: See :doc:`troubleshooting` section on "Track File Issues"
- **Units and variable specifications**: See :doc:`troubleshooting` section on "Units and Variable Specifications"
- **Memory or performance issues**: See :doc:`troubleshooting` section on "Memory and Performance Issues"
- **Configuration errors**: See :doc:`troubleshooting` section on "Configuration Issues"
- **CDS API problems**: See :doc:`troubleshooting` section on "CDS API Issues"

For a comprehensive list of common issues and solutions, visit the :doc:`troubleshooting` page.

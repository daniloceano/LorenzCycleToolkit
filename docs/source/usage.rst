Usage
=====

Fixed Framework
---------------

**Prerequisites**
- NetCDF file with necessary variables (u, v, omega, air temperature, geopotential or geopotential height).
- `inputs/box_limits` file defining the domain.
- `inputs/namelist` file specifying variable names and units.

**Command**
::

   python lorenzcycletoolkit.py path/to/infile.nc -r -f

Moving Framework
----------------

**Pre-defined Domain**
- Follow initial steps for the Fixed Framework.
- Use a `inputs/track_file` to define the system's center over time.

**Command**
::

   python lorenzcycletoolkit.py path/to/infile.nc -r -t

Interactive Domain Selection
----------------------------
- For interactive domain selection, run::

   python lorenzcycletoolkit.py path/to/infile.nc -r -c

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

Troubleshooting and FAQ
=======================

This page addresses common issues and frequently asked questions when using the LorenzCycleToolkit.

Track File Issues
-----------------

Track file has more timesteps than input data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Your track file extends beyond the temporal coverage of your NetCDF data file.

**Example**:

.. code-block:: text

   Track file: 2020-01-01 00:00 to 2020-01-10 23:00
   Input data:  2020-01-01 00:00 to 2020-01-08 18:00

**Error message**:

.. code-block:: text

   ❌ Track final timestamp (2020-01-10 23:00:00) is later than data final timestamp 
   (2020-01-08 18:00:00). Please adjust the track file or re-download the data.

**Solutions**:

1. **Trim your track file** to match the data coverage:

   .. code-block:: python
   
      import pandas as pd
      
      # Load track file
      track = pd.read_csv('inputs/track', sep=';', parse_dates=['time'], index_col='time')
      
      # Trim to match data period
      track_trimmed = track['2020-01-01':'2020-01-08 18:00']
      
      # Save trimmed track
      track_trimmed.to_csv('inputs/track_trimmed', sep=';')

2. **Download additional data** to cover the full track period (if using CDS API):

   .. code-block:: bash
   
      python lorenzcycletoolkit.py extended_data.nc -t -r --cdsapi --trackfile inputs/track

3. **Use a pre-downloaded file** that covers the full period

Track file has fewer timesteps than input data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Your track file covers a shorter period than your input data.

**Example**:

.. code-block:: text

   Track file: 2020-01-03 00:00 to 2020-01-05 18:00
   Input data:  2020-01-01 00:00 to 2020-01-10 23:00

**Behavior**: This is **perfectly fine**! The toolkit will automatically use only the timesteps from the input data that match the track file period.

**What happens**:

.. code-block:: python

   # The toolkit automatically does this:
   data = data.sel(time=track.index.values)

Your analysis will run from 2020-01-03 00:00 to 2020-01-05 18:00, using only the relevant portion of the input data.

**Use case**: This is useful when you have a large dataset but only want to analyze a specific event or time period.

Track file has timesteps not present in data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Some track timesteps don't match any data timesteps (e.g., track has 1-hour intervals but data has 3-hour intervals).

**Example**:

.. code-block:: text

   Track timesteps: 00:00, 01:00, 02:00, 03:00, 04:00, ...
   Data timesteps:  00:00, 03:00, 06:00, 09:00, ...

**Error message**:

.. code-block:: text

   KeyError: 'Some timesteps from track not found in data'

**Solutions**:

1. **Resample your track file** to match data temporal resolution:

   .. code-block:: python
   
      import pandas as pd
      
      # Load track
      track = pd.read_csv('inputs/track', sep=';', parse_dates=['time'], index_col='time')
      
      # Resample to 3-hour intervals (matching data)
      track_resampled = track.resample('3H').first()  # or .interpolate() for interpolation
      
      # Remove any NaN rows
      track_resampled = track_resampled.dropna()
      
      # Save resampled track
      track_resampled.to_csv('inputs/track_resampled', sep=';')

2. **If using CDS API**, the toolkit automatically handles this by resampling the track:

   .. code-block:: bash
   
      # The toolkit will resample your track to match --time-resolution
      python lorenzcycletoolkit.py data.nc -t -r --cdsapi --time-resolution 3 --trackfile inputs/track

Track temporal resolution differs from data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Track has higher temporal resolution than the data (e.g., 1-hour track vs 6-hour data).

**Behavior**: The toolkit checks this and raises an error if track resolution is higher than data resolution.

**Error message**:

.. code-block:: text

   ❌ Track time step is higher than data time step. Please resample the track 
   to match the data time step.

**Why this matters**: You cannot analyze a system every hour if your atmospheric data is only available every 6 hours.

**Solution**: Resample your track file to match or exceed the data temporal resolution:

.. code-block:: python

   import pandas as pd
   
   track = pd.read_csv('inputs/track', sep=';', parse_dates=['time'], index_col='time')
   
   # If data is 6-hourly, resample track to 6-hourly
   track_6h = track.resample('6H').first()  # Takes the first value in each 6h window
   # OR use interpolation if you want intermediate positions
   track_6h = track.resample('6H').interpolate()
   
   track_6h = track_6h.dropna()
   track_6h.to_csv('inputs/track_6h', sep=';')

Units and Variable Specifications
----------------------------------

Do variables need to match units specified in namelist?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Short answer**: No, but the namelist must specify the **actual units in your file**.

**How it works**:

1. The namelist tells the toolkit what units are **currently** in your NetCDF file
2. The toolkit reads the data with those units
3. The toolkit **automatically converts** everything to SI units internally using MetPy

**Example**:

If your ERA5 file has temperature in Kelvin:

.. code-block:: text

   ;standard_name;Variable;Units
   Air Temperature;air_temperature;t;K

The toolkit will:

- Read temperature as Kelvin
- Keep it in Kelvin for calculations (Kelvin is already SI)

If your file had temperature in Celsius (unusual):

.. code-block:: text

   ;standard_name;Variable;Units
   Air Temperature;air_temperature;t;degC

The toolkit would:

- Read temperature as Celsius
- Convert to Kelvin internally
- Perform all calculations in SI units

**Important**: The "Units" column in the namelist should match what's **actually in your file**, not what you want. The toolkit handles conversions automatically.

Pressure levels must be in Pa or hPa?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Answer**: Either is fine! The toolkit automatically converts to Pascals (Pa) internally.

**Example from the code**:

.. code-block:: python

   # The toolkit does this automatically:
   levels_Pa = (data[LevelIndexer] * units(str(data[LevelIndexer].units))).metpy.convert_units("Pa")

Your file can have:

- ``pressure_level`` in hPa: [1000, 925, 850, 700, 500, ...]
- ``pressure_level`` in Pa: [100000, 92500, 85000, 70000, 50000, ...]

Both will work correctly. Just make sure the units attribute in your NetCDF file is set correctly.

Why does my analysis fail with "units" error?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Error mentioning units or unit conversion.

**Common causes**:

1. **Missing units attribute in NetCDF file**:

   .. code-block:: python
   
      import xarray as xr
      
      # Check if variables have units
      ds = xr.open_dataset('your_file.nc')
      print(ds['t'].attrs)  # Should show 'units': 'K'
      print(ds['pressure_level'].attrs)  # Should show 'units': 'hPa' or 'Pa'

   **Solution**: Add units to your NetCDF file:

   .. code-block:: python
   
      ds['t'].attrs['units'] = 'K'
      ds['pressure_level'].attrs['units'] = 'hPa'
      ds.to_netcdf('your_file_fixed.nc')

2. **Non-standard unit names**: Use standard CF conventions:

   - Temperature: ``K`` or ``degC`` (not ``Kelvin`` or ``celsius``)
   - Pressure: ``Pa`` or ``hPa`` (not ``mb`` or ``millibars``)
   - Winds: ``m/s`` or ``m s-1`` (not ``m/sec``)

3. **Geopotential vs Geopotential Height**:

   - Geopotential: units should be ``m**2 s**-2`` or ``m2 s-2``
   - Geopotential Height: units should be ``m``
   
   The toolkit handles both, but you must specify correctly in your namelist.

Data Format and Loading Issues
-------------------------------

"Could not open file" error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Error message**:

.. code-block:: text

   ❌ Could not open file. Check if path, namelist file, and file format (.nc) are correct.

**Checklist**:

1. **File path is correct**:

   .. code-block:: bash
   
      ls -lh path/to/your/file.nc  # Verify file exists

2. **File is valid NetCDF**:

   .. code-block:: bash
   
      ncdump -h your_file.nc  # Should show header without errors

3. **File is not corrupted**:

   .. code-block:: python
   
      import xarray as xr
      ds = xr.open_dataset('your_file.nc')  # Should open without errors

4. **You have read permissions**:

   .. code-block:: bash
   
      chmod 644 your_file.nc  # Add read permissions if needed

"Namelist does not match the data" error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Error message**:

.. code-block:: text

   ❌ The variable list does not match the data. Check if the 'namelist' text file is correct.

Followed by a detailed list of available coordinates and variables in your dataset.

**Solution process**:

1. **Check the error output** - it shows all coordinates and variables in your file:

   .. code-block:: text
   
      📊 DATASET INFORMATION
      ======================================================================
      
      🗺️  Available Coordinates:
         • latitude (degrees_north) - latitude
         • longitude (degrees_east) - longitude
         • level (hPa) - pressure_level
         • time - valid_time
      
      📋 Available Variables:
         • t (K) - Temperature [air_temperature]
         • u (m/s) - U component of wind [eastward_wind]
         • v (m/s) - V component of wind [northward_wind]
         • w (Pa/s) - Vertical velocity [lagrangian_tendency_of_air_pressure]
         • z (m**2/s**2) - Geopotential [geopotential]

2. **Update your namelist** to use the actual variable names:

   .. code-block:: text
   
      ;standard_name;Variable;Units
      Air Temperature;air_temperature;t;K
      Geopotential;geopotential;z;m**2/s**2
      Omega Velocity;omega;w;Pa/s
      Eastward Wind Component;eastward_wind;u;m/s
      Northward Wind Component;northward_wind;v;m/s
      Longitude;;longitude
      Latitude;;latitude
      Time;;time
      Vertical Level;;level

3. **Use a preset if available**:

   .. code-block:: bash
   
      # For ERA5 data
      cp inputs/namelist_ERA5-cdsapi inputs/namelist
      
      # For NCEP-R2 data
      cp inputs/namelist_NCEP-R2 inputs/namelist

Memory and Performance Issues
------------------------------

"MemoryError" or system becomes unresponsive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Your dataset is too large for available RAM.

**Immediate solutions**:

1. **Reduce spatial domain**: Pre-process your data to a smaller region:

   .. code-block:: python
   
      import xarray as xr
      
      ds = xr.open_dataset('large_file.nc')
      
      # Subset to smaller region
      ds_subset = ds.sel(
          latitude=slice(-50, -10),
          longitude=slice(-70, -30)
      )
      
      ds_subset.to_netcdf('smaller_file.nc')

2. **Reduce temporal coverage**: Analyze shorter periods:

   .. code-block:: python
   
      ds_subset = ds.sel(time=slice('2020-01-01', '2020-01-07'))

3. **Reduce vertical levels**: Keep only tropospheric levels:

   .. code-block:: python
   
      ds_subset = ds.sel(level=slice(1000, 100))  # 1000 to 100 hPa

4. **Increase temporal resolution**: Use 6-hour instead of 3-hour data:

   .. code-block:: python
   
      ds_subset = ds.sel(time=ds.time.dt.hour.isin([0, 6, 12, 18]))

5. **Remove unnecessary variables**:

   .. code-block:: python
   
      # Keep only required variables
      required_vars = ['t', 'u', 'v', 'w', 'z', 'latitude', 'longitude', 'level', 'time']
      ds_subset = ds[required_vars]

**See also**: :doc:`usage` section on "System Requirements and Performance"

Analysis takes too long
~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Processing is slower than expected.

**Optimization tips**:

1. **Use appropriate temporal resolution**: 6-hour data processes faster than 1-hour
2. **Reduce domain size**: Analyze only the region of interest
3. **Disable plots for initial runs**: Remove ``-p`` flag during testing
4. **Use fixed framework instead of interactive**: ``-f`` is much faster than ``-c``
5. **Check your data preprocessing**: Pre-subset data before running the toolkit

Configuration Issues
--------------------

"No such file or directory: 'inputs/namelist'"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: You didn't create the namelist file.

**Solution**: Copy one of the presets:

.. code-block:: bash

   # For ERA5 data
   cp inputs/namelist_ERA5-cdsapi inputs/namelist
   
   # For NCEP Reanalysis 2
   cp inputs/namelist_NCEP-R2 inputs/namelist
   
   # For NCEP Reanalysis 1
   cp inputs/namelist_NCEP-R1 inputs/namelist

**See also**: :doc:`configuration` for detailed instructions.

"No such file or directory: 'inputs/box_limits'" (fixed framework)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Using ``-f`` flag but haven't created box_limits file.

**Solution**: Create the file with your domain bounds:

.. code-block:: bash

   cat > inputs/box_limits << EOF
   min_lon;-60
   max_lon;-30
   min_lat;-45
   max_lat;-20
   EOF

"No such file or directory: 'inputs/track'" (moving framework)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Using ``-t`` flag but haven't specified a track file.

**Solution**: Either:

1. **Create a track file** at ``inputs/track``
2. **Use the --trackfile flag** to specify a different path:

   .. code-block:: bash
   
      python lorenzcycletoolkit.py data.nc -r -t --trackfile path/to/your/track.csv

CDS API Issues
--------------

CDS API authentication fails
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Error**: Connection or authentication errors when using ``--cdsapi``.

**Solution checklist**:

1. **Check credentials file exists**:

   .. code-block:: bash
   
      cat ~/.cdsapirc

2. **Verify format**:

   .. code-block:: text
   
      url: https://cds.climate.copernicus.eu/api
      key: YOUR-UID:YOUR-API-KEY

3. **Verify credentials are correct**: Log in to https://cds.climate.copernicus.eu/user and check your UID and API key

4. **Check permissions**:

   .. code-block:: bash
   
      chmod 600 ~/.cdsapirc

5. **Accept Terms and Conditions**: You must accept the Copernicus terms at https://cds.climate.copernicus.eu/ before API access works

CDS API download is very slow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: Downloads take hours or timeout.

**Causes and solutions**:

1. **CDS servers are busy**: Try during off-peak hours (evenings/weekends in Europe)
2. **Domain is too large**: Reduce spatial coverage or temporal resolution
3. **Requesting too many levels**: Consider requesting only tropospheric levels
4. **Use coarser temporal resolution**:

   .. code-block:: bash
   
      # Use 6-hour instead of 3-hour
      python lorenzcycletoolkit.py file.nc -t -r --cdsapi --time-resolution 6

Still Having Issues?
--------------------

If your problem isn't covered here:

1. **Enable verbose logging** to see detailed information:

   .. code-block:: bash
   
      python lorenzcycletoolkit.py your_file.nc -r -f -v

2. **Check the log file** in ``LEC_Results/*/log.*`` for detailed error messages

3. **Review the documentation**:
   
   - :doc:`configuration` - Setup requirements
   - :doc:`usage` - Command-line options and data requirements
   - :doc:`examples_and_tutorials` - Step-by-step examples

4. **Seek support**:
   
   - GitHub Issues: https://github.com/daniloceano/LorenzCycleToolkit/issues
   - Email: danilo.oceano@gmail.com

When reporting issues, please include:

- Full error message and traceback
- Your command line
- Relevant parts of the log file
- Description of your input data (source, resolution, coverage)
- LorenzCycleToolkit version: check with ``git describe --tags``

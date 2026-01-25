CDS API - Complete Example
===========================

This guide provides a complete walkthrough for using the CDS API feature to automatically download ERA5 data.

Step-by-Step Example
---------------------

**1. Setup CDS API Credentials**

First, register at the `Copernicus Climate Data Store <https://cds.climate.copernicus.eu/user/register>`_ and configure your API key.

Create or edit ``~/.cdsapirc``:

.. code-block:: bash

   cat > ~/.cdsapirc << EOF
   url: https://cds.climate.copernicus.eu/api
   key: YOUR-UID:YOUR-API-KEY
   EOF

Replace ``YOUR-UID:YOUR-API-KEY`` with your actual credentials from your CDS profile.

**2. Prepare Your Track File**

Create a track file with your system's trajectory. The file must use semicolons as delimiters and include at least ``time``, ``Lat``, and ``Lon`` columns:

.. code-block:: text

   time;Lat;Lon
   1979-02-05 00:00:00;-35.5;-45.2
   1979-02-05 06:00:00;-36.0;-44.8
   1979-02-05 12:00:00;-36.5;-44.3

Save this file (e.g., as ``inputs/track_19790205.csv``).

**Important**: Your track can start and end at any hour of the day. The toolkit will automatically download only the necessary hours:

- Track from 18:00 on day 1 → Downloads start at 18:00 (not 00:00)
- Track ending at 15:00 on last day → Downloads stop at 15:00 (not 23:00)
- Middle days are always downloaded completely (00:00 to 23:00)

**3. Configure Namelist for ERA5**

Copy the ERA5-specific namelist:

.. code-block:: bash

   cp inputs/namelist_ERA5-cdsapi inputs/namelist

**4. Run the Toolkit with CDS API**

Use the following command structure:

.. code-block:: bash

   python lorenzcycletoolkit.py OUTPUT_FILENAME.nc -t -r -p --cdsapi --trackfile PATH/TO/TRACKFILE.csv

Optionally, specify the temporal resolution (default is 3 hours):

.. code-block:: bash

   python lorenzcycletoolkit.py OUTPUT_FILENAME.nc -t -r -p --cdsapi --time-resolution 6 --trackfile PATH/TO/TRACKFILE.csv

**Complete Example**:

.. code-block:: bash

   python lorenzcycletoolkit.py 19790205_ERA5.nc -t -r -p --cdsapi --trackfile inputs/track_19790205.csv

**Important**: 
- The first argument (``19790205_ERA5.nc``) is the output filename where ERA5 data will be saved
- This file doesn't need to exist beforehand—it will be created by the download process
- If it already exists, the download is skipped

What Happens Behind the Scenes
-------------------------------

When you run the command with ``--cdsapi``:

1. **Track Analysis**: The program reads your track file to determine:
   
   - Start date and time: earliest timestamp in track
   - End date and time: latest timestamp in track
   - Spatial bounds: min/max latitude and longitude

2. **Smart Time Range Detection**: The program intelligently determines which hours to download:
   
   - **First day**: Downloads only from the track's start time onwards (not from 00Z)
   - **Last day**: Downloads only until the track's end time (not until 23Z)
   - **Middle days**: Downloads complete days (00Z to 23Z)
   - All times are rounded to match your ``--time-resolution`` setting
   
   Example: If your track runs from June 6 18:00 to June 8 15:00 with 3-hour resolution:
   
   - June 6: Downloads 18:00, 21:00 only
   - June 7: Downloads 00:00, 03:00, 06:00, ..., 21:00 (full day)
   - June 8: Downloads 00:00, 03:00, 06:00, 09:00, 12:00, 15:00 only

3. **Automatic Buffering**: A 15° buffer is added to all spatial boundaries to ensure adequate coverage around your system

4. **Day-by-Day Download**: Data is downloaded separately for each day, showing progress for each file

5. **Data Request**: The program requests from CDS:
   
   - **Variables**: u, v, t (temperature), w (omega), z (geopotential)
   - **Pressure levels**: All 37 standard levels (1000 to 1 hPa)
   - **Temporal resolution**: Configurable via ``--time-resolution`` flag (default: 3 hours)
   - **Spatial coverage**: Buffered domain around your track

6. **Concatenation**: All daily files are merged into a single NetCDF file

7. **Processing**: The toolkit continues with LEC calculations using the downloaded data

Expected Output
---------------

During execution, you'll see detailed progress information:

.. code-block:: text

   INFO - Starting LEC analysis
   INFO - Will download 3 day(s) of data: 1992-06-06 to 1992-06-08
   INFO - Starting data download from CDS API...
   INFO - Expected files to download: 3
   INFO -   - Day 1: 1992-06-06
   INFO -   - Day 2: 1992-06-07
   INFO -   - Day 3: 1992-06-08
   INFO - Downloading day 1/3: 1992-06-06
   DEBUG -   Times: 18:00 to 21:00 (2 timesteps)
   INFO -   ✓ Downloaded successfully: 45.23 MB
   INFO - Downloading day 2/3: 1992-06-07
   DEBUG -   Times: 00:00 to 21:00 (8 timesteps)
   INFO -   ✓ Downloaded successfully: 178.45 MB
   INFO - Downloading day 3/3: 1992-06-08
   DEBUG -   Times: 00:00 to 18:00 (7 timesteps)
   INFO -   ✓ Downloaded successfully: 156.78 MB
   INFO - ============================================================
   INFO - All daily files downloaded successfully!
   INFO - Total files: 3
   INFO - ============================================================
   INFO - Concatenating daily files into single dataset...
   INFO - Merging datasets along time dimension...
   INFO - Saving final dataset to: output.nc
   INFO - ============================================================
   INFO - ✓ DATA DOWNLOAD COMPLETED SUCCESSFULLY!
   INFO - Final file: output.nc
   INFO - File size: 380.46 MB
   INFO - Time coverage: 1992-06-06 18:00:00 to 1992-06-08 18:00:00
   INFO - ============================================================

Key features of the new logging:

- **Progress tracking**: See which day is currently being downloaded
- **Time range info**: Know exactly which hours are being requested
- **Success indicators**: Checkmarks (✓) show successful downloads
- **File sizes**: Monitor download sizes for each day
- **Clear summaries**: Final summary shows total coverage and file size

The download typically takes 5-15 minutes depending on:
- Data volume (duration and spatial extent)
- CDS server load
- Your internet connection speed

**Note**: Use the ``-v`` (verbose) flag to see detailed timing information for each day.

Common Use Cases
----------------

**Case 1: Single Cyclone Analysis**

.. code-block:: bash

   python lorenzcycletoolkit.py cyclone_catarina_2004.nc -t -r -p --cdsapi --trackfile tracks/catarina_track.csv

**Case 2: Multiple Time Periods**

For analyzing the same region at different times, download once and reuse:

.. code-block:: bash

   # First run - downloads data
   python lorenzcycletoolkit.py region_jan2020.nc -t -r -p --cdsapi --trackfile tracks/jan2020.csv
   
   # Subsequent runs - reuses existing file (no download)
   python lorenzcycletoolkit.py region_jan2020.nc -t -r -p --trackfile tracks/jan2020.csv

**Case 3: Different Analysis Parameters**

Once data is downloaded, experiment with different flags without re-downloading:

.. code-block:: bash

   # With plots and residuals
   python lorenzcycletoolkit.py mydata.nc -t -r -p --trackfile tracks/mytrack.csv
   
   # Without residuals
   python lorenzcycletoolkit.py mydata.nc -t -p --trackfile tracks/mytrack.csv

Troubleshooting
---------------

**Problem**: ``Authentication failed``

**Solution**: Check your ``.cdsapirc`` file:

.. code-block:: bash

   cat ~/.cdsapirc

Ensure it contains valid credentials and the URL is correct.

**Problem**: ``Request too large`` or ``cost limits exceeded``

**Solution**: 
- Use a larger time resolution: ``--time-resolution 6`` or ``--time-resolution 12``
- Reduce the spatial extent of your track
- Use shorter time periods
- Check track file doesn't span too many days

**Example**:

.. code-block:: bash

   # If 3-hour resolution fails, try 6 hours
   python lorenzcycletoolkit.py output.nc -t -r --cdsapi --time-resolution 6 --trackfile track.csv

**Problem**: Download hangs or is very slow

**Solution**:
- Use ``-v`` flag for verbose logging to monitor progress
- Check CDS service status at https://cds.climate.copernicus.eu
- Try during off-peak hours (CDS is busiest during European business hours)

**Problem**: ``FileNotFoundError`` for track file

**Solution**: Use absolute paths:

.. code-block:: bash

   python lorenzcycletoolkit.py output.nc -t -r --cdsapi --trackfile /full/path/to/track.csv

Advanced: Custom Temporal Resolution
-------------------------------------

The ``--time-resolution`` flag allows you to control the temporal resolution of downloaded ERA5 data. This is crucial for managing API cost limits.

**Available Resolutions**:

- **1 hour**: Highest resolution, but likely to exceed cost limits for large domains or long periods
- **3 hours** (default): Good balance between resolution and data volume
- **6 hours**: Recommended for large domains or long time periods
- **12 hours**: Suitable for extended analyses or very large domains
- **24 hours**: Daily data

**Examples**:

.. code-block:: bash

   # 3-hour resolution (default)
   python lorenzcycletoolkit.py output.nc -t -r --cdsapi --trackfile track.csv
   
   # 6-hour resolution
   python lorenzcycletoolkit.py output.nc -t -r --cdsapi --time-resolution 6 --trackfile track.csv
   
   # 12-hour resolution for large domain
   python lorenzcycletoolkit.py output.nc -t -r --cdsapi --time-resolution 12 --trackfile track.csv

**Cost Limit Guidelines**:

- Small domain (< 20° x 20°), short period (< 7 days): 3-hour resolution usually works
- Medium domain (20-40°), medium period (7-14 days): use 6-hour resolution
- Large domain (> 40°), long period (> 14 days): use 12-hour resolution

**Note**: After downloading with a specific resolution, the toolkit will interpolate your track data to match the downloaded time steps during processing.

Notes and Limitations
---------------------

- **Data source**: Only ERA5 reanalysis is supported via CDS API
- **Spatial resolution**: ERA5 native resolution (0.25° × 0.25°)
- **Pressure levels**: All 37 standard levels are downloaded (cannot be customized)
- **Variables**: Fixed set required for LEC calculations (u, v, t, omega, geopotential)
- **File format**: NetCDF output only
- **Authentication**: Requires valid CDS account and API key
- **Rate limits**: Subject to CDS API rate limiting during high-demand periods

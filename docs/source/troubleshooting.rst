Troubleshooting and FAQ
=======================

This page is the single place for solutions to common problems with the LorenzCycleToolkit. Each section is
organised by symptom; the last sections answer frequent questions and explain how to report a problem.

.. contents:: On this page
   :local:
   :depth: 2

.. _ts-installation:

1. Installation and Configuration
---------------------------------

``ModuleNotFoundError`` when running the toolkit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause**: the environment with the toolkit's dependencies is not active.

**Solution**: activate the environment created during :doc:`installation` (for example
``conda activate lorenz``) and run the command again from the repository root.

"No such file or directory: 'inputs/namelist'"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: You didn't create the namelist file.

**Solution**: When ``inputs/namelist`` is absent the toolkit uses ``inputs/namelist.default``; this message
means neither exists in the directory you run from. Run from the repository root, or copy one of the presets:

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

.. _ts-input-data:

2. Input Data and Tracks
------------------------

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

Why does my analysis fail with "units" error?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Track file extends beyond the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: the first track time is earlier than the first data time, or the last track time is later than the
last data time.

**Error message** (example):

.. code-block:: text

   ❌ Track final timestamp is later than data final timestamp!
   Track final timestamp (2020-01-10 23:00:00) is later than data final timestamp (2020-01-08 18:00:00).
   Please adjust the track file or re-download the data.

**Solutions**: provide data that cover the whole track (with ``--cdsapi`` or ``--ncar``, use a new output file
name: an existing file is not extended), or shorten the track file yourself and keep a record of what you changed.

Track file covers a shorter period than the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Behavior**: this is expected. Only the data times of the track period are used; the rest of the file is ignored.

Track times missing from the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: some track times are not data times, and the data are **not** simply sampled less often than the track
(for example a gap in the data, or a track at 01, 04, 07 UTC with data at 00, 03, 06 UTC).

**Behavior**: the run stops with an error that names the missing times and says whether they lie before, after or
inside the data period. Missing data times are never filled and the track is never interpolated.

**Solutions**: provide data with these times, or adjust the track file yourself (and keep a record of what you
changed).

.. _ts-track-alignment:

Track sampled more often than the data (automatic alignment)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: The track is sampled more often than the data, for example a 3-hourly track with 6-hourly ERA5
data (``--cdsapi`` or ``--ncar`` with ``--time-resolution 6``, or a local file).

**Behavior**: The toolkit keeps the track positions at the data times, with their original coordinates, and
leaves the other positions out. It is a selection by exact time, not an interpolation, and the track file is
not modified. The same rule applies to local files, ``--cdsapi`` and ``--ncar``. A warning is written once to
the terminal and to the run log, with the counts of the run:

.. code-block:: text

   ⚠️  WARNING — Temporal resolution mismatch
   Original track interval:        3 h
   Meteorological data interval:   6 h
   Action: track timestamps aligned with meteorological data.

   Original positions: 55
   Selected positions: 27 (1979-12-29T06:00:00 to 1980-01-04T18:00:00)
   Excluded positions: 28 (first 1979-12-29T03:00:00, last 1980-01-04T21:00:00)

   No spatial interpolation was performed.
   Original track file was not modified.

``lec_metadata.json`` records the alignment under ``track_time_alignment``: ``adapted``, ``rule``, the track, data
and selected intervals, the counts of original, selected and excluded positions, the selected times, and
``spatial_interpolation: false`` and ``trackfile_modified: false``. The moving framework's track output
(``<name>_track_trackfile``) holds only the positions used. With ``--cdsapi`` and ``--ncar`` the track is also checked against the requested hours **before**
the download, so a track that the request cannot serve is refused before any transfer, with the message
"request not made" and the reason.

**Do not interpolate the track** to work around this situation: the toolkit already uses the exact track
positions at the data times. Interpolating would introduce positions that are not in the track.

**When it is still an error**: the alignment is made only when it is unambiguous. The run stops with an
explanation when:

- the track or the data times are duplicated or not in increasing order;
- the track interval is irregular, or the data interval is irregular inside the track period (a gap);
- the data interval is not a whole multiple of the track interval (e.g. 4-hourly track, 6-hourly data);
- the data times do not fall on track times (aligning would need an interpolation);
- track times on the data grid are missing from the data (a gap or truncation, not a subsampled time);
- fewer than two positions would remain.

When the data are sampled **more** often than the track (e.g. 6-hourly track, 3-hourly data), every track
position is used, as before.

.. _ts-remote-acquisition:

3. Remote Acquisition
---------------------

The methods themselves (accounts, coverage, limits) are described in :doc:`remote_data_acquisition`.

CDS API: authentication fails or HTTP 403
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Checklist**:

1. **The credentials file exists and is private** — check without printing its content:

   .. code-block:: bash

      ls -l ~/.cdsapirc        # should exist, mode -rw------- (600)
      chmod 600 ~/.cdsapirc

2. **Its format is** (open it in an editor; do not paste the key into terminals, issues or logs):

   .. code-block:: text

      url: https://cds.climate.copernicus.eu/api
      key: <your-personal-access-token>

3. **The key is current**: compare it with the one shown in your CDS profile after logging in.
4. **The ERA5 licence is accepted** in the CDS web interface; without it requests are refused with HTTP 403.

CDS API: "Request too large" or cost limits exceeded
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 37 levels are always requested, so reduce the hours or the area:

- use a larger ``--time-resolution`` (6 or 12 hours);
- shorten the track period or reduce its extent.

.. code-block:: bash

   python lorenzcycletoolkit.py case.nc -t -r --cdsapi --time-resolution 6 --trackfile inputs/track_case

A track sampled more often than the requested hours is aligned automatically (see :ref:`ts-track-alignment`).

CDS API: download hangs or is very slow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the time is usually spent in the CDS queue, which varies from minutes to hours. Use ``-v`` to follow each
day, check the CDS service status page, and reduce the volume as above. **Stopping the client does not cancel the
job**: it may stay queued or running and use the account's quota; check the account for pending requests before
submitting again.

``--ncar``: the date is not available
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GDEX publishes ERA5 months after the CDS. The last day of the track is requested first, so an unpublished date
fails at once, before anything is downloaded, with a message suggesting ``--cdsapi``. Use ``--cdsapi`` for recent
cases.

``--ncar``: timeouts, HTTP 503 or 504
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each request is attempted up to 4 times. If all attempts fail, run the same command again later: the responses
already verified (in ``<file>.parts/``) are reused and only the missing ones are requested. In the D-Bench
tests such failures occurred on some domains of 1.9–3.9 GiB; their cause has not been established.

``--ncar``: the existing file is not reused
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An existing file is never overwritten, and it is reused only if its provenance record matches it and it covers the
whole request with the same variable names and numeric type. The message gives the reason. Use another file name
for a new download, or run without ``--ncar`` to read the file as a plain local input.

``--ncar``: "another retrieval is writing" the file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A lock file ``<file>.lock`` prevents two downloads of the same file (for instance two runs started together).
Wait for the first one, or use different file names. A killed download does not block the next one; on systems
where the lock is left behind, the message says how to remove it.

"request not made" before a ``--cdsapi`` or ``--ncar`` download
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The track was checked against the hours of ``--time-resolution`` before downloading and the request could not serve
it (for example track times that do not fall on the requested hours). Nothing was downloaded. Choose a
``--time-resolution`` whose hours fall on the track times; see :ref:`ts-track-alignment`.

.. _ts-performance:

4. Performance and Memory
-------------------------

"MemoryError" or system becomes unresponsive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: the dataset is too large for the available memory.

**Solutions**: pre-process the data to what the analysis needs (see "Input Data Requirements" in :doc:`usage`):

1. **Smaller domain**:

   .. code-block:: python

      import xarray as xr

      ds = xr.open_dataset('large_file.nc')
      ds_subset = ds.sel(latitude=slice(-50, -10), longitude=slice(-70, -30))
      ds_subset.to_netcdf('smaller_file.nc')

2. **Shorter period**: ``ds.sel(time=slice('2020-01-01', '2020-01-07'))``
3. **Fewer times**, e.g. 6-hourly instead of 3-hourly: ``ds.sel(time=ds.time.dt.hour.isin([0, 6, 12, 18]))``
4. **Only the five required variables** and their coordinates.
5. **Fewer levels**, e.g. ``ds.sel(level=slice(1000, 100))``. This changes the integrated layer and therefore the
   results; see :ref:`rda-science` before relying on reduced level sets.

**See also**: "System Requirements and Performance" in :doc:`usage`.

Analysis takes too long
~~~~~~~~~~~~~~~~~~~~~~~

1. Use 6-hourly instead of hourly data when the study allows it.
2. Reduce the domain to the region of interest.
3. Omit ``-p`` for initial runs.
4. Use ``-f`` or ``-t`` instead of ``-c`` when an interactive choice is not needed.
5. Prefer an uncompressed, already subset NetCDF file: compressed files with large blocks can be read much more
   slowly (an example is given in :ref:`rda-comparison`).

.. _ts-results-figures:

5. Results and Figures
----------------------

The run ends with "PARTIAL output" and exit status 3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: With ``-p``, the log ends with "Analysis complete with PARTIAL output" and the command exits with status 3.

**Meaning**: The computation is complete and the CSV results are written, but at least one requested figure was not generated. ``run_status.json`` in the results directory lists, for each figure, its state and the cause (``reason`` and ``error``). The other figures were still attempted. Use ``--allow-partial-plots`` if a workflow must accept partial figures with exit status 0; the run is still reported as partial.

Maps fail when working offline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: ``box_limits_map``, ``domain_maps`` or ``track_map`` fail with a ``URLError`` (for example "nodename nor servname provided").

**Cause**: The maps use Natural Earth shapefiles (coastline, land, ocean and land borders at 110m, 50m or 10m depending on the map extent; countries, states and cities at 50m). Cartopy downloads a missing file the first time it is needed, which fails without network access.

**Solution**: Set ``CARTOPY_DATA_DIR`` to a directory holding the files listed in ``tests/data/naturalearth_manifest.tsv`` (layout ``shapefiles/natural_earth/{physical,cultural}/``), or run once with network access so that Cartopy stores them. Use ``--offline-maps`` to forbid any download attempt: the maps are then drawn from local data only, and a missing file is reported as ``NaturalEarthOfflineError``. The computation does not depend on the maps: without them the results are the same and the run is reported as partial.

The life-cycle periods and the Lorenz Phase Space are unavailable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: In the moving framework, ``periods`` is ``unavailable`` ("usable time steps < minimum ...") and ``lec_periods`` and ``lps`` are unavailable because they depend on it.

**Cause**: The cyclophaser filter and smoothing windows are derived from the length of the series; cyclophaser 1.7.3 rejects series shorter than 8 time steps with the high-resolution or processed-vorticity options (3 with the low-resolution ones). Such series, like the 5 hourly steps of the ERA5 sample, are not passed to cyclophaser. The CSV results are not affected and the run is complete (exit status 0).

**Note**: Longer series are accepted by cyclophaser, but acceptance does not make the phases physically meaningful; ``run_status.json`` records this, and for fewer than 20 steps it also notes that cyclophaser replaces the whole band-pass filtered series with its low-pass series.

The periods fail with CyclophaserCompatibilityError
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause**: The minimum lengths are known only for cyclophaser 1.7.3, the options of the toolkit and regular time steps of 1, 3 or 6 hours. With another version, other options or another time step the periods are reported as failed without calling cyclophaser, rather than applying limits that may be wrong. Use the tested cyclophaser version, or ask for the limits to be established for your configuration.

.. _ts-faq:

6. Frequently Asked Questions
-----------------------------

Which acquisition method should I use?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It depends on the purpose: recent dates, account, transfer volume, tolerance to failures, and whether the study
depends on magnitudes or vertical structure. See :ref:`rda-conclusions`.

Do I need all 37 levels?
~~~~~~~~~~~~~~~~~~~~~~~~

Not necessarily. In the case assessed, 13 or 17 levels gave the same time evolution, signs and timing of Ca, Ce and
Ck, but magnitudes changed (mainly Ca, about 5 % in the time mean). For magnitudes, budgets or vertical structure use
the same full level set throughout. See :ref:`rda-science`.

Will CDS and NCAR files give identical results?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not always to the last digit. The fields agree within tight tolerances, but small temperature differences can move
tendency and residual terms beyond 10\ :sup:`-4` of each column. Do not mix sources in one comparison that needs
that level of agreement. See :ref:`rda-precision`.

Is the track ever interpolated?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No. When the data are sampled less often than the track, the positions at the data times are used unchanged; in any
other mismatch the run stops. See :ref:`ts-track-alignment`.

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

.. _ts-reporting:

7. Reporting a Problem
----------------------

If your problem isn't covered here:

1. **Enable verbose logging**:

   .. code-block:: bash

      python lorenzcycletoolkit.py your_file.nc -r -f -v

2. **Check the log file** in ``LEC_Results/*/log.*`` and ``run_status.json`` in the same directory.

3. **Review the documentation**: :doc:`configuration`, :doc:`usage`, :doc:`remote_data_acquisition`,
   :doc:`examples_and_tutorials`.

4. **Seek support**:

   - GitHub Issues: https://github.com/daniloceano/LorenzCycleToolkit/issues
   - Email: danilo.oceano@gmail.com

When reporting issues, please include:

- Full error message and traceback
- Your command line
- Relevant parts of the log file
- Description of your input data (source, resolution, coverage)
- LorenzCycleToolkit version: check with ``git describe --tags``

**Never include credentials** (``~/.cdsapirc`` content, API keys or tokens) in a report, log excerpt or screenshot;
remove them before sharing.

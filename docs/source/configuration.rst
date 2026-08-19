Configuration
=============

The LorenzCycleToolkit requires specific configuration files for both fixed and moving frameworks. Below are examples and descriptions of the required configuration files.

Namelist File
-------------

**Important**: You must create a ``namelist`` file in the ``inputs/`` directory before running the toolkit. This file is not tracked by version control, as it is customized by each user for their specific dataset.

The ``namelist`` file specifies the variable names and units used in the input NetCDF file. The toolkit provides several preset namelist files that you can use as templates:

- ``inputs/namelist_ERA5-cdsapi`` - For ERA5 data from CDS API
- ``inputs/namelist_NCEP-R1`` - For NCEP Reanalysis 1 data
- ``inputs/namelist_NCEP-R2`` - For NCEP Reanalysis 2 data
- ``inputs/namelist_ERA5`` - For ERA5 data (other sources)
- ``inputs/namelist_MPAS-A`` - For MPAS-A model data

**Creating your namelist file**:

1. Choose the preset that best matches your data source
2. Copy it to ``inputs/namelist``:

   .. code-block:: bash

      cp inputs/namelist_NCEP-R2 inputs/namelist

3. Edit if needed to match your specific dataset variable names

Here is an example ``namelist`` file for the NCEP-R2 dataset:

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

Box Limits File (Optional)
--------------------------

The `box_limits` file defines the spatial domain for the fixed framework. It is a CSV file with the following format:

.. code-block:: text

    min_lon;-60
    max_lon;-30
    min_lat;-42.5
    max_lat;-17.5

Track File (Optional)
---------------------

For the moving framework, a `track_file` can be used to define the system's center over time. Optionally, length and width columns can be added to adjust the domain size. Here is an example `track_file`:

.. code-block:: text

    time;lon;lat;length;width
    2021-01-01 00:00:00;-50;0;10;10
    2021-01-01 06:00:00;-49;1;10;10
    2021-01-01 12:00:00;-48;2;10;10

**Note**: The track file uses semicolon (``;``) as the delimiter and must include a ``time`` column that can be parsed as datetime.

Namelist for ERA5 Data (CDS API)
---------------------------------

When using the ``--cdsapi`` flag to automatically download ERA5 data, the toolkit **automatically** uses the ERA5-compatible namelist (``inputs/namelist_ERA5-cdsapi``). You do not need to manually create or configure the namelist file when using this feature.

The ``inputs/namelist_ERA5-cdsapi`` preset contains:

.. code-block:: text

    ;standard_name;Variable;Units
    Air Temperature;air_temperature;t;K
    Geopotential;geopotential;z;m**2/s**2
    Omega Velocity;omega;w;Pa/s
    Eastward Wind Component;eastward_wind;u;m/s
    Northward Wind Component;northward_wind;v;m/s
    Longitude;;longitude
    Latitude;;latitude
    Time;;valid_time
    Vertical Level;;pressure_level

If you have ERA5 data from other sources (not downloaded via ``--cdsapi``), you should manually create your namelist:

.. code-block:: bash

   cp inputs/namelist_ERA5-cdsapi inputs/namelist

Next Steps
----------

With your configuration files ready, you can:

- See :doc:`usage` for command-line options and framework selection
- Follow :doc:`examples_and_tutorials` for complete workflows with sample data
- Check :doc:`cdsapi_example` if you need to automatically download ERA5 data
- Review :doc:`results` to understand the output structure


Configuration
=============

The LorenzCycleToolkit requires specific configuration files for both fixed and moving frameworks. Below are examples and descriptions of the required configuration files.

Namelist File
-------------

The `namelist` file specifies the variable names and units used in the input NetCDF file. Here is an example `namelist` file for the NCEP-R2 dataset:

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

    time, lon, lat, length, width
    2021-01-01 00:00:00, -50, 0, 10, 10
    2021-01-01 06:00:00, -49, 1, 10, 10
    2021-01-01 12:00:00, -48, 2, 10, 10


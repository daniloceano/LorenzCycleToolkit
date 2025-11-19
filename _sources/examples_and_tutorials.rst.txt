Examples and Tutorials
======================

The preferred method for using the LorenzCycleToolkit is through command line arguments executed from the top-level directory of the project. Below are examples demonstrating how to configure and run the toolkit using different frameworks.

Fixed (Eulerian) Framework Example
----------------------------------
1. Prepare the `inputs/box_limits` file.  
   The `box_limits` file should define the spatial domain as a CSV file with the following format:

   .. code-block:: text

       min_lon;-60
       max_lon;-30
       min_lat;-42.5
       max_lat;-17.5

2. Prepare the `inputs/namelist` file.  
   The `namelist` file specifies the variable names and units used in the input NetCDF file. Here is an example for the National Center for Environmental Protection Reanalysis 2 (NCEP-R2) dataset:

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

3. From the top-level directory of the project, run the following command::

       python lorenzcycletoolkit.py samples/testdata_NCEP-R2.nc -r -f

   This will execute the LorenzCycleToolkit using the fixed (Eulerian) framework.

Moving (Semi-Lagrangian) Framework Example  
------------------------------------------

1. Prepare the `inputs/track_file`.  
   The `track_file` should define the system's center over time. Optionally, add length and width columns to adjust the domain size. Example:

   .. code-block:: text

       time;Lat;Lon
       2005-08-08-0000;-22.5;-45.0
       2005-08-08-0600;-22.6;-44.9
       2005-08-08-1200;-22.7;-44.8
       2005-08-08-1800;-22.8;-44.7
       2005-08-09-0000;-22.9;-44.6

   A default track file for testing purposes is also provided and can be used by running the following command:

   .. code-block:: shell

       cp inputs/track_testdata_NCEP-R2 inputs/track

   The program will always use the `inputs/track` file for processing. Therefore, if you have a custom track file, you need to overwrite the `inputs/track` file with the desired track data before running the program.

2. From the top-level directory of the project, run the following command::

       python lorenzcycletoolkit.py samples/testdata__NCEP-R2.nc -r -t

   This will execute the LorenzCycleToolkit using the moving (Semi-Lagrangian) framework.

Interactive Domain Selection Example
------------------------------------
1. From the top-level directory of the project, run the following command::

       python lorenzcycletoolkit.py samples/testdata__NCEP-R2.nc -r -c

   This command enables interactive map-based domain selection. First, define a slice of global data to enhance visualization and processing speed. Then, select the computational area for each timestep interactively.

   Each timestep prompts the user to select a computational area:

   .. image:: https://user-images.githubusercontent.com/56005607/214922008-5b7c094f-c160-4415-a528-07cc58730827.png
      :width: 350

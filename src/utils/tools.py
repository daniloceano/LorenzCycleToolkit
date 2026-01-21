# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    tools.py                                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:33:03 by daniloceano       #+#    #+#              #
#    Updated: 2026/01/21 10:19:30 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import argparse
import logging
import os
from datetime import timedelta

import dask
import numpy as np
import pandas as pd
import xarray as xr
from metpy.units import units

from .select_area import slice_domain


def initialize_logging(results_subdirectory, args):
    """
    Initializes the logging configuration for the application.

    Args:
        results_subdirectory (str): Directory path to save the log file.
        args (object): The argparse object containing the command line arguments.
    """

    verbose = True if args.verbosity else False

    # Set root logger to higher severity level (INFO or ERROR)
    root_log_level = logging.ERROR if not verbose else logging.INFO
    logging.basicConfig(
        level=root_log_level, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Create a separate logger for the application
    app_logger = logging.getLogger("lorenzcycletoolkit")
    app_log_level = logging.DEBUG if verbose else logging.INFO
    app_logger.setLevel(app_log_level)
    # Prevent the logger from propagating messages to the root logger
    app_logger.propagate = False

    # Create file handler for saving logs
    log_file_name = f'log.{os.path.basename(args.infile).split(".")[0]}'
    log_file = os.path.join(results_subdirectory, log_file_name)
    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setLevel(app_log_level)
    file_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(file_formatter)
    app_logger.addHandler(file_handler)

    # Create a console handler for app logger
    console_handler = logging.StreamHandler()
    console_handler.setLevel(app_log_level)
    console_handler.setFormatter(file_formatter)
    app_logger.addHandler(console_handler)

    return app_logger


def convert_longitude_range(df: xr.Dataset, lon_indexer: str) -> xr.Dataset:
    """
    Convert longitude range from 0-360 degrees to -180 to 180 degrees.

    This function modifies the longitude coordinates of the provided xarray Dataset
    and sorts the data based on these updated longitudes.

    Args:
        df (xr.Dataset): The dataset containing longitude coordinates.
        lon_indexer (str): The name of the longitude coordinate in the dataset.

    Returns:
        xr.Dataset: The dataset with updated longitude coordinates.
    """
    df.coords[lon_indexer] = (df.coords[lon_indexer] + 180) % 360 - 180
    df = df.sortby(df[lon_indexer])
    return df


def find_extremum_coordinates(
    data: xr.DataArray, lat: xr.DataArray, lon: xr.DataArray, variable: str
) -> tuple:
    """
    Finds the indices of extremum values for a given variable in a dataset.

    Args:
        data (xr.DataArray): DataArray containing the data for calculation.
        lat (xr.DataArray): DataArray containing the latitudes.
        lon (xr.DataArray): DataArray containing the longitudes.
        variable (str): Name of the variable to find extremum indices.

    Returns:
        tuple: Tuple containing latitude and longitude of the extremum value.
    """
    lat_values = lat.values
    lon_values = lon.values
    min_lat = lat_values.min()

    data = np.array(data)  # Ensure data is a numpy array

    if variable == "min_max_zeta":
        index = np.unravel_index(
            data.argmin() if min_lat < 0 else data.argmax(), data.shape
        )
    elif variable in ["min_hgt", "max_wind"]:
        index = np.unravel_index(
            data.argmin() if variable == "min_hgt" else data.argmax(), data.shape
        )
    else:
        logging.error(f"Invalid variable specified: {variable}")
        raise ValueError(f"Invalid variable specified: {variable}")

    return lat_values[index[0]], lon_values[index[1]]


def get_cdsapi_data(
    args: argparse.Namespace, track: pd.DataFrame, app_logger: logging.Logger
) -> str:
    """
    Retrieve ERA5 pressure level data from the CDS API based on track information.

    This function downloads ERA5 reanalysis data day-by-day to avoid issues with
    large file downloads, then concatenates the daily files into a single output file.
    Temporary daily files are automatically deleted after concatenation.

    Args:
        args (argparse.Namespace): Command-line arguments containing the output file path.
        track (pd.DataFrame): DataFrame with 'Lat', 'Lon' columns and datetime index.
        app_logger (logging.Logger): Logger for the application.

    Returns:
        str: Path to the downloaded NetCDF file.

    Raises:
        FileNotFoundError: If temporary files are not created during download.
        Exception: For other errors during data retrieval or concatenation.
    """
    import math
    import tempfile

    import cdsapi

    # Extract bounding box (lat/lon limits) from track
    min_lat, max_lat = float(track["Lat"].min()), float(track["Lat"].max())
    min_lon, max_lon = float(track["Lon"].min()), float(track["Lon"].max())

    # Apply a 15-degree buffer and round to nearest integer
    buffered_max_lat = math.ceil(max_lat + 15)
    buffered_min_lon = math.floor(min_lon - 15)
    buffered_min_lat = math.floor(min_lat - 15)
    buffered_max_lon = math.ceil(max_lon + 15)

    # Define the area for the request [North, West, South, East]
    area = [buffered_max_lat, buffered_min_lon, buffered_min_lat, buffered_max_lon]

    # Pressure levels in hPa
    pressure_levels = [
        "1", "2", "3", "5", "7", "10", "20", "30", "50", "70",
        "100", "125", "150", "175", "200", "225", "250", "300",
        "350", "400", "450", "500", "550", "600", "650", "700",
        "750", "775", "800", "825", "850", "875", "900", "925",
        "950", "975", "1000",
    ]

    # Variables to retrieve
    variables = [
        "u_component_of_wind",
        "v_component_of_wind",
        "temperature",
        "vertical_velocity",
        "geopotential",
    ]

    # Convert track index to DatetimeIndex and find the date range
    track_datetime_index = pd.DatetimeIndex(track.index)
    first_timestamp = track_datetime_index.min()
    last_timestamp = track_datetime_index.max()

    # Get all unique dates needed for download
    dates = pd.date_range(
        start=first_timestamp.date(),
        end=last_timestamp.date(),
        freq='D'
    )

    # Use time_resolution from args (default is 3 hours)
    time_step = args.time_resolution
    app_logger.debug(f"Using time resolution from args: {time_step} hour(s)")

    # Log track file bounds and requested data bounds
    app_logger.debug(
        f"Track File Limits: lon range: [{min_lon:.2f}, {max_lon:.2f}], "
        f"lat range: [{min_lat:.2f}, {max_lat:.2f}]"
    )
    app_logger.debug(
        f"Buffered Data Bounds: lon range: [{buffered_min_lon}, {buffered_max_lon}], "
        f"lat range: [{buffered_min_lat}, {buffered_max_lat}]"
    )
    app_logger.debug(
        f"Track period: {first_timestamp} to {last_timestamp}"
    )
    app_logger.info(
        f"Will download {len(dates)} day(s) of data: {dates[0].date()} to {dates[-1].date()}"
    )

    # Create temporary directory for daily files
    temp_dir = tempfile.mkdtemp(prefix="cdsapi_daily_")
    app_logger.debug(f"Created temporary directory: {temp_dir}")
    
    # Download data day by day
    app_logger.info(f"Starting data download from CDS API...")
    app_logger.info(f"Expected files to download: {len(dates)}")
    for i, date in enumerate(dates, 1):
        app_logger.info(f"  - Day {i}: {date.date()}")
    
    daily_files = []
    client = cdsapi.Client(timeout=600, retry_max=500)
    
    try:
        for i, date in enumerate(dates, 1):
            year = date.strftime('%Y')
            month = date.strftime('%m')
            day = date.strftime('%d')
            
            # Determine which hours to download for this day
            # For the first day, start from the first timestamp hour
            # For the last day, end at the last timestamp hour
            # For middle days, download all hours
            if date.date() == first_timestamp.date():
                # First day: start from first timestamp hour
                start_hour = first_timestamp.hour
                # Round down to nearest time_step
                start_hour = (start_hour // time_step) * time_step
                if date.date() == last_timestamp.date():
                    # Same day: both start and end on same day
                    end_hour = last_timestamp.hour
                    # Round up to nearest time_step, ensuring at least one timestep
                    end_hour = ((end_hour + time_step - 1) // time_step) * time_step
                    if end_hour <= start_hour:
                        end_hour = start_hour + time_step
                else:
                    # First day but not last: download until end of day
                    end_hour = 24
            elif date.date() == last_timestamp.date():
                # Last day (not first): start from beginning of day
                start_hour = 0
                end_hour = last_timestamp.hour
                # Round up to nearest time_step
                end_hour = ((end_hour + time_step - 1) // time_step) * time_step
                # If last timestamp is at midnight (hour 0), download that hour
                if end_hour == 0:
                    end_hour = time_step
            else:
                # Middle day: download all hours
                start_hour = 0
                end_hour = 24
            
            # Generate times list for this specific day
            times = [f"{hour:02d}:00" for hour in range(start_hour, end_hour, time_step)]
            
            if not times:
                app_logger.warning(f"No times to download for {year}-{month}-{day}, skipping...")
                continue
            
            # Create temporary file for this day
            temp_file = os.path.join(temp_dir, f"era5_{year}{month}{day}.nc")
            daily_files.append(temp_file)
            
            app_logger.info(f"Downloading day {i}/{len(dates)}: {year}-{month}-{day}")
            app_logger.debug(f"  Times: {times[0]} to {times[-1]} ({len(times)} timesteps)")
            
            # Prepare request for single day
            request = {
                "product_type": "reanalysis",
                "format": "netcdf",
                "variable": variables,
                "pressure_level": pressure_levels,
                "year": year,
                "month": month,
                "day": day,
                "time": times,
                "area": area,
            }
            
            # Download data for this day
            try:
                client.retrieve(
                    "reanalysis-era5-pressure-levels",
                    request,
                    temp_file,
                )
                
                # Verify file was created
                if not os.path.exists(temp_file):
                    raise FileNotFoundError(
                        f"Daily file not created at expected path: {temp_file}"
                    )
                
                file_size = os.path.getsize(temp_file)
                app_logger.info(f"  ✓ Downloaded successfully: {file_size / (1024**2):.2f} MB")
                
            except Exception as e:
                app_logger.error(f"✗ Error downloading data for {year}-{month}-{day}: {e}")
                raise
        
        app_logger.info("=" * 60)
        app_logger.info("All daily files downloaded successfully!")
        app_logger.info(f"Total files: {len(daily_files)}")
        app_logger.info("=" * 60)
        
        # Concatenate daily files
        app_logger.info("Concatenating daily files into single dataset...")
        datasets = []
        for idx, daily_file in enumerate(daily_files, 1):
            app_logger.debug(f"  Loading file {idx}/{len(daily_files)}: {os.path.basename(daily_file)}")
            ds = xr.open_dataset(daily_file)
            datasets.append(ds)
        
        # Concatenate along time dimension
        app_logger.info("Merging datasets along time dimension...")
        combined_ds = xr.concat(datasets, dim='valid_time')
        
        # Save concatenated dataset
        app_logger.info(f"Saving final dataset to: {args.infile}")
        combined_ds.to_netcdf(args.infile)
        combined_ds.close()
        
        # Close all datasets
        for ds in datasets:
            ds.close()
        
        # Verify the final file was created
        if not os.path.exists(args.infile):
            raise FileNotFoundError(
                f"CDS API file not created at expected path: {args.infile}"
            )
        
        final_file_size = os.path.getsize(args.infile)
        app_logger.info("=" * 60)
        app_logger.info("✓ DATA DOWNLOAD COMPLETED SUCCESSFULLY!")
        app_logger.info(f"Final file: {args.infile}")
        app_logger.info(f"File size: {final_file_size / (1024**2):.2f} MB")
        app_logger.info(f"Time coverage: {first_timestamp} to {last_timestamp}")
        app_logger.info("=" * 60)
        
    finally:
        # Clean up temporary files and directory
        app_logger.debug("Cleaning up temporary files...")
        for temp_file in daily_files:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                    app_logger.debug(f"Deleted temporary file: {temp_file}")
            except Exception as e:
                app_logger.warning(f"Could not delete temporary file {temp_file}: {e}")
        
        try:
            os.rmdir(temp_dir)
            app_logger.debug(f"Deleted temporary directory: {temp_dir}")
        except Exception as e:
            app_logger.warning(f"Could not delete temporary directory {temp_dir}: {e}")
    
    return args.infile


def get_data(args: argparse.Namespace, app_logger: logging.Logger) -> xr.Dataset:
    """
    Opens a NetCDF file and extracts variables specified in a CSV file.

    Args:
        args (argparse.Namespace): Command-line arguments.
        variable_list_df (pd.DataFrame): DataFrame containing the variables to extract.
        app_logger (logging.Logger): Logger for the application.

    Returns:
        xr.Dataset: Dataset containing extracted variables.

    Raises:
        FileNotFoundError: If CSV file or NetCDF file is not found.
        Exception: For other errors occurring during file opening.
    """
    infile = args.infile

    if args.cdsapi:
        if not os.path.exists(infile):
            app_logger.debug(
                "CDS API data not found. Attempting to retrieve data from CDS API..."
            )
            track = pd.read_csv(
                args.trackfile, parse_dates=[0], delimiter=";", index_col="time"
            )
            infile = get_cdsapi_data(args, track, app_logger)
            app_logger.debug(f"CDS API data: {infile}")
        else:
            app_logger.info("CDS API data found.")

    app_logger.debug("Opening input data... ")
    try:
        with dask.config.set(array={"slicing": {"split_large_chunks": True}}):
            data = xr.open_dataset(infile)
    except FileNotFoundError:
        app_logger.error(
            "Could not open file. Check if path, namelist file, and file format (.nc) are correct."
        )
        raise
    except Exception as e:
        app_logger.exception("An exception occurred: {}".format(e))
        raise
    app_logger.debug("Ok.")

    return data


def process_data(
    data: xr.Dataset,
    args: argparse.Namespace,
    variable_list_df: pd.DataFrame,
    app_logger: logging.Logger,
) -> xr.Dataset:
    """
    Process the given data and return a modified dataset.

    Parameters:
    - data: A dataset containing the data to be processed (type: xr.Dataset).
    - args: An argparse.Namespace object containing the command line arguments (type: argparse.Namespace).
    - variable_list_df: A DataFrame containing a list of variables (type: pd.DataFrame).
    - app_logger: A logger object for logging debug messages (type: logging.Logger).

    Returns:
    - data: A modified dataset after processing (type: xr.Dataset).
    """
    # Select only data matching the track dates
    if args.track:
        app_logger.debug("Selecting only data matching the track dates... ")
        track_file = args.trackfile
        track = pd.read_csv(
            track_file, parse_dates=[0], delimiter=";", index_col="time"
        )
        TimeIndexer = variable_list_df.loc["Time"]["Variable"]
        # If using CDS API, resample track to data time step
        if args.cdsapi:
            time_delta = int(
                (data[TimeIndexer][1] - data[TimeIndexer][0]) / np.timedelta64(1, "h")
            )
            track = track[track.index.hour % time_delta == 0]
        data = data.sel({TimeIndexer: track.index.values})

    if (
        data[variable_list_df.loc["Longitude"]["Variable"]].min() < -180
        or data[variable_list_df.loc["Longitude"]["Variable"]].max() > 180
    ):
        data = convert_longitude_range(
            data, variable_list_df.loc["Longitude"]["Variable"]
        )

    LonIndexer = variable_list_df.loc["Longitude"]["Variable"]
    LatIndexer = variable_list_df.loc["Latitude"]["Variable"]
    LevelIndexer = variable_list_df.loc["Vertical Level"]["Variable"]

    app_logger.debug("Assigning geospatial coordinates in radians... ")
    data = data.assign_coords({"rlats": np.deg2rad(data[LatIndexer])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data[LatIndexer]))})
    data = data.assign_coords({"rlons": np.deg2rad(data[LonIndexer])})
    app_logger.debug("Ok.")

    # Drop unnecessary dimensions
    if 'expver' in data.coords:
        data = data.drop('expver')
    if 'number' in data.coords:
        data = data.drop('number')

    levels_Pa = (
        data[LevelIndexer] * units(str(data[LevelIndexer].units))
    ).metpy.convert_units("Pa")
    data = data.assign_coords({LevelIndexer: levels_Pa})

    data = (
        data.sortby(LonIndexer)
        .sortby(LevelIndexer, ascending=True)
        .sortby(LatIndexer, ascending=True)
    )

    lowest_level = float(data[LevelIndexer].max())
    data = data.sel({LevelIndexer: slice(1000, lowest_level)})

    if args.mpas:
        data = data.drop_dims("standard_height")

    app_logger.debug("Data opened successfully.")
    return data


def prepare_data(
    args, varlist: str = "inputs/namelist", app_logger: logging.Logger = None
) -> xr.Dataset:
    """
    Prepare the data for further analysis.

    Parameters:
        args (object): The arguments for the function.
        varlist (str): The file path to the variable list file (namelist).
        app_logger (logging.Logger): The logger for the application.

    Returns:
        method (str): The method used for the analysis: fixed, track or choose.
        xr.Dataset: The prepared dataset for analysis.
    """

    app_logger.debug(f"Variables specified by the user in: {varlist}")
    app_logger.debug(f"Attempting to read {varlist} file...")
    try:
        variable_list_df = pd.read_csv(varlist, sep=";", index_col=0, header=0)
    except FileNotFoundError:
        app_logger.error("The 'namelist' text file could not be found.")
        raise
    except pd.errors.EmptyDataError:
        app_logger.error("The 'namelist' text file is empty.")
        raise
    app_logger.debug("List of variables found:\n" + str(variable_list_df))

    data = get_data(args, app_logger)

    # Check if variable_list_df matches the data
    if not set(variable_list_df["Variable"]).issubset(set(data.variables)):
        app_logger.error(
            "The variable list does not match the data. Check if the 'namelist' text file is correct."
        )
        raise ValueError("'namelist' text file does not match the data.")

    processed_data = process_data(data, args, variable_list_df, app_logger)
    sliced_data = slice_domain(processed_data, args, variable_list_df)
    return sliced_data

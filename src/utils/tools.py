# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    tools.py                                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:33:03 by daniloceano       #+#    #+#              #
#    Updated: 2024/01/18 08:21:32 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import dask
import logging
import xarray as xr
import pandas as pd
import numpy as np
import argparse
from metpy.units import units
from .select_area import slice_domain

def initialize_logging(results_subdirectory, verbose=False):
    """
    Initializes the logging configuration for the application.

    Args:
        results_subdirectory (str): Directory path to save the log file.
        verbose (bool): Flag to set logging level to DEBUG for detailed logging.
    """
    # Set root logger to higher severity level (INFO or ERROR)
    root_log_level = logging.ERROR if not verbose else logging.INFO
    logging.basicConfig(level=root_log_level, format='%(asctime)s - %(levelname)s - %(message)s')

    # Create a separate logger for the application
    app_logger = logging.getLogger('lorenz_cycle')
    app_log_level = logging.DEBUG if verbose else logging.INFO
    app_logger.setLevel(app_log_level)
    app_logger.propagate = False  # Prevent the logger from propagating messages to the root logger

    # Create file handler for saving logs
    log_file = os.path.join(results_subdirectory, 'log.txt')
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(app_log_level)
    file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
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

def find_extremum_coordinates(data: xr.DataArray, lat: xr.DataArray, lon: xr.DataArray, variable: str) -> tuple:
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

    if variable == 'min_max_zeta':
        index = np.unravel_index(data.argmin() if min_lat < 0 else data.argmax(), data.shape)
    elif variable in ['min_hgt', 'max_wind']:
        index = np.unravel_index(data.argmin() if variable == 'min_hgt' else data.argmax(), data.shape)
    else:
        logging.error(f"Invalid variable specified: {variable}")
        raise ValueError(f"Invalid variable specified: {variable}")

    return lat_values[index[0]], lon_values[index[1]]

def get_cdsapi_data(args: argparse.Namespace, track: pd.DataFrame, app_logger: logging.Logger) -> xr.Dataset:
    import cdsapi

    # Extract bounding box (lat/lon limits) from track
    min_lat, max_lat = track['Lat'].min(), track['Lat'].max()
    min_lon, max_lon = track['Lon'].min(), track['Lon'].max()

    pressure_levels = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70',
                       '100', '125', '150', '175', '200', '225', '250', '300', '350',
                       '400', '450', '500', '550', '600', '650', '700', '750', '775',
                       '800', '825', '850', '875', '900', '925', '950', '975', '1000']
    
    variables = ["u_component_of_wind", "v_component_of_wind", "temperature",
                 "vertical_velocity", "geopotential"]
    

    # Convert unique dates to string format for the request
    dates = track.index.strftime('%Y%m%d').unique().tolist()
    time_range = f"{dates[0]}/{dates[-1]}"
    area = f"{max_lat+15}/{min_lon-15}/{min_lat-15}/{max_lon+15}"
    time_step = str(int((track.index[1] - track.index[0]).total_seconds() / 3600))
    app_logger.debug(f"Requesting data for area: {area}, time range: {time_range} and time step: {time_step}...")
    
    # Load ERA5 data
    app_logger.info("Retrieving data from CDS API...")
    c = cdsapi.Client()
    c.retrieve(
        "reanalysis-era5-pressure-levels",
        {
            "product_type":"reanalysis",
            "format": "netcdf",
            "pressure_level":pressure_levels,
            "date": time_range,
            "area": area,
            'time':f'00/to/23/by/{time_step}',
            "variable":variables,
        }, args.infile # save file as passed in arguments
    )

    if not os.path.exists(args.infile):
        raise FileNotFoundError("CDS API file not created.")
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
            app_logger.debug("CDS API data not found. Attempting to retrieve data from CDS API...")
            track = pd.read_csv('inputs/track', parse_dates=[0], delimiter=';', index_col='time')
            infile = get_cdsapi_data(args, track, app_logger)
            app_logger.debug(f"CDS API data: {infile}")
        else:
            app_logger.info("CDS API data found.")

    app_logger.debug("Opening input data... ")
    try:
        with dask.config.set(array={'slicing': {'split_large_chunks': True}}):
            data = xr.open_dataset(infile)
    except FileNotFoundError:
        app_logger.error("Could not open file. Check if path, fvars file, and file format (.nc) are correct.")
        raise
    except Exception as e:
        app_logger.exception("An exception occurred: {}".format(e))
        raise
    app_logger.debug("Ok.")

    return data

def process_data(data: xr.Dataset, args: argparse.Namespace, variable_list_df: pd.DataFrame, app_logger: logging.Logger) -> xr.Dataset:
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
        track = pd.read_csv('inputs/track', parse_dates=[0], delimiter=';', index_col='time')
        data = data.sel(time=track.index.values)

    else:
        data = convert_longitude_range(data, variable_list_df.loc['Longitude']['Variable'])

    LonIndexer = variable_list_df.loc["Longitude"]["Variable"]
    LatIndexer = variable_list_df.loc["Latitude"]["Variable"]
    LevelIndexer = variable_list_df.loc["Vertical Level"]["Variable"]

    app_logger.debug("Assigning geospatial coordinates in radians... ")
    data = data.assign_coords({"rlats": np.deg2rad(data[LatIndexer])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data[LatIndexer]))})
    data = data.assign_coords({"rlons": np.deg2rad(data[LonIndexer])})
    app_logger.debug("Ok.")

    levels_Pa = (data[LevelIndexer] * units(str(data[LevelIndexer].units))).metpy.convert_units("Pa")
    data = data.assign_coords({LevelIndexer: levels_Pa})
    
    data = data.sortby(LonIndexer).sortby(LevelIndexer, ascending=True).sortby(LatIndexer, ascending=True)

    lowest_level = float(data[LevelIndexer].max())
    data = data.sel({LevelIndexer: slice(1000, lowest_level)})

    if args.mpas:
        data = data.drop_dims('standard_height')

    app_logger.debug("Data opened successfully.")
    return data

def prepare_data(args, varlist: str = 'inputs/fvars', app_logger: logging.Logger = None) -> xr.Dataset:
    """
    Prepare the data for further analysis.

    Parameters:
        args (object): The arguments for the function.
        varlist (str): The file path to the variable list file (fvars).
        app_logger (logging.Logger): The logger for the application.

    Returns:
        method (str): The method used for the analysis: fixed, track or choose.
        xr.Dataset: The prepared dataset for analysis.
    """

    app_logger.debug(f"Variables specified by the user in: {varlist}")
    app_logger.debug(f"Attempting to read {varlist} file...")
    try:
        variable_list_df = pd.read_csv(varlist, sep=';', index_col=0, header=0)
    except FileNotFoundError:
        app_logger.error("The 'fvar' text file could not be found.")
        raise
    except pd.errors.EmptyDataError:
        app_logger.error("The 'fvar' text file is empty.")
        raise
    app_logger.debug("List of variables found:\n" + str(variable_list_df))


    data = get_data(args, app_logger)
    processed_data = process_data(data, args, variable_list_df, app_logger)
    sliced_data = slice_domain(processed_data, args, variable_list_df)
    return sliced_data
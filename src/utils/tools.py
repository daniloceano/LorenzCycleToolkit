# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    tools.py                                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:33:03 by daniloceano       #+#    #+#              #
#    Updated: 2024/01/16 10:22:44 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import dask
import logging
import xarray as xr
import pandas as pd
import numpy as np
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

def get_climate_data():
    import climetlab as cml

    dataset_name = "era5-pressure-levels"
    variables = ["temperature", "u_component_of_wind", "v_component_of_wind", "vertical_velocity", "geopotential"]  
    date_range = ["2020-01-01", "2020-01-31"] 
    pressure_levels = [500, 700, 850]  
    area = [50, -40, 20, 60] 

    # Retrieve the dataset
    dataset = cml.load_dataset(
        dataset_name,
        variable=variables,
        pressure_level=pressure_levels,
        date_range=date_range,
        area=area
    )

def get_data(infile: str, varlist: str, climet: bool = False) -> xr.Dataset:
    """
    Opens a NetCDF file and extracts variables specified in a CSV file.

    Args:
        infile (str): Path to the NetCDF file.
        varlist (str): Path to the CSV file listing variables to extract.
        climet (bool): Whether or not to use climetlab for retrieving data.

    Returns:
        xr.Dataset: Dataset containing extracted variables.

    Raises:
        FileNotFoundError: If CSV file or NetCDF file is not found.
        Exception: For other errors occurring during file opening.
    """
    logging.debug(f"Variables specified by the user in: {varlist}")
    logging.debug(f"Attempting to read {varlist} file...")

    try:
        variable_list_df = pd.read_csv(varlist, sep=';', index_col=0, header=0)
    except FileNotFoundError:
        logging.error("The 'fvar' text file could not be found.")
        raise
    except pd.errors.EmptyDataError:
        logging.error("The 'fvar' text file is empty.")
        raise
    logging.debug("List of variables found:\n" + str(variable_list_df))

    LonIndexer = variable_list_df.loc["Longitude"]["Variable"]
    LatIndexer = variable_list_df.loc["Latitude"]["Variable"]
    LevelIndexer = variable_list_df.loc["Vertical Level"]["Variable"]

    if climet:
        data = get_climate_data()

    logging.debug("Opening input data... ")
    try:
        with dask.config.set(array={'slicing': {'split_large_chunks': True}}):
            data = convert_longitude_range(
                xr.open_dataset(infile),
                variable_list_df.loc['Longitude']['Variable']
            )
    except FileNotFoundError:
        logging.error("Could not open file. Check if path, fvars file, and file format (.nc) are correct.")
        raise
    except Exception as e:
        logging.exception("An exception occurred: {}".format(e))
        raise
    logging.debug("Ok.")

    logging.debug("Assigning geospatial coordinates in radians... ")
    data = data.assign_coords({"rlats": np.deg2rad(data[LatIndexer])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data[LatIndexer]))})
    data = data.assign_coords({"rlons": np.deg2rad(data[LonIndexer])})
    logging.debug("Ok.")

    levels_Pa = (data[LevelIndexer] * units(str(data[LevelIndexer].units))).metpy.convert_units("Pa")
    data = data.assign_coords({LevelIndexer: levels_Pa})
    
    data = data.sortby(LonIndexer).sortby(LevelIndexer, ascending=True).sortby(LatIndexer, ascending=True)

    lowest_level = float(data[LevelIndexer].max())
    data = data.sel({LevelIndexer: slice(1000, lowest_level)})

    logging.debug("Data opened successfully.")
    return data

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

def prepare_data(args, fvars: str = 'inputs/fvars', climet: bool = False) -> xr.Dataset:
    """
    Prepare the data for further analysis.

    Parameters:
        args (object): The arguments for the function.
        fvars (str): The file path to the fvars.
        climet (bool): Whether or not to use climetlab for retrieving data.

    Returns:
        method (str): The method used for the analysis: fixed, track or choose.
        xr.Dataset: The prepared dataset for analysis.
    """
    data = get_data(args.infile, fvars, climet)
    if args.mpas:
        data = data.drop_dims('standard_height')
    
    return slice_domain(data, args, fvars)
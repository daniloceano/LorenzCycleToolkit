# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    preprocessing.py                                   :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2025/02/26 10:45:00 by daniloceano       #+#    #+#              #
#    Updated: 2025/02/26 10:45:00 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
Data preprocessing and preparation functions for the Lorenz Energy Cycle Toolkit.

This module handles data loading, preprocessing, and preparation for analysis,
including timestep validation, unit conversion, and coordinate transformation.
"""

import argparse
import logging
import os

import dask
import numpy as np
import pandas as pd
import xarray as xr
from metpy.units import units

from .select_area import slice_domain
from .tools import convert_longitude_range, get_cdsapi_data
from .validation import validate_track_file, validate_namelist_file


def get_data(args: argparse.Namespace, app_logger: logging.Logger) -> xr.Dataset:
    """
    Opens a NetCDF file and extracts variables specified in a CSV file.

    Args:
        args (argparse.Namespace): Command-line arguments.
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
                "🌐 CDS API data not found. Attempting to retrieve data from CDS API..."
            )
            # Validate track file format and detect delimiter
            delimiter, _ = validate_track_file(args.trackfile, app_logger)
            track = pd.read_csv(
                args.trackfile, 
                delimiter=delimiter, 
                index_col="time",
                parse_dates=["time"],
                date_format="%Y-%m-%d-%H%M"
            )
            infile = get_cdsapi_data(args, track, app_logger)
            app_logger.debug(f"✅ CDS API data ready: {infile}")
        else:
            app_logger.info("✅ CDS API data already exists, skipping download.")

    app_logger.debug("📂 Opening input data... ")
    try:
        with dask.config.set(array={"slicing": {"split_large_chunks": True}}):
            data = xr.open_dataset(infile)
    except FileNotFoundError:
        app_logger.error("❌ Input file not found!")
        app_logger.error("\n" + "="*70)
        app_logger.error("📁 FILE NOT FOUND ERROR")
        app_logger.error("="*70)
        app_logger.error(f"Looking for: {os.path.abspath(infile)}")
        app_logger.error(f"Current directory: {os.getcwd()}")
        app_logger.error("\n💡 User Solutions:")
        app_logger.error("   1. Check if the file path is correct")
        app_logger.error("   2. Use absolute path instead of relative path")
        app_logger.error("   3. Verify the file exists: ls -lh <path>")
        app_logger.error("\n🔧 Developer Info:")
        app_logger.error(f"   infile argument: {infile}")
        app_logger.error(f"   os.path.exists(infile): {os.path.exists(infile)}")
        app_logger.error("="*70 + "\n")
        raise FileNotFoundError(
            f"Input file not found: {os.path.abspath(infile)}. "
            f"Check path and file existence."
        )
    except OSError as e:
        app_logger.error("❌ Error opening file - OS/File system issue!")
        app_logger.error("\n" + "="*70)
        app_logger.error("🚫 FILE SYSTEM ERROR")
        app_logger.error("="*70)
        app_logger.error(f"File: {os.path.abspath(infile)}")
        app_logger.error(f"Error: {type(e).__name__}: {e}")
        app_logger.error("\n💡 User Solutions:")
        app_logger.error("   1. Check file permissions: ls -l <path>")
        app_logger.error("   2. Verify the file is not corrupted")
        app_logger.error("   3. Check if file format is NetCDF (.nc)")
        app_logger.error("   4. Try: ncdump -h <file> to test file integrity")
        app_logger.error("\n🔧 Developer Info:")
        app_logger.error(f"   Exception type: {type(e).__name__}")
        app_logger.error(f"   Exception message: {str(e)}")
        app_logger.error(f"   File exists: {os.path.exists(infile)}")
        if os.path.exists(infile):
            app_logger.error(f"   File size: {os.path.getsize(infile)} bytes")
            app_logger.error(f"   Is readable: {os.access(infile, os.R_OK)}")
        app_logger.error("="*70 + "\n")
        raise OSError(
            f"Cannot open file {os.path.abspath(infile)}: {e}. "
            f"Check file permissions and format."
        )
    except Exception as e:
        app_logger.error("❌ Unexpected error opening NetCDF file!")
        app_logger.error("\n" + "="*70)
        app_logger.error("⚠️  UNEXPECTED ERROR")
        app_logger.error("="*70)
        app_logger.error(f"File: {os.path.abspath(infile)}")
        app_logger.error(f"Error type: {type(e).__name__}")
        app_logger.error(f"Error message: {str(e)}")
        app_logger.error("\n💡 Possible causes:")
        app_logger.error("   1. File is not a valid NetCDF format")
        app_logger.error("   2. File is corrupted or incomplete")
        app_logger.error("   3. Missing required libraries (netCDF4, xarray)")
        app_logger.error("   4. Memory issues with large files")
        app_logger.error("\n🔧 Developer Info:")
        app_logger.error(f"   Exception: {type(e).__name__}")
        app_logger.error(f"   Message: {str(e)}")
        if hasattr(e, '__traceback__'):
            import traceback
            app_logger.error("   Traceback:")
            for line in traceback.format_tb(e.__traceback__):
                app_logger.error(f"     {line.strip()}")
        app_logger.error("="*70 + "\n")
        raise Exception(
            f"Failed to open NetCDF file {os.path.abspath(infile)}: "
            f"{type(e).__name__}: {e}"
        )
    app_logger.debug("✅ Data opened successfully.")

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
        app_logger.debug("📅 Selecting only data matching the track dates... ")
        track_file = args.trackfile
        
        # Validate track file format and detect delimiter
        delimiter, _ = validate_track_file(track_file, app_logger)
        
        # Read track file with custom date parser for YYYY-MM-DD-HHMM format
        track = pd.read_csv(
            track_file, 
            delimiter=delimiter, 
            index_col="time",
            parse_dates=["time"],
            date_format="%Y-%m-%d-%H%M"
        )
        TimeIndexer = variable_list_df.loc["Time"]["Variable"]
        # Check input data and track time steps
        data_time_delta = int(
            (data[TimeIndexer][1] - data[TimeIndexer][0]) / np.timedelta64(1, "h")
        )
        track_time_delta = int(
            (track.index[1] - track.index[0]) / np.timedelta64(1, "h")
        )
        
        # Log timestep information
        app_logger.debug(f"📊 Data time resolution: {data_time_delta} hour(s)")
        app_logger.debug(f"📊 Track time resolution: {track_time_delta} hour(s)")
        app_logger.debug(f"📅 Data time range: {data[TimeIndexer][0].values} to {data[TimeIndexer][-1].values}")
        app_logger.debug(f"📍 Track time range: {track.index[0]} to {track.index[-1]}")
        
        # If data dt is higher than track dt, raise error (can't select timesteps that don't exist in data)
        if data_time_delta > track_time_delta:
            app_logger.error("❌ Data time step is higher than track time step!")
            app_logger.error("\n" + "="*70)
            app_logger.error("⏱️  TIME RESOLUTION MISMATCH")
            app_logger.error("="*70)
            app_logger.error(f"Data time resolution:  {data_time_delta} hour(s)")
            app_logger.error(f"Track time resolution: {track_time_delta} hour(s)")
            app_logger.error("")
            app_logger.error("📅 Data timestamps:")
            app_logger.error(f"   First: {data[TimeIndexer][0].values}")
            app_logger.error(f"   Last:  {data[TimeIndexer][-1].values}")
            app_logger.error("")
            app_logger.error("📍 Track timestamps:")
            app_logger.error(f"   First: {track.index[0]}")
            app_logger.error(f"   Last:  {track.index[-1]}")
            app_logger.error("")
            app_logger.error("⚠️  Problem: The track requests timesteps that don't exist in the data.")
            app_logger.error(f"   Track needs data every {track_time_delta} hour(s), but data only has timesteps every {data_time_delta} hour(s).")
            app_logger.error("")
            app_logger.error("💡 Solutions:")
            app_logger.error(f"   1. Resample the track file to match or exceed the data time step ({data_time_delta}h)")
            app_logger.error(f"   2. Re-download data with higher temporal resolution (≤{track_time_delta}h)")
            app_logger.error("="*70 + "\n")
            raise ValueError(
                f"Data time step ({data_time_delta}h) is higher than track time step ({track_time_delta}h). "
                f"Cannot select track timesteps that don't exist in data. "
                f"Please resample the track or re-download data with higher temporal resolution."
            )
        # Check data and track initial and final timestamps
        if track.index[0] < data[TimeIndexer][0].values:
            app_logger.error("❌ Track initial timestamp is earlier than data initial timestamp!")
            app_logger.error("\n" + "="*70)
            app_logger.error("📅 TIMESTAMP MISMATCH - Track starts too early")
            app_logger.error("="*70)
            app_logger.error("Data timestamps:")
            app_logger.error(f"   First: {data[TimeIndexer][0].values}")
            app_logger.error(f"   Last:  {data[TimeIndexer][-1].values}")
            app_logger.error("")
            app_logger.error("Track timestamps:")
            app_logger.error(f"   First: {track.index[0]} ← TOO EARLY!")
            app_logger.error(f"   Last:  {track.index[-1]}")
            app_logger.error("")
            app_logger.error("💡 Solution: Adjust the track file to start at or after the data start time.")
            app_logger.error("="*70 + "\n")
            raise ValueError(
                f"Track initial timestamp ({track.index[0]}) is earlier than data initial timestamp ({data[TimeIndexer][0].values}). "
                f"Please adjust the track file."
            ) 
        if track.index[-1] > data[TimeIndexer][-1].values:
            app_logger.error("❌ Track final timestamp is later than data final timestamp!")
            app_logger.error("\n" + "="*70)
            app_logger.error("📅 TIMESTAMP MISMATCH - Track extends beyond data")
            app_logger.error("="*70)
            app_logger.error("Data timestamps:")
            app_logger.error(f"   First: {data[TimeIndexer][0].values}")
            app_logger.error(f"   Last:  {data[TimeIndexer][-1].values} ← DATA ENDS HERE")
            app_logger.error("")
            app_logger.error("Track timestamps:")
            app_logger.error(f"   First: {track.index[0]}")
            app_logger.error(f"   Last:  {track.index[-1]} ← TOO LATE!")
            app_logger.error("")
            app_logger.error("💡 Solution: Either adjust the track file to end at or before the data end time,")
            app_logger.error("   or re-download the data to cover the full track period.")
            app_logger.error("="*70 + "\n")
            raise ValueError(
                f"Track final timestamp ({track.index[-1]}) is later than data final timestamp ({data[TimeIndexer][-1].values}). "
                f"Please adjust the track file or re-download the data."
            )
        # If using CDS API, resample track to data time step
        if args.cdsapi:
            time_delta = int(
                (data[TimeIndexer][1] - data[TimeIndexer][0]) / np.timedelta64(1, "h")
            )
            track = track[track.index.hour % time_delta == 0]
        data = data.sel({TimeIndexer: track.index.values})
        app_logger.debug("✅ Track dates selection complete.")
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

    app_logger.debug("🌐 Assigning geospatial coordinates in radians... ")
    data = data.assign_coords({"rlats": np.deg2rad(data[LatIndexer])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data[LatIndexer]))})
    data = data.assign_coords({"rlons": np.deg2rad(data[LonIndexer])})
    app_logger.debug("✅ Geospatial coordinates assigned.")

    # Drop unnecessary dimensions
    if 'expver' in data.coords:
        data = data.drop('expver')
    if 'number' in data.coords:
        data = data.drop('number')

    # Convert pressure levels to Pa
    app_logger.debug(f"🔄 Converting vertical levels to Pascals...")
    try:
        if not hasattr(data[LevelIndexer], 'units'):
            app_logger.warning(
                f"⚠️  Vertical level coordinate '{LevelIndexer}' has no units attribute. "
                f"Assuming hPa (hectopascals)."
            )
            # Assume hPa if no units
            levels_Pa = data[LevelIndexer] * units('hPa')
            levels_Pa = levels_Pa.metpy.convert_units("Pa")
        else:
            levels_Pa = (
                data[LevelIndexer] * units(str(data[LevelIndexer].units))
            ).metpy.convert_units("Pa")
        data = data.assign_coords({LevelIndexer: levels_Pa})
        app_logger.debug(f"✅ Vertical levels converted to Pa (range: {float(levels_Pa.min()):.0f} - {float(levels_Pa.max()):.0f} Pa)")
    except ValueError as e:
        app_logger.error("❌ Error converting vertical level units!")
        app_logger.error("\n" + "="*70)
        app_logger.error("🔄 UNIT CONVERSION ERROR")
        app_logger.error("="*70)
        app_logger.error(f"Vertical level coordinate: {LevelIndexer}")
        if hasattr(data[LevelIndexer], 'units'):
            app_logger.error(f"Current units: {data[LevelIndexer].units}")
        else:
            app_logger.error("Current units: NOT SPECIFIED")
        app_logger.error(f"Target units: Pa (Pascals)")
        app_logger.error(f"Error: {e}")
        app_logger.error("\n💡 User Solutions:")
        app_logger.error("   1. Check if vertical level units in your data are pressure units")
        app_logger.error("   2. Expected units: Pa, hPa, mb, mbar")
        app_logger.error("   3. If units are missing, add them to the file:")
        app_logger.error("      ncatted -a units,<level_var>,o,c,'hPa' input.nc")
        app_logger.error("\n🔧 Developer Info:")
        app_logger.error(f"   Level variable: {LevelIndexer}")
        app_logger.error(f"   Level values: {data[LevelIndexer].values}")
        app_logger.error(f"   Has units attr: {hasattr(data[LevelIndexer], 'units')}")
        if hasattr(data[LevelIndexer], 'units'):
            app_logger.error(f"   Units value: {data[LevelIndexer].units}")
        app_logger.error(f"   MetPy error: {str(e)}")
        app_logger.error("="*70 + "\n")
        raise ValueError(
            f"Cannot convert vertical level units to Pa. "
            f"Check if '{LevelIndexer}' has valid pressure units."
        )
    except Exception as e:
        app_logger.error("❌ Unexpected error in unit conversion!")
        app_logger.error("\n" + "="*70)
        app_logger.error("⚠️  UNIT CONVERSION FAILURE")
        app_logger.error("="*70)
        app_logger.error(f"Error type: {type(e).__name__}")
        app_logger.error(f"Error message: {e}")
        app_logger.error("\n🔧 Developer Info:")
        app_logger.error(f"   Level variable: {LevelIndexer}")
        app_logger.error(f"   Exception: {type(e).__name__}: {str(e)}")
        app_logger.error("="*70 + "\n")
        raise

    data = (
        data.sortby(LonIndexer)
        .sortby(LevelIndexer, ascending=True)
        .sortby(LatIndexer, ascending=True)
    )

    lowest_level = float(data[LevelIndexer].max())
    data = data.sel({LevelIndexer: slice(1000, lowest_level)})

    if args.mpas:
        data = data.drop_dims("standard_height")

    app_logger.debug("✅ Data processing complete.")
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
        xr.Dataset: The prepared dataset for analysis.
    """
    # Automatically use ERA5-cdsapi namelist when --cdsapi flag is set
    if args.cdsapi:
        varlist = "inputs/namelist_ERA5-cdsapi"
        app_logger.info("🌐 CDS API mode detected: automatically using ERA5-compatible namelist")
        app_logger.debug(f"📋 Using namelist: {varlist}")

    # Validate and load namelist using validation module
    variable_list_df = validate_namelist_file(varlist, app_logger)

    # Load data
    data = get_data(args, app_logger)

    # Import validation functions for coordinate and variable checks
    from .validation import validate_variable_match, validate_required_coordinates
    
    # Validate variable match
    validate_variable_match(variable_list_df, data, varlist, app_logger)
    
    # Validate required coordinates
    validate_required_coordinates(variable_list_df, data, varlist, app_logger)

    # Process and slice data
    processed_data = process_data(data, args, variable_list_df, app_logger)
    sliced_data = slice_domain(processed_data, args, variable_list_df)
    
    return sliced_data

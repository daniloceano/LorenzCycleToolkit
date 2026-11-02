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

from .ncar_gdex import get_ncar_data
from .run_metadata import record_top_cut
from .run_timing import RunTiming
from .select_area import slice_domain
from .tools import convert_longitude_range, get_cdsapi_data
from .track_alignment import (TrackAlignmentError, align_track_times, planned_request_times,
                              record_alignment, warn_adapted)
from .validation import validate_track_file, validate_namelist_file


def check_planned_times(track, args, source, app_logger):
    """
    Before a ``--cdsapi``/``--ncar`` download, check the track against the times the request will retrieve.

    The requested hours are known before the download (``--time-resolution``), so a track that they
    cannot serve is refused here instead of after the transfer. This does not replace the check
    against the times of the file received (:func:`process_data`).

    Args:
        track (pd.DataFrame): Track with a datetime index.
        args (argparse.Namespace): Command-line arguments (``time_resolution``).
        source (str): Name of the data source in the messages.
        app_logger (logging.Logger): Logger for the application.

    Raises:
        TrackAlignmentError: If the requested times cannot serve the track.
    """
    planned = planned_request_times(track.index, args.time_resolution)
    label = f"the {source} request (--time-resolution {args.time_resolution})"
    try:
        alignment = align_track_times(track.index, planned, label)
        missing = alignment.selected[~alignment.selected.isin(planned)]
        if len(missing):
            raise TrackAlignmentError(
                f"{len(missing)} track time(s) are not among the times of {label}, first {missing[0]}: "
                f"the request would not retrieve them. Use a --time-resolution whose times fall on the track times."
            )
    except TrackAlignmentError as error:
        app_logger.error(f"❌ {source} request not made: {error}")
        raise
    if alignment.adapted:
        app_logger.info(
            f"⏱️ {source} request every {alignment.data_step_h:g} h for a track every {alignment.track_step_h:g} h: "
            f"after the download the track will be aligned with the data times ({len(alignment.selected)} of "
            f"{alignment.n_original} positions, checked against the file received)."
        )


def get_data(args: argparse.Namespace, app_logger: logging.Logger, varlist: str = None,
             timing=None) -> xr.Dataset:
    """
    Opens a NetCDF file and extracts variables specified in a CSV file.

    Args:
        args (argparse.Namespace): Command-line arguments.
        app_logger (logging.Logger): Logger for the application.
        varlist (str): Namelist of the run; with --ncar the downloaded file
            uses its variable names (default: inputs/namelist).
        timing (RunTiming): If given, receives the acquisition time (the
            call to get_cdsapi_data or get_ncar_data) and the origin of the
            input file. Timing only: the data are not changed.

    Returns:
        xr.Dataset: Dataset containing extracted variables.

    Raises:
        FileNotFoundError: If CSV file or NetCDF file is not found.
        Exception: For other errors occurring during file opening.
    """
    infile = args.infile
    timing = RunTiming() if timing is None else timing
    timing.set_source("local", "local")

    if args.cdsapi:
        timing.set_source("cds", "reused")
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
            check_planned_times(track, args, "CDS API", app_logger)
            timing.set_source("cds", None)  # "downloaded" only once the call succeeds
            with timing.phase("acquisition"):
                infile = get_cdsapi_data(args, track, app_logger)
            timing.set_source("cds", "downloaded")
            app_logger.debug(f"✅ CDS API data ready: {infile}")
        else:
            if getattr(args, "cds_max_seconds", None) is not None:
                from .cds_guard import CdsGuardError
                raise CdsGuardError(
                    f"guarded CDS mode makes a new acquisition, but {infile} already exists; it is not reused. "
                    f"Use a new destination."
                )
            app_logger.info("✅ CDS API data already exists, skipping download.")
    elif getattr(args, "ncar", False):
        # Same track reading as the CDS API route; an existing file is reused
        # only if its provenance record verifies it
        delimiter, _ = validate_track_file(args.trackfile, app_logger)
        track = pd.read_csv(
            args.trackfile,
            delimiter=delimiter,
            index_col="time",
            parse_dates=["time"],
            date_format="%Y-%m-%d-%H%M"
        )
        if not os.path.exists(infile):
            check_planned_times(track, args, "NCAR GDEX", app_logger)
        timing.set_source("ncar", None)  # set once the call has returned
        with timing.phase("acquisition"):
            infile = get_ncar_data(args, track, app_logger, namelist=varlist or "inputs/namelist",
                                   on_retrieved=timing.record_ncar_download)
        if timing.input_data is None:
            timing.set_source("ncar", "reused")  # returned without retrieving: an existing file was verified

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
        data_times = pd.DatetimeIndex(data[TimeIndexer].values)

        # Log timestep information
        app_logger.debug(f"📊 Data time steps: {sorted(set(np.diff(data_times.values) / np.timedelta64(1, 'h')))} hour(s)")
        app_logger.debug(f"📊 Track time steps: {sorted(set(np.diff(track.index.values) / np.timedelta64(1, 'h')))} hour(s)")
        app_logger.debug(f"📅 Data time range: {data[TimeIndexer][0].values} to {data[TimeIndexer][-1].values}")
        app_logger.debug(f"📍 Track time range: {track.index[0]} to {track.index[-1]}")

        # Align the track with the data times: when the data are sampled less
        # often than the track (e.g. 6-hourly data, 3-hourly track), keep the
        # original track positions at the data times (track_alignment rules)
        try:
            alignment = align_track_times(track.index, data_times, "the meteorological data")
        except TrackAlignmentError as error:
            app_logger.error("❌ Track times cannot be aligned with the data times!")
            app_logger.error("\n" + "="*70)
            app_logger.error("⏱️  TIME RESOLUTION MISMATCH")
            app_logger.error("="*70)
            app_logger.error(str(error))
            app_logger.error(f"   Data:  {data_times[0]} to {data_times[-1]} ({len(data_times)} times)")
            app_logger.error(f"   Track: {track.index[0]} to {track.index[-1]} ({len(track)} positions)")
            app_logger.error("="*70 + "\n")
            raise
        if alignment.adapted:
            track = track.loc[alignment.selected]
            warn_adapted(alignment, app_logger)
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
        missing = track.index[~track.index.isin(data_times)]
        if len(missing):
            app_logger.error(f"❌ {len(missing)} track time(s) are missing from the data: "
                             f"{', '.join(str(t) for t in missing[:6])}{' ...' if len(missing) > 6 else ''}")
            raise ValueError(
                f"{len(missing)} track time(s) are not in the data (first: {missing[0]}). Missing data times "
                f"are not filled: adjust the track file or provide data with these times."
            )
        data = data.sel({TimeIndexer: track.index.values})
        record_alignment(data, alignment, track_file)
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

    input_grid = [float(p) for p in data[LevelIndexer].values]
    lowest_level = float(data[LevelIndexer].max())
    data = data.sel({LevelIndexer: slice(1000, lowest_level)})
    record_top_cut(data, input_grid, LevelIndexer, app_logger)

    if args.mpas:
        data = data.drop_dims("standard_height")

    app_logger.debug("✅ Data processing complete.")
    return data


def prepare_data(
    args, varlist: str = "inputs/namelist", app_logger: logging.Logger = None, timing=None
) -> xr.Dataset:
    """
    Prepare the data for further analysis.

    Parameters:
        args (object): The arguments for the function.
        varlist (str): The file path to the variable list file (namelist).
        app_logger (logging.Logger): The logger for the application.
        timing (RunTiming): If given, receives the acquisition time and the
            origin of the input file (see :func:`get_data`).

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
    data = get_data(args, app_logger, varlist, timing=timing)

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

# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    tools.py                                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:33:03 by daniloceano       #+#    #+#              #
#    Updated: 2026/02/26 09:29:00 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import argparse
import logging
import os
import re
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


def validate_track_file(
    track_file: str, app_logger: logging.Logger
) -> tuple[str, bool]:
    """
    Validate the track file format and detect the delimiter.
    
    Expected format:
    - Delimiter: ';' (semicolon) - standard, or ',' (comma) - alternative
    - Required columns: 'time', 'Lat', 'Lon'
    - Date format: YYYY-MM-DD-HHMM (e.g., 2005-08-08-0000)
    - First line must be header with column names
    
    Args:
        track_file (str): Path to the track file.
        app_logger (logging.Logger): Logger for the application.
    
    Returns:
        tuple: (delimiter, has_warnings) where delimiter is the detected separator
               and has_warnings indicates if any format issues were found.
    
    Raises:
        ValueError: If the track file format is invalid or missing required columns.
    """
    app_logger.debug(f"🔍 Validating track file format: {track_file}")
    
    if not os.path.exists(track_file):
        app_logger.error(f"❌ Track file not found: {track_file}")
        raise FileNotFoundError(f"Track file not found: {track_file}")
    
    # Read first few lines to detect format
    with open(track_file, 'r') as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()
    
    # Detect delimiter by checking header
    delimiter = None
    has_warnings = False
    
    if ';' in first_line:
        delimiter = ';'
    elif ',' in first_line:
        delimiter = ','
        app_logger.warning(
            "⚠️  Track file uses ',' as delimiter instead of the standard ';'"
        )
        app_logger.warning(
            "    The file will be read correctly, but consider converting to ';' separator."
        )
        has_warnings = True
    else:
        delimiter = ';'  # Default fallback
        app_logger.error(
            f"❌ Could not detect delimiter in track file header: {first_line}"
        )
        raise ValueError(
            f"Invalid track file format. Header should contain ';' or ',' separators.\n"
            f"Found: {first_line}"
        )
    
    # Parse header to check columns
    header_columns = [col.strip() for col in first_line.split(delimiter)]
    
    app_logger.debug(f"📋 Detected columns: {header_columns}")
    app_logger.debug(f"🔧 Detected delimiter: '{delimiter}'")
    
    # Check required columns
    required_columns = ['time', 'Lat', 'Lon']
    missing_columns = [col for col in required_columns if col not in header_columns]
    
    if missing_columns:
        app_logger.error("❌ Track file is missing required columns!")
        app_logger.error(f"   Required columns: {required_columns}")
        app_logger.error(f"   Found columns: {header_columns}")
        app_logger.error(f"   Missing: {missing_columns}")
        app_logger.error("\n" + "="*70)
        app_logger.error("📝 EXPECTED TRACK FILE FORMAT:")
        app_logger.error("="*70)
        app_logger.error("time;Lat;Lon")
        app_logger.error("2005-08-08-0000;-22.5;-45")
        app_logger.error("2005-08-08-0600;-22.5;-45")
        app_logger.error("2005-08-08-1200;-22.5;-45")
        app_logger.error("...")
        app_logger.error("="*70)
        app_logger.error("Required:")
        app_logger.error("  • Delimiter: ';' (semicolon)")
        app_logger.error("  • Columns: time, Lat, Lon (case-sensitive)")
        app_logger.error("  • Date format: YYYY-MM-DD-HHMM")
        app_logger.error("  • Optional: additional columns (e.g., SLP_hPa)")
        app_logger.error("="*70 + "\n")
        raise ValueError(
            f"Track file missing required columns: {missing_columns}\n"
            f"Expected: {required_columns}\n"
            f"Found: {header_columns}"
        )
    
    # Validate date format in second line
    if second_line:
        date_str = second_line.split(delimiter)[0].strip()
        
        # Check if date matches expected format YYYY-MM-DD-HHMM
        date_pattern = r'^\d{4}-\d{2}-\d{2}-\d{4}$'
        
        if not re.match(date_pattern, date_str):
            app_logger.error("❌ Track file has invalid date format!")
            app_logger.error(f"   Found: '{date_str}'")
            app_logger.error(f"   Expected format: YYYY-MM-DD-HHMM (e.g., 2005-08-08-0000)")
            app_logger.error("\n" + "="*70)
            app_logger.error("📅 DATE FORMAT EXAMPLES:")
            app_logger.error("="*70)
            app_logger.error("✓ Correct: 2005-08-08-0000 (year-month-day-hourminute)")
            app_logger.error("✓ Correct: 2021-06-26-1800")
            app_logger.error("✗ Wrong: 2005-08-08 00:00 (space and colon)")
            app_logger.error("✗ Wrong: 2005/08/08-0000 (forward slashes)")
            app_logger.error("✗ Wrong: 08-08-2005-0000 (day-month-year)")
            app_logger.error("="*70 + "\n")
            raise ValueError(
                f"Invalid date format in track file: '{date_str}'\n"
                f"Expected: YYYY-MM-DD-HHMM (e.g., 2005-08-08-0000)"
            )
        
        app_logger.debug(f"✓ Date format validated: {date_str}")
    
    # Additional validation: try to read the file with detected delimiter
    try:
        test_df = pd.read_csv(track_file, delimiter=delimiter, nrows=2)
        app_logger.debug(f"✓ File successfully parsed with delimiter '{delimiter}'")
    except Exception as e:
        app_logger.error(f"❌ Error parsing track file: {e}")
        raise ValueError(f"Track file could not be parsed: {e}")
    
    if has_warnings:
        app_logger.info("⚠️  Track file format has minor issues but will be processed.")
    else:
        app_logger.debug("✅ Track file format validation passed.")
    
    return delimiter, has_warnings


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
    app_logger.debug(f"⏱️ Using time resolution from args: {time_step} hour(s)")

    # Log track file bounds and requested data bounds
    app_logger.debug(
        f"📍 Track File Limits: lon range: [{min_lon:.2f}, {max_lon:.2f}], "
        f"lat range: [{min_lat:.2f}, {max_lat:.2f}]"
    )
    app_logger.debug(
        f"🗺️ Buffered Data Bounds: lon range: [{buffered_min_lon}, {buffered_max_lon}], "
        f"lat range: [{buffered_min_lat}, {buffered_max_lat}]"
    )
    app_logger.debug(
        f"📅 Track period: {first_timestamp} to {last_timestamp}"
    )
    app_logger.info(
        f"📥 Will download {len(dates)} day(s) of data: {dates[0].date()} to {dates[-1].date()}"
    )

    # Create temporary directory for daily files
    temp_dir = tempfile.mkdtemp(prefix="cdsapi_daily_")
    app_logger.debug(f"📁 Created temporary directory: {temp_dir}")
    
    # Download data day by day
    app_logger.info(f"🌐 Starting ERA5 data download from CDS API...")
    app_logger.info(f"📦 Expected files to download: {len(dates)}")
    for i, date in enumerate(dates, 1):
        app_logger.info(f"  📅 Day {i}: {date.date()}")
    
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
                    # Round up to include the last timestamp - go one step beyond
                    end_hour = ((end_hour // time_step) + 1) * time_step
                    if end_hour <= start_hour:
                        end_hour = start_hour + time_step
                    # Cap at 24 hours (end of day)
                    if end_hour > 24:
                        end_hour = 24
                else:
                    # First day but not last: download until end of day
                    end_hour = 24
            elif date.date() == last_timestamp.date():
                # Last day (not first): start from beginning of day
                start_hour = 0
                end_hour = last_timestamp.hour
                # Round up to nearest time_step to INCLUDE the last timestamp
                # We need to go one step further to include the last hour in the range
                end_hour = ((end_hour // time_step) + 1) * time_step
                # If last timestamp is at midnight (hour 0), download that hour
                if end_hour == 0:
                    end_hour = time_step
                # Cap at 24 hours (end of day)
                if end_hour > 24:
                    end_hour = 24
            else:
                # Middle day: download all hours
                start_hour = 0
                end_hour = 24 
            
            # Generate times list for this specific day
            times = [f"{hour:02d}:00" for hour in range(start_hour, end_hour, time_step)]
            
            if not times:
                app_logger.warning(f"⚠️ No times to download for {year}-{month}-{day}, skipping...")
                continue
            
            # Create temporary file for this day
            temp_file = os.path.join(temp_dir, f"era5_{year}{month}{day}.nc")
            daily_files.append(temp_file)
            
            app_logger.info(f"⬇️ Downloading day {i}/{len(dates)}: {year}-{month}-{day}")
            app_logger.debug(f"  ⏰ Times: {times[0]} to {times[-1]} ({len(times)} timesteps)")
            
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
        app_logger.info("✅ All daily files downloaded successfully!")
        app_logger.info(f"📦 Total files: {len(daily_files)}")
        app_logger.info("=" * 60)
        
        # Concatenate daily files
        app_logger.info("🔗 Concatenating daily files into single dataset...")
        datasets = []
        for idx, daily_file in enumerate(daily_files, 1):
            app_logger.debug(f"  📂 Loading file {idx}/{len(daily_files)}: {os.path.basename(daily_file)}")
            ds = xr.open_dataset(daily_file)
            datasets.append(ds)
        
        # Concatenate along time dimension
        app_logger.info("🔀 Merging datasets along time dimension...")
        combined_ds = xr.concat(datasets, dim='valid_time')
        
        # Save concatenated dataset
        app_logger.info(f"💾 Saving final dataset to: {args.infile}")
        combined_ds.to_netcdf(args.infile)
        combined_ds.close()
        
        # Close all datasets
        for ds in datasets:
            ds.close()
        
        # Verify the final file was created
        if not os.path.exists(args.infile):
            raise FileNotFoundError(
                f"❌ CDS API file not created at expected path: {args.infile}"
            )
        
        final_file_size = os.path.getsize(args.infile)
        app_logger.info("=" * 60)
        app_logger.info("🎉 DATA DOWNLOAD COMPLETED SUCCESSFULLY!")
        app_logger.info(f"📄 Final file: {args.infile}")
        app_logger.info(f"📊 File size: {final_file_size / (1024**2):.2f} MB")
        app_logger.info(f"📅 Time coverage: {first_timestamp} to {last_timestamp}")
        app_logger.info("=" * 60)
        
    finally:
        # Clean up temporary files and directory
        app_logger.debug("🧹 Cleaning up temporary files...")
        for temp_file in daily_files:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                    app_logger.debug(f"🗑️ Deleted temporary file: {temp_file}")
            except Exception as e:
                app_logger.warning(f"⚠️ Could not delete temporary file {temp_file}: {e}")
        
        try:
            os.rmdir(temp_dir)
            app_logger.debug(f"🗑️ Deleted temporary directory: {temp_dir}")
        except Exception as e:
            app_logger.warning(f"⚠️ Could not delete temporary directory {temp_dir}: {e}")
    
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
        app_logger.error(
            "❌ Could not open file. Check if path, namelist file, and file format (.nc) are correct."
        )
        raise
    except Exception as e:
        app_logger.exception("❌ An exception occurred: {}".format(e))
        raise
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
        method (str): The method used for the analysis: fixed, track or choose.
        xr.Dataset: The prepared dataset for analysis.
    """
    # Automatically use ERA5-cdsapi namelist when --cdsapi flag is set
    if args.cdsapi:
        varlist = "inputs/namelist_ERA5-cdsapi"
        app_logger.info("🌐 CDS API mode detected: automatically using ERA5-compatible namelist")
        app_logger.debug(f"📋 Using namelist: {varlist}")

    app_logger.debug(f"📝 Variables specified by the user in: {varlist}")
    app_logger.debug(f"📂 Attempting to read {varlist} file...")
    try:
        variable_list_df = pd.read_csv(varlist, sep=";", index_col=0, header=0)
    except FileNotFoundError:
        app_logger.error("❌ The 'namelist' text file could not be found.")
        raise
    except pd.errors.EmptyDataError:
        app_logger.error("❌ The 'namelist' text file is empty.")
        raise
    app_logger.debug("✅ Variable list loaded:\n" + str(variable_list_df))

    data = get_data(args, app_logger)

    # Check if variable_list_df matches the data
    if not set(variable_list_df["Variable"]).issubset(set(data.variables)):
        app_logger.error(
            "❌ The variable list does not match the data. Check if the 'namelist' text file is correct."
        )
        
        # Display dataset information to help user configure the namelist
        app_logger.error("\n" + "="*70)
        app_logger.error("📊 DATASET INFORMATION")
        app_logger.error("="*70)
        
        # List coordinates
        app_logger.error("\n🗺️  Available Coordinates:")
        for coord_name in data.coords:
            coord = data.coords[coord_name]
            units_str = f" ({coord.units})" if hasattr(coord, 'units') else ""
            long_name = f" - {coord.long_name}" if hasattr(coord, 'long_name') else ""
            app_logger.error(f"   • {coord_name}{units_str}{long_name}")
        
        # List variables
        app_logger.error("\n📋 Available Variables:")
        for var_name in data.data_vars:
            var = data[var_name]
            units_str = f" ({var.units})" if hasattr(var, 'units') else ""
            long_name = f" - {var.long_name}" if hasattr(var, 'long_name') else ""
            standard_name = f" [{var.standard_name}]" if hasattr(var, 'standard_name') else ""
            app_logger.error(f"   • {var_name}{units_str}{long_name}{standard_name}")
        
        app_logger.error("\n" + "="*70)
        app_logger.error("💡 TIP: Update your namelist file with the correct variable names")
        app_logger.error("    from the list above, or use one of the preset namelists:")
        app_logger.error("    - inputs/namelist_ERA5-cdsapi (for ERA5 data)")
        app_logger.error("    - inputs/namelist_NCEP-R1 (for NCEP Reanalysis 1)")
        app_logger.error("    - inputs/namelist_NCEP-R2 (for NCEP Reanalysis 2)")
        app_logger.error("="*70 + "\n")
        
        raise ValueError("'namelist' text file does not match the data.")

    processed_data = process_data(data, args, variable_list_df, app_logger)
    sliced_data = slice_domain(processed_data, args, variable_list_df)
    return sliced_data

# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    validation.py                                      :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2026/02/26 10:40:00 by daniloceano       #+#    #+#              #
#    Updated: 2026/02/26 10:40:00 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
Validation functions for the Lor enz Energy Cycle Toolkit.

This module contains functions to validate input files, namelists, coordinates,
and other configuration parameters before processing.
"""

import logging
import os
import re

import pandas as pd
import xarray as xr


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
        FileNotFoundError: If the track file does not exist.
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


def validate_namelist_file(
    varlist: str, app_logger: logging.Logger
) -> pd.DataFrame:
    """
    Validate and load the namelist file.
    
    Args:
        varlist (str): Path to the namelist file.
        app_logger (logging.Logger): Logger for the application.
    
    Returns:
        pd.DataFrame: The loaded variable list DataFrame.
    
    Raises:
        FileNotFoundError: If the namelist file is not found.
        pd.errors.EmptyDataError: If the namelist file is empty.
        Exception: For other parsing errors.
    """
    app_logger.debug(f"📝 Variables specified by the user in: {varlist}")
    app_logger.debug(f"📂 Attempting to read {varlist} file...")
    
    try:
        variable_list_df = pd.read_csv(varlist, sep=";", index_col=0, header=0)
    except FileNotFoundError:
        app_logger.error("❌ The 'namelist' file could not be found!")
        app_logger.error("\n" + "="*70)
        app_logger.error("📄 NAMELIST FILE NOT FOUND")
        app_logger.error("="*70)
        app_logger.error(f"Looking for: {os.path.abspath(varlist)}")
        app_logger.error(f"Current directory: {os.getcwd()}")
        app_logger.error("\n💡 User Solutions:")
        app_logger.error("   1. Create a namelist file from one of the presets:")
        app_logger.error("      cp inputs/namelist_ERA5-cdsapi inputs/namelist")
        app_logger.error("   2. Available preset namelists:")
        app_logger.error("      - inputs/namelist_ERA5-cdsapi (for ERA5 data)")
        app_logger.error("      - inputs/namelist_NCEP-R1 (for NCEP Reanalysis 1)")
        app_logger.error("      - inputs/namelist_NCEP-R2 (for NCEP Reanalysis 2)")
        app_logger.error("   3. Customize the namelist to match your dataset")
        app_logger.error("\n🔗 Documentation:")
        app_logger.error("   See Configuration guide: docs/source/configuration.rst")
        app_logger.error("="*70 + "\n")
        raise FileNotFoundError(
            f"Namelist file not found: {os.path.abspath(varlist)}. "
            f"Please create one from the preset namelists in inputs/ directory."
        )
    except pd.errors.EmptyDataError:
        app_logger.error("❌ The 'namelist' file is empty!")
        app_logger.error("\n" + "="*70)
        app_logger.error("📄 EMPTY NAMELIST FILE")
        app_logger.error("="*70)
        app_logger.error(f"File: {os.path.abspath(varlist)}")
        app_logger.error("\n💡 Solution:")
        app_logger.error("   The namelist file exists but contains no data.")
        app_logger.error("   Copy from one of the preset namelists:")
        app_logger.error("   cp inputs/namelist_ERA5-cdsapi inputs/namelist")
        app_logger.error("="*70 + "\n")
        raise pd.errors.EmptyDataError(
            f"Namelist file is empty: {os.path.abspath(varlist)}"
        )
    except Exception as e:
        app_logger.error(f"❌ Error reading namelist file!")
        app_logger.error("\n" + "="*70)
        app_logger.error("📄 NAMELIST PARSING ERROR")
        app_logger.error("="*70)
        app_logger.error(f"File: {os.path.abspath(varlist)}")
        app_logger.error(f"Error: {type(e).__name__}: {e}")
        app_logger.error("\n💡 Possible causes:")
        app_logger.error("   1. File format is incorrect (should be CSV with ';' delimiter)")
        app_logger.error("   2. File encoding issues")
        app_logger.error("   3. Malformed CSV structure")
        app_logger.error("\n🔧 Developer Info:")
        app_logger.error(f"   Exception: {type(e).__name__}")
        app_logger.error(f"   Message: {str(e)}")
        app_logger.error("="*70 + "\n")
        raise
    
    app_logger.debug("✅ Variable list loaded:\n" + str(variable_list_df))
    return variable_list_df


def validate_variable_match(
    variable_list_df: pd.DataFrame, data: xr.Dataset, varlist: str, app_logger: logging.Logger
) -> None:
    """
    Validate that variables in the namelist exist in the dataset.
    
    Args:
        variable_list_df (pd.DataFrame): The variable list from namelist.
        data (xr.Dataset): The loaded dataset.
        varlist (str): Path to the namelist file (for error messages).
        app_logger (logging.Logger): Logger for the application.
    
    Raises:
        ValueError: If variables don't match between namelist and dataset.
    """
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


def validate_required_coordinates(
    variable_list_df: pd.DataFrame, data: xr.Dataset, varlist: str, app_logger: logging.Logger
) -> None:
    """
    Validate that all required coordinates exist in the dataset.
    
    Args:
        variable_list_df (pd.DataFrame): The variable list from namelist.
        data (xr.Dataset): The loaded dataset.
        varlist (str): Path to the namelist file (for error messages).
        app_logger (logging.Logger): Logger for the application.
    
    Raises:
        ValueError: If required coordinates are missing.
    """
    app_logger.debug("🔍 Validating required coordinates...")
    required_coords = ["Longitude", "Latitude", "Vertical Level", "Time"]
    missing_coords = []
    
    for coord_type in required_coords:
        coord_name = variable_list_df.loc[coord_type]["Variable"]
        if coord_name not in data.coords and coord_name not in data.variables:
            missing_coords.append((coord_type, coord_name))
    
    if missing_coords:
        app_logger.error("❌ Required coordinates are missing in the dataset!")
        app_logger.error("\n" + "="*70)
        app_logger.error("🗺️  MISSING COORDINATES")
        app_logger.error("="*70)
        app_logger.error("Required coordinates from namelist that are missing:")
        for coord_type, coord_name in missing_coords:
            app_logger.error(f"   • {coord_type}: '{coord_name}'")
        
        app_logger.error("\nAvailable coordinates in dataset:")
        for coord in data.coords:
            app_logger.error(f"   • {coord}")
        
        app_logger.error("\n💡 User Solutions:")
        app_logger.error("   1. Update the namelist to use the correct coordinate names")
        app_logger.error("   2. Check if your data has these dimensions:")
        app_logger.error("      • Longitude (e.g., 'lon', 'longitude', 'x')")
        app_logger.error("      • Latitude (e.g., 'lat', 'latitude', 'y')")
        app_logger.error("      • Vertical Level (e.g., 'level', 'pressure_level', 'lev')")
        app_logger.error("      • Time (e.g., 'time', 'valid_time', 't')")
        
        app_logger.error("\n🔧 Developer Info:")
        app_logger.error(f"   Dataset dimensions: {list(data.dims)}")
        app_logger.error(f"   Dataset coordinates: {list(data.coords)}")
        app_logger.error(f"   Namelist file: {varlist}")
        app_logger.error("="*70 + "\n")
        
        raise ValueError(
            f"Missing required coordinates: {[c[0] for c in missing_coords]}. "
            f"Update namelist to match your dataset coordinate names."
        )
    
    app_logger.debug("✅ All required coordinates validated.")

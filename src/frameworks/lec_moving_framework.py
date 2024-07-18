# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lec_moving_framework.py                            :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:32:55 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/18 00:36:10 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import argparse
import logging
import os

import numpy as np
import pandas as pd
import xarray as xr
from metpy.calc import vorticity, wind_speed
from metpy.constants import g
from metpy.units import units

from ..analysis.boundary_terms import BoundaryTerms
from ..analysis.conversion_terms import ConversionTerms
from ..analysis.energy_contents import EnergyContents
from ..analysis.generation_and_dissipation_terms import \
    GenerationDissipationTerms
from ..utils.box_data import BoxData
from ..utils.calc_budget_and_residual import calc_budget_diff, calc_residuals
from ..utils.select_area import draw_box_map, plot_domain_attributes
from ..utils.tools import find_extremum_coordinates, initialize_logging


def create_terms_dict(args):
    """
    Create a dictionary for storing computation results of various meteorological terms.

    Args:
        args (argparse.Namespace): Script arguments.

    Returns:
        dict: Dictionary with keys for each meteorological term initialized to empty lists.
    """
    energy_terms = ["Az", "Ae", "Kz", "Ke"]
    conversion_terms = ["Cz", "Ca", "Ck", "Ce"]
    boundary_terms = ["BAz", "BAe", "BKz", "BKe", "BΦZ", "BΦE"]
    generation_dissipation_terms = (
        ["Gz", "Ge", "Dz", "De"] if not args.residuals else ["Gz", "Ge"]
    )

    terms = (
        energy_terms + conversion_terms + boundary_terms + generation_dissipation_terms
    )
    return {term: [] for term in terms}


def handle_track_file(
    data, times, LonIndexer, LatIndexer, TimeIndexer, args, app_logger
):
    """
    Handles the track file by validating its time and spatial limits against the provided dataset.

    Args:
        data (xr.Dataset): A Xarray Dataset containing the data to compute the energy cycle.
        times (pd.DatetimeIndex): The time series of the dataset.
        LonIndexer (str): The name of the longitude coordinate in the dataset.
        LatIndexer (str): The name of the latitude coordinate in the dataset.
        args (argparse.Namespace): Arguments provided to the script.
        app_logger (logging.Logger): Logger for logging messages.

    Returns:
        pd.DataFrame: DataFrame containing the track information if the track file is valid.

    Raises:
        FileNotFoundError: If the track file is not found.
        ValueError: If the time or spatial limits of the track file do not match the dataset.
    """

    trackfile = args.trackfile
    app_logger.debug(f"Using track file: {trackfile}")

    track = pd.read_csv(trackfile, parse_dates=[0], delimiter=";", index_col="time")

    # Remove timezone information from track index if it's tz-aware
    if track.index.tz is not None:
        track.index = track.index.tz_localize(None)

    if args.cdsapi:
        time_delta = int(
            (data[TimeIndexer][1] - data[TimeIndexer][0]) / np.timedelta64(1, "h")
        )
        track = track[track.index.hour % time_delta == 0]

    try:
        data_lon_max, data_lon_min = float(data[LonIndexer].max()), float(
            data[LonIndexer].min()
        )
        data_lat_max, data_lat_min = float(data[LatIndexer].max()), float(
            data[LatIndexer].min()
        )

        app_logger.debug(
            f"Data spatial limits --> lon_min: {data_lon_min:.2f}, lon_max: {data_lon_max:.2f}, "
            f"lat_min: {data_lat_min:.2f}, lat_max: {data_lat_max:.2f}"
        )
        app_logger.debug(
            f"Track spatial limits --> lon_min: {track['Lon'].min():.2f}, lon_max: {track['Lon'].max():.2f}, "
            f"lat_min: {track['Lat'].min():.2f}, lat_max: {track['Lat'].max():.2f}"
        )

        if track.index[0] < times.min() or track.index[-1] > times.max():
            app_logger.error("Track time limits do not match with data time limits.")
            raise ValueError("Track time limits do not match with data time limits.")

        # Check longitude limits
        if track["Lon"].max() > data_lon_max:
            app_logger.error(
                f"Track file longitude max limit ({track['Lon'].max():.2f}) exceeds data max "
                f"longitude limit ({data_lon_max:.2f})."
            )
            raise ValueError(
                f"Track file longitude max limit ({track['Lon'].max():.2f}) exceeds data max "
                f"longitude limit ({data_lon_max:.2f})."
            )
        if track["Lon"].min() < data_lon_min:
            app_logger.error(
                f"Track file longitude min limit ({track['Lon'].min():.2f}) is below data min "
                f"longitude limit ({data_lon_min:.2f})."
            )
            raise ValueError(
                f"Track file longitude min limit ({track['Lon'].min():.2f}) is below data min "
                f"longitude limit ({data_lon_min:.2f})."
            )

        # Check latitude limits
        if track["Lat"].max() > data_lat_max:
            app_logger.error(
                f"Track file latitude max limit ({track['Lat'].max():.2f}) exceeds data max "
                f"latitude limit ({data_lat_max:.2f})."
            )
            raise ValueError(
                f"Track file latitude max limit ({track['Lat'].max():.2f}) exceeds data max "
                f"latitude limit ({data_lat_max:.2f})."
            )
        if track["Lat"].min() < data_lat_min:
            app_logger.error(
                f"Track file latitude min limit ({track['Lat'].min():.2f}) is below data min "
                f"latitude limit ({data_lat_min:.2f})."
            )
            raise ValueError(
                f"Track file latitude min limit ({track['Lat'].min():.2f}) is below data min "
                f"latitude limit ({data_lat_min:.2f})."
            )

        return track

    except FileNotFoundError:
        app_logger.error(f"Track file {trackfile} not found.")
        raise


def extract_wind_and_height_components(idata, variable_list_df, args):
    """
    Extract wind components and geopotential height from the dataset for a specific time.

    Args:
        idata (xr.Dataset): The dataset at a specific time.
        variable_list_df (pd.DataFrame): DataFrame containing variable mappings.
        args (argparse.Namespace): Script arguments.

    Returns:
        tuple: Tuple containing u-component, v-component, and geopotential height DataArrays.
    """
    u_var = variable_list_df.loc["Eastward Wind Component"]["Variable"]
    v_var = variable_list_df.loc["Northward Wind Component"]["Variable"]
    if "Geopotential Height" in variable_list_df.index:
        geopotential_var = variable_list_df.loc["Geopotential Height"]["Variable"]
    else:
        geopotential_var = variable_list_df.loc["Geopotential"]["Variable"]
    iu = idata[u_var].compute() * units(
        variable_list_df.loc["Eastward Wind Component"]["Units"]
    ).to("m/s")
    iv = idata[v_var].compute() * units(
        variable_list_df.loc["Northward Wind Component"]["Units"]
    ).to("m/s")
    ihgt = idata[geopotential_var].compute() * units(
        variable_list_df[variable_list_df["Variable"] == geopotential_var][
            "Units"
        ].iloc[0]
    )

    if "Geopotential" in variable_list_df.index:
        ihgt = (ihgt / g).metpy.convert_units("gpm")

    return iu, iv, ihgt


def get_limits(args, t, data850, track=None):
    """
    Determine the central position of the system at a specific time and the box dimensions.
    This function calculates the central latitude and longitude, as well as the dimensions
    of the box (width and length) based on the tracking data or user input.

    Args:
        args (argparse.Namespace): Script arguments.
        t (pd.Timestamp): Current time.
        track (pd.DataFrame): DataFrame containing the tracking information.
        data850 (dict): Dictionary containing meteorological data at 850 hPa.

    Returns:
        dict: Dictionary containing the central position (latitude, longitude) and limits (width, length) of the system.
    """
    # Initialize default values
    central_lat, central_lon, width, length = None, None, None, None

    if args.track:
        # Find the closest time index in the track DataFrame
        closest_index = np.argmin(np.abs(track.index - t))
        track_row = track.iloc[closest_index]
        central_lat, central_lon = track_row["Lat"], track_row["Lon"]

        # Use provided width and length if available, otherwise default values
        width = track_row.get("width", 15)  # Default width if not specified
        length = track_row.get("length", 15)  # Default length if not specified

    elif args.choose:
        # For the 'choose' option, interactively determine the box limits
        (
            izeta_850,
            ihgt_850,
        ) = (
            data850["izeta_850"],
            data850["ihgt_850"],
        )
        iu_850, iv_850 = data850["iu_850"], data850["iv_850"]
        lat, lon = data850["lat"], data850["lon"]
        limits = draw_box_map(iu_850, iv_850, izeta_850, ihgt_850, lat, lon, t)
        central_lat, central_lon = (limits["max_lat"] + limits["min_lat"]) / 2, (
            limits["max_lon"] + limits["min_lon"]
        ) / 2
        width, length = (
            limits["max_lon"] - limits["min_lon"],
            limits["max_lat"] - limits["min_lat"],
        )

    # Calculate the bounding box limits
    min_lon, max_lon = central_lon - (width / 2), central_lon + (width / 2)
    min_lat, max_lat = central_lat - (length / 2), central_lat + (length / 2)

    # Format time for output
    datestr = pd.to_datetime(t).strftime("%Y-%m-%d-%H%M")

    limits = {
        "datestr": datestr,
        "central_lat": central_lat,
        "central_lon": central_lon,
        "length": length,
        "width": width,
        "min_lon": min_lon,
        "max_lon": max_lon,
        "min_lat": min_lat,
        "max_lat": max_lat,
    }

    return limits


def get_position(
    track, limits, izeta_850, ihgt_850, iwspd_850, LatIndexer, LonIndexer, args
):
    """
    Retrieves or calculates the values of 'min_max_zeta_850', 'min_hgt_850',
    and 'max_wind_850' based on the given track and data slices.

    Parameters:
        track (DataFrame): The track file containing 'min_max_zeta_850', 'min_hgt_850', and 'max_wind_850' columns.
        track_itime (int): The index of the track file to retrieve the values from.
        central_lat (float): The central latitude.
        central_lon (float): The central longitude.
        min_lat (float): The minimum latitude.
        izeta_850 (DataArray): The data slice for 'izeta_850'.
        izeta_850_slice (DataArray): The data slice for 'izeta_850' with the required latitude and longitude values.
        ight_850_slice (DataArray): The data slice for 'ight_850'.
        iwspd_850_slice (DataArray): The data slice for 'iwspd_850'.

    Returns:
        Tuple[float, float, float]: A tuple containing the values of 'min_max_zeta_850', 'min_hgt_850', and 'max_wind_850'.
    """
    # Track timestep closest to the model timestep, just in case
    # the track file has a poorer temporal resolution
    track_itime = track.loc[limits["datestr"]].name if args.track else None

    limits["central_lat"]
    limits["central_lon"]
    min_lat = limits["min_lat"]
    min_lon = limits["min_lon"]
    max_lat = limits["max_lat"]
    max_lon = limits["max_lon"]

    # Slice data for defined box
    izeta_850_slice = izeta_850.sel(
        {LatIndexer: slice(min_lat, max_lat), LonIndexer: slice(min_lon, max_lon)}
    )
    ihgt_850_slice = ihgt_850.sel(
        {LatIndexer: slice(min_lat, max_lat), LonIndexer: slice(min_lon, max_lon)}
    )
    iwspd_850_slice = iwspd_850.sel(
        {LatIndexer: slice(min_lat, max_lat), LonIndexer: slice(min_lon, max_lon)}
    )

    # Check if 'min_max_zeta_850' is present and valid in the track DataFrame
    if track is not None and "min_max_zeta_850" in track.columns:
        min_max_zeta = float(track.loc[track_itime, "min_max_zeta_850"])
    else:
        # Fallback logic if 'min_max_zeta_850' is not available
        if track is not None and args.zeta and "min_max_zeta_850" not in track.columns:
            # When args.zeta, whe trust the vorticity from the trackfile. If
            # 'min_max_zeta_850' is not available, we need to calculate it
            min_max_zeta_unformatted = izeta_850.sel(
                latitude=limits["central_lat"],
                longitude=limits["central_lon"],
                method="nearest",
            )
        elif track is not None and args.zeta and "min_max_zeta_850" in track.columns:
            # When args.zeta, whe trust the vorticity from the trackfile. If
            # 'min_max_zeta_850' is available, we use it
            min_max_zeta_unformatted = track.loc[track_itime, "min_max_zeta_850"]
        else:
            # When not args.zeta, we need to find the min vorticity in the box
            # and calculate it
            izeta_850_slice = izeta_850.sel(
                {
                    LatIndexer: slice(limits["min_lat"], limits["max_lat"]),
                    LonIndexer: slice(limits["min_lon"], limits["max_lon"]),
                }
            )
            min_max_zeta_unformatted = izeta_850_slice.min()

        if limits["min_lat"] < 0:
            min_max_zeta = float(np.nanmin(min_max_zeta_unformatted))
        else:
            min_max_zeta = float(np.nanmax(min_max_zeta_unformatted))

    # Check if 'min_hgt_850' is present and valid in the track DataFrame
    if (
        track is not None
        and "min_hgt_850" in track.columns
        and not pd.isna(track.loc[track_itime, "min_hgt_850"])
    ):
        min_hgt = float(track.loc[track_itime, "min_hgt_850"])
    else:
        # Fallback logic if 'min_hgt_850' is not available
        ihgt_850_slice = ihgt_850.sel(
            {
                LatIndexer: slice(limits["min_lat"], limits["max_lat"]),
                LonIndexer: slice(limits["min_lon"], limits["max_lon"]),
            }
        )
        min_hgt = float(ihgt_850_slice.min())

    # Check if 'max_wind_850' is present and valid in the track DataFrame
    if (
        track is not None
        and "max_wind_850" in track.columns
        and not pd.isna(track.loc[track_itime, "max_wind_850"])
    ):
        max_wind = float(track.loc[track_itime, "max_wind_850"])
    else:
        # Fallback logic if 'max_wind_850' is not available
        iwspd_850_slice = iwspd_850.sel(
            {
                LatIndexer: slice(limits["min_lat"], limits["max_lat"]),
                LonIndexer: slice(limits["min_lon"], limits["max_lon"]),
            }
        )
        max_wind = float(iwspd_850_slice.max())

    # Find extremum coordinates
    lat_slice, lon_slice = izeta_850_slice[LatIndexer], izeta_850_slice[LonIndexer]
    min_max_zeta_lat, min_max_zeta_lon = find_extremum_coordinates(
        izeta_850_slice, lat_slice, lon_slice, "min_max_zeta"
    )
    min_hgt_lat, min_hgt_lon = find_extremum_coordinates(
        ihgt_850_slice, lat_slice, lon_slice, "min_hgt"
    )
    max_wind_lat, max_wind_lon = find_extremum_coordinates(
        iwspd_850_slice, lat_slice, lon_slice, "max_wind"
    )

    position = {
        "min_max_zeta_850": {
            "latitude": min_max_zeta_lat,
            "longitude": min_max_zeta_lon,
            "value": min_max_zeta,
        },
        "min_hgt_850": {
            "latitude": min_hgt_lat,
            "longitude": min_hgt_lon,
            "value": min_hgt,
        },
        "max_wind_850": {
            "latitude": max_wind_lat,
            "longitude": max_wind_lon,
            "value": max_wind,
        },
    }

    return position


def flatten_position(position):
    """Flatten the position dictionary to include latitude, longitude, and value as separate keys."""
    flattened = {}
    for key, value in position.items():
        flattened[f"{key}_lat"] = value["latitude"]
        flattened[f"{key}_lon"] = value["longitude"]
        flattened[f"{key}"] = value["value"]
    return flattened


def compute_and_store_terms(box_obj, terms_dict, app_logger):
    """
    Compute various meteorological terms using the provided BoxData object
    and store the results in the provided dictionary, with specific error handling.

    Args:
        box_obj (BoxData): BoxData object containing the data for computations.
        terms_dict (dict): Dictionary to store the computed terms.

    Returns:
        dict: Updated dictionary with the computed terms.
    """
    # Energy Contents
    app_logger.info("Computing Energy Contents...")
    try:
        ec_obj = EnergyContents(box_obj, method="moving", app_logger=app_logger)
        terms_dict["Az"].append(ec_obj.calc_az())
        terms_dict["Ae"].append(ec_obj.calc_ae())
        terms_dict["Kz"].append(ec_obj.calc_kz())
        terms_dict["Ke"].append(ec_obj.calc_ke())
    except Exception as e:
        app_logger.exception(f"Error in computing Energy Contents: {e}")
        raise

    # Conversion Terms
    app_logger.info("Computing Conversion Terms...")
    try:
        ct_obj = ConversionTerms(box_obj, method="moving", app_logger=app_logger)
        terms_dict["Cz"].append(ct_obj.calc_cz())
        terms_dict["Ca"].append(ct_obj.calc_ca())
        terms_dict["Ck"].append(ct_obj.calc_ck())
        terms_dict["Ce"].append(ct_obj.calc_ce())
    except Exception as e:
        app_logger.exception(f"Error in computing Conversion Terms: {e}")
        raise

    # Boundary Terms
    app_logger.info("Computing Boundary Terms...")
    try:
        bt_obj = BoundaryTerms(box_obj, method="moving", app_logger=app_logger)
        terms_dict["BAz"].append(bt_obj.calc_baz())
        terms_dict["BAe"].append(bt_obj.calc_bae())
        terms_dict["BKz"].append(bt_obj.calc_bkz())
        terms_dict["BKe"].append(bt_obj.calc_bke())
        terms_dict["BΦZ"].append(bt_obj.calc_boz())
        terms_dict["BΦE"].append(bt_obj.calc_boe())
    except Exception as e:
        app_logger.exception(f"Error in computing Boundary Terms: {e}")
        raise

    # Generation/Dissipation Terms
    app_logger.info("Computing Generation/Dissipation Terms...")
    try:
        gdt_obj = GenerationDissipationTerms(
            box_obj, method="moving", app_logger=app_logger
        )
        terms_dict["Gz"].append(gdt_obj.calc_gz())
        terms_dict["Ge"].append(gdt_obj.calc_ge())
        if not box_obj.args.residuals:
            terms_dict["Dz"].append(gdt_obj.calc_dz())
            terms_dict["De"].append(gdt_obj.calc_de())
    except Exception as e:
        app_logger.exception(f"Error in computing Generation/Dissipation Terms: {e}")
        raise

    return terms_dict


def finalize_results(
    times, terms_dict, args, results_subdirectory, out_track, app_logger
):
    """
    Process and finalize results.
    This includes creating a DataFrame with results and saving it to a CSV file.
    Also saves the system track as a CSV file for replicability.

    Args:
        times (array-like): Array of times for the computation.
        terms_dict (dict): Dictionary containing computed terms.
        args (argparse.Namespace): Script arguments.
        results_subdirectory (str): Directory to save results.
    """
    # Convert the terms_dict to DataFrame
    df = pd.DataFrame(terms_dict, index=pd.to_datetime(times), dtype=float)

    # Estimating budget terms (∂X/∂t) using finite differences
    app_logger.info("Estimating budget terms...")
    df = calc_budget_diff(df, times, app_logger)

    # Computing residuals, if required
    if args.residuals:
        app_logger.info("Computing residuals...")
        df = calc_residuals(df, app_logger)

    # Constructing output filename
    method = "track" if args.track else "choose"
    infile_name = os.path.basename(args.infile).split(".nc")[0]
    results_filename = "".join(f"{infile_name}_{method}_results")
    results_file = os.path.join(results_subdirectory, f"{results_filename}.csv")

    # Saving the DataFrame to a CSV file
    df.to_csv(results_file)
    app_logger.info(f"Results saved to {results_file}")

    # Save system position as a csv file for replicability
    out_track = out_track.rename(
        columns={"datestr": "time", "central_lat": "Lat", "central_lon": "Lon"}
    )
    out_trackfile_name = "".join(f"{infile_name}_{method}_trackfile")
    output_trackfile = os.path.join(results_subdirectory, out_trackfile_name)
    out_track.to_csv(output_trackfile, index=False, sep=";")
    app_logger.info(f"System track saved to {output_trackfile}")

    return results_file, df


def lec_moving(
    data: xr.Dataset,
    variable_list_df: pd.DataFrame,
    dTdt: xr.Dataset,
    results_subdirectory: str,
    figures_directory: str,
    results_subdirectory_vertical_levels: str,
    app_logger: logging.Logger,
    args: argparse.Namespace,
):
    """
    Computes the Lorenz Energy Cycle using a moving (semi-lagrangian) framework.

    Args:
        data (xr.Dataset): A Xarray Dataset containing the data to compute the energy cycle.
        variable_list_df (pd.DataFrame): DataFrame with variable mappings.
        dTdt (xr.Dataset): Dataset containing temperature gradient data.
        results_subdirectory (str): Directory for saving results.
        figures_directory (str): Directory for saving figures.
        args (argparse.Namespace): Arguments provided to the script.

    Returns:
        None
    """
    app_logger.info("Computing energetics using moving framework")

    # Indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = (
        variable_list_df.loc["Longitude"]["Variable"],
        variable_list_df.loc["Latitude"]["Variable"],
        variable_list_df.loc["Time"]["Variable"],
        variable_list_df.loc["Vertical Level"]["Variable"],
    )

    PressureData = data[VerticalCoordIndexer]

    # Create csv files for storing vertical results
    for term in [
        "Az",
        "Ae",
        "Kz",
        "Ke",
        "Ge",
        "Gz",
        "Cz",
        "Cz_1",
        "Cz_2",
        "Ca",
        "Ca_1",
        "Ca_2",
        "Ce",
        "Ce_1",
        "Ce_2",
        "Ck",
        "Ck_1",
        "Ck_2",
        "Ck_3",
        "Ck_4",
        "Ck_5",
    ]:
        columns = [TimeName] + [float(i) for i in PressureData.values]
        df = pd.DataFrame(columns=columns)
        file_name = term + "_" + VerticalCoordIndexer + ".csv"
        file_path = os.path.join(results_subdirectory_vertical_levels, file_name)
        df.to_csv(file_path, index=None)
        app_logger.info(f"{file_path} created (but still empty)")

    # Get the time array
    times = pd.to_datetime(data[TimeName].values)

    # Get the system track and check for errors
    if args.track:
        track = handle_track_file(
            data, times, LonIndexer, LatIndexer, TimeName, args, app_logger
        )

    # Dictionary for saving system position and attributes
    results_keys = [
        "datestr",
        "central_lat",
        "central_lon",
        "length",
        "width",
        "min_max_zeta_850",
        "min_hgt_850",
        "max_wind_850",
    ]
    out_track = pd.DataFrame(columns=results_keys)

    # Dictionary for storing results
    terms_dict = create_terms_dict(args)

    # Iterating over times
    for t in times:
        app_logger.info(f"Processing data at time: {t}...")
        try:
            idata, idTdt = data.sel({TimeName: t}), dTdt.sel({TimeName: t})
            if idata[TimeName].shape != ():
                idata, idTdt = data.isel({TimeName: 1}), dTdt.isel({TimeName: 1})
        except KeyError as e:
            app_logger.error(f"Time indexing error: {e}")
            continue

        # Wind components and geopotential height
        iu, iv, ihgt = extract_wind_and_height_components(idata, variable_list_df, args)

        # Extract 850 hPa data
        iu_850, iv_850, ihgt_850 = (
            iu.sel({VerticalCoordIndexer: 85000}),
            iv.sel({VerticalCoordIndexer: 85000}),
            ihgt.sel({VerticalCoordIndexer: 85000}),
        )

        # Compute wind speed and vorticity at 850 hPa
        iwspd_850, izeta_850 = (
            wind_speed(iu_850, iv_850),
            vorticity(iu_850, iv_850).metpy.dequantify(),
        )
        lat, lon = idata[LatIndexer], idata[LonIndexer]
        # data850 = construct_data850(izeta_850, ihgt_850, iwspd_850, lat, lon)
        data850 = {
            "izeta_850": izeta_850,
            "ihgt_850": ihgt_850,
            "iwspd_850": iwspd_850,
            "iu_850": iu_850,
            "iv_850": iv_850,
            "lat": lat,
            "lon": lon,
        }

        # Get box attributes for current time
        limits = get_limits(args, t, data850, track if args.track else None)
        app_logger.info(
            f"central lat: {limits['central_lat']}, central lon: {limits['central_lon']}, "
            f"size: {limits['length']} x {limits['width']}, "
            f"lon range: {limits['min_lon']} to {limits['max_lon']}, "
            f"lat range: {limits['min_lat']} to {limits['max_lat']}"
        )

        # Get position of  850 hPaextreme values for current time
        position = get_position(
            track if args.track else None,
            limits,
            izeta_850,
            ihgt_850,
            iwspd_850,
            LatIndexer,
            LonIndexer,
            args,
        )
        app_logger.info(
            f"850 hPa diagnostics --> "
            f"min/max ζ: {position['min_max_zeta_850']['value']:.2e}, "
            f"min geopotential height: {position['min_hgt_850']['value']:.0f}, "
            f"max wind speed: {position['max_wind_850']['value']:.2f}"
        )
        app_logger.info(
            f"850 hPa positions (lat/lon) --> "
            f"min/max ζ: {position['min_max_zeta_850']['latitude']:.2f}, {position['min_max_zeta_850']['longitude']:.2f}, "
            f"min geopotential height: {position['min_hgt_850']['latitude']:.2f}, {position['min_hgt_850']['longitude']:.2f}, "
            f"max wind speed: {position['max_wind_850']['latitude']:.2f}, {position['max_wind_850']['longitude']:.2f}"
        )

        # Store results
        flattened_position = flatten_position(position)
        limits_and_position = {**limits, **flattened_position}
        new_entry = pd.DataFrame(
            [limits_and_position], index=[t.strftime("%Y-%m-%d-%H%M")]
        )
        if out_track.empty:
            out_track = new_entry
        else:
            out_track = pd.concat([out_track, new_entry], ignore_index=True)

        # Save figure with the domain box, extreme values, vorticity and
        # geopotential height
        if args.plots:
            plot_domain_attributes(data850, limits, position, figures_directory)

        # Create box object
        app_logger.info("Creating box object...")
        try:
            box_obj = BoxData(
                data=idata.compute(),
                variable_list_df=variable_list_df,
                args=args,
                western_limit=limits["min_lon"],
                eastern_limit=limits["max_lon"],
                southern_limit=limits["min_lat"],
                northern_limit=limits["max_lat"],
                results_subdirectory=results_subdirectory,
                results_subdirectory_vertical_levels=results_subdirectory_vertical_levels,
                dTdt=idTdt,
            )
        except Exception as e:
            app_logger.exception(f"Error creating BoxData object: {e}")
            raise
        app_logger.info("Ok.")

        # Compute and store various meteorological terms
        terms_dict = compute_and_store_terms(box_obj, terms_dict, app_logger)

        # Log that the processing for this time is done
        app_logger.info("Done.\n")

    # Finalize and process results
    results_file, df = finalize_results(
        times, terms_dict, args, results_subdirectory, out_track, app_logger
    )

    if args.cdsapi:
        app_logger.info("Deleting files from CDS...")
        os.remove(args.infile)
        app_logger.info("Done.")

    if args.plots:
        from glob import glob

        from ..plots.map_track import map_track
        from ..plots.plot_boxplot import boxplot_terms
        from ..plots.plot_hovmoller import plot_hovmoller
        from ..plots.plot_LEC import plot_lorenzcycletoolkit
        from ..plots.plot_LPS import plot_LPS
        from ..plots.plot_periods import plot_periods
        from ..plots.timeseries_terms import plot_timeseries
        from ..plots.timeseries_zeta_and_Z import plot_min_zeta_hgt

        app_logger.info("Generating plots..")
        figures_directory = os.path.join(results_subdirectory, "Figures")
        track_file = glob(os.path.join(results_subdirectory, "*trackfile"))[0]

        # Basic diagnostics plots
        map_track(results_file, track_file, figures_directory, app_logger)
        plot_min_zeta_hgt(track_file, figures_directory, app_logger)
        plot_timeseries(results_file, figures_directory, app_logger)
        plot_hovmoller(results_file, figures_directory, app_logger)
        boxplot_terms(results_file, results_subdirectory, figures_directory, app_logger)

        # Determine periods: if vorticity is already processed, don't filter
        # the series
        processed_vorticity = (
            True if args.zeta and "min_max_zeta_850" in track.columns else False
        )
        plot_periods(
            out_track,
            times,
            lat,
            results_subdirectory,
            figures_directory,
            app_logger,
            processed_vorticity=processed_vorticity,
        )

        # Use found periods for other plots
        plot_lorenzcycletoolkit(
            results_file,
            figures_directory,
            os.path.join(results_subdirectory, "periods.csv"),
            app_logger,
        )
        plot_LPS(df, args.infile, results_subdirectory, figures_directory, app_logger)

        app_logger.info("Done.")


if __name__ == "__main__":
    from ..utils.tools import prepare_data

    args = argparse.Namespace(
        infile="samples/Reg1-Representative_NCEP-R2.nc",
        residuals=True,
        fixed=False,
        geopotential=False,
        track=True,
        choose=False,
        zeta=False,
        mpas=False,
        plots=True,
        outname=None,
        verbosity=False,
    )

    app_logger = initialize_logging()

    varlist = "../inputs/namelist_NCEP-R2"
    variable_list_df = pd.read_csv(varlist, sep=";", index_col=0, header=0)

    data, method = prepare_data(args, varlist, app_logger)

    dTdt = data[variable_list_df.loc["Air Temperature"]["Variable"]].differentiate(
        variable_list_df.loc["Time"]["Variable"], datetime_unit="s"
    ) * units("K/s")

    resuts_directory = "../LEC_Results/"
    results_subdirectory = os.path.join(
        resuts_directory,
        "".join(args.infile.split("/")[-1].split(".nc")) + "_" + method,
    )
    figures_directory = os.path.join(results_subdirectory, "Figures")
    os.makedirs(figures_directory, exist_ok=True)
    os.makedirs(results_subdirectory, exist_ok=True)

    lec_moving(
        data, variable_list_df, dTdt, results_subdirectory, figures_directory, args
    )

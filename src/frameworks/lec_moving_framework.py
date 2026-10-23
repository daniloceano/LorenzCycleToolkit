# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lec_moving_framework.py                            :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:32:55 by daniloceano       #+#    #+#              #
#    Updated: 2026/01/25 17:40:51 by daniloceano      ###   ########.fr        #
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
from ..analysis.mass_continuity import MassContinuity
from ..utils.box_data import BoxData
from ..utils.calc_budget_and_residual import calc_budget_diff, calc_residuals
from ..utils.output_names import results_stem
from ..utils.run_metadata import RunMetadata, input_extent, top_cut
from ..utils.run_status import RunStatus
from ..utils.select_area import draw_box_map, plot_domain_attributes
from ..utils.tools import find_extremum_coordinates, initialize_logging
from ..utils.vertical_harmonization import harmonize_levels


def create_terms_dict(args):
    """
    Create a dictionary for storing computation results of various meteorological terms.

    Args:
        args (argparse.Namespace): Script arguments.

    Returns:
        dict: Dictionary with keys for each meteorological term initialized to empty lists.
    """
    energy_terms = ["Az", "Ae", "Kz", "Ke"]
    conversion_terms = ["Cz", "Ca", "Ck", "Ce", "C_overturning"]
    diagnostic_terms = ["M"]
    boundary_terms = ["BAz", "BAe", "BKz", "BKe", "BΦZ", "BΦE"]
    generation_dissipation_terms = (
        ["Gz", "Ge", "Dz", "De"] if not args.residuals else ["Gz", "Ge"]
    )

    terms = (
        energy_terms
        + conversion_terms
        + boundary_terms
        + diagnostic_terms
        + generation_dissipation_terms
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
            app_logger.error("❌ Track time limits do not match with data time limits.")
            raise ValueError("Track time limits do not match with data time limits.")

        # Check longitude limits
        if track["Lon"].max() > data_lon_max:
            app_logger.error(
                f"❌ Track file longitude max limit ({track['Lon'].max():.2f}) exceeds data max "
                f"longitude limit ({data_lon_max:.2f})."
            )
            raise ValueError(
                f"Track file longitude max limit ({track['Lon'].max():.2f}) exceeds data max "
                f"longitude limit ({data_lon_max:.2f})."
            )
        if track["Lon"].min() < data_lon_min:
            app_logger.error(
                f"❌ Track file longitude min limit ({track['Lon'].min():.2f}) is below data min "
                f"longitude limit ({data_lon_min:.2f})."
            )
            raise ValueError(
                f"Track file longitude min limit ({track['Lon'].min():.2f}) is below data min "
                f"longitude limit ({data_lon_min:.2f})."
            )

        # Check latitude limits
        if track["Lat"].max() > data_lat_max:
            app_logger.error(
                f"❌ Track file latitude max limit ({track['Lat'].max():.2f}) exceeds data max "
                f"latitude limit ({data_lat_max:.2f})."
            )
            raise ValueError(
                f"Track file latitude max limit ({track['Lat'].max():.2f}) exceeds data max "
                f"latitude limit ({data_lat_max:.2f})."
            )
        if track["Lat"].min() < data_lat_min:
            app_logger.error(
                f"❌ Track file latitude min limit ({track['Lat'].min():.2f}) is below data min "
                f"latitude limit ({data_lat_min:.2f})."
            )
            raise ValueError(
                f"Track file latitude min limit ({track['Lat'].min():.2f}) is below data min "
                f"latitude limit ({data_lat_min:.2f})."
            )

        return track

    except FileNotFoundError:
        app_logger.error(f"❌ Track file {trackfile} not found.")
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


def select_time(data, dTdt, TimeName, t):
    """
    Select one time step of the data and of the temperature tendency.

    Args:
        data (xr.Dataset): The data of the run.
        dTdt (xr.DataArray): Temperature tendency.
        TimeName (str): Name of the time coordinate.
        t (pd.Timestamp): Time step.

    Returns:
        tuple: The data and the temperature tendency at ``t``.

    Raises:
        KeyError: If ``t`` is not a time of the data.
    """
    idata, idTdt = data.sel({TimeName: t}), dTdt.sel({TimeName: t})
    if idata[TimeName].shape != ():
        idata, idTdt = data.isel({TimeName: 1}), dTdt.isel({TimeName: 1})
    return idata, idTdt


def fields_at_850(idata, variable_list_df, args, VerticalCoordIndexer, LatIndexer, LonIndexer):
    """
    Wind, geopotential height, wind speed and vorticity at 850 hPa for one time step.

    They locate the system and draw the domain maps; the energetics are
    computed over the levels of the vertical harmonization.

    Returns:
        dict: ``izeta_850``, ``ihgt_850``, ``iwspd_850``, ``iu_850``,
        ``iv_850``, ``lat`` and ``lon``.
    """
    iu, iv, ihgt = extract_wind_and_height_components(idata, variable_list_df, args)
    iu_850, iv_850, ihgt_850 = (
        iu.sel({VerticalCoordIndexer: 85000}),
        iv.sel({VerticalCoordIndexer: 85000}),
        ihgt.sel({VerticalCoordIndexer: 85000}),
    )
    iwspd_850, izeta_850 = (
        wind_speed(iu_850, iv_850),
        vorticity(iu_850, iv_850).metpy.dequantify(),
    )
    return {
        "izeta_850": izeta_850,
        "ihgt_850": ihgt_850,
        "iwspd_850": iwspd_850,
        "iu_850": iu_850,
        "iv_850": iv_850,
        "lat": idata[LatIndexer],
        "lon": idata[LonIndexer],
    }


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

            # Get mininum vorticity if in the southern hemisphere
            if limits["min_lat"] < 0:
                min_max_zeta_unformatted = izeta_850_slice.min()
            # Get maximum vorticity if in the northern hemisphere
            else:
                min_max_zeta_unformatted = izeta_850_slice.max()

        # Get minimum vorticity if in the southern hemisphere
        if limits["min_lat"] < 0:
            min_max_zeta = float(np.nanmin(min_max_zeta_unformatted))
        # Get maximum vorticity if in the northern hemisphere
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
    app_logger.info("⚡ Computing Energy Contents...")
    try:
        ec_obj = EnergyContents(box_obj, method="moving", app_logger=app_logger)
        terms_dict["Az"].append(ec_obj.calc_az())
        terms_dict["Ae"].append(ec_obj.calc_ae())
        terms_dict["Kz"].append(ec_obj.calc_kz())
        terms_dict["Ke"].append(ec_obj.calc_ke())
    except Exception as e:
        app_logger.exception(f"❌ Error in computing Energy Contents: {e}")
        raise

    # Conversion Terms
    app_logger.info("🔄 Computing Conversion Terms...")
    try:
        ct_obj = ConversionTerms(box_obj, method="moving", app_logger=app_logger)
        terms_dict["Cz"].append(ct_obj.calc_cz())
        terms_dict["Ca"].append(ct_obj.calc_ca())
        terms_dict["Ck"].append(ct_obj.calc_ck())
        terms_dict["Ce"].append(ct_obj.calc_ce())
        terms_dict["C_overturning"].append(ct_obj.calc_c_overturning())
    except Exception as e:
        app_logger.exception(f"❌ Error in computing Conversion Terms: {e}")
        raise

    # Boundary Terms
    app_logger.info("🏁 Computing Boundary Terms...")
    try:
        bt_obj = BoundaryTerms(box_obj, method="moving", app_logger=app_logger)
        terms_dict["BAz"].append(bt_obj.calc_baz())
        terms_dict["BAe"].append(bt_obj.calc_bae())
        terms_dict["BKz"].append(bt_obj.calc_bkz())
        terms_dict["BKe"].append(bt_obj.calc_bke())
        terms_dict["BΦZ"].append(bt_obj.calc_boz())
        terms_dict["BΦE"].append(bt_obj.calc_boe())
    except Exception as e:
        app_logger.exception(f"❌ Error in computing Boundary Terms: {e}")
        raise

    # Mass-continuity diagnostic (profile plus column-integrated residual).
    app_logger.info("⚖️ Computing mass-continuity residual...")
    try:
        terms_dict["M"].append(
            MassContinuity(box_obj, method="moving", app_logger=app_logger)
            .calc_mass_residual()
        )
    except Exception as e:
        app_logger.exception(f"❌ Error in computing mass continuity: {e}")
        raise

    # Generation/Dissipation Terms
    app_logger.info("🔥 Computing Generation/Dissipation Terms...")
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
        app_logger.exception(f"❌ Error in computing Generation/Dissipation Terms: {e}")
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
    app_logger.info("📈 Estimating budget terms...")
    df = calc_budget_diff(df, times, app_logger)

    # Computing residuals, if required
    if args.residuals:
        app_logger.info("🧮 Computing residuals...")
        df = calc_residuals(df, app_logger)

    # Constructing output filename
    method = "track" if args.track else "choose"
    results_file = os.path.join(results_subdirectory, f"{results_stem(args, method)}_results.csv")

    # Saving the DataFrame to a CSV file
    df.to_csv(results_file)
    app_logger.info(f"💾 Results saved to {results_file}")

    # Save system position as a csv file for replicability
    out_track = out_track.rename(
        columns={"datestr": "time", "central_lat": "Lat", "central_lon": "Lon"}
    )
    output_trackfile = os.path.join(results_subdirectory, f"{results_stem(args, method)}_trackfile")
    out_track.to_csv(output_trackfile, index=False, sep=";")
    app_logger.info(f"📍 System track saved to {output_trackfile}")

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
    run_status: RunStatus = None,
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
        run_status (RunStatus): Status record of the run (outputs and figures).

    Returns:
        None
    """
    app_logger.info("🌀 Computing energetics using moving framework...")
    if run_status is None:
        run_status = RunStatus(args, results_subdirectory, app_logger)

    # Indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = (
        variable_list_df.loc["Longitude"]["Variable"],
        variable_list_df.loc["Latitude"]["Variable"],
        variable_list_df.loc["Time"]["Variable"],
        variable_list_df.loc["Vertical Level"]["Variable"],
    )

    # Get the time array
    times = pd.to_datetime(data[TimeName].values)

    # Get the system track and check for errors
    if args.track:
        track = handle_track_file(
            data, times, LonIndexer, LatIndexer, TimeName, args, app_logger
        )

    # The box of every time step is decided before any term is computed, so
    # that all time steps share one vertical control volume (DEC-10 B, in every
    # framework): a level incomplete inside the box of any time step is excluded
    # from all of them. With --choose every box is drawn first.
    mode = "track" if args.track else "choose"
    if args.choose:
        app_logger.info(
            f"🖱️ --choose mode: draw the box of each of the {len(times)} time steps; the "
            "computation starts after the last one, so that the levels can be decided for "
            "all of them together."
        )
    boxes = {}
    for t in times:
        if args.choose:
            idata, _ = select_time(data, dTdt, TimeName, t)
            data850 = fields_at_850(idata, variable_list_df, args, VerticalCoordIndexer, LatIndexer, LonIndexer)
        else:
            data850 = None
        boxes[t] = get_limits(args, t, data850, track if args.track else None)
    harmonization = harmonize_levels(
        [
            (t, select_time(data, dTdt, TimeName, t)[0],
             (boxes[t]["min_lon"], boxes[t]["max_lon"], boxes[t]["min_lat"], boxes[t]["max_lat"]))
            for t in times
        ],
        variable_list_df,
        app_logger,
        mode,
    )

    # Create csv files for storing vertical results. As in the fixed framework,
    # the headers list the levels actually used.
    level_files = []
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
        "C_overturning",
        "M",
        "Ck",
        "Ck_1",
        "Ck_2",
        "Ck_3",
        "Ck_4",
        "Ck_5",
    ]:
        columns = [TimeName] + harmonization.used_levels
        df = pd.DataFrame(columns=columns)
        file_name = term + "_" + VerticalCoordIndexer + ".csv"
        file_path = os.path.join(results_subdirectory_vertical_levels, file_name)
        df.to_csv(file_path, index=None)
        level_files.append(file_path)
        app_logger.info(f"{file_path} created (but still empty)")

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

    # What each time step actually used (domain_selection.csv, lec_metadata.json)
    metadata = RunMetadata(args, mode, variable_list_df, input_extent(data), top_cut(data))
    metadata.record_harmonization(harmonization)

    # Iterating over times
    for t in times:
        app_logger.info(f"⏰ Processing data at time: {t}...")
        try:
            idata, idTdt = select_time(data, dTdt, TimeName, t)
        except KeyError as e:
            app_logger.error(f"❌ Time indexing error: {e}")
            continue

        # Wind, geopotential height, vorticity and wind speed at 850 hPa
        data850 = fields_at_850(idata, variable_list_df, args, VerticalCoordIndexer, LatIndexer, LonIndexer)
        izeta_850, ihgt_850, iwspd_850 = data850["izeta_850"], data850["ihgt_850"], data850["iwspd_850"]
        lat = data850["lat"]

        # Box of the current time (decided before the loop)
        limits = boxes[t]
        app_logger.info(
            f"🗺️ Box: center=({limits['central_lat']:.2f}, {limits['central_lon']:.2f}), "
            f"size={limits['length']}°x{limits['width']}°"
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
            f"🌡️ 850 hPa: ζ={position['min_max_zeta_850']['value']:.2e}, "
            f"hgt={position['min_hgt_850']['value']:.0f}m, "
            f"wspd={position['max_wind_850']['value']:.1f}m/s"
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
        # geopotential height. Only this call is guarded: a cartographic
        # failure (e.g. Natural Earth data unavailable offline) is recorded
        # and the computation of this time step goes on (E0-05).
        if args.plots:
            run_status.run_product(
                "domain_maps", plot_domain_attributes, data850, limits, position, figures_directory,
                expected=[f"Figures/boxes/box_{limits['datestr']}.png"], step=limits["datestr"],
            )

        # Create box object over the levels common to all time steps
        app_logger.info("📦 Creating box object...")
        try:
            box_obj = BoxData(
                data=harmonization.restrict(idata).compute(),
                variable_list_df=variable_list_df,
                args=args,
                western_limit=limits["min_lon"],
                eastern_limit=limits["max_lon"],
                southern_limit=limits["min_lat"],
                northern_limit=limits["max_lat"],
                results_subdirectory=results_subdirectory,
                results_subdirectory_vertical_levels=results_subdirectory_vertical_levels,
                dTdt=harmonization.restrict(idTdt),
            )
        except Exception as e:
            app_logger.exception(f"❌ Error creating BoxData object: {e}")
            raise
        if box_obj.dropped_levels:
            raise RuntimeError(
                f"BoxData dropped {box_obj.dropped_levels} at {t} although the vertical "
                "harmonization found these levels complete in every box"
            )
        app_logger.info("✅ Box object created.")
        metadata.record_box(
            box_obj,
            {"west": limits["min_lon"], "east": limits["max_lon"],
             "south": limits["min_lat"], "north": limits["max_lat"]},
            t.strftime("%Y-%m-%d %H:%M:%S"),
            idata[VerticalCoordIndexer].values,
            app_logger,
        )

        # Compute and store various meteorological terms
        terms_dict = compute_and_store_terms(box_obj, terms_dict, app_logger)

        # Log that the processing for this time is done
        app_logger.info("✅ Timestep complete.\n")

    if args.plots:
        run_status.close_product("domain_maps")

    # Finalize and process results
    results_file, df = finalize_results(
        times, terms_dict, args, results_subdirectory, out_track, app_logger
    )
    method = "track" if args.track else "choose"
    track_file = os.path.join(results_subdirectory, f"{results_stem(args, method)}_trackfile")
    domain_file = metadata.write_domain_selection(results_subdirectory)
    metadata_file = metadata.write_metadata(
        results_subdirectory, results_file, level_files, df.columns, df.index.values, trackfile=track_file
    )
    app_logger.info(f"🧾 Metadata saved to {metadata_file} and {domain_file}")
    run_status.record_scientific_outputs(
        {
            "results_csv": [results_file],
            "trackfile": [track_file],
            "vertical_level_csvs": level_files,
            "lec_metadata": [metadata_file],
            "domain_selection": [domain_file],
        }
    )

    if args.cdsapi:
        app_logger.info("🗑️ Cleaning up CDS API temporary files...")
        os.remove(args.infile)
        app_logger.info("✅ Cleanup complete.")

    if args.plots:
        app_logger.info("🎨 Generating plots...")
        figures_directory = os.path.join(results_subdirectory, "Figures")
        generate_moving_figures(
            run_status, args, results_file, track_file, df, out_track, times, lat,
            track if args.track else None, results_subdirectory, figures_directory, app_logger,
        )
        if run_status.figures_ok() and not run_status.unavailable_products():
            app_logger.info("🎉 All plots generated successfully!")
        elif run_status.figures_ok():
            app_logger.info(
                "🎨 Plots generated; unavailable for these data: "
                + ", ".join(run_status.unavailable_products())
                + f" (see {run_status.path})"
            )
        else:
            app_logger.warning(
                "⚠️ Plots finished with failures: "
                + ", ".join(run_status.failed_figures())
                + f" (see {run_status.path})"
            )


def generate_moving_figures(
    run_status, args, results_file, track_file, df, out_track, times, lat, track,
    results_subdirectory, figures_directory, app_logger,
):
    """
    Generate the figures of the moving framework after the results are saved.

    Each figure is recorded in run_status.json and counts as generated only if
    its call returns and writes its expected files. A failure is logged and
    the remaining figures are still attempted. The Lorenz cycle diagrams of
    the periods and the Lorenz Phase Space need periods.csv: when the periods
    were not determined they are reported as failed because of that.
    """
    figures = (
        "track_map", "zeta_hgt_timeseries", "timeseries", "hovmoller", "boxplots",
        "periods", "lec_cycle", "lec_periods", "lps",
    )
    try:
        from ..plots.map_track import map_track
        from ..plots.plot_boxplot import boxplot_terms
        from ..plots.plot_hovmoller import plot_hovmoller
        from ..plots.plot_LEC import (plot_lorenzcycletoolkit,
                                      plot_lorenzcycletoolkit_periods)
        from ..plots.plot_LPS import plot_LPS
        from ..plots.plot_periods import cyclophaser_inputs, plot_periods
        from ..utils.cyclophaser_compat import CyclophaserCompatibilityError
        from ..utils.cyclophaser_compat import assess as assess_cyclophaser
        from ..plots.timeseries_terms import plot_timeseries
        from ..plots.timeseries_zeta_and_Z import plot_min_zeta_hgt
        from ..plots.utils import TERM_DETAILS
        from .lec_fixed_framework import lec_daily_figures
    except Exception as exc:
        run_status.fail_products(figures, exc, "plotting modules could not be imported")
        return

    # Basic diagnostics plots
    run_status.run_product(
        "track_map", map_track, results_file, track_file, figures_directory, app_logger,
        expected=["Figures/track.png"],
    )
    run_status.run_product(
        "zeta_hgt_timeseries", plot_min_zeta_hgt, track_file, figures_directory, app_logger,
        expected=["Figures/timeseries-min_zeta_hgt.png"],
    )
    run_status.run_product(
        "timeseries", plot_timeseries, results_file, figures_directory, app_logger,
        expected=[f"Figures/timeseries/timeseires_{k}.png" for k in TERM_DETAILS],
    )
    run_status.run_product(
        "hovmoller", plot_hovmoller, results_file, figures_directory, app_logger,
        expected=["Figures/hovmollers/hovmoller_energy.png"],
    )
    run_status.run_product(
        "boxplots", boxplot_terms, results_file, results_subdirectory, figures_directory, app_logger,
        expected=[
            "Figures/boxplots/boxplot_terms_energy.png",
            "Figures/boxplots/boxplot_vertical_levels_energy.png",
            "Figures/boxplots/vertical_timeseries_energy.png",
        ],
    )

    # Determine periods: if vorticity is already processed, don't filter
    # the series. Before calling cyclophaser, check that the series is long
    # enough for the tested version and options (E0-06): below the minimum
    # the periods are unavailable and cyclophaser is not called.
    periods_file = os.path.join(results_subdirectory, "periods.csv")
    periods_expected = ["periods.csv", "Figures/Periods/periods.png", "Figures/Periods/periods_steps.png"]
    try:
        processed_vorticity = (
            True if args.zeta and "min_max_zeta_850" in track.columns else False
        )
        vorticity_data, branch, options, _ = cyclophaser_inputs(out_track, lat, processed_vorticity)
        assessment = assess_cyclophaser(vorticity_data, times, branch, options)
    except CyclophaserCompatibilityError as exc:
        run_status.fail_products(["periods"], exc, "cyclophaser compatibility check failed")
    except Exception as exc:
        run_status.fail_products(["periods"], exc, "cyclophaser inputs could not be prepared")
    else:
        if not assessment.applicable:
            run_status.mark_unavailable("periods", assessment.reason, expected=periods_expected)
        else:
            run_status.run_product(
                "periods", plot_periods, out_track, times, lat, results_subdirectory,
                figures_directory, app_logger, processed_vorticity=processed_vorticity,
                expected=periods_expected,
            )
            run_status.add_notes("periods", assessment.notes)

    # Daily Lorenz cycle diagrams (do not need the periods)
    run_status.run_product(
        "lec_cycle", plot_lorenzcycletoolkit, results_file, figures_directory, app_logger=app_logger,
        expected=lec_daily_figures(df.index),
    )

    # Use found periods for other plots
    dt_hours = int((df.index[1] - df.index[0]).total_seconds() / 3600) if len(df.index) > 1 else None
    lps_figures = [
        f"Figures/LPS/LPS_{kind}{zoom}.png"
        for kind in (f"{dt_hours}h", "periods", "1d") for zoom in ("", "_zoom")
    ]
    if not run_status.succeeded("periods"):
        run_status.dependency_failed("lec_periods", "periods")
        run_status.dependency_failed("lps", "periods", expected=lps_figures)
        return

    try:
        expected_periods = period_figures(periods_file, df.index)
    except Exception as exc:
        run_status.fail_products(["lec_periods"], exc, "periods.csv could not be read")
    else:
        run_status.run_product(
            "lec_periods", plot_lorenzcycletoolkit_periods, results_file, figures_directory, periods_file,
            app_logger, expected=expected_periods,
        )
    run_status.run_product(
        "lps", plot_LPS, df, args.infile, results_subdirectory, figures_directory, app_logger,
        expected=lps_figures,
    )


def period_figures(periods_file, index):
    """
    Lorenz cycle diagrams expected for the periods in periods.csv: one per
    period that contains at least one time of the results (periods without
    data are skipped, with a warning, by plot_period_means).
    """
    periods = pd.read_csv(periods_file, parse_dates=["start", "end"], index_col=0)
    times = pd.DatetimeIndex(index)
    return [
        f"Figures/LEC/LEC_{name}.png"
        for name, row in periods.iterrows()
        if ((times >= row["start"]) & (times <= row["end"])).any()
    ]


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

# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lec_moving_framework.py                            :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:32:55 by daniloceano       #+#    #+#              #
#    Updated: 2024/01/03 00:47:39 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import logging
import argparse
import pandas as pd
import xarray as xr
import numpy as np

from pathlib import Path

from metpy.units import units
from metpy.calc import vorticity, wind_speed
from metpy.constants import g

from ..utils.select_area import draw_box_map, plot_domain_attributes
from ..utils.tools import find_extremum_coordinates, initialize_logging
from ..utils.calc_budget_and_residual import calc_budget_diff, calc_residuals
from ..utils.box_data import BoxData
from ..analysis.energy_contents import EnergyContents
from ..analysis.conversion_terms import ConversionTerms
from ..analysis.boundary_terms import BoundaryTerms
from ..analysis.generation_and_dissipation_terms import GenerationDissipationTerms

def create_terms_dict(args):
    """
    Create a dictionary for storing computation results of various meteorological terms.

    Args:
        args (argparse.Namespace): Script arguments.

    Returns:
        dict: Dictionary with keys for each meteorological term initialized to empty lists.
    """
    energy_terms = ['Az', 'Ae', 'Kz', 'Ke']
    conversion_terms = ['Cz', 'Ca', 'Ck', 'Ce']
    boundary_terms = ['BAz', 'BAe', 'BKz', 'BKe', 'BΦZ', 'BΦE']
    generation_dissipation_terms = ['Gz', 'Ge', 'Dz', 'De'] if not args.residuals else ['Gz', 'Ge']

    terms = energy_terms + conversion_terms + boundary_terms + generation_dissipation_terms
    return {term: [] for term in terms}

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
    u_var = variable_list_df.loc['Eastward Wind Component']['Variable']
    v_var = variable_list_df.loc['Northward Wind Component']['Variable']
    geopotential_var = variable_list_df.loc['Geopotential Height']['Variable'] if not args.geopotential else variable_list_df.loc['Geopotential']['Variable']

    iu = (idata[u_var].compute() * units(variable_list_df.loc['Eastward Wind Component']['Units']).to('m/s'))
    iv = (idata[v_var].compute() * units(variable_list_df.loc['Northward Wind Component']['Units']).to('m/s'))
    ihgt = (idata[geopotential_var].compute() * units(
        variable_list_df[variable_list_df['Variable'] == geopotential_var]['Units'].iloc[0])
        )
    
    if args.geopotential:
        ihgt = (ihgt / g).metpy.convert_units('gpm')

    return iu, iv, ihgt

def get_limits(args, t, data850, iu_850, iv_850, track=None):
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
        central_lat, central_lon = track_row['Lat'], track_row['Lon']

        # Use provided width and length if available, otherwise default values
        width = track_row.get('width', 15)  # Default width if not specified
        length = track_row.get('length', 15)  # Default length if not specified

    elif args.choose:
        # For the 'choose' option, interactively determine the box limits
        izeta_850, ihgt_850, = data850['min_max_zeta']['data'], data850['min_hgt']['data']
        limits = draw_box_map(iu_850, iv_850, izeta_850, ihgt_850, data850['lat'], data850['lon'], t)
        central_lat, central_lon = (limits['max_lat'] + limits['min_lat']) / 2, (limits['max_lon'] + limits['min_lon']) / 2
        width, length = limits['max_lon'] - limits['min_lon'], limits['max_lat'] - limits['min_lat']

    # Calculate the bounding box limits
    min_lon, max_lon = central_lon - (width / 2), central_lon + (width / 2)
    min_lat, max_lat = central_lat - (length / 2), central_lat + (length / 2)

    # Format time for output
    datestr = pd.to_datetime(t).strftime('%Y-%m-%d-%H%M')

    limits = {
        'datestr': datestr,
        'central_lat': central_lat,
        'central_lon': central_lon,
        'length': length,
        'width': width,
        'min_lon': min_lon,
        'max_lon': max_lon,
        'min_lat': min_lat,
        'max_lat': max_lat
    }

    return limits

def get_position(track, limits, izeta_850, ihgt_850, iwspd_850, LatIndexer, LonIndexer, args):
    """
    Retrieves or calculates the values of 'min_max_zeta_850', 'min_hgt_850', and 'max_wind_850' based on the given track and data slices.

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
    track_itime = track.iloc[-1]

    central_lat = limits['central_lat']
    central_lon = limits['central_lon']
    min_lat = limits['min_lat']
    min_lon = limits['min_lon']
    max_lat = limits['max_lat']
    max_lon = limits['max_lon']

    # Slice data for defined box
    izeta_850_slice = izeta_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
    ihgt_850_slice = ihgt_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
    iwspd_850_slice = iwspd_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})

    try:
        min_max_zeta = float(track.loc[track_itime]['min_max_zeta_850'])
    except KeyError:
        if args.zeta:
            min_max_zeta_unformatted = izeta_850.sel(latitude=central_lat, longitude=central_lon, method='nearest')
        else:
            min_max_zeta_unformatted = izeta_850_slice.min()
        if min_lat < 0:
            min_max_zeta = float(np.nanmin(min_max_zeta_unformatted))
        else:
            min_max_zeta = float(np.nanmax(min_max_zeta_unformatted))

    try:
        min_hgt = float(track.loc[track_itime]['min_hgt_850'])
    except KeyError:
        min_hgt = float(ihgt_850_slice.min())

    try:
        max_wind = float(track.loc[track_itime]['max_wind_850'])
    except KeyError:
        max_wind = float(iwspd_850_slice.max())

    position = {
        'min_max_zeta_850': min_max_zeta,
        'min_hgt_850': min_hgt,
        'max_wind_850': max_wind
    }

    return position

def construct_data850(izeta_850, ihgt_850, iwspd_850, lat, lon):
    """
    Construct a dictionary with key meteorological features at 850 hPa and their corresponding
    latitude, longitude, and data.

    Args:
        izeta_850 (xr.DataArray): Vorticity data at 850 hPa.
        ight_850 (xr.DataArray): Geopotential height data at 850 hPa.
        iwspd_850 (xr.DataArray): Wind speed data at 850 hPa.
        lat (xr.DataArray): Latitude data array.
        lon (xr.DataArray): Longitude data array.

    Returns:
        dict: Dictionary containing meteorological features at 850 hPa with their coordinates and data.
    """
    # Find extremum coordinates for each feature
    min_max_zeta_lat, min_max_zeta_lon = find_extremum_coordinates(izeta_850, lat, lon, 'min_max_zeta')
    min_hgt_lat, min_hgt_lon = find_extremum_coordinates(ihgt_850, lat, lon, 'min_hgt')
    max_wind_lat, max_wind_lon = find_extremum_coordinates(iwspd_850, lat, lon, 'max_wind')

    data850 = {
        'min_max_zeta': {
            'latitude': min_max_zeta_lat,
            'longitude': min_max_zeta_lon,
            'data': izeta_850
        },
        'min_hgt': {
            'latitude': min_hgt_lat,
            'longitude': min_hgt_lon,
            'data': ihgt_850
        },
        'max_wind': {
            'latitude': max_wind_lat,
            'longitude': max_wind_lon,
            'data': iwspd_850
        },
        'lat': lat,
        'lon': lon,
    }

    return data850 

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
        ec_obj = EnergyContents(box_obj, method='moving', app_logger=app_logger)
        terms_dict['Az'].append(ec_obj.calc_az())
        terms_dict['Ae'].append(ec_obj.calc_ae())
        terms_dict['Kz'].append(ec_obj.calc_kz())
        terms_dict['Ke'].append(ec_obj.calc_ke())
    except Exception as e:
        app_logger.exception(f"Error in computing Energy Contents: {e}")
        raise

    # Conversion Terms
    app_logger.info("Computing Conversion Terms...")
    try:
        ct_obj = ConversionTerms(box_obj, method='moving', app_logger=app_logger)
        terms_dict['Cz'].append(ct_obj.calc_cz())
        terms_dict['Ca'].append(ct_obj.calc_ca())
        terms_dict['Ck'].append(ct_obj.calc_ck())
        terms_dict['Ce'].append(ct_obj.calc_ce())
    except Exception as e:
        app_logger.exception(f"Error in computing Conversion Terms: {e}")
        raise

    # Boundary Terms
    app_logger.info("Computing Boundary Terms...")
    try:
        bt_obj = BoundaryTerms(box_obj, method='moving', app_logger=app_logger)
        terms_dict['BAz'].append(bt_obj.calc_baz())
        terms_dict['BAe'].append(bt_obj.calc_bae())
        terms_dict['BKz'].append(bt_obj.calc_bkz())
        terms_dict['BKe'].append(bt_obj.calc_bke())
        terms_dict['BΦZ'].append(bt_obj.calc_boz())
        terms_dict['BΦE'].append(bt_obj.calc_boe())
    except Exception as e:
        app_logger.exception(f"Error in computing Boundary Terms: {e}")
        raise

    # Generation/Dissipation Terms
    app_logger.info("Computing Generation/Dissipation Terms...")
    try:
        gdt_obj = GenerationDissipationTerms(box_obj, method='moving', app_logger=app_logger)
        terms_dict['Gz'].append(gdt_obj.calc_gz())
        terms_dict['Ge'].append(gdt_obj.calc_ge())
        if not box_obj.args.residuals:
            terms_dict['Dz'].append(gdt_obj.calc_dz())
            terms_dict['De'].append(gdt_obj.calc_de())
    except Exception as e:
        app_logger.exception(f"Error in computing Generation/Dissipation Terms: {e}")
        raise

    return terms_dict

def finalize_results(times, terms_dict, args, results_subdirectory, out_track, app_logger):
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
    df = pd.DataFrame(terms_dict, index=pd.to_datetime(times), dtype = float)

    # Estimating budget terms (∂X/∂t) using finite differences
    app_logger.info("Estimating budget terms...")
    df = calc_budget_diff(df, times, app_logger)
    
    # Computing residuals, if required
    if args.residuals:
        app_logger.info("Computing residuals...")
        df = calc_residuals(df, app_logger)

    # Constructing output filename
    method = 'track' if args.track else 'choose'
    infile_name = os.path.basename(args.infile).split('.nc')[0]
    results_filename = ''.join(f'{infile_name}_{method}_results')
    results_file = os.path.join(results_subdirectory, f'{results_filename}.csv')

    # Saving the DataFrame to a CSV file
    df.to_csv(results_file)
    app_logger.info(f'Results saved to {results_file}')

    # Save system position as a csv file for replicability
    out_track = out_track.rename(columns={'datestr':'time','central_lat':'Lat','central_lon':'Lon'})
    output_trackfile = os.path.join(results_subdirectory, results_filename+'_trackfile')
    out_track.to_csv(output_trackfile, index=False, sep=";")
    app_logger.info(f'System track saved to {output_trackfile}')

    return results_file, df
    

def lec_moving(data: xr.Dataset, variable_list_df: pd.DataFrame, dTdt: xr.Dataset,
               results_subdirectory: str, figures_directory: str,
               app_logger: logging.Logger, args: argparse.Namespace):
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
    app_logger.info('Computing energetics using moving framework')

    # Indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = (
        variable_list_df.loc['Longitude']['Variable'],
        variable_list_df.loc['Latitude']['Variable'],
        variable_list_df.loc['Time']['Variable'],
        variable_list_df.loc['Vertical Level']['Variable']
    )
    
    PressureData = data[VerticalCoordIndexer]
    
    # Create csv files for storing vertical results
    for term in ['Az', 'Ae', 'Kz', 'Ke', 'Cz', 'Ca', 'Ck', 'Ce', 'Ge', 'Gz']:
        columns = [TimeName] + [float(i) for i in PressureData.values]
        df = pd.DataFrame(columns=columns)
        file_name = term + '_' + VerticalCoordIndexer + '.csv'
        file_path = os.path.join(results_subdirectory, file_name)
        df.to_csv(file_path, index=None) 
        app_logger.info(f'{file_path} created (but still empty)')

    # Track file handling
    if args.track:
        trackfile = 'inputs/track'
        try:
            track = pd.read_csv(trackfile, parse_dates=[0], delimiter=';', index_col='time')
        except FileNotFoundError:
            app_logger.error(f"Track file {trackfile} not found.")
            raise

    # Dictionary for saving system position and attributes
    results_keys = ['datestr', 'central_lat', 'central_lon', 'length', 'width',
                    'min_max_zeta_850', 'min_hgt_850', 'max_wind_850']
    out_track = pd.DataFrame(columns=results_keys)

    # Dictionary for storing results
    terms_dict = create_terms_dict(args)

    # Slice the time array to match the track file
    times = pd.to_datetime(data[TimeName].values)
    if args.track:
        times = times[(times >= track.index[0]) & (times <= track.index[-1])]
        if len(times) == 0:
            app_logger.error("Mismatch between trackfile and data times.")
            raise ValueError("Mismatch between trackfile and data times.")
        
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
            ihgt.sel({VerticalCoordIndexer: 85000})
        )

        # Compute wind speed and vorticity at 850 hPa
        iwspd_850, izeta_850 = wind_speed(iu_850, iv_850), vorticity(iu_850, iv_850).metpy.dequantify()
        lat, lon = idata[LatIndexer], idata[LonIndexer]
        data850 = construct_data850(izeta_850, ihgt_850, iwspd_850, lat, lon)

        # Get box attributes for current time
        limits = get_limits(args, t, data850, iu_850, iv_850, track if args.track else None)
        app_logger.info(
            f"central lat: {limits['central_lat']}, central lon: {limits['central_lon']}, "
            f"size: {limits['length']} x {limits['width']}, "
            f"lon range: {limits['min_lon']} to {limits['max_lon']}, "
            f"lat range: {limits['min_lat']} to {limits['max_lat']}"
            )


        # Get position of  850 hPaextreme values for current time
        position = get_position(track, limits, izeta_850, ihgt_850, iwspd_850, LatIndexer, LonIndexer, args)
        app_logger.info(
            f"850 hPa diagnostics --> "
            f"min/max ζ: {position['min_max_zeta_850']:.2e}, "
            f"min geopotential height: {position['min_hgt_850']:.0f}, "
            f"max wind speed: {position['max_wind_850']:.4f}"
        )

        # Store results
        limits_and_position = {**limits, **position}
        new_entry = pd.DataFrame([limits_and_position], index=[t.strftime('%Y-%m-%d-%H%M')])
        if out_track.empty:
            out_track = new_entry
        else:
            out_track = pd.concat([out_track, new_entry], ignore_index=True)

        # Save figure with the domain box, extreme values, vorticity and geopotential height
        plot_domain_attributes(data850, limits, figures_directory)

        # Create box object
        app_logger.info("Creating box object...")
        try:
            box_obj = BoxData(
                data=idata.compute(),
                variable_list_df=variable_list_df,
                args=args,
                western_limit=limits['min_lon'],
                eastern_limit=limits['max_lon'],
                southern_limit=limits['min_lat'],
                northern_limit=limits['max_lat'],
                output_dir=results_subdirectory,
                dTdt=idTdt
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
    results_file, df  = finalize_results(times, terms_dict, args, results_subdirectory, out_track, app_logger)

    if args.plots:
        from ..plots.timeseries_terms import plot_timeseries
        from ..plots.timeseries_zeta_and_Z import plot_min_zeta_hgt
        from ..plots.map_box_limits import plot_box_limits
        from ..plots.plot_boxplot import boxplot_terms
        from ..plots.plot_periods import plot_periods
        from ..plots.plot_LPS import plot_LPS

        app_logger.info('Generating plots..')
        figures_directory = os.path.join(results_subdirectory, 'Figures')
        plot_timeseries(results_file, figures_directory, app_logger)
        boxplot_terms(results_file, results_subdirectory, figures_directory, app_logger)
        plot_periods(out_track, times, lat, figures_directory, results_subdirectory, app_logger)
        plot_LPS(df, args.infile, results_subdirectory, figures_directory, app_logger)


        app_logger.info('Done.')

if __name__ == '__main__':
    from ..utils.tools import prepare_data

    args = argparse.Namespace(
        infile='samples/Reg1-Representative_NCEP-R2.nc',
        residuals=True,
        fixed=False,
        geopotential=False,
        track=True,
        choose=False,
        zeta=False,
        mpas=False,
        plots=True,
        outname=None,
        verbosity=False
    )

    initialize_logging()

    varlist = "../inputs/fvars_NCEP-R2"
    variable_list_df = pd.read_csv(varlist, sep=';', index_col=0, header=0)

    data, method = prepare_data(args, varlist)

    dTdt =  data[variable_list_df.loc['Air Temperature']['Variable']].differentiate(
                variable_list_df.loc['Time']['Variable'],datetime_unit='s') * units('K/s')
    
    resuts_directory = "../LEC_Results/"
    results_subdirectory = os.path.join(
        resuts_directory, "".join(args.infile.split('/')[-1].split('.nc')) + '_' + method)
    figures_directory = os.path.join(results_subdirectory, 'Figures')
    os.makedirs(figures_directory, exist_ok=True)
    os.makedirs(results_subdirectory, exist_ok=True)

    lec_moving(data, variable_list_df, dTdt, results_subdirectory, figures_directory, args)
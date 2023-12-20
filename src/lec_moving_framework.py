# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lec_moving_framework.py                            :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:32:55 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/20 15:23:59 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import logging
import argparse
import pandas as pd
import xarray as xr
import numpy as np

from metpy.units import units
from metpy.calc import vorticity, wind_speed
from metpy.constants import g

from cyclophaser import determine_periods
from select_area import draw_box_map, plot_domain_attributes
from tools import find_extremum_coordinates, initialize_logging
from EnergyContents import EnergyContents
from ConversionTerms import ConversionTerms
from BoundaryTerms import BoundaryTerms
from GenerationDissipationTerms import GenerationDissipationTerms
from BoxData import BoxData
from BudgetResidual import calc_budget_diff, calc_residuals

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
    ihgt = (idata[geopotential_var].compute() * units(variable_list_df.loc[geopotential_var]['Units']))
    
    if args.geopotential:
        ihgt = (ihgt / g).metpy.convert_units('gpm')

    return iu, iv, ihgt

def get_position_and_limits(args, t, track, iu_850, iv_850, izeta_850, ihgt_850, lat, lon):
    """
    Determine the central position of the system at a specific time and the box dimensions.
    This function calculates the central latitude and longitude, as well as the dimensions
    of the box (width and length) based on the tracking data or user input.

    Args:
        args (argparse.Namespace): Script arguments.
        t (pd.Timestamp): Current time.
        track (pd.DataFrame): DataFrame containing the tracking information.
        iu_850, iv_850, izeta_850, ihgt_850 (xr.DataArray): Data arrays for specific calculations.
        lat, lon (xr.DataArray): Latitude and longitude data arrays.

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
        # This part will depend on how you implement the interactive selection
        limits = draw_box_map(iu_850, iv_850, izeta_850, ihgt_850, lat, lon, t)
        central_lat, central_lon = (limits['max_lat'] + limits['min_lat']) / 2, (limits['max_lon'] + limits['min_lon']) / 2
        width, length = limits['max_lon'] - limits['min_lon'], limits['max_lat'] - limits['min_lat']

    # Calculate the bounding box limits
    min_lon, max_lon = central_lon - (width / 2), central_lon + (width / 2)
    min_lat, max_lat = central_lat - (length / 2), central_lat + (length / 2)

    # Format time for output
    datestr = pd.to_datetime(t).strftime('%Y-%m-%d-%H%M')

    position = {
        'datestr': datestr,
        'central_lat': central_lat,
        'central_lon': central_lon,
        'length': length,
        'width': width
    }

    limits = {
        'min_lon': min_lon,
        'max_lon': max_lon,
        'min_lat': min_lat,
        'max_lat': max_lat
    }

    return position, limits

def compute_and_store_terms(box_obj, terms_dict):
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
    try:
        ec_obj = EnergyContents(box_obj, method='moving')
        terms_dict['Az'].append(ec_obj.calc_az())
        terms_dict['Ae'].append(ec_obj.calc_ae())
        terms_dict['Kz'].append(ec_obj.calc_kz())
        terms_dict['Ke'].append(ec_obj.calc_ke())
    except Exception as e:
        logging.exception(f"Error in computing Energy Contents: {e}")

    # Conversion Terms
    try:
        ct_obj = ConversionTerms(box_obj, method='moving')
        terms_dict['Cz'].append(ct_obj.calc_cz())
        terms_dict['Ca'].append(ct_obj.calc_ca())
        terms_dict['Ck'].append(ct_obj.calc_ck())
        terms_dict['Ce'].append(ct_obj.calc_ce())
    except Exception as e:
        logging.exception(f"Error in computing Conversion Terms: {e}")

    # Boundary Terms
    try:
        bt_obj = BoundaryTerms(box_obj, method='moving')
        terms_dict['BAz'].append(bt_obj.calc_baz())
        terms_dict['BAe'].append(bt_obj.calc_bae())
        terms_dict['BKz'].append(bt_obj.calc_bkz())
        terms_dict['BKe'].append(bt_obj.calc_bke())
        terms_dict['BΦZ'].append(bt_obj.calc_boz())
        terms_dict['BΦE'].append(bt_obj.calc_boe())
    except Exception as e:
        logging.exception(f"Error in computing Boundary Terms: {e}")

    # Generation/Dissipation Terms
    try:
        gdt_obj = GenerationDissipationTerms(box_obj, method='moving')
        terms_dict['Gz'].append(gdt_obj.calc_gz())
        terms_dict['Ge'].append(gdt_obj.calc_ge())
        if not box_obj.args.residuals:
            terms_dict['Dz'].append(gdt_obj.calc_dz())
            terms_dict['De'].append(gdt_obj.calc_de())
    except Exception as e:
        logging.exception(f"Error in computing Generation/Dissipation Terms: {e}")

    return terms_dict

def finalize_results(times, terms_dict, args, results_sub_directory):
    """
    Process and finalize results. 
    This includes creating a DataFrame with results and saving it to a CSV file.
    Also saves the system track as a CSV file for replicability.

    Args:
        times (array-like): Array of times for the computation.
        terms_dict (dict): Dictionary containing computed terms.
        args (argparse.Namespace): Script arguments.
        results_sub_directory (str): Directory to save results.
    """
    # Convert the terms_dict to DataFrame
    df = pd.DataFrame(terms_dict, index=pd.to_datetime(times), orient ='columns',dtype = float)

    # Estimating budget terms (∂X/∂t) using finite differences
    df = calc_budget_diff(df, times, args)
    
    # Computing residuals, if required
    if args.residuals:
        df = calc_residuals(df, args)

    # Constructing output filename
    method = 'track' if args.track else 'choose'
    outfile_name = args.outname if args.outname else ''.join(args.infile.split('/')[-1].split('.nc')) + f'_{method}'
    outfile_path = os.path.join(results_sub_directory, f'{outfile_name}.csv')

    # Saving the DataFrame to a CSV file
    df.to_csv(outfile_path)
    logging.info(f'Results saved to {outfile_path}')

    # Save system position as a csv file for replicability
    out_track = out_track.rename(columns={'datestr':'time','central_lat':'Lat','central_lon':'Lon'})
    output_trackfile =  results_sub_directory+outfile_name+'_track'
    out_track.to_csv(output_trackfile, index=False, sep=";")

    return outfile_path, df

def plot_results(outfile_path, results_sub_directory, args):
    """
    Plot the results using the specified plotting scripts or functions.

    Args:
        outfile_path (str): Path to the file containing the results.
        results_sub_directory (str): Directory to save plots.
        args (argparse.Namespace): Script arguments.
    """
    # Example plotting commands, replace with actual plotting functions or scripts
    plot_flag = ' -r' if args.residuals else ''
    os.system(f"python ../plots/plot_timeseries.py {outfile_path}{plot_flag}")
    os.system(f"python ../plots/plot_vertical.py {results_sub_directory}")
    os.system(f"python ../plots/plot_boxplot.py {results_sub_directory}{plot_flag}")
    os.system(f"python ../plots/plot_LEC.py {outfile_path}")
    os.system(f"python ../plots/plot_LPS.py {outfile_path}")
    os.system(f"python ../plots/plot_track.py {outfile_path}")

def LEC_moving(data: xr.Dataset, variable_list_df: pd.DataFrame, dTdt: xr.Dataset,
               results_sub_directory: str, figures_directory: str,
               args: argparse.Namespace):
    """
    Computes the Lorenz Energy Cycle using a moving (semi-lagrangian) framework.
    
    Args:
        data (xr.Dataset): A Xarray Dataset containing the data to compute the energy cycle.
        variable_list_df (pd.DataFrame): DataFrame with variable mappings.
        dTdt (xr.Dataset): Dataset containing temperature gradient data.
        results_sub_directory (str): Directory for saving results.
        figures_directory (str): Directory for saving figures.
        args (argparse.Namespace): Arguments provided to the script.

    Returns:
        None
    """
    logging.info('Computing energetics using moving framework')

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
        file_path = os.path.join(results_sub_directory, file_name)
        df.to_csv(file_path, index=None) 
        logging.info(f'{file_path} created (but still empty)')

    # Track file handling
    if args.track:
        trackfile = '../inputs/track'
        try:
            track = pd.read_csv(trackfile, parse_dates=[0], delimiter=';', index_col='time')
        except FileNotFoundError:
            logging.error(f"Track file {trackfile} not found.")
            raise

    # Dictionary for saving system position and attributes
    results_keys = ['datestr', 'central_lat', 'central_lon', 'length', 'width',
                    'min_max_zeta_850', 'min_hgt_850', 'max_wind_850']
    out_track = pd.DataFrame(columns=results_keys)

    # Dictionary for storing results
    terms_dict = {term: [] for term in 
                  ['Az', 'Ae', 'Kz', 'Ke', 'Cz', 'Ca', 'Ck', 'Ce', 'Ge', 'Gz', 
                   'BAz', 'BAe', 'BKz', 'BKe', 'BΦZ', 'BΦE', 'Dz', 'De']}

    # Slice the time array to match the track file
    times = pd.to_datetime(data[TimeName].values)
    if args.track:
        times = times[(times >= track.index[0]) & (times <= track.index[-1])]
        if len(times) == 0:
            logging.error("Mismatch between trackfile and data times.")
            raise ValueError("Mismatch between trackfile and data times.")

    # Iterating over times
    for t in times:
        try:
            idata, idTdt = data.sel({TimeName: t}), dTdt.sel({TimeName: t})
            if idata[TimeName].shape != ():
                idata, idTdt = data.isel({TimeName: 1}), dTdt.isel({TimeName: 1})
        except KeyError as e:
            logging.error(f"Time indexing error: {e}")
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

        # Get current time attributes
        position, limits = get_position_and_limits(args, t, track, iu_850, iv_850, izeta_850, ihgt_850, lat, lon)
        out_track = out_track.append(position, ignore_index=True)
        plot_domain_attributes(iu_850, iv_850, izeta_850, ihgt_850, lat, lon, t, figures_directory)

        # Create box object
        try:
            box_obj = BoxData(
                data=idata.compute(),
                variable_list_df=variable_list_df,
                args=args,
                western_limit=limits['min_lon'],
                eastern_limit=limits['max_lon'],
                southern_limit=limits['min_lat'],
                northern_limit=limits['max_lat'],
                output_dir=results_sub_directory,
                dTdt=idTdt
            )
        except Exception as e:
            logging.exception(f"Error creating BoxData object: {e}")
            continue  # Skip to next iteration

        # Compute and store various meteorological terms
        terms_dict = compute_and_store_terms(box_obj, terms_dict)

    # Finalize and process results
    df, outfile_path = finalize_results(times, terms_dict, args, results_sub_directory)

    if args.plots:
        plot_results(outfile_path, results_sub_directory, args)
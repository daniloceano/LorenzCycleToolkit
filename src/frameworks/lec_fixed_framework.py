# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lec_fixed_framework.py                             :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:32:59 by daniloceano       #+#    #+#              #
#    Updated: 2024/05/03 16:11:03 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import pandas as pd
import logging
import argparse
import xarray as xr
from pathlib import Path
from ..analysis.energy_contents import EnergyContents
from ..analysis.conversion_terms import ConversionTerms
from ..analysis.boundary_terms import BoundaryTerms
from ..analysis.generation_and_dissipation_terms import GenerationDissipationTerms
from ..utils.box_data import BoxData
from ..utils.calc_budget_and_residual import calc_budget_diff, calc_residuals

def lec_fixed(data: xr.Dataset, variable_list_df: pd.DataFrame, results_subdirectory: str,
              app_logger: logging.Logger, args: argparse.Namespace):
    """
    Computes the Lorenz Energy Cycle (LEC) using a fixed framework.

    Args:
        data (xr.Dataset): Dataset containing the atmospheric data for LEC computation.
        variable_list_df (pd.DataFrame): DataFrame with variable mappings used in the LEC analysis.
        results_subdirectory (str): Directory path to save the results.
        args (argparse.Namespace): Arguments provided to the script, including configurations 
                                   for the LEC computation.

    Raises:
        ValueError: If the bounding box limits are invalid.
        Exception: General exceptions for processing errors.

    Note:
        This function computes various aspects of the LEC and saves the results in CSV format.
        It also triggers plotting scripts if specified in the arguments.
    """
    logging.info('--- Computing energetics using fixed framework ---')
    
    box_limits_file = args.box_limits
    dfbox = pd.read_csv(box_limits_file, header=None, delimiter=';', index_col=0)
    min_lon, max_lon = dfbox.loc['min_lon'].iloc[0], dfbox.loc['max_lon'].iloc[0]
    min_lat, max_lat = dfbox.loc['min_lat'].iloc[0], dfbox.loc['max_lat'].iloc[0]
    
    if min_lon > max_lon:
        error_message = f"Error in box_limits: min_lon ({min_lon}) is greater than max_lon ({max_lon})"
        app_logger.error(error_message)
        raise ValueError(error_message)

    if min_lat > max_lat:
        error_message = f"Error in box_limits: min_lat ({min_lat}) is greater than max_lat ({max_lat})"
        app_logger.error(error_message)
        raise ValueError(error_message)
    
    app_logger.debug('Loading data into memory..')
    data = data.compute()
    app_logger.debug('Data loaded into memory.')

    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = (
        variable_list_df.loc['Longitude']['Variable'],
        variable_list_df.loc['Latitude']['Variable'],
        variable_list_df.loc['Time']['Variable'],
        variable_list_df.loc['Vertical Level']['Variable']
    )
    
    pres = data[VerticalCoordIndexer] * data[VerticalCoordIndexer].metpy.units
    app_logger.info(f'Bounding box: min_lon={min_lon}, max_lon={max_lon}, min_lat={min_lat}, max_lat={max_lat}')
    
    for term in ['Az', 'Ae', 'Kz', 'Ke', 'Cz', 'Ca', 'Ck', 'Ce', 'Ge', 'Gz']:
        columns = [TimeName] + [float(i)/100 for i in pres.values]
        output_path = Path(results_subdirectory, f'{term}_{VerticalCoordIndexer}.csv')
        pd.DataFrame(columns=columns).to_csv(output_path, index=None)

    try:
        box_obj = BoxData(data, variable_list_df, min_lon, max_lon, min_lat, max_lat, args, results_subdirectory)
    except Exception as e:
        app_logger.exception("An exception occurred while creating BoxData object")
        raise

    try:
        ec_obj = EnergyContents(box_obj, 'fixed', app_logger)
        energy_list = [ec_obj.calc_az(), ec_obj.calc_ae(), ec_obj.calc_kz(), ec_obj.calc_ke()]
    except Exception as e:
        app_logger.exception("An exception occurred while creating EnergyContents object")
        raise
    app_logger.info('Computed energy contents')

    try:
        ct_obj = ConversionTerms(box_obj, 'fixed', app_logger)
        conversion_list = [ct_obj.calc_cz(), ct_obj.calc_ca(), ct_obj.calc_ck(), ct_obj.calc_ce()]
    except Exception as e:
        app_logger.exception("An exception occurred while creating ConversionTerms object")
        raise
    app_logger.info('Computed conversion terms')

    try:
        bt_obj = BoundaryTerms(box_obj, 'fixed', app_logger)
        boundary_list = [bt_obj.calc_baz(), bt_obj.calc_bae(), bt_obj.calc_bkz(), bt_obj.calc_bke(), bt_obj.calc_boz(), bt_obj.calc_boe()]
    except Exception as e:
        app_logger.exception("An exception occurred while creating BoundaryTerms object")
        raise
    app_logger.info('Computed boundary terms')

    try:
        gdt_obj = GenerationDissipationTerms(box_obj, 'fixed', app_logger)
        gen_diss_list = [gdt_obj.calc_gz(), gdt_obj.calc_ge()] if args.residuals else [gdt_obj.calc_gz(), gdt_obj.calc_ge(), gdt_obj.calc_dz(), gdt_obj.calc_de()]
    except Exception as e:
        app_logger.exception("An exception occurred while creating GenerationDissipationTerms object")
        raise
    app_logger.info('Computed generation dissipation terms')

    dates = data[TimeName].values
    df = pd.DataFrame(index=dates.astype('datetime64'))
    for i, col in enumerate(['Az', 'Ae', 'Kz', 'Ke']):
        df[col] = energy_list[i]
    for i, col in enumerate(['Cz', 'Ca', 'Ck', 'Ce']):
        df[col] = conversion_list[i]
    for i, col in enumerate(['BAz', 'BAe', 'BKz', 'BKe', 'Gz', 'Ge', 'Dz', 'De'][:len(gen_diss_list) + 4]):
        df[col] = boundary_list[i] if i < 4 else gen_diss_list[i-4]

    df = calc_budget_diff(df, dates, app_logger)
    df = calc_residuals(df, app_logger)
    app_logger.info('Computed budget and residuals')

    if args.outname:
        results_filename = args.outname
    else:
        infile_name = os.path.basename(args.infile).split('.nc')[0]
        results_filename = ''.join(f'{infile_name}_fixed_results')
    results_file = Path(results_subdirectory, f'{results_filename}.csv')
    df.to_csv(results_file)
    app_logger.info(f'Results saved to {results_file}')

    if args.plots:
        from ..plots.timeseries_terms import plot_timeseries
        from ..plots.map_box_limits import plot_box_limits
        from ..plots.plot_boxplot import boxplot_terms
        from ..plots.plot_LEC import plot_lorenz_cycle
        from ..plots.plot_hovmoller import plot_hovmoller
        
        app_logger.info('Generating plots..')
        figures_directory = os.path.join(results_subdirectory, 'Figures')
        plot_timeseries(results_file, figures_directory, app_logger)
        plot_box_limits(box_limits_file, figures_directory, app_logger)
        boxplot_terms(results_file, results_subdirectory, figures_directory, app_logger)
        plot_lorenz_cycle(results_file, figures_directory, app_logger=app_logger)
        plot_hovmoller(results_file, figures_directory, app_logger)
# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lec_fixed_framework.py                             :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:32:59 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/19 22:56:16 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os 
import pandas as pd
import logging
import argparse
import xarray as xr
from EnergyContents import EnergyContents
from ConversionTerms import ConversionTerms
from BoundaryTerms import BoundaryTerms
from GenerationDissipationTerms import GenerationDissipationTerms
from BoxData import BoxData
from BudgetResidual import calc_budget_diff, calc_residuals

def lec_fixed(data: xr.Dataset, variable_list_df: pd.DataFrame, results_subdirectory: str, args: argparse.Namespace):
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

    logging.info('Computing energetics using fixed framework')
    
    logging.info('Opening data into memory..')
    data = data.compute()
    logging.info('Data loaded into memory.')

    dfbox = pd.read_csv('../inputs/box_limits', header=None, delimiter=';', index_col=0)
    min_lon, max_lon = dfbox.loc['min_lon'].iloc[0], dfbox.loc['max_lon'].iloc[0]
    min_lat, max_lat = dfbox.loc['min_lat'].iloc[0], dfbox.loc['max_lat'].iloc[0]
    
    if min_lon > max_lon:
        raise ValueError('Error in box_limits: min_lon > max_lon')
    if min_lat > max_lat:
        raise ValueError('Error in box_limits: min_lat > max_lat')

    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = (
        variable_list_df.loc['Longitude']['Variable'],
        variable_list_df.loc['Latitude']['Variable'],
        variable_list_df.loc['Time']['Variable'],
        variable_list_df.loc['Vertical Level']['Variable']
    )
    
    pres = data[VerticalCoordIndexer] * data[VerticalCoordIndexer].metpy.units
    
    logging.info('Parameters specified for the bounding box: min_lon={}, max_lon={}, min_lat={}, max_lat={}'.format(min_lon, max_lon, min_lat, max_lat))
    
    # Create CSV files for storing vertical results
    for term in ['Az', 'Ae', 'Kz', 'Ke', 'Cz', 'Ca', 'Ck', 'Ce', 'Ge', 'Gz']:
        columns = [TimeName] + [float(i)/100 for i in pres.values]
        pd.DataFrame(columns=columns).to_csv(os.path.join(results_subdirectory, f'{term}_{VerticalCoordIndexer}.csv'), index=None)
    
    try:
        box_obj = BoxData(data=data, variable_list_df=variable_list_df, args=args, western_limit=min_lon, eastern_limit=max_lon, southern_limit=min_lat, northern_limit=max_lat, output_dir=results_subdirectory)
        ec_obj = EnergyContents(box_obj, method='fixed')
        energy_list = [ec_obj.calc_az(), ec_obj.calc_ae(), ec_obj.calc_kz(), ec_obj.calc_ke()]

        ct_obj = ConversionTerms(box_obj, method='fixed')
        conversion_list = [ct_obj.calc_cz(), ct_obj.calc_ca(), ct_obj.calc_ck(), ct_obj.calc_ce()]

        bt_obj = BoundaryTerms(box_obj, method='fixed')
        boundary_list = [bt_obj.calc_baz(), bt_obj.calc_bae(), bt_obj.calc_bkz(), bt_obj.calc_bke(), bt_obj.calc_boz(), bt_obj.calc_boe()]

        gdt_obj = GenerationDissipationTerms(box_obj, method='fixed')
        gen_diss_list = [gdt_obj.calc_gz(), gdt_obj.calc_ge()] if args.residuals else [gdt_obj.calc_gz(), gdt_obj.calc_ge(), gdt_obj.calc_dz(), gdt_obj.calc_de()]
    except Exception as e:
        logging.exception("An exception occurred: {}".format(e))
        raise

    # Organizing results in a DataFrame
    dates = data[TimeName].values
    df = pd.DataFrame({'Date': dates.astype('datetime64[D]'), 'Hour': pd.to_datetime(dates).hour})
    for i, col in enumerate(['Az', 'Ae', 'Kz', 'Ke']):
        df[col] = energy_list[i]
    for i, col in enumerate(['Cz', 'Ca', 'Ck', 'Ce']):
        df[col] = conversion_list[i]
    for i, col in enumerate(['BAz', 'BAe', 'BKz', 'BKe', 'Gz', 'Ge', 'Dz', 'De'][:len(gen_diss_list) + 4]):
        df[col] = boundary_list[i] if i < 4 else gen_diss_list[i-4]

    df = calc_budget_diff(df, dates, args)
    df = calc_residuals(df, args)

    outfile_name = args.outname if args.outname else ''.join(args.infile.split('/')[-1].split('.nc')) + '_fixed' 
    outfile = os.path.join(results_subdirectory, f'{outfile_name}.csv')
    df.to_csv(outfile)
    logging.info(f'Results saved to {outfile}')

    # Generate plots
    plot_flag = ' -r' if args.residuals else ' '
    if args.plots:
        plot_scripts = ["plot_timeseries.py", "plot_vertical.py", "plot_boxplot.py", "plot_LEC.py", "plot_LPS.py"]
        for script in plot_scripts:
            os.system(f"python ../plots/{script} {outfile}{plot_flag}")
        os.system(f"python ../plots/plot_area.py {min_lon} {max_lon} {min_lat} {max_lat} {results_subdirectory}")

if __name__ == '__main__':
    from tools import initialize_logging, prepare_data

    # Mocking the argparse.Namespace for testing
    args = argparse.Namespace(
        infile='../samples/Reg1-Representative_NCEP-R2.nc',
        residuals=True,
        fixed=True,
        geopotential=False,
        track=False,
        choose=False,
        zeta=False,
        mpas=False,
        plots=False,
        outname=None,
        verbosity=False
    )

    initialize_logging()

    varlist = "../inputs/fvars_NCEP-R2"
    variable_list_df = pd.read_csv(varlist, sep=';', index_col=0, header=0)

    # Assuming create_arg_parser and prepare_data are modified to work with the mocked args
    data, _ = prepare_data(args, varlist)

    results_subdirectory = "../LEC_results"
    os.makedirs(results_subdirectory, exist_ok=True)

    lec_fixed(data, variable_list_df, results_subdirectory, args)
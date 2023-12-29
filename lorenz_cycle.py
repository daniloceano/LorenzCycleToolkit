# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lorenz_cycle.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/04/20 10:05:52 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/29 14:18:26 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main script for the Lorenz Energy Cycle (LEC) Program.

This script handles user inputs to open the data, specifies the bounding box for the calculations,
and validates user inputs for processing the Lorenz Energy Cycle using either a fixed or moving framework.

Created by Danilo Couto de Souza, Universidade de São Paulo (USP), Instituto de Astronomia, Ciências Atmosféricas e Geociências, São Paulo - Brazil
Contact: danilo.oceano@gmail.com
"""

import os
import sys
import time
import logging
import argparse
import warnings
import pandas as pd
from metpy.units import units
from src.utils.tools import prepare_data, initialize_logging
from src.frameworks.lec_fixed_framework import lec_fixed
from src.frameworks.lec_moving_framework import lec_moving

# Suppress specific RuntimeWarning from xarray
warnings.filterwarnings("ignore", category=RuntimeWarning, module='xarray.backends.plugins')

def create_arg_parser():
    """
    Create an argument parser for the Lorenz Energy Cycle (LEC) program.

    Returns:
        argparse.ArgumentParser: The argument parser object.
    """
    parser = argparse.ArgumentParser(description="Lorenz Energy Cycle (LEC) program.")
    parser.add_argument("infile", help="Input .nc file with temperature, geopotential, and wind components.")
    parser.add_argument("-r", "--residuals", action='store_true', help="Compute the Dissipation and Generation terms as residuals.")
    parser.add_argument("-g", "--geopotential", action='store_true', help="Use geopotential data instead of geopotential height.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fixed", action='store_true', help="Compute the energetics for a fixed domain specified by the 'box_lims' file.")
    group.add_argument("-t", "--track", action='store_true', help="Define the box using a track file specified by the 'track' file.")
    group.add_argument("-c", "--choose", action='store_true', help="Choose the domain for each time step by clicking on the screen.")
    parser.add_argument("-z", "--zeta", action='store_true', help="Use this flag if the track file was created using vorticity.")
    parser.add_argument("-m", "--mpas", action='store_true', help="for MPAS-A data processed with MPAS-BR routines")
    parser.add_argument("-p", "--plots", action='store_true', help="wether or not to make plots.")
    parser.add_argument("-o", "--outname", type=str, help="Choose a name for saving results.")
    parser.add_argument("-v", "--verbosity", action='store_true', help="Increase output verbosity.")
    return parser

def setup_results_directory(args, method):
    """
    Create and setup the results directory for storing output files and figures.

    Args:
        args (object): The arguments object containing input file information.
        method (str): The method used for processing the input file.

    Returns:
        tuple: A tuple containing the path to the results directory and the path to the figures directory.
    """
    resuts_directory = "../LEC_Results/"
    results_subdirectory = os.path.join(
        resuts_directory, "".join(args.infile.split('/')[-1].split('.nc')) + '_' + method)
    figures_directory = os.path.join(results_subdirectory, 'Figures')
    os.makedirs(figures_directory, exist_ok=True)
    os.makedirs(results_subdirectory, exist_ok=True)
    os.makedirs(results_subdirectory, exist_ok=True)

    return results_subdirectory, figures_directory

def run_lec_analysis(data, args, results_subdirectory, figures_directory, app_logger):
    """
    Runs the LEc analysis on the given data using the specified method and arguments.

    Args:
        data (pd.DataFrame): The input data for the analysis.
        method (str): The method to be used for the analysis.
        args (dict): Additional arguments for the analysis.

    Returns:
        None

    Raises:
        None
    """
    start_time = time.time()

    variable_list_df = pd.read_csv("inputs/fvars", sep=';', index_col=0, header=0)

    if args.fixed:
        lec_fixed(data, variable_list_df, results_subdirectory, app_logger, args)
        app_logger.info("--- %s seconds running fixed framework ---" % (time.time() - start_time))

    if args.track or args.choose:
        dTdt =  data[variable_list_df.loc['Air Temperature']['Variable']].differentiate(
                variable_list_df.loc['Time']['Variable'],datetime_unit='s') * units('K/s')
        
        lec_moving(data, variable_list_df, dTdt, results_subdirectory, figures_directory, app_logger, args)
        app_logger.info("--- %s seconds running moving framework ---" % (time.time() - start_time))

def main():
    """
    Executes the main function.

    This function is the entry point of the program. It creates an argument parser, parses the command line arguments, initializes logging, prepares the data, and runs the LEC analysis.

    Parameters:
        None.

    Returns:
        None.
    """
    
    # Parse command line arguments
    parser = create_arg_parser()
    args = parser.parse_args()

    # # Example usage for debugging
    # print("WARNING: USING EXAMPLE ARGUMENTS")
    # args = parser.parse_args(['samples/Reg1-Representative_NCEP-R2.nc', '-f', '-r', '-p'])

    # Set method
    if args.fixed:
        method = 'fixed'
    elif args.track:
        method = 'track'
    elif args.choose:
        method = 'choose'

    # Setup results directory
    results_subdirectory, figures_directory = setup_results_directory(args, method)

    # Initialize logging
    app_logger = initialize_logging(results_subdirectory, args.verbosity)
    app_logger.info("Starting LEC analysis")
    app_logger.info(f"Command line arguments: {args}")

    # Prepare data
    data = prepare_data(args)
    
    # Run LEC analysis
    run_lec_analysis(data, args, results_subdirectory, figures_directory, app_logger)

if __name__ == "__main__":
    main()
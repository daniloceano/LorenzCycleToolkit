# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lorenzcycletoolkit.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/04/20 10:05:52 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/14 12:46:56 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main script for the Lorenz Energy Cycle (LEC) Program.

This script handles user inputs to open the data, specifies the bounding box for the calculations,
and validates user inputs for processing the Lorenz Energy Cycle using either a fixed or moving framework.

Created by:

Danilo Couto de Souza
Universidade de São Paulo (USP) Instituto de Astronomia, Ciências Atmosféricas e Geociências,
São Paulo - Brazil

Contact: danilo.oceano@gmail.com
"""

import argparse
import os
import sys
import time
import warnings

import pandas as pd
from metpy.units import units

from src.frameworks.lec_fixed_framework import lec_fixed
from src.frameworks.lec_moving_framework import lec_moving
from src.utils.tools import initialize_logging, prepare_data

# Suppress specific RuntimeWarning from xarray
warnings.filterwarnings(
    "ignore", category=RuntimeWarning, module="xarray.backends.plugins"
)


def create_arg_parser():
    """
    Create an argument parser for the Lorenz Energy Cycle (LEC) program.

    Returns:
        argparse.ArgumentParser: The argument parser object.
    """
    parser = argparse.ArgumentParser(description="Lorenz Energy Cycle (LEC) program.")
    parser.add_argument(
        "infile",
        help="Input .nc file with temperature, geopotential/geopotential height, and wind components data.",
    )
    parser.add_argument(
        "-r",
        "--residuals",
        action="store_true",
        help="Compute the Dissipation and Generation terms as residuals.",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-f",
        "--fixed",
        action="store_true",
        help="Compute the energetics for a fixed domain specified by the 'box_limits' file.",
    )
    group.add_argument(
        "-t",
        "--track",
        action="store_true",
        help="Define the domain using a track file.",
    )
    group.add_argument(
        "-c",
        "--choose",
        action="store_true",
        help="Interactively select the domain for each time step.",
    )
    parser.add_argument(
        "-z",
        "--zeta",
        action="store_true",
        help="Use the vorticity from the track file instead of computing it at 850 hPa.",
    )
    parser.add_argument(
        "-m",
        "--mpas",
        action="store_true",
        help="Specify this flag if working with MPAS-A data processed with MPAS-BR routines.",
    )
    parser.add_argument("-p", "--plots", action="store_true", help="Generate plots.")
    parser.add_argument(
        "-v", "--verbosity", action="store_true", help="Logger level set to debug mode."
    )
    parser.add_argument(
        "--cdsapi",
        action="store_true",
        help="Use CDS API for downloading data (experimental).",
    )
    parser.add_argument(
        "--trackfile",
        type=str,
        default="inputs/track",
        help="Specify a custom track file. Default is 'inputs/track'.",
    )
    parser.add_argument(
        "--box_limits",
        type=str,
        default="inputs/box_limits",
        help="Specify a custom box limits file. Default is 'inputs/box_limits'.",
    )
    parser.add_argument(
        "-o", "--outname", type=str, help="Specify an output name for the results."
    )
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
    resuts_directory = "./LEC_Results/"
    results_subdirectory = os.path.join(
        resuts_directory,
        "".join(args.infile.split("/")[-1].split(".nc")) + "_" + method,
    )
    results_subdirectory_vertical_levels = os.path.join(
        results_subdirectory, "results_vertical_levels")
    figures_directory = os.path.join(results_subdirectory, "Figures")
    os.makedirs(figures_directory, exist_ok=True)
    os.makedirs(results_subdirectory, exist_ok=True)
    os.makedirs(results_subdirectory_vertical_levels, exist_ok=True)

    return results_subdirectory, figures_directory, results_subdirectory_vertical_levels


def run_lec_analysis(data, args, results_subdirectory, figures_directory, results_subdirectory_vertical_levels, app_logger):
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

    variable_list_df = pd.read_csv("inputs/namelist", sep=";", index_col=0, header=0)

    if args.fixed:
        lec_fixed(data, variable_list_df, results_subdirectory, results_subdirectory_vertical_levels, app_logger, args)
        app_logger.info(
            "--- %s seconds running fixed framework ---" % (time.time() - start_time)
        )

    if args.track or args.choose:
        dTdt = data[variable_list_df.loc["Air Temperature"]["Variable"]].differentiate(
            variable_list_df.loc["Time"]["Variable"], datetime_unit="s"
        ) * units("K/s")

        lec_moving(
            data,
            variable_list_df,
            dTdt,
            results_subdirectory,
            figures_directory,
            results_subdirectory_vertical_levels,
            app_logger,
            args,
        )
        app_logger.info(
            "--- %s seconds running moving framework ---" % (time.time() - start_time)
        )


def main():
    """
    Executes the main function.

    This function is the entry point of the program. It creates an argument parser, parses the command line arguments,
    initializes logging, prepares the data, and runs the LEC analysis.

    Parameters:
        None.

    Returns:
        None.
    """

    # Parse command line arguments
    if len(sys.argv) > 1:
        parser = create_arg_parser()
        args = parser.parse_args()

    else:
        # Example usage for debugging
        parser = create_arg_parser()
        print(
            "----------------------------------------------------------------------------"
        )
        print("WARNING: USING EXAMPLE ARGUMENTS")
        import shutil

        shutil.copy("inputs/namelist_NCEP-R2", "inputs/namelist")
        shutil.copy("inputs/track_testdata_NCEP-R2", "inputs/track")
        args = parser.parse_args(
            ["samples/testdata_NCEP-R2.nc", "-t", "-r", "-p", "-v"]
        )
        # args = parser.parse_args(['20070536_ERA5_sliced.nc', '-t', '-r', '-p', '-v', '-g', '--cdsapi', '--trackfile',
        # 'inputs/track_20070536.csv'])
        print(
            "----------------------------------------------------------------------------"
        )

    # Set method
    if args.fixed:
        method = "fixed"
    elif args.track:
        method = "track"
    elif args.choose:
        method = "choose"

    # Setup results directory
    results_subdirectory, figures_directory, results_subdirectory_vertical_levels = setup_results_directory(args, method)

    # Initialize logging
    app_logger = initialize_logging(results_subdirectory, args)
    app_logger.info("Starting LEC analysis")
    app_logger.info(f"Command line arguments: {args}")

    # Prepare data
    data = prepare_data(args, "inputs/namelist", app_logger)

    # Run LEC analysis
    run_lec_analysis(data, args, results_subdirectory, figures_directory, results_subdirectory_vertical_levels, app_logger)


if __name__ == "__main__":
    main()

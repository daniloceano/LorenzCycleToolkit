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
from src.utils.input_files import resolve_input_file
from src.utils.offline_maps import offline_maps
from src.utils.output_names import results_directory, validate_outname
from src.utils.preprocessing import prepare_data
from src.utils.run_status import EXIT_PARTIAL, RunStatus, snapshot_directory
from src.utils.tools import initialize_logging

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
        "--allow-partial-plots",
        action="store_true",
        help=(
            "Exit with status 0 when the computation is complete but some requested "
            f"figures could not be generated (default: exit status {EXIT_PARTIAL}). "
            "run_status.json still reports the run as partial."
        ),
    )
    parser.add_argument(
        "--offline-maps",
        action="store_true",
        help=(
            "Draw maps from local Natural Earth data only (CARTOPY_DATA_DIR or Cartopy's "
            "data directory) and never download them; a map whose data are missing is "
            "reported as failed (partial output)."
        ),
    )
    parser.add_argument(
        "-v", "--verbosity", action="store_true", help="Logger level set to debug mode."
    )
    parser.add_argument(
        "--cdsapi",
        action="store_true",
        help="Use CDS API for downloading data (experimental).",
    )
    parser.add_argument(
        "--time-resolution",
        type=int,
        default=3,
        help="Temporal resolution in hours for CDS API data download (default: 3). Use 3, 6, or 12 to avoid cost limits.",
    )
    parser.add_argument(
        "--ncar",
        action="store_true",
        help=(
            "Download the ERA5 data from the NCAR GDEX copy of ERA5 (dataset d633000, no account) instead of "
            "the CDS API when the input file does not exist (experimental). Same area, 37 pressure levels and "
            "hours as --cdsapi (track file and --time-resolution). A provenance record is written next to the "
            "file, and the file is reused only if it matches that record."
        ),
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
        "-o",
        "--outname",
        type=str,
        help=(
            "Name of the results: LEC_Results/NAME_<method>/ with NAME_<method>_results.csv, "
            "NAME_<method>_trackfile (track) and log.NAME. Default: the input file name."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Reuse a results directory that already contains files. Without it the "
            "program stops (exit status 2) before writing anything. Files are rewritten "
            "in place; files of an earlier run that this run does not produce are kept "
            "and listed in run_status.json (nothing is deleted)."
        ),
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
    results_subdirectory = results_directory(args, method)
    results_subdirectory_vertical_levels = os.path.join(
        results_subdirectory, "results_vertical_levels")
    figures_directory = os.path.join(results_subdirectory, "Figures")
    os.makedirs(figures_directory, exist_ok=True)
    os.makedirs(results_subdirectory, exist_ok=True)
    os.makedirs(results_subdirectory_vertical_levels, exist_ok=True)

    return results_subdirectory, figures_directory, results_subdirectory_vertical_levels


def run_lec_analysis(data, args, results_subdirectory, figures_directory, results_subdirectory_vertical_levels, app_logger, run_status=None):
    """
    Runs the LEc analysis on the given data using the specified method and arguments.

    Args:
        data (pd.DataFrame): The input data for the analysis.
        method (str): The method to be used for the analysis.
        args (dict): Additional arguments for the analysis.
        run_status (RunStatus): Status record of the run; a new one is
            created (and written to run_status.json) if not given.

    Returns:
        RunStatus: The status of the run, also written to run_status.json.

    Raises:
        Any exception of the scientific computation, after recording it in
        run_status.json. Figure failures are recorded, not raised.
    """
    start_time = time.time()

    if run_status is None:
        run_status = RunStatus(args, results_subdirectory, app_logger)
        run_status.write()

    with run_status.guard(), offline_maps(getattr(args, "offline_maps", False)):
        variable_list_df = pd.read_csv(resolve_input_file("inputs/namelist"), sep=";", index_col=0, header=0)

        if args.fixed:
            framework = "Fixed"
            lec_fixed(data, variable_list_df, results_subdirectory, results_subdirectory_vertical_levels, app_logger, args,
                      run_status=run_status)

        if args.track or args.choose:
            framework = "Moving"
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
                run_status=run_status,
            )

        run_status.finish()

    message = run_status.final_message(framework, time.time() - start_time)
    if run_status.status == "success":
        app_logger.info(message)
    else:
        app_logger.warning(message)
    return run_status


def main():
    """
    Executes the main function.

    This function is the entry point of the program. It creates an argument parser, parses the command line arguments,
    initializes logging, prepares the data, and runs the LEC analysis.

    Every run writes run_status.json in its results directory (first as
    "running", then with the final state).

    Parameters:
        None.

    Returns:
        int: Exit status. 0 when every requested product was generated;
        3 when the computation is complete but a requested figure failed
        (0 with --allow-partial-plots). Errors of the computation propagate
        (exit status 1); argument errors exit with status 2.
    """

    # Parse command line arguments. With no arguments argparse prints the
    # usage and exits with status 2; nothing is downloaded or written.
    parser = create_arg_parser()
    args = parser.parse_args()

    if args.cdsapi and args.ncar:
        parser.error("--cdsapi and --ncar are alternative data sources; use one of them")

    # Set method
    if args.fixed:
        method = "fixed"
    elif args.track:
        method = "track"
    elif args.choose:
        method = "choose"

    if args.outname is not None:
        try:
            validate_outname(args.outname)
        except ValueError as exc:
            parser.error(f"argument -o/--outname: {exc}")

    # Never overwrite the results of an earlier run unless asked to: stop
    # before anything is written (no directory, log or status file changes).
    previous_files = snapshot_directory(results_directory(args, method))
    if previous_files and not args.overwrite:
        listed = sorted(previous_files)
        shown = ", ".join(listed[:5]) + (f", ... ({len(listed)} files)" if len(listed) > 5 else "")
        parser.error(
            f"the results directory {results_directory(args, method)} already contains "
            f"files of an earlier run ({shown}). Use --overwrite to reuse it, or "
            f"-o/--outname to write elsewhere."
        )

    # Setup results directory
    results_subdirectory, figures_directory, results_subdirectory_vertical_levels = setup_results_directory(args, method)

    # Initialize logging
    app_logger = initialize_logging(results_subdirectory, args)
    app_logger.info("🚀 Starting LEC analysis")
    app_logger.info(f"⚙️ Command line arguments: {args}")

    # Status of the run, written before any computation so that an
    # interrupted run never leaves a complete-looking status behind
    run_status = RunStatus(args, results_subdirectory, app_logger,
                           previous_files=previous_files if args.overwrite else None)
    if previous_files:
        app_logger.warning(
            f"⚠️ --overwrite: reusing {results_subdirectory}, which already holds "
            f"{len(previous_files)} files; they are rewritten in place and none is deleted. "
            f"Files this run does not produce are listed in run_status.json (previous_results)."
        )
    run_status.write()

    # Prepare data
    with run_status.guard():
        data = prepare_data(args, resolve_input_file("inputs/namelist"), app_logger)

    # Run LEC analysis
    run_lec_analysis(data, args, results_subdirectory, figures_directory, results_subdirectory_vertical_levels, app_logger,
                     run_status=run_status)

    return run_status.exit_code


if __name__ == "__main__":
    sys.exit(main())

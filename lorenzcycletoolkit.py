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
from src.utils.input_cleanup import cleanup_input, downloaded_input
from src.utils.input_files import resolve_input_file
from src.utils.offline_maps import offline_maps
from src.utils.output_names import results_directory, validate_outname
from src.utils.preprocessing import namelist_for_run, prepare_data
from src.utils.run_status import EXIT_PARTIAL, RunStatus, snapshot_directory
from src.utils.run_timing import RunTiming
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
        help="Download the ERA5 data with the CDS API when the input file does not exist (37 pressure levels).",
    )
    parser.add_argument(
        "--delete-downloaded-input",
        action="store_true",
        help=(
            "With --cdsapi, --ncar, --wekeo, --arco or --destine: delete the input NetCDF file after the run, only if this run downloaded it, "
            "the run succeeded and the file is unchanged. Without this flag no input file is ever deleted."
        ),
    )
    parser.add_argument(
        "--time-resolution",
        type=int,
        default=3,
        help=("Temporal resolution, in whole hours, of the data downloaded by any data source (--cdsapi, --ncar, "
              "--wekeo, --arco, --destine; default: 3): the same hours of each day for all of them. Restrictions: "
              "--arco only multiples of 6 (6-hourly product: 6, 12, 18, 24), refused otherwise before anything is "
              "read; --cdsapi: 3, 6 or 12 keep the request within the CDS cost limits; --ncar, --wekeo and "
              "--destine (hourly products): any whole number of hours. Not used without a data source."),
    )
    parser.add_argument(
        "--ncar",
        action="store_true",
        help=(
            "Download the ERA5 data from the NCAR GDEX copy of ERA5 (dataset d633000) instead of "
            "the CDS API when the input file does not exist. Same area, 37 pressure levels and "
            "hours as --cdsapi (track file and --time-resolution). The access method is chosen with "
            "--ncar-method (default: kerchunk). A provenance record is written next to the "
            "file, and the file is reused only if it matches that record."
        ),
    )
    parser.add_argument(
        "--wekeo",
        action="store_true",
        help=(
            "Download the ERA5 data from WEkEO (Harmonised Data Access, dataset "
            "EO:ECMWF:DAT:REANALYSIS_ERA5_PRESSURE_LEVELS) instead of the CDS API when the input file does not "
            "exist: same area, 37 pressure levels, hours and file contents as --cdsapi, one request per month. "
            "Needs the optional hda package (pip install \"lorenzcycletoolkit[wekeo]\") and a WEkEO account in "
            "~/.hdarc (or HDA_USER and HDA_PASSWORD). Volume limit: --max-gib (default 20 GiB). Orders are "
            "recorded in <file>.wekeo-orders.json before submission and resumed by a new run; a provenance "
            "record is written next to the file."
        ),
    )
    parser.add_argument(
        "--arco",
        action="store_true",
        help=(
            "Read the ERA5 data from the CDS ARCO Zarr store (Analysis Ready, Cloud Optimised) instead of the CDS "
            "API when the input file does not exist: same area, hours and file contents as --cdsapi, BUT only "
            "13 pressure levels (1000 to 50 hPa), every 6 h (--time-resolution must be a multiple of 6), and a "
            "transfer several times larger than the data kept (whole chunks). Terms of the cycle may differ from "
            "routes with 37 levels (in the LZV case, Ca -5 %% on average). Needs the optional packages "
            "(pip install \"lorenzcycletoolkit[arco]\") and the CDS API key of --cdsapi (CDSAPI_KEY, or "
            "CDSAPI_RC / ~/.cdsapirc). Volume limit: --max-gib (default 20 GiB), counted on the bytes "
            "transferred. A provenance record is written next to the file."
        ),
    )
    parser.add_argument(
        "--destine",
        action="store_true",
        help=(
            "Read the ERA5 data from the DestinE Earth Data Hub Zarr store (the Zarr v3 store, updated monthly; "
            "see --destine-store) instead of the CDS API when the input file does not exist: same area, hours "
            "and file contents as --cdsapi, BUT the product's own 19 pressure levels (1000 to 1 hPa, not "
            "interpolated) and values stored with reduced precision (10 bits of mantissa): in the 2b comparison, "
            "18 of 240 terms of the cycle differed between 1 %% and 2.7 %% from the CDS (boundary, generation and "
            "residual terms); energy contents and conversions did not change. Every chunk is one HTTP request against the monthly quota of the account "
            "(500 000): the requests are estimated and limited before any read (--max-requests, default 10000), "
            "and the transfer is many times the data kept, whole chunks of 60 days x 15 x 15 degrees in the v3 "
            "store (--max-gib, default 20 GiB, is usually too small: the refusal says the value needed). Needs "
            "the optional packages (pip install \"lorenzcycletoolkit[destine]\") and the Earth "
            "Data Hub API key in ~/.netrc (machine data.earthdatahub.destine.eu, password <key>) or "
            "DESTINE_EDH_TOKEN. A provenance record is written next to the file."
        ),
    )
    parser.add_argument(
        "--destine-store",
        choices=("v3", "v2"),
        default=None,
        help=("With --destine: Earth Data Hub store (default: v3). v3: Zarr v3 store, updated monthly, chunks of "
              "60 days x 1 level x 15 x 15 degrees. v2: Zarr v2 store of the 2026-11 "
              "benchmark, no longer updated since 2026-07-06 and possibly withdrawn without notice, chunks of "
              "5 days x 1 level x 64 x 64 degrees; fewer requests and bytes for short cases. Same values in both "
              "(identical bit for bit in the 2026-11 test)."),
    )
    parser.add_argument(
        "--max-requests",
        type=int,
        default=None,
        help=("With --destine: largest number of HTTP requests accepted for the download (default: 10000). The "
              "estimate is shown and checked before any chunk is read and, if it passes the limit, the run stops "
              "saying the --max-requests needed; every attempt counts."),
    )
    parser.add_argument(
        "--ncar-method",
        choices=("kerchunk", "subset", "ncss"),
        default=None,
        help=(
            "With --ncar: access method (default: kerchunk). kerchunk: byte ranges of the GDEX files, no "
            "account; its reference files are kept in --kerchunk-refs and downloaded once when missing. "
            "subset: GDEX Get a Subset, needs a GDEX token (GDEX_API_TOKEN or ~/.gdex_api_token). "
            "ncss (experimental): NetCDF Subset Service; no success in the 2026-11 benchmark."
        ),
    )
    parser.add_argument(
        "--kerchunk-refs",
        type=str,
        default=None,
        help=("With --ncar and the kerchunk method: directory of the GDEX kerchunk reference files "
              "(default: inputs/kerchunk_refs)."),
    )
    parser.add_argument(
        "--max-gib",
        type=float,
        default=None,
        help=("With --ncar, --wekeo, --arco or --destine: largest volume accepted for the download, in GiB (default: 20 GiB), for any "
              "--ncar-method. The estimate is checked before anything is downloaded and, if it passes the "
              "limit, the run stops saying the --max-gib needed; bytes of failed attempts and retries count "
              "too. Same meaning as in python -m src.utils.ncar_gdex."),
    )
    guard = parser.add_argument_group(
        "guarded CDS mode (opt-in, for controlled runs such as D-Bench; see src/utils/cds_guard.py)")
    guard.add_argument(
        "--cds-max-seconds", type=float, default=None,
        help=("With --cdsapi: global deadline in seconds of the acquisition (queue, transfers and file assembly); "
              "enables the guarded mode, which also needs --cds-ledger, an explicit credential (CDSAPI_URL and "
              "CDSAPI_KEY, or CDSAPI_RC) and a new destination file. Stopping the client never cancels a remote request."),
    )
    guard.add_argument("--cds-ledger", type=str, default=None,
                       help="Private ledger (JSON lines, mode 600) written before each CDS submission.")
    guard.add_argument("--cds-retries", type=int, default=None,
                       help="HTTP tries of the CDS client in guarded mode (default 3).")
    guard.add_argument("--cds-timeout", type=float, default=None,
                       help="Network timeout in seconds of the CDS client in guarded mode (default 120).")
    guard.add_argument("--cds-case-id", type=str, default=None,
                       help="Case identifier in the ledger (default: the input file name).")
    guard.add_argument("--cds-delete-on-stop", action="store_true", default=None,
                       help=("After a deadline or an interruption, send a best-effort delete of the job just "
                             "submitted (that job only) and record whether the API confirmed it; new submissions "
                             "stay blocked either way."))
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


def run_lec_analysis(data, args, results_subdirectory, figures_directory, results_subdirectory_vertical_levels, app_logger, run_status=None,
                     downloaded=None, namelist=None):
    """
    Runs the LEc analysis on the given data using the specified method and arguments.

    Args:
        data (pd.DataFrame): The input data for the analysis.
        method (str): The method to be used for the analysis.
        args (dict): Additional arguments for the analysis.
        run_status (RunStatus): Status record of the run; a new one is
            created (and written to run_status.json) if not given.
        downloaded (tuple): Fingerprint of the input file when this run
            downloaded it (see src.utils.input_cleanup), else None; only
            such a file can be deleted, with --delete-downloaded-input.
        namelist (str): Namelist with which the data were prepared (see
            src.utils.preprocessing.namelist_for_run); inputs/namelist when
            not given.

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
        variable_list_df = pd.read_csv(namelist or resolve_input_file("inputs/namelist"), sep=";", index_col=0, header=0)

        # Timing of the framework only (the namelist read above is preparation)
        with run_status.timing.phase("lec"):
            if args.fixed:
                framework = "Fixed"
                lec_fixed(data, variable_list_df, results_subdirectory, results_subdirectory_vertical_levels, app_logger,
                          args, run_status=run_status)

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

    # input files are kept unless --delete-downloaded-input allows deleting a file this run downloaded
    cleanup_input(args, run_status, downloaded, data=data, logger=app_logger)
    # timing first: the closing message stays the last line of the log
    app_logger.info(run_status.timing.summary())
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
    timing = RunTiming()  # monotonic start of the whole command (total_seconds)
    parser = create_arg_parser()
    args = parser.parse_args()

    if sum(bool(x) for x in (args.cdsapi, args.ncar, args.wekeo, args.arco, args.destine)) > 1:
        parser.error("--cdsapi, --ncar, --wekeo, --arco and --destine are alternative data sources; use one of them")
    if args.arco:
        from src.utils.arco_zarr import ArcoRequestError, check_time_resolution
        try:
            check_time_resolution(args.time_resolution)  # before anything is read
        except ArcoRequestError as error:
            parser.error(str(error))
    if args.ncar_method is not None and not args.ncar:
        parser.error("--ncar-method needs --ncar")
    if args.max_gib is not None and not (args.ncar or args.wekeo or args.arco or args.destine):
        parser.error("--max-gib needs --ncar or --wekeo (or --arco, --destine)")
    if args.max_requests is not None and not args.destine:
        parser.error("--max-requests needs --destine")
    if args.destine_store is not None and not args.destine:
        parser.error("--destine-store needs --destine")
    if args.max_requests is not None and not args.max_requests > 0:
        parser.error("--max-requests must be greater than 0")
    if args.max_gib is not None and not args.max_gib > 0:
        parser.error("--max-gib must be greater than 0")
    if args.kerchunk_refs is not None and not (args.ncar and (args.ncar_method or "kerchunk") == "kerchunk"):
        parser.error("--kerchunk-refs needs --ncar with the kerchunk method")
    if args.delete_downloaded_input and not (args.cdsapi or args.ncar or args.wekeo or args.arco or args.destine):
        parser.error("--delete-downloaded-input needs --cdsapi, --ncar, --wekeo or --arco (or --destine) "
                     "(local input files are never deleted)")
    guard_options = [name for name in ("cds_ledger", "cds_retries", "cds_timeout", "cds_case_id", "cds_delete_on_stop")
                     if getattr(args, name) is not None]
    if args.cds_max_seconds is not None or guard_options:
        if not args.cdsapi:
            parser.error("the guarded CDS options (--cds-*) need --cdsapi")
        if args.cds_max_seconds is None:
            parser.error("--cds-ledger, --cds-retries, --cds-timeout, --cds-case-id and --cds-delete-on-stop "
                         "need --cds-max-seconds")
        if args.cds_max_seconds <= 0 or not args.cds_ledger:
            parser.error("--cds-max-seconds must be positive and needs --cds-ledger")
        if args.cds_retries is not None and args.cds_retries < 1:
            parser.error("--cds-retries must be at least 1")

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
                           previous_files=previous_files if args.overwrite else None, timing=timing)
    if previous_files:
        app_logger.warning(
            f"⚠️ --overwrite: reusing {results_subdirectory}, which already holds "
            f"{len(previous_files)} files; they are rewritten in place and none is deleted. "
            f"Files this run does not produce are listed in run_status.json (previous_results)."
        )
    run_status.write()

    # Prepare data
    existed_before = os.path.exists(args.infile)
    namelist = namelist_for_run(args)  # the same file prepares the data and feeds the LEC
    with run_status.guard(), timing.phase("preparation"):
        data = prepare_data(args, namelist, app_logger, timing=timing)

    # Run LEC analysis
    run_lec_analysis(data, args, results_subdirectory, figures_directory, results_subdirectory_vertical_levels, app_logger,
                     run_status=run_status, downloaded=downloaded_input(args, existed_before, timing),
                     namelist=namelist)

    return run_status.exit_code


if __name__ == "__main__":
    sys.exit(main())

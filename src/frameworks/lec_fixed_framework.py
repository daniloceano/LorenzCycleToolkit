# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lec_fixed_framework.py                             :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:32:59 by daniloceano       #+#    #+#              #
#    Updated: 2026/02/26 10:53:21 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import argparse
import logging
import os
from pathlib import Path

import pandas as pd
import xarray as xr

from ..analysis.boundary_terms import BoundaryTerms
from ..analysis.conversion_terms import ConversionTerms
from ..analysis.energy_contents import EnergyContents
from ..analysis.generation_and_dissipation_terms import \
    GenerationDissipationTerms
from ..utils.box_data import BoxData
from ..utils.calc_budget_and_residual import calc_budget_diff, calc_residuals


def lec_fixed(
    data: xr.Dataset,
    variable_list_df: pd.DataFrame,
    results_subdirectory: str,
    results_subdirectory_vertical_levels: str,
    app_logger: logging.Logger,
    args: argparse.Namespace,
):
    """
    Computes the Lorenz Energy Cycle (LEC) using a fixed framework.

    Args:
        data (xr.Dataset): Dataset containing the atmospheric data for LEC computation.
        variable_list_df (pd.DataFrame): DataFrame with variable mappings used in the LEC analysis.
        results_subdirectory (str): Directory path to save the results.
        results_subdirectory_vertical_levels (str): Directory path to save the vertical level results.
        args (argparse.Namespace): Arguments provided to the script, including configurations
                                   for the LEC computation.

    Raises:
        ValueError: If the bounding box limits are invalid.
        Exception: General exceptions for processing errors.

    Note:
        This function computes various aspects of the LEC and saves the results in CSV format.
        It also triggers plotting scripts if specified in the arguments.
    """
    logging.info("📊 Computing energetics using fixed framework...")

    box_limits_file = args.box_limits
    
    # Read and validate box_limits file
    try:
        dfbox = pd.read_csv(box_limits_file, header=None, delimiter=";", index_col=0)
    except FileNotFoundError:
        app_logger.error("❌ Box limits file not found!")
        app_logger.error("\n" + "="*70)
        app_logger.error("📁 BOX LIMITS FILE NOT FOUND")
        app_logger.error("="*70)
        app_logger.error(f"Looking for: {os.path.abspath(box_limits_file)}")
        app_logger.error(f"Current directory: {os.getcwd()}")
        app_logger.error("\n💡 User Solutions:")
        app_logger.error("   1. Create a box_limits file with the domain boundaries")
        app_logger.error("   2. Use the default: inputs/box_limits")
        app_logger.error("   3. Specify custom file with: --box_limits <path>")
        app_logger.error("\n📝 Expected format:")
        app_logger.error("   min_lon;-60")
        app_logger.error("   max_lon;-30")
        app_logger.error("   min_lat;-50")
        app_logger.error("   max_lat;-20")
        app_logger.error("\n🔧 Developer Info:")
        app_logger.error(f"   args.box_limits: {box_limits_file}")
        app_logger.error("="*70 + "\n")
        raise FileNotFoundError(
            f"Box limits file not found: {os.path.abspath(box_limits_file)}. "
            f"Create one or use --box_limits to specify path."
        )
    except pd.errors.EmptyDataError:
        app_logger.error("❌ Box limits file is empty!")
        app_logger.error(f"File: {os.path.abspath(box_limits_file)}")
        raise pd.errors.EmptyDataError(f"Box limits file is empty: {box_limits_file}")
    except Exception as e:
        app_logger.error(f"❌ Error reading box_limits file!")
        app_logger.error(f"File: {os.path.abspath(box_limits_file)}")
        app_logger.error(f"Error: {type(e).__name__}: {e}")
        app_logger.error("\n💡 Check file format (should be CSV with ';' delimiter)")
        raise
    
    # Validate required fields exist
    required_fields = ["min_lon", "max_lon", "min_lat", "max_lat"]
    missing_fields = [field for field in required_fields if field not in dfbox.index]
    
    if missing_fields:
        app_logger.error("❌ Box limits file is missing required fields!")
        app_logger.error("\n" + "="*70)
        app_logger.error("📋 MISSING BOX LIMITS FIELDS")
        app_logger.error("="*70)
        app_logger.error(f"File: {box_limits_file}")
        app_logger.error(f"Missing fields: {missing_fields}")
        app_logger.error(f"Found fields: {list(dfbox.index)}")
        app_logger.error("\n📝 Required format:")
        app_logger.error("   min_lon;<value>")
        app_logger.error("   max_lon;<value>")
        app_logger.error("   min_lat;<value>")
        app_logger.error("   max_lat;<value>")
        app_logger.error("="*70 + "\n")
        raise ValueError(
            f"Box limits file missing required fields: {missing_fields}. "
            f"Found: {list(dfbox.index)}"
        )
    
    min_lon, max_lon = dfbox.loc["min_lon"].iloc[0], dfbox.loc["max_lon"].iloc[0]
    min_lat, max_lat = dfbox.loc["min_lat"].iloc[0], dfbox.loc["max_lat"].iloc[0]

    if min_lon > max_lon:
        app_logger.error("❌ Invalid box limits: min_lon > max_lon!")
        app_logger.error("\n" + "="*70)
        app_logger.error("🗺️  INVALID BOX LIMITS - LONGITUDE")
        app_logger.error("="*70)
        app_logger.error(f"min_lon: {min_lon}")
        app_logger.error(f"max_lon: {max_lon}")
        app_logger.error("\n💡 Solution:")
        app_logger.error("   Ensure min_lon < max_lon in your box_limits file")
        app_logger.error(f"   File: {box_limits_file}")
        app_logger.error("="*70 + "\n")
        raise ValueError(
            f"Invalid box_limits: min_lon ({min_lon}) > max_lon ({max_lon}). "
            f"Check {box_limits_file}"
        )

    if min_lat > max_lat:
        app_logger.error("❌ Invalid box limits: min_lat > max_lat!")
        app_logger.error("\n" + "="*70)
        app_logger.error("🗺️  INVALID BOX LIMITS - LATITUDE")
        app_logger.error("="*70)
        app_logger.error(f"min_lat: {min_lat}")
        app_logger.error(f"max_lat: {max_lat}")
        app_logger.error("\n💡 Solution:")
        app_logger.error("   Ensure min_lat < max_lat in your box_limits file")
        app_logger.error(f"   File: {box_limits_file}")
        app_logger.error("="*70 + "\n")
        raise ValueError(
            f"Invalid box_limits: min_lat ({min_lat}) > max_lat ({max_lat}). "
            f"Check {box_limits_file}"
        )

    app_logger.debug("💾 Loading data into memory..")
    data = data.compute()
    app_logger.debug("✅ Data loaded into memory.")

    _, _, TimeName, VerticalCoordIndexer = (
        variable_list_df.loc["Longitude"]["Variable"],
        variable_list_df.loc["Latitude"]["Variable"],
        variable_list_df.loc["Time"]["Variable"],
        variable_list_df.loc["Vertical Level"]["Variable"],
    )

    PressureData = data[VerticalCoordIndexer] * data[VerticalCoordIndexer].metpy.units
    app_logger.info(
        f"🗺️ Bounding box: lon=[{min_lon}, {max_lon}], lat=[{min_lat}, {max_lat}]"
    )

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
        output_path = Path(results_subdirectory_vertical_levels, f"{term}_{VerticalCoordIndexer}.csv")
        pd.DataFrame(columns=columns).to_csv(output_path, index=None)

    try:
        box_obj = BoxData(
            data,
            variable_list_df,
            min_lon,
            max_lon,
            min_lat,
            max_lat,
            args,
            results_subdirectory,
            results_subdirectory_vertical_levels,
        )
    except Exception:
        app_logger.exception("❌ An exception occurred while creating BoxData object")
        raise

    try:
        ec_obj = EnergyContents(box_obj, "fixed", app_logger)
        energy_list = [
            ec_obj.calc_az(),
            ec_obj.calc_ae(),
            ec_obj.calc_kz(),
            ec_obj.calc_ke(),
        ]
    except Exception:
        app_logger.exception(
            "❌ An exception occurred while computing EnergyContents"
        )
        raise
    app_logger.info("⚡ Computed energy contents (Az, Ae, Kz, Ke)")

    try:
        ct_obj = ConversionTerms(box_obj, "fixed", app_logger)
        conversion_list = [
            ct_obj.calc_cz(),
            ct_obj.calc_ca(),
            ct_obj.calc_ck(),
            ct_obj.calc_ce(),
        ]
    except Exception:
        app_logger.exception(
            "❌ An exception occurred while computing ConversionTerms"
        )
        raise
    app_logger.info("🔄 Computed conversion terms (Cz, Ca, Ck, Ce)")

    try:
        bt_obj = BoundaryTerms(box_obj, "fixed", app_logger)
        boundary_list = [
            bt_obj.calc_baz(),
            bt_obj.calc_bae(),
            bt_obj.calc_bkz(),
            bt_obj.calc_bke(),
            bt_obj.calc_boz(),
            bt_obj.calc_boe(),
        ]
    except Exception:
        app_logger.exception(
            "❌ An exception occurred while computing BoundaryTerms"
        )
        raise
    app_logger.info("🏁 Computed boundary terms (BAz, BAe, BKz, BKe, BΦZ, BΦE)")

    try:
        gdt_obj = GenerationDissipationTerms(box_obj, "fixed", app_logger)
        gen_diss_list = (
            [gdt_obj.calc_gz(), gdt_obj.calc_ge()]
            if args.residuals
            else [
                gdt_obj.calc_gz(),
                gdt_obj.calc_ge(),
                gdt_obj.calc_dz(),
                gdt_obj.calc_de(),
            ]
        )
    except Exception:
        app_logger.exception(
            "❌ An exception occurred while computing GenerationDissipationTerms"
        )
        raise
    app_logger.info("🔥 Computed generation/dissipation terms (Gz, Ge, Dz, De)")

    dates = data[TimeName].values
    df = pd.DataFrame(index=dates.astype("datetime64"))
    for i, col in enumerate(["Az", "Ae", "Kz", "Ke"]):
        df[col] = energy_list[i]
    for i, col in enumerate(["Cz", "Ca", "Ck", "Ce"]):
        df[col] = conversion_list[i]
    for i, col in enumerate(
        ["BAz", "BAe", "BKz", "BKe", "Gz", "Ge", "Dz", "De"][: len(gen_diss_list) + 4]
    ):
        df[col] = boundary_list[i] if i < 4 else gen_diss_list[i - 4]

    df = calc_budget_diff(df, dates, app_logger)
    df = calc_residuals(df, app_logger)
    app_logger.info("📈 Computed budget and residuals")

    if args.outname:
        results_filename = args.outname
    else:
        infile_name = os.path.basename(args.infile).split(".nc")[0]
        results_filename = "".join(f"{infile_name}_fixed_results")
    results_file = Path(results_subdirectory, f"{results_filename}.csv")
    df.to_csv(results_file)
    app_logger.info(f"💾 Results saved to {results_file}")

    if args.plots:
        from ..plots.map_box_limits import plot_box_limits
        from ..plots.plot_boxplot import boxplot_terms
        from ..plots.plot_hovmoller import plot_hovmoller
        from ..plots.plot_LEC import plot_lorenzcycletoolkit
        from ..plots.timeseries_terms import plot_timeseries

        app_logger.info("🎨 Generating plots...")
        figures_directory = os.path.join(results_subdirectory, "Figures")

        # Plot time series
        try:
            plot_timeseries(results_file, figures_directory, app_logger)
            app_logger.info("  ✅ Time series plot generated")
        except Exception as e:
            app_logger.error(f"  ❌ Error generating time series plot: {e}")

        # Plot box limits
        try:
            plot_box_limits(box_limits_file, figures_directory, app_logger)
            app_logger.info("  ✅ Box limits plot generated")
        except Exception as e:
            app_logger.error(f"  ❌ Error generating box limits plot: {e}")

        # Plot boxplot terms
        try:
            boxplot_terms(results_file, results_subdirectory, figures_directory, app_logger)
            app_logger.info("  ✅ Boxplot terms generated")
        except Exception as e:
            app_logger.error(f"  ❌ Error generating boxplot terms: {e}")

        # Plot Hovmöller diagram
        try:
            plot_hovmoller(results_file, figures_directory, app_logger)
            app_logger.info("  ✅ Hovmöller plot generated")
        except Exception as e:
            app_logger.error(f"  ❌ Error generating Hovmöller plot: {e}")

        # Plot Lorenz cycle toolkit
        try:
            plot_lorenzcycletoolkit(results_file, figures_directory, app_logger=app_logger)
            app_logger.info("  ✅ Lorenz cycle plot generated")
        except Exception as e:
            app_logger.error(f"  ❌ Error generating Lorenz cycle plot: {e}")

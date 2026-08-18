# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lec_fixed_framework.py                             :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:32:59 by daniloceano       #+#    #+#              #
#    Updated: 2026/01/25 17:40:26 by daniloceano      ###   ########.fr        #
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
from ..analysis.mass_continuity import MassContinuity
from ..utils.box_data import BoxData
from ..utils.calc_budget_and_residual import calc_budget_diff, calc_residuals
from ..utils.longitude import canonicalize_box


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
    if not os.path.exists(box_limits_file) and os.path.exists(
        f"{box_limits_file}.default"
    ):
        box_limits_file = f"{box_limits_file}.default"
    app_logger.info(f"📐 Fixed-domain limits read from: {box_limits_file}")
    dfbox = pd.read_csv(box_limits_file, header=None, delimiter=";", index_col=0)
    min_lon, max_lon = float(dfbox.loc["min_lon"].iloc[0]), float(
        dfbox.loc["max_lon"].iloc[0]
    )
    min_lat, max_lat = float(dfbox.loc["min_lat"].iloc[0]), float(
        dfbox.loc["max_lat"].iloc[0]
    )

    if min_lat > max_lat:
        error_message = f"❌ Error in box_limits: min_lat ({min_lat}) is greater than max_lat ({max_lat})"
        app_logger.error(error_message)
        raise ValueError(error_message)

    # express the requested longitudes in the dataset's own
    # convention. A wrapped interval raises an explicit, actionable error
    # rather than being rejected as a malformed input or, worse, silently
    # producing a partial slice.
    LonIndexerForBox = variable_list_df.loc["Longitude"]["Variable"]
    try:
        min_lon, max_lon = canonicalize_box(
            min_lon,
            max_lon,
            data[LonIndexerForBox].values,
            app_logger=app_logger,
            context=f"box_limits ({box_limits_file})",
        )
    except ValueError as exc:
        app_logger.error(f"❌ {exc}")
        raise

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

    # The per-level CSV headers must list the levels ACTUALLY used, which are
    # only known after BoxData has fixed the vertical control volume
    #. Creating them earlier from the full level set silently
    # misaligned every archived profile whenever a level was excluded.
    used_levels = [float(i) for i in box_obj.PressureData.metpy.dequantify().values]
    app_logger.info(
        f"🧾 Vertical-level CSVs will carry {len(used_levels)} levels: "
        f"{used_levels[0]:.0f} Pa to {used_levels[-1]:.0f} Pa"
    )
    for term in [
        "Az", "Ae", "Kz", "Ke", "Ge", "Gz",
        "Cz", "Cz_1", "Cz_2", "Ca", "Ca_1", "Ca_2",
        "Ce", "Ce_1", "Ce_2", "C_sobreposicao", "M",
        "Ck", "Ck_1", "Ck_2", "Ck_3", "Ck_4", "Ck_5",
    ]:
        columns = [TimeName] + used_levels
        output_path = Path(
            results_subdirectory_vertical_levels, f"{term}_{VerticalCoordIndexer}.csv"
        )
        pd.DataFrame(columns=columns).to_csv(output_path, index=None)

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
            ct_obj.calc_c_sobreposicao(),
        ]
    except Exception:
        app_logger.exception(
            "❌ An exception occurred while computing ConversionTerms"
        )
        raise
    app_logger.info(
        "🔄 Computed conversion terms (Cz, Ca, Ck, Ce, C_sobreposicao)"
    )

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
        mass_residual = MassContinuity(
            box_obj, "fixed", app_logger
        ).calc_mass_residual()
    except Exception:
        app_logger.exception(
            "❌ An exception occurred while computing mass continuity"
        )
        raise
    app_logger.info("⚖️ Computed mass-continuity residual M")

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
    for i, col in enumerate(["Cz", "Ca", "Ck", "Ce", "C_sobreposicao"]):
        df[col] = conversion_list[i]
    # all six boundary diagnostics are computed, so all six are
    # exported. Previously BΦZ and BΦE were computed and then discarded, which
    # also made the fixed and moving frameworks emit different column sets.
    # They are exported here now that their formulations have been re-derived
    # They do NOT enter the residuals.
    for i, col in enumerate(["BAz", "BAe", "BKz", "BKe", "BΦZ", "BΦE"]):
        df[col] = boundary_list[i]
    df["M"] = mass_residual
    for i, col in enumerate(["Gz", "Ge", "Dz", "De"][: len(gen_diss_list)]):
        df[col] = gen_diss_list[i]

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

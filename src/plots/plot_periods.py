# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    plot_periods.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2024/01/02 23:38:49 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/14 10:10:32 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import logging
import os

import numpy as np
import pandas as pd
import xarray as xr


def plot_periods(
    out_track: pd.DataFrame,
    times: pd.DatetimeIndex,
    lat: xr.DataArray,
    results_subdirectory: str,
    figures_directory: str,
    app_logger: logging.Logger,
    processed_vorticity: bool = False,
):
    """
    Calls the cyclophaser function to determine the life cycle periods.

    Parameters:
        out_track (pd.DataFrame): The output track data.
        times (pd.DatetimeIndex): The time values.
        lat (xr.DataArray): The latitude values.
        app_logger: The application logger.

    Returns:
        The result of the determine_periods function.
    """
    from cyclophaser import determine_periods

    app_logger.info("Calling Cyclophaser for determining life cycle periods...")

    periods_figure_directory = os.path.join(figures_directory, "Periods")
    os.makedirs(periods_figure_directory, exist_ok=True)
    vorticity_data = out_track["min_max_zeta_850"].to_list()

    options_high_res = {
        "use_filter": "auto",
        "replace_endpoints_with_lowpass": 24,
        "use_smoothing": "auto",
        "use_smoothing_twice": False,
        "savgol_polynomial": 3,
        "cutoff_low": 168,
        "cutoff_high": 48,
    }

    options_low_res = {
        "use_filter": len(vorticity_data) // 6 if len(vorticity_data) // 6 > 3 else 4,
        "replace_endpoints_with_lowpass": 24,
        "use_smoothing": (
            len(vorticity_data) // 8 | 1 if len(vorticity_data) // 8 | 1 > 3 else 4
        ),
        "use_smoothing_twice": False,
        "savgol_polynomial": 3,
        "cutoff_low": 168,
        "cutoff_high": 12,
    }

    # Select options based on data resolution
    dx = float(np.abs(lat[1] - lat[0]).values)
    if dx < 2:
        options = options_high_res
        app_logger.info(
            f"Using high resolution options for cyclophaser (dx = {dx}): %s", options
        )
    else:
        options = options_low_res
        app_logger.info(
            f"Using low resolution options for cyclophaser (dx = {dx}): %s", options
        )

    # If vorticity is given and the information is on trackfile, it means that
    # it is already processed
    if processed_vorticity:
        options = {
            "use_filter": False,
            "use_smoothing": "auto",
            "use_smoothing_twice": "auto",
        }
        app_logger.info(
            "Vorticity already processed, using default smoothing options for cyclophaser: %s",
            options,
        )

    periods_args = {
        "plot": os.path.join(periods_figure_directory, "periods"),
        "plot_steps": os.path.join(periods_figure_directory, "periods_steps"),
        "export_dict": os.path.join(results_subdirectory, "periods"),
        "process_vorticity_args": options,
    }

    try:
        determine_periods(vorticity_data, x=times, **periods_args)
    except Exception as e:
        app_logger.error(f"Error calling cyclophaser: {e}")
        return

# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    plot_LPS.py                                        :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/06/14 16:32:27 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/14 19:36:12 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from lorenz_phase_space.phase_diagrams import Visualizer


def create_and_save_plot(
    dataframe: pd.DataFrame,
    title: str,
    zoom: bool,
    plot_type: str,
    figures_dir: Path,
    app_logger: logging.Logger,
):

    x_min, x_max = dataframe["Ck"].min(), dataframe["Ck"].max()
    y_min, y_max = dataframe["Ca"].min(), dataframe["Ca"].max()
    color_min, color_max = dataframe["Ge"].min(), dataframe["Ge"].max()
    size_min, size_max = dataframe["Ke"].min(), dataframe["Ke"].max()

    # Initialize the Visualizer with or without zoom
    lps = Visualizer(
        LPS_type="mixed",
        zoom=zoom,
        x_limits=[x_min, x_max],
        y_limits=[y_min, y_max],
        color_limits=[color_min, color_max],
        marker_limits=[size_min, size_max],
    )

    # Plot the data
    lps.plot_data(
        x_axis=dataframe["Ck"].values,
        y_axis=dataframe["Ca"].values,
        marker_color=dataframe["Ge"].values,
        marker_size=dataframe["Ke"].values,
    )

    plt.text(
        0,
        1,
        title,
        va="bottom",
        ha="left",
        fontsize=14,
        fontweight="bold",
        color="gray",
        transform=plt.gca().transAxes,
    )

    # Construct filename and save the plot
    suffix = f"_{plot_type}_zoom" if zoom else f"_{plot_type}"
    file_path = os.path.join(figures_dir, f"LPS{suffix}.png")
    plt.savefig(file_path, dpi=300)
    plt.close()  # Close the plot to free memory
    app_logger.info(f"LPS plot saved to {file_path}")


def plot_LPS(dataframe, infile, results_subdirectory, figures_directory, app_logger):
    app_logger.info("Plotting Lorenz Phase Space...")
    infile_name = Path(infile).stem

    figures_LPS_directory = f"{figures_directory}/LPS"
    os.makedirs(figures_LPS_directory, exist_ok=True)

    # Calculate the time difference between two consecutive timestamps
    time_difference = dataframe.index[1] - dataframe.index[0]

    # Convert the time difference to total seconds and then to hours
    dt_hours = int(time_difference.total_seconds() / 3600)

    periods_file = f"{results_subdirectory}/periods.csv"
    try:
        periods_df = pd.read_csv(
            periods_file, parse_dates=["start", "end"], index_col=0
        )
    except FileNotFoundError:
        app_logger.error("Periods file not found.")
        return
    except Exception as e:
        app_logger.error(f"Error while reading periods file: {e}")
        return

    # Initialize an empty DataFrame to store period means
    period_means_df = pd.DataFrame()

    # Iterate through each period and calculate means
    for period_name, row in periods_df.iterrows():
        start, end = row["start"], row["end"]
        df_period = dataframe.loc[start:end]

        # Check if the period DataFrame is not empty
        if not df_period.empty:
            # Calculate mean for the period
            period_mean = df_period.mean()
            # Add the mean to the period_means_df DataFrame
            period_means_df = pd.concat(
                [period_means_df, pd.DataFrame(period_mean).transpose()],
                ignore_index=True,
            )
            period_means_df.index = period_means_df.index.map(
                lambda x: period_name if x == period_means_df.index[-1] else x
            )
        else:
            app_logger.warning(f"No data available for the period: {period_name}")

    try:
        title, datasource = infile_name.split("_")
    except ValueError:
        title, datasource = infile_name, "unknown"

    figure_title = (
        f"System: {title}\n"
        f"Datasource: {datasource}\n"
        f"Start: {dataframe.index[0].strftime('%Y-%m-%d %HZ')}\n"
        f"End: {dataframe.index[-1].strftime('%Y-%m-%d %HZ')}"
    )

    for zoom in [False, True]:
        create_and_save_plot(
            dataframe,
            figure_title,
            zoom,
            f"{dt_hours}h",
            figures_LPS_directory,
            app_logger,
        )
        create_and_save_plot(
            period_means_df,
            figure_title,
            zoom,
            "periods",
            figures_LPS_directory,
            app_logger,
        )
        df_daily = dataframe.resample("1D").mean()
        create_and_save_plot(
            df_daily, figure_title, zoom, "1d", figures_LPS_directory, app_logger
        )

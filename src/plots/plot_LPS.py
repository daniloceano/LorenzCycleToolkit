# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    plot_LPS.py                                        :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/06/14 16:32:27 by daniloceano       #+#    #+#              #
#    Updated: 2024/01/03 00:46:23 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import pandas as pd
from pathlib import Path
import logging
from lorenz_phase_space.LPS import LorenzPhaseSpace as LPS

def create_and_save_plot(
        dataframe: pd.DataFrame,
        title: str,
        datasource: str,
        zoom: bool,
        plot_type: str,
        figures_dir: Path,
        app_logger: logging.Logger,
        periods: bool = False
    ):
    suffix = f"_{plot_type}_zoom" if zoom else f"_{plot_type}"
    fig, _ = create_plot(dataframe, zoom, title, datasource, periods)
    file_path = os.path.join(figures_dir, f"LPS{suffix}.png")
    fig.savefig(file_path, dpi=300)
    app_logger.info(f"LPS plot saved to {file_path}")


def create_plot(dataframe, zoom, title, datasource, periods=False):
    x_axis, y_axis = dataframe['Ck'], dataframe['Ca']
    marker_color, marker_size = dataframe['Ge'], dataframe['Ke']

    if not periods:
        start, end = map(lambda x: pd.to_datetime(x).strftime('%Y-%m-%d %H:%M'), [dataframe.index[0], dataframe.index[-1]])
    else:
        start, end = None, None
        title, datasource = None, None

    try:
        lps_mixed = LPS(x_axis, y_axis, marker_color, marker_size, zoom=zoom, title=title, datasource=datasource, start=start, end=end)
    except Exception as e:
        raise
    
    return lps_mixed.plot()

def plot_LPS(dataframe, infile, results_subdirectory, figures_directory, app_logger):
    app_logger.info("Plotting Lorenz Phase Space...")
    infile_name = Path(infile).stem
    title, datasource = infile_name.split('_')

    figures_LPS_directory = f"{figures_directory}/LPS"
    os.makedirs(figures_LPS_directory, exist_ok=True)

    # Calculate the time difference between two consecutive timestamps
    time_difference = (dataframe.index[1] - dataframe.index[0])

    # Convert the time difference to total seconds and then to hours
    dt_hours = int(time_difference.total_seconds() / 3600)

    periods_file = f"{results_subdirectory}/periods.csv"
    try:
        periods_df = pd.read_csv(periods_file, parse_dates=['start', 'end'], index_col=0)
    except FileNotFoundError:
        app_logger.error(f"Periods file not found.")
        raise
    except Exception as e:
        app_logger.error(f"Error while reading periods file: {e}")
        raise

    # Initialize an empty DataFrame to store period means
    period_means_df = pd.DataFrame()

    # Iterate through each period and calculate means
    for period_name, row in periods_df.iterrows():
        start, end = row['start'], row['end']
        df_period = dataframe.loc[start:end]

        # Check if the period DataFrame is not empty
        if not df_period.empty:
            # Calculate mean for the period
            period_mean = df_period.mean()
            # Add the mean to the period_means_df DataFrame
            period_means_df = pd.concat([period_means_df, pd.DataFrame(period_mean).transpose()], ignore_index=True)
            period_means_df.index = period_means_df.index.map(lambda x: period_name if x == period_means_df.index[-1] else x)
        else:
            app_logger.warning(f"No data available for the period: {period_name}")

    for zoom in [False, True]:
        create_and_save_plot(dataframe, title, datasource, zoom, f"{dt_hours}h", figures_LPS_directory, app_logger)
        create_and_save_plot(period_means_df, title, datasource, zoom, "periods", figures_LPS_directory, app_logger, periods=True)
        df_daily = dataframe.resample('1D').mean()
        create_and_save_plot(df_daily, title, datasource, zoom, "1d", figures_LPS_directory, app_logger)
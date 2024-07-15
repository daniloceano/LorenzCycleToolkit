# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    timeseries_zeta_and_Z.py                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/28 14:48:27 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/29 13:54:38 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

COLORS = ["#BF3D3B", "#023e8a"]
LINEWIDTH = 3


def plot_min_zeta_hgt(track_file, figures_directory, app_logger=False):
    """
    Plots minimum vorticity and geopotential height from a tracking file.

    Parameters:
    track_file (str): Path to the CSV file containing track data.
    figures_directory (str): Directory path to save the plot.
    """
    try:
        track = pd.read_csv(
            track_file, parse_dates=[0], delimiter=";", index_col="time"
        )
    except FileNotFoundError:
        app_logger.info("Error: Track file was not found.") if app_logger else None
        return
    except pd.errors.ParserError:
        (
            app_logger.info("Error: Parsing CSV file failed. Check the file format.")
            if app_logger
            else None
        )
        return

    plt.close("all")
    fig, ax1 = plt.subplots(figsize=(10, 8))

    ax1.plot(
        track.index,
        track["min_max_zeta_850"],
        "--",
        c=COLORS[0],
        lw=LINEWIDTH,
        label=r"$\zeta$",
    )
    ax2 = ax1.twinx()
    ax2.plot(
        track.index, track["min_hgt_850"], "-", c=COLORS[1], lw=LINEWIDTH, label=r"$Z$"
    )

    ax1.set_ylabel("Vorticity [s⁻¹]", fontsize=16, color=COLORS[0])
    ax2.set_ylabel("Geopotential height [m]", fontsize=16, color=COLORS[1])
    ax1.tick_params(axis="x", labelrotation=20, labelsize=14)
    ax1.tick_params(axis="y", colors=COLORS[0], labelsize=14)
    ax2.tick_params(axis="y", colors=COLORS[1], labelsize=14)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H"))

    fig.legend(
        loc="upper right",
        bbox_to_anchor=(1, 1),
        bbox_transform=ax1.transAxes,
        fontsize=16,
    )
    ax1.grid(c="gray", linewidth=0.25, linestyle="dashdot", axis="x")
    fig.suptitle(
        r"Central Vorticity ($\zeta$) and Geopotential Height (Z) at 850 hPa",
        fontsize=18,
    )

    # Check if the directory exists and create it if not
    os.makedirs(figures_directory, exist_ok=True)
    plt.savefig(
        os.path.join(figures_directory, "timeseries-min_zeta_hgt.png"),
        bbox_inches="tight",
    )


if __name__ == "__main__":
    track_file = "samples/sample_track.csv"
    figures_directory = "samples/Figures/"
    os.makedirs(figures_directory, exist_ok=True)
    plot_min_zeta_hgt(track_file, figures_directory)

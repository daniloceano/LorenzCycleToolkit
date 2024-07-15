# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    map_track.py                                       :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/07/05 19:48:17 by daniloceano       #+#    #+#              #
#    Updated: 2024/01/04 10:40:24 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com

"""

import os
from bisect import bisect_left

import cartopy.crs as ccrs
import cmocean
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

import src.plots.utils as utils
from src.plots.utils import (map_borders, read_results, read_track,
                             setup_gridlines, setup_map)

COLORS = ["#BF3D3B", "#3B95BF"]


def set_plot_legend(fig, labels, title, msizes):
    """
    Set the plot legend with custom scatter plot elements.

    Parameters:
    fig (matplotlib.figure.Figure): The figure object to which the legend is added.
    labelx (float): The x-coordinate for the legend placement.
    labels (list of str): Labels for each category in the legend.
    title (str): Title of the legend.
    msizes (list of int): Marker sizes corresponding to each label.
    """
    legend_elements = [
        plt.scatter([], [], c=utils.TEXT_COLOR, s=size, label=label)
        for label, size in zip(labels, msizes)
    ]
    leg = fig.legend(
        handles=legend_elements,
        ncol=1,
        frameon=False,
        fontsize=utils.LEGEND_FONT_SIZE,
        handlelength=0.3,
        handleheight=4,
        borderpad=1.5,
        scatteryoffsets=[0.1],
        framealpha=1,
        handletextpad=1.5,
        title=title,
        scatterpoints=1,
        loc="center left",
        bbox_to_anchor=(0.8, 0.2),
        labelcolor=utils.TEXT_COLOR,
    )
    leg._legend_box.align = "center"
    plt.setp(
        leg.get_title(),
        color=utils.TEXT_COLOR,
        ha="center",
        fontsize=utils.AXIS_LABEL_FONT_SIZE,
    )


def map_track(results_file_path, trackfile, figures_subdirectory, app_logger=False):
    """
    Generate a map plot of track data with scatter plot elements.

    Parameters:
    results_file_path (str): Path to the results CSV file.
    trackfile (str): Path to the track CSV file.
    figures_subdirectory (str): Directory path to save the generated plot.
    """
    app_logger.info("Plotting track...") if app_logger else print("Plotting track...")

    track = read_track(trackfile, app_logger)
    df = read_results(results_file_path, app_logger)

    # Extract track latitudes and longitudes
    try:
        track_lon = track["Lon"]
        track_lat = track["Lat"]
    except KeyError:
        (
            app_logger.error("Error: Lon and Lat columns not found in track CSV file.")
            if app_logger
            else print("Error: Lon and Lat columns not found in track CSV file.")
        )
        return

    plt.close("all")
    fig, ax = setup_figure(track_lon, track_lat)
    plot_track_data(fig, ax, df, track)

    figure_file = os.path.join(figures_subdirectory, "track.png")
    plt.savefig(figure_file, bbox_inches="tight")
    (
        app_logger.info(f"Track plot saved to {figure_file}")
        if app_logger
        else print(f"Track plot saved to {figure_file}")
    )


def setup_figure(track_lon, track_lat):
    min_lon, max_lon, min_lat, max_lat = (
        track_lon.min(),
        track_lon.max(),
        track_lat.min(),
        track_lat.max(),
    )
    lenx, leny = max_lon - min_lon, max_lat - min_lat
    fig_size = (12, 9) if lenx <= leny else (14, 9)
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_axes(
        [0.05, 0.05, 0.9, 0.9], projection=ccrs.PlateCarree(), frameon=True
    )
    ax.set_extent(
        [min_lon - 10, max_lon + 10, max_lat + 10, min_lat - 10], crs=ccrs.PlateCarree()
    )
    setup_map(ax)
    setup_gridlines(ax)
    map_borders(ax)
    return fig, ax


def plot_track_data(fig, ax, df, track, app_logger=False):
    # Define sizes and intervals for scatter plot markers
    num_quantiles = 5  # Number of quantiles (or size categories)
    msizes = np.linspace(100, 1000, num_quantiles)

    # Extract Ke and Ae values
    try:
        Ke = df["Ke"]
        Ae = df["Ae"]
    except KeyError:
        (
            app_logger.error("Error: Ke and Ae columns not found in results CSV file.")
            if app_logger
            else print("Error: Ke and Ae columns not found in results CSV file.")
        )
        return

    # Extract track latitudes and longitudes
    track_lon, track_lat = track["Lon"], track["Lat"]

    # Scale down Ke and Ae by 10^5
    Ke_scaled = Ke / 1e5
    Ae_scaled = Ae / 1e5

    # Calculate quantile values for scaled Ke
    quantile_values = Ke_scaled.quantile(
        [i / (num_quantiles - 1) for i in np.arange(1, num_quantiles)]
    ).tolist()

    # Assign sizes based on quantiles
    def assign_Ke_size(value):
        index = bisect_left(quantile_values, value)
        # Ensure index is within the range of msizes
        return msizes[min(index, len(msizes) - 1)]

    df["sizes"] = Ke_scaled.apply(assign_Ke_size)

    # Create labels for each quantile interval with rounded values
    # Label for values below the first quantile
    labels = [f"< {quantile_values[0]:.2f}"]
    for i in range(len(quantile_values) - 1):
        labels.append(f"{quantile_values[i]:.2f} - {quantile_values[i + 1]:.2f}")
    # Label for values above the last quantile
    labels.append(f"> {quantile_values[-1]:.2f}")

    # Add label for the max range
    labels[-1] = f"> {quantile_values[-1]:.2f}"
    title = r"Eddy Kinect Energy" + "\n" r"(Ke) [x $10^{5}$ J m$^{-2}$]"
    set_plot_legend(fig, labels, title, msizes)

    # Plotting the track and scatter plot
    ax.plot(track_lon, track_lat, c=utils.TEXT_COLOR)
    norm = colors.TwoSlopeNorm(
        vmin=min(Ae_scaled), vcenter=np.mean(Ae_scaled), vmax=max(Ae_scaled)
    )

    # Adjust DataFrame if necessary
    if len(track) < len(df):
        df = df[df["Datetime"].isin(track.index)]
    elif len(df) < len(track):
        track = track[track.index.isin(df["Datetime"])]

    # Scatter plot for tracking data
    dots = ax.scatter(
        track_lon,
        track_lat,
        c=Ae_scaled,
        cmap=cmocean.cm.matter,
        s=df["sizes"],
        zorder=100,
        edgecolors=utils.MARKER_EDGE_COLOR,
        norm=norm,
    )

    # Colorbar setup
    cax = fig.add_axes(
        [
            ax.get_position().x1 + 0.02,
            ax.get_position().y0 + 0.38,
            0.02,
            ax.get_position().height / 1.74,
        ]
    )
    cbar = plt.colorbar(dots, extend="both", cax=cax)
    cbar.ax.set_ylabel(
        "Eddy Potential Energy (Ae) [x $10^{5}$ " + r" $J\,m^{-2}$" + "]",
        rotation=270,
        fontsize=utils.AXIS_LABEL_FONT_SIZE,
        verticalalignment="bottom",
        color=utils.TEXT_COLOR,
        labelpad=15,
    )

    # Mark the beginning and end of the system
    start = [track_lon.iloc[0], track_lat.iloc[0]]
    end = [track_lon.iloc[-1], track_lat.iloc[-1]]
    ax.text(
        *start,
        "A",
        zorder=101,
        fontsize=24,
        horizontalalignment="center",
        verticalalignment="center",
    )
    ax.text(
        *end,
        "Z",
        zorder=101,
        fontsize=24,
        horizontalalignment="center",
        verticalalignment="center",
    )


if __name__ == "__main__":
    results_file_path = "samples/sample_results.csv"
    track_file = "samples/sample_track.csv"
    figures_subdirectory = "samples/Figures/"
    os.makedirs(figures_subdirectory, exist_ok=True)
    map_track(results_file_path, track_file, figures_subdirectory)

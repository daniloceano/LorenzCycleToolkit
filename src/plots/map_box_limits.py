# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    map_box_limits.py                                  :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/06/14 14:58:51 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/29 14:02:22 by daniloceano      ###   ########.fr        #
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

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon

from src.plots.utils import (map_borders, read_box_limits, setup_gridlines,
                             setup_map)


def plot_box_limits(box_limits_file, figures_directory, app_logger=False):
    """
    Plots the box limits on a map.

    Parameters:
        df (pandas.DataFrame): DataFrame containing the box limits.
        figures_directory (str): Directory to save the figure.
        app_logger: Optional logger for logging messages.
    """
    (
        app_logger.info("Plotting box limits on map...")
        if app_logger
        else print("Plotting box limits on map...")
    )
    df = read_box_limits(box_limits_file)

    min_lon, max_lon = df.loc["min_lon"].item(), df.loc["max_lon"].item()
    min_lat, max_lat = df.loc["min_lat"].item(), df.loc["max_lat"].item()

    plt.close("all")
    fig, ax = plt.subplots(
        figsize=(8, 8), subplot_kw={"projection": ccrs.PlateCarree()}
    )
    ax.set_extent(
        [min_lon - 20, max_lon + 20, max_lat + 20, min_lat - 20], crs=ccrs.PlateCarree()
    )
    map_borders(ax)
    setup_map(ax)
    setup_gridlines(ax)

    # Plot selected domain
    pgon = Polygon(
        [
            (min_lon, min_lat),
            (min_lon, max_lat),
            (max_lon, max_lat),
            (max_lon, min_lat),
            (min_lon, min_lat),
        ]
    )
    ax.add_geometries(
        [pgon],
        crs=ccrs.PlateCarree(),
        facecolor="none",
        edgecolor="#BF3D3B",
        linewidth=3,
        alpha=1,
        zorder=3,
    )

    fig.savefig(os.path.join(figures_directory, "box_limits.png"))
    if app_logger:
        app_logger.info("Figure saved in directory: {}".format(figures_directory))
    else:
        print("Figure saved in directory: {}".format(figures_directory))


if __name__ == "__main__":
    box_limits_file = "samples/sample_box_limits"
    figures_directory = "samples/Figures"
    os.makedirs(figures_directory, exist_ok=True)
    plot_box_limits(box_limits_file, figures_directory)

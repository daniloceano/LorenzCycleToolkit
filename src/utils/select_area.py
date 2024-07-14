#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 12:33:00 2023

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import logging
import os
import time

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cartopy.feature import BORDERS, COASTLINE, NaturalEarthFeature
from metpy.calc import vorticity
from shapely.geometry.polygon import Polygon

nclicks = 2

# define all CRS
crs_longlat = ccrs.PlateCarree()

# Transformation function


def coordXform(orig_crs, target_crs, x, y):
    return target_crs.transform_points(orig_crs, x, y)


def tellme(s):
    plt.title(s, fontsize=16)
    plt.draw()


def draw_box(ax, limits, crs):
    max_lon, min_lon = limits["max_lon"], limits["min_lon"]
    max_lat, min_lat = limits["max_lat"], limits["min_lat"]
    pgon = Polygon(
        (
            (min_lon, min_lat),
            (min_lon, max_lat),
            (max_lon, max_lat),
            (max_lon, min_lat),
            (min_lon, min_lat),
        )
    )
    ax.add_geometries(
        [pgon], crs=crs, facecolor="None", edgecolor="k", linewidth=3, alpha=1, zorder=3
    )


def plot_zeta(ax, zeta, lat, lon, hgt=None):
    if np.abs(zeta.min()) < np.abs(zeta.max()):
        norm = colors.TwoSlopeNorm(vmin=-zeta.max(), vcenter=0, vmax=zeta.max())
    else:
        norm = colors.TwoSlopeNorm(vmin=zeta.min(), vcenter=0, vmax=-zeta.min())
    cmap = cmo.balance
    # plot contours
    cf1 = ax.contourf(
        lon, lat, zeta, cmap=cmap, norm=norm, levels=51, transform=crs_longlat
    )
    plt.colorbar(cf1, pad=0.12, orientation="vertical", shrink=0.5)
    if hgt is not None:
        cs = ax.contour(
            lon,
            lat,
            hgt,
            levels=11,
            colors="#344e41",
            linestyles="dashed",
            linewidths=1.3,
            transform=crs_longlat,
        )
        ax.clabel(cs, cs.levels, inline=True, fontsize=10)


def map_decorators(ax):
    ax.coastlines()
    gl = ax.gridlines(
        draw_labels=True,
        zorder=2,
        linestyle="dashed",
        alpha=0.7,
        linewidth=0.5,
        color="#383838",
    )
    gl.xlabel_style = {"size": 14, "color": "#383838"}
    gl.ylabel_style = {"size": 14, "color": "#383838"}
    gl.top_labels = None
    gl.right_labels = None


def plot_min_zeta(ax, zeta, lat, lon, limits):
    max_lon, min_lon = limits["max_lon"], limits["min_lon"]
    max_lat, min_lat = limits["max_lat"], limits["min_lat"]

    # Plot mininum zeta point whithin box
    izeta = zeta.sel(
        {lon.dims[0]: slice(min_lon, max_lon), lat.dims[0]: slice(min_lat, max_lat)}
    )
    if min_lat < 0:
        min_max_zeta = izeta.min()
    else:
        min_max_zeta = izeta.max()

    min_max_zeta_loc = izeta.where(izeta == min_max_zeta, drop=True).squeeze()

    # sometimes there are multiple minimuns
    if min_max_zeta_loc.shape:
        if len(min_max_zeta_loc.shape) > 1:
            for points in min_max_zeta_loc:
                for point in points:
                    ax.scatter(
                        point[lon.dims[0]],
                        point[lat.dims[0]],
                        marker="o",
                        facecolors="none",
                        linewidth=3,
                        edgecolor="k",
                        s=200,
                    )
        else:
            for point in min_max_zeta_loc:
                ax.scatter(
                    point[lon.dims[0]],
                    point[lat.dims[0]],
                    marker="o",
                    facecolors="none",
                    linewidth=3,
                    edgecolor="k",
                    s=200,
                )
    else:
        ax.scatter(
            min_max_zeta_loc[lon.dims[0]],
            min_max_zeta_loc[lat.dims[0]],
            marker="o",
            facecolors="none",
            linewidth=3,
            edgecolor="k",
            s=200,
        )


def initial_domain(zeta, lat, lon):
    plt.close("all")
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=crs_longlat)
    fig.add_axes(ax)
    ax.set_global()
    plot_zeta(ax, zeta.compute(), lat, lon)
    map_decorators(ax)
    plt.subplots_adjust(bottom=0, top=1.2)

    while True:
        pts = []
        while len(pts) < nclicks:
            tellme("Select an initial spatial domain")
            pts = np.asarray(
                plt.ginput(nclicks, timeout=15, mouse_stop="MouseButton.MIDDLE")
            )
            if len(pts) < nclicks:
                tellme("Too few points, starting over")
                time.sleep(1)  # Wait a second

        xmin, xmax = min([pts[0, 0], pts[1, 0]]), max([pts[0, 0], pts[1, 0]])
        ymin, ymax = min([pts[0, 1], pts[1, 1]]), max([pts[0, 1], pts[1, 1]])
        xs, ys = np.array((xmin, xmax)), np.array((ymin, ymax))
        lls = coordXform(crs_longlat, crs_longlat, xs, ys)
        min_lon, min_lat = lls[0, 0], lls[0, 1]
        max_lon, max_lat = lls[1, 0], lls[1, 1]

        limits = {
            "min_lon": min_lon,
            "max_lon": max_lon,
            "min_lat": min_lat,
            "max_lat": max_lat,
        }
        draw_box(ax, limits, crs_longlat)

        tellme("Happy? Key press any keyboard key for yes, mouse click for no")
        if plt.waitforbuttonpress():
            break

    return limits


def draw_box_map(u, v, zeta, hgt, lat, lon, timestr):
    plt.close("all")
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=crs_longlat)
    fig.add_axes(ax)

    plot_zeta(ax, zeta, lat, lon, hgt)
    ax.streamplot(
        lon.values,
        lat.values,
        u.values,
        v.values,
        color="#2A1D21",
        transform=crs_longlat,
    )
    map_decorators(ax)

    while True:
        pts = []
        while len(pts) < nclicks:
            time_formated = timestr.strftime("%Y-%m-%d %HZ")
            tellme(f"Select box corners \nModel time step: {time_formated}")
            pts = np.asarray(
                plt.ginput(nclicks, timeout=15, mouse_stop="MouseButton.MIDDLE")
            )
            if len(pts) < nclicks:
                tellme("Too few points, starting over")
                time.sleep(1)  # Wait a second

        xmin, xmax = min([pts[0, 0], pts[1, 0]]), max([pts[0, 0], pts[1, 0]])
        ymin, ymax = min([pts[0, 1], pts[1, 1]]), max([pts[0, 1], pts[1, 1]])
        xs, ys = np.array((xmin, xmax)), np.array((ymin, ymax))
        lls = coordXform(crs_longlat, crs_longlat, xs, ys)
        min_lon, min_lat = lls[0, 0], lls[0, 1]
        max_lon, max_lat = lls[1, 0], lls[1, 1]

        limits = {
            "min_lon": min_lon,
            "max_lon": max_lon,
            "min_lat": min_lat,
            "max_lat": max_lat,
        }
        draw_box(ax, limits, crs_longlat)
        plot_min_zeta(ax, zeta, lat, lon, limits)

        tellme("Happy? Key press any keyboard key for yes, mouse click for no")

        if plt.waitforbuttonpress():
            break

    return limits


def slice_domain(NetCDF_data, args, variable_list_df):
    """
    Slice a NetCDF dataset based on given limits and return the sliced dataset.

    Parameters:
        NetCDF_data (xarray.Dataset): The NetCDF dataset to be sliced.
        args (argparse.Namespace): Command-line arguments.
        variable_list_df (pandas.DataFrame): The DataFrame containing the variable list.

    Returns:
        xarray.Dataset: The sliced NetCDF dataset.
    """
    # Data indexers
    LonIndexer = variable_list_df.loc["Longitude"]["Variable"]
    LatIndexer = variable_list_df.loc["Latitude"]["Variable"]
    TimeIndexer = variable_list_df.loc["Time"]["Variable"]
    LevelIndexer = variable_list_df.loc["Vertical Level"]["Variable"]

    if args.fixed:
        dfbox = pd.read_csv(
            "inputs/box_limits", header=None, delimiter=";", index_col=0
        )
        WesternLimit = float(
            NetCDF_data[LonIndexer].sel(
                {LonIndexer: float(dfbox.loc["min_lon"].values[0])}, method="nearest"
            )
        )
        EasternLimit = float(
            NetCDF_data[LonIndexer].sel(
                {LonIndexer: float(dfbox.loc["max_lon"].values[0])}, method="nearest"
            )
        )
        SouthernLimit = float(
            NetCDF_data[LatIndexer].sel(
                {LatIndexer: float(dfbox.loc["min_lat"].values[0])}, method="nearest"
            )
        )
        NorthernLimit = float(
            NetCDF_data[LatIndexer].sel(
                {LatIndexer: float(dfbox.loc["max_lat"].values[0])}, method="nearest"
            )
        )

    elif args.track:
        dx, dy = (
            NetCDF_data[LonIndexer][1] - NetCDF_data[LonIndexer][0],
            NetCDF_data[LatIndexer][1] - NetCDF_data[LatIndexer][0],
        )

        trackfile = args.trackfile if args.trackfile else "inputs/track"
        track = pd.read_csv(trackfile, parse_dates=[0], delimiter=";", index_col="time")
        if "width" in track.columns:
            max_width, max_length = track["width"].max(), track["length"].max()
        else:
            max_width, max_length = 15, 15

        WesternLimit = track["Lon"].min() - (max_width / 2) - dx
        EasternLimit = track["Lon"].max() + (max_width / 2) + dx
        SouthernLimit = track["Lat"].min() - (max_length / 2) - dy
        NorthernLimit = track["Lat"].max() + (max_length / 2) + dy

    elif args.choose:
        iu_850 = NetCDF_data.isel({TimeIndexer: 0}).sel(
            {LevelIndexer: 8500}, method="nearest"
        )[variable_list_df.loc["Eastward Wind Component"]["Variable"]]
        iv_850 = NetCDF_data.isel({TimeIndexer: 0}).sel(
            {LevelIndexer: 8500}, method="nearest"
        )[variable_list_df.loc["Northward Wind Component"]["Variable"]]
        zeta = vorticity(iu_850, iv_850).metpy.dequantify()

        lat, lon = iu_850[LatIndexer], iu_850[LonIndexer]
        domain_limits = initial_domain(zeta, lat, lon)
        WesternLimit = domain_limits["min_lon"]
        EasternLimit = domain_limits["max_lon"]
        SouthernLimit = domain_limits["min_lat"]
        NorthernLimit = domain_limits["max_lat"]

    NetCDF_data = NetCDF_data.sel(
        **{
            LatIndexer: slice(SouthernLimit, NorthernLimit),
            LonIndexer: slice(WesternLimit, EasternLimit),
        }
    )

    return NetCDF_data


def plot_domain_attributes(data850, limits, position, figures_directory):

    time = limits["datestr"]
    central_lon, central_lat = limits["central_lon"], limits["central_lat"]
    width, length = limits["width"], limits["length"]
    min_lon = central_lon - (width / 2)
    max_lon = central_lon + (width / 2)
    min_lat = central_lat - (length / 2)
    max_lat = central_lat + (length / 2)

    # Create figure
    plt.close("all")
    fig, ax = plt.subplots(
        figsize=(8, 8.5), subplot_kw=dict(projection=ccrs.PlateCarree())
    )

    # Set map extent and features
    ax.set_extent(
        [min_lon - 20, max_lon + 20, max_lat + 20, min_lat - 20], crs=ccrs.PlateCarree()
    )

    # Plot selected domain
    # Create a sample polygon, `pgon`
    pgon = Polygon(
        (
            (min_lon, min_lat),
            (min_lon, max_lat),
            (max_lon, max_lat),
            (max_lon, min_lat),
            (min_lon, min_lat),
        )
    )
    ax.add_geometries(
        [pgon],
        crs=ccrs.PlateCarree(),
        facecolor="None",
        edgecolor="gray",
        linewidth=3,
        alpha=1,
        zorder=3,
    )

    # Plot central point, mininum vorticity, minimum hgt and maximum wind
    ax.scatter(central_lon, central_lat, marker="o", c="#31332e", s=100, zorder=4)
    ax.scatter(
        position["min_max_zeta_850"]["longitude"],
        position["min_max_zeta_850"]["latitude"],
        marker="s",
        c="#31332e",
        s=100,
        zorder=4,
        label="min/max zeta",
    )
    ax.scatter(
        position["min_hgt_850"]["longitude"],
        position["min_hgt_850"]["latitude"],
        marker="x",
        c="#31332e",
        s=100,
        zorder=4,
        label="min hgt",
    )
    ax.scatter(
        position["max_wind_850"]["longitude"],
        position["max_wind_850"]["latitude"],
        marker="^",
        c="#31332e",
        s=100,
        zorder=4,
        label="max wind",
    )

    plt.legend(loc="upper left", frameon=True, fontsize=14, bbox_to_anchor=(1.1, 1.2))

    # Add gridlines
    gl = ax.gridlines(draw_labels=True, zorder=2)
    gl.xlabel_style = {"size": 16}
    gl.ylabel_style = {"size": 16}

    # Add title
    plt.title("Box defined for computations\n", fontsize=22)

    plot_zeta(
        ax, data850["izeta_850"], data850["lat"], data850["lon"], data850["izeta_850"]
    )
    ax.add_feature(COASTLINE, edgecolor="#283618", linewidth=1)
    ax.add_feature(BORDERS, edgecolor="#283618", linewidth=1)
    _ = ax.add_feature(
        cfeature.NaturalEarthFeature(
            "physical", "land", "50m", edgecolor="face", facecolor="None"
        )
    )

    states = NaturalEarthFeature(
        category="cultural",
        scale="50m",
        facecolor="none",
        name="admin_1_states_provinces_lines",
    )
    _ = ax.add_feature(states, edgecolor="#283618", linewidth=1)

    cities = NaturalEarthFeature(
        category="cultural", scale="50m", facecolor="none", name="populated_places"
    )
    _ = ax.add_feature(cities, edgecolor="#283618", linewidth=1)

    # Save figure
    boxes_directory = os.path.join(figures_directory, "boxes")
    os.makedirs(boxes_directory, exist_ok=True)
    filename = os.path.join(boxes_directory, f"box_{time}.png")
    plt.savefig(filename)
    logging.debug(f"\nCreated figure with box defined for computations at {filename}")

# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    utils.py                                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/29 10:02:20 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/18 09:24:41 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import re
from glob import glob

import cartopy
import cartopy.feature as cfeature
import pandas as pd
from cartopy.feature import NaturalEarthFeature

# Global constants for plotting
TERM_DETAILS = {
    "energy": {"terms": ["Az", "Ae", "Kz", "Ke"], "label": "Energy", "unit": "J·m⁻²"},
    "conversion": {
        "terms": ["Cz", "Ca", "Ck", "Ce"],
        "label": "Conversion",
        "unit": "W·m⁻²",
    },
    "boundary": {
        "terms": ["BAz", "BAe", "BKz", "BKe"],
        "label": "Transport across boundaries",
        "unit": "W·m⁻²",
    },
    "budget_diff": {
        "terms": [
            "∂Az/∂t (finite diff.)",
            "∂Ae/∂t (finite diff.)",
            "∂Kz/∂t (finite diff.)",
            "∂Ke/∂t (finite diff.)",
        ],
        "label": "Energy budgets (estimated using finite diffs.)",
        "unit": "W·m⁻²",
    },
    "residuals": {
        "terms": ["RGz", "RKz", "RGe", "RKe"],
        "label": "Residuals",
        "unit": "W·m⁻²",
    },
    "generation_dissipation": {
        "terms": ["Gz", "Ge", "Dz", "De"],
        "label": "Generation/Dissipation",
        "unit": "W·m⁻²",
    },
    "comparing_generation": {
        "terms": ["RGz", "RGe", "Gz", "Ge"],
        "label": "Comparing Generation",
        "unit": "W·m⁻²",
    },
    "comparing_dissipation": {
        "terms": ["RKz", "Dz", "RKe", "De"],
        "label": "Comparing Dissipation",
        "unit": "W·m⁻²",
    },
}

COLORS = ["#3B95BF", "#87BF4B", "#BFAB37", "#BF3D3B", "#873e23", "#A13BF0"]
MARKERS = ["s", "o", "^", "v", "<", ">"]
MARKER_COLORS = ["#59c0f0", "#b0fa61", "#f0d643", "#f75452", "#f07243", "#bc6ff7"]
LINESTYLE = "-"
LINEWIDTH = 3
TEXT_COLOR = "#383838"
MARKER_EDGE_COLOR = "grey"
LEGEND_FONT_SIZE = 10
AXIS_LABEL_FONT_SIZE = 12
TITLE_FONT_SIZE = 18


def read_results(results_file_path, app_logger=False):
    """
    Read the results from a CSV file and return a pandas DataFrame.

    Parameters:
        results_file_path (str): The path to the CSV file containing the results.

    Returns:
        pandas.DataFrame or None: The DataFrame containing the results, or None if an error occurs.
    """
    try:
        df = pd.read_csv(results_file_path, index_col=[0])
        df["Datetime"] = pd.to_datetime(df.index)
        df.set_index("Datetime", inplace=True)
        return df
    except FileNotFoundError:
        (
            app_logger.error(f"Error: File {results_file_path} was not found.")
            if app_logger
            else print(f"Error: File {results_file_path} was not found.")
        )
        return None
    except pd.errors.ParserError:
        (
            app_logger.error(f"Error: {results_file_path} is not a valid CSV file.")
            if app_logger
            else print(f"Error: {results_file_path} is not a valid CSV file.")
        )
        return None


def read_track(trackfile, app_logger=False):
    """
    Reads a track file and returns a pandas DataFrame.

    Parameters:
        trackfile (str): The path to the track file.

    Returns:
        pandas.DataFrame or None: The track data as a DataFrame, or None if an error occurs.
    """
    try:
        return pd.read_csv(trackfile, parse_dates=[0], delimiter=";", index_col="time")
    except FileNotFoundError:
        (
            app_logger.error(f"Error: File {trackfile} was not found.")
            if app_logger
            else print(f"Error: File {trackfile} was not found.")
        )
        return None
    except Exception as e:
        (
            app_logger.error(f"Error reading {trackfile}: {e}")
            if app_logger
            else print(f"Error reading {trackfile}: {e}")
        )
        return None


def read_box_limits(box_limits_file, app_logger=None):
    """
    Reads the box limits from a CSV file.

    Parameters:
        box_limits_file (str): Path to the CSV file containing box limits.
        app_logger: Optional logger for logging messages.

    Returns:
        pandas.DataFrame: DataFrame with box limits or None if an error occurs.
    """
    try:
        df = pd.read_csv(box_limits_file, sep=";", index_col=0, header=None)
        df = df.loc[["min_lon", "max_lon", "min_lat", "max_lat"]]
        return df
    except (FileNotFoundError, pd.errors.ParserError, KeyError) as e:
        if app_logger:
            app_logger.error(f"Error reading box limits file: {e}")
        else:
            print(f"Error reading box limits file: {e}")
        return None


def get_data_vertical_levels(results_subdirectory):
    """
    Retrieves data from CSV files in a given directory and organizes it based on the term found in each file.

    Parameters:
        results_directory (str): The directory containing the CSV files.

    Returns:
        data (dict): A dictionary where the keys are terms found in the CSV file names and the values are pandas DataFrames
        containing the data from the corresponding files.
    """
    # Set results subdirectory
    results_subdirectory_vertical_levels = os.path.join(
        results_subdirectory, "results_vertical_levels"
    )

    # Get list of CSV files but remove files for split terms
    files = [
        file
        for file in glob(os.path.join(results_subdirectory_vertical_levels, "*.csv"))
        if not re.search(r"_[0-9]", os.path.basename(file))
    ]

    data = {}
    for file in files:
        term = os.path.splitext(os.path.basename(file))[0].split("_")[0]
        try:
            data[term] = pd.read_csv(file, header=0, index_col=0, parse_dates=True)
            data[term].columns = [float(x) for x in data[term].columns]
        except pd.errors.ParserError as e:
            print(f"Error reading {file}: {e}")

    return data


def setup_gridlines(ax):
    gl = ax.gridlines(
        draw_labels=True,
        zorder=2,
        linestyle="-",
        alpha=0.8,
        color=TEXT_COLOR,
        linewidth=0.25,
    )
    gl.xlabel_style = {"size": 14, "color": TEXT_COLOR}
    gl.ylabel_style = {"size": 14, "color": TEXT_COLOR}
    gl.bottom_labels = None
    gl.right_labels = None


def setup_map(ax):
    ax.coastlines(zorder=1)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN, facecolor="lightblue")


def map_borders(ax):
    # Add land feature (no facecolor)
    ax.add_feature(
        cfeature.NaturalEarthFeature(
            "physical", "land", "50m", edgecolor="face", facecolor="none"
        )
    )

    # Add state borders
    states = NaturalEarthFeature(
        category="cultural",
        scale="50m",
        facecolor="none",
        name="admin_1_states_provinces_lines",
    )
    ax.add_feature(states, edgecolor="#283618", linewidth=1)

    # Add populated places (cities)
    cities = NaturalEarthFeature(
        category="cultural", scale="50m", facecolor="none", name="populated_places"
    )
    ax.add_feature(cities, edgecolor="#283618", linewidth=1)

    # Add country borders
    countries = NaturalEarthFeature(
        category="cultural", scale="50m", facecolor="none", name="admin_0_countries"
    )
    ax.add_feature(countries, edgecolor="black", linewidth=1)

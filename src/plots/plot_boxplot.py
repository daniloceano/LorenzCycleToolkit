#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 14:22:40 2022

This script reads CSV files with energy and conversion terms from the Lorenz
Energy Cycle and make boxplots.


Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import os
from datetime import datetime

import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

import src.plots.utils as utils
from src.plots.utils import get_data_vertical_levels, read_results


def boxplot_time(dict_vertical, figures_subdirectory, app_logger):
    """
    Generates a boxplot for each term in the given dictionary and saves the resulting figures in the specified subdirectory.

    Parameters:
    - dict_vertical: A dictionary containing vertical data for different terms.
    - figures_subdirectory: The subdirectory where the generated figures will be saved.
    - app_logger: (Optional) The application logger to log information about the saved figures.

    Returns:
    None
    """
    term_details = utils.TERM_DETAILS
    all_terms = {
        key: [term for term in details["terms"] if term in dict_vertical]
        for key, details in term_details.items()
        if key in ["energy", "conversion", "generation_dissipation"]
    }

    dummy_data = dict_vertical[next(iter(all_terms["energy"]), None)]  # Safe access
    t_id = dummy_data.columns[0]  # Time indexer
    times = dummy_data[t_id].index

    for term_type, terms in all_terms.items():
        if not terms:
            continue  # Skip if no terms to plot

        plt.close("all")
        fig_layout = (10, 10) if len(terms) == 4 else (10, 5)
        fig = plt.figure(figsize=fig_layout)
        n_rows, n_cols = (2, 2) if len(terms) == 4 else (1, 2)

        gs = gridspec.GridSpec(nrows=n_rows, ncols=n_cols, hspace=0.25, wspace=0.3)

        for i, term in enumerate(terms):
            ax = fig.add_subplot(gs[i])
            if i % n_cols == 0:
                label = f'{term_details[term_type]["label"]} ({term_details[term_type]["unit"]})'
                ax.set_ylabel(label, fontsize=18)
            for t in times:
                data = dict_vertical[term].loc[t].values
                time_step = mdates.date2num(
                    datetime.fromisoformat(t.strftime("%Y-%m-%d %H:%M"))
                )
                bplot = ax.boxplot(data, positions=[time_step], patch_artist=True)
                bplot["boxes"][-1].set_facecolor(utils.COLORS[i % len(utils.COLORS)])
                bplot["boxes"][-1].set_alpha(0.85)
                for median in bplot["medians"]:
                    median.set_color("k")

            locator = mdates.AutoDateLocator(minticks=5, maxticks=12)
            formatter = mdates.ConciseDateFormatter(locator)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(formatter)
            ax.tick_params(axis="x", labelrotation=20)
            ax.xaxis.set_tick_params(labelsize=16)
            ax.yaxis.set_tick_params(labelsize=16)
            ax.set_title(term, fontsize=20)

        file_name = f"{figures_subdirectory}/vertical_timeseries_{term_type}.png"
        plt.savefig(file_name, bbox_inches="tight")
        (
            app_logger.info(f"Figure saved in directory: {figures_subdirectory}")
            if app_logger
            else print(f"Figure saved in directory: {figures_subdirectory}")
        )
        plt.close("all")


def boxplot_vertical(dict_vertical, figures_subdirectory, app_logger):
    # Define the labels lists based on utils.TERM_DETAILS
    term_details = utils.TERM_DETAILS
    all_terms = {
        key: [term for term in details["terms"] if term in dict_vertical]
        for key, details in term_details.items()
        if key in ["energy", "conversion", "generation_dissipation"]
    }

    linecolors = utils.COLORS

    for term_type, terms in all_terms.items():
        if not terms:
            continue  # Skip if no terms to plot

        plt.close("all")
        fig_layout = (10, 10) if len(terms) == 4 else (10, 5)
        fig = plt.figure(figsize=fig_layout)
        n_rows, n_cols = (2, 2) if len(terms) == 4 else (1, 2)

        wspace = 0.3 if term_type != "generation_dissipation" else 0.4

        gs = gridspec.GridSpec(nrows=n_rows, ncols=n_cols, hspace=0.15, wspace=wspace)

        label = (
            f'{term_details[term_type]["label"]} ({term_details[term_type]["unit"]})'
        )
        for i, term in enumerate(terms):
            ax = fig.add_subplot(gs[i])

            if i % 2 == 0:
                ax.set_ylabel(label, fontsize=18)

            # Assuming second column onward are levels
            levs = dict_vertical[term].columns[1:]
            for j, lev in enumerate(levs):
                bplot = ax.boxplot(
                    dict_vertical[term][lev].values,
                    positions=[j / 3],
                    tick_labels=[lev],
                    patch_artist=True,
                )
                bplot["boxes"][-1].set_facecolor(linecolors[i % len(linecolors)])
                bplot["boxes"][-1].set_alpha(0.85)
                for median in bplot["medians"]:
                    median.set_color("k")

            ax.tick_params(axis="x", labelrotation=60)
            if i < len(terms) - 2:  # Assuming last row should have x labels
                ax.axes.xaxis.set_ticklabels([])
            ax.xaxis.set_tick_params(labelsize=16)
            ax.yaxis.set_tick_params(labelsize=16)
            ax.set_title(term, fontsize=20)

        outfile = f"{figures_subdirectory}/boxplot_vertical_levels_{term_type}.png"
        plt.savefig(outfile, bbox_inches="tight")
        (
            app_logger.info(f"Created {outfile}")
            if app_logger
            else print(f"Created {outfile}")
        )
        plt.close("all")


def plot_boxplot_terms(df, figures_subdirectory, app_logger):
    """
    Generates a boxplot for each term in the term_list from the given DataFrame and saves the resulting figure.

    Parameters:
    - df: The DataFrame containing the data to be plotted.
    - term_list: List of terms to plot.
    - figures_subdirectory: The directory where the figure will be saved.
    - app_logger: Application logger for logging messages.
    """
    # Define the labels lists based on utils.TERM_DETAILS
    term_details = utils.TERM_DETAILS
    all_terms = {
        key: [term for term in details["terms"] if term in df.columns]
        for key, details in term_details.items()
        if key in term_details.keys()
    }

    for term_type, term_list in all_terms.items():
        if not term_list:
            continue  # Skip if no terms to plot

        # Determine the figure size based on the number of terms
        fig_layout = (8, 8) if len(term_list) == 4 else (6, 8)
        plt.figure(figsize=fig_layout)

        for i, term in enumerate(term_list):
            bplot = plt.boxplot(
                df[term],
                positions=[i],
                vert=True,
                patch_artist=True,
                notch=True,
                tick_labels=[term],
            )
            box_color = utils.COLORS[i % len(utils.COLORS)]
            for box in bplot["boxes"]:
                box.set_facecolor(box_color)
                box.set_alpha(0.85)
            for median in bplot["medians"]:
                median.set_color("k")

            # Horizontal line for 0 if not an energy term
            if term not in utils.TERM_DETAILS["energy"]["terms"]:
                plt.axhline(
                    y=0, color="k", linestyle="-", linewidth=1, zorder=1, alpha=0.8
                )

            label = utils.TERM_DETAILS[term_type]["label"]
            plt.ylabel(
                f'{label} ({utils.TERM_DETAILS[term_type]["unit"]})', fontsize=14
            )
            plt.tick_params(axis="x", labelsize=12)
            plt.tick_params(axis="y", labelsize=12)

        # Set x-ticks rotation if budget difference terms are included
        if any(
            term in utils.TERM_DETAILS["budget_diff"]["terms"] for term in term_list
        ):
            plt.xticks(rotation=25)

        plt.tight_layout()

        # Saving figure
        fname = f"boxplot_terms_{term_type}"
        plt.savefig(f"{figures_subdirectory}/{fname}.png")
        (
            app_logger.info(f"{fname}.png created")
            if app_logger
            else print(f"{fname}.png created")
        )
        plt.close("all")


def boxplot_terms(results_file, results_subdirectory, figures_directory, app_logger=False):
    (
        app_logger.info("Plotting boxplots...")
        if app_logger
        else print("Plotting boxplots...")
    )

    df = read_results(results_file, app_logger)

    dict_vertical = get_data_vertical_levels(results_subdirectory)

    figures_subdirectory = os.path.join(figures_directory, "boxplots")
    os.makedirs(figures_subdirectory, exist_ok=True)

    (
        app_logger.info("Creating boxplot for each time")
        if app_logger
        else print("Creating boxplot for each time")
    )
    boxplot_time(dict_vertical, figures_subdirectory, app_logger)

    (
        app_logger.info("Creating boxplot for each vertical level")
        if app_logger
        else print("Creating boxplot for each vertical level")
    )
    boxplot_vertical(dict_vertical, figures_subdirectory, app_logger)

    (
        app_logger.info("Creating boxplot for each term")
        if app_logger
        else print("Creating boxplot for each term")
    )
    plot_boxplot_terms(df, figures_subdirectory, app_logger)


if __name__ == "__main__":
    results_file = "samples/sample_results.csv"
    results_directory = "samples/"
    figures_directory = "samples/Figures/"
    boxplot_terms(results_file, results_directory, figures_directory, app_logger=False)

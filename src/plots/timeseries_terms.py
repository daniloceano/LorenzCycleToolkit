# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    timeseries_terms.py                                :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/05/22 17:56:57 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/29 14:04:06 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
This script reads an CSV file with energy and conversion terms from the Lorenz
Energy Cycle (as input from user) and plot a timeseries for each.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""
import os

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

import src.plots.utils as utils
from src.plots.utils import read_results


def _plotter(df, term_list, figures_directory):
    """Plots a timeseries for given terms."""

    plt.close("all")
    plt.figure(figsize=(10, 8))
    ax = plt.gca()

    terms = utils.TERM_DETAILS[term_list]["terms"]
    title = utils.TERM_DETAILS[term_list]["label"]
    y_label = f'[{utils.TERM_DETAILS[term_list]["unit"]}]'
    times = df.index

    for i, term in enumerate(terms):
        if term in df.columns:
            plt.plot(
                times,
                df[term],
                c=utils.COLORS[i],
                label=term,
                lw=utils.LINEWIDTH,
                marker=utils.MARKERS[i],
                markersize=8,
                markerfacecolor=utils.COLORS[i],
                markeredgecolor=utils.MARKER_EDGE_COLOR,
                markeredgewidth=0.5,
                linestyle=utils.LINESTYLE,
            )

    plt.grid(c="gray", linewidth=0.15, linestyle="-")
    plt.legend(prop={"size": 18})
    plt.title(title, fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel(y_label, fontsize=18)
    plt.xlim(times[0], times[-1])
    plt.tick_params(axis="x", labelrotation=20)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H"))

    plt.savefig(os.path.join(figures_directory, f"timeseires_{term_list}.png"))


def plot_timeseries(results_file, figures_directory, app_logger=False):
    """Reads data from an results_file and plots timeseries for LEC terms."""
    (
        app_logger.info("Plotting timeseries...")
        if app_logger
        else print("Plotting timeseries...")
    )
    df = read_results(results_file, app_logger)

    figures_subdirectory = os.path.join(figures_directory, "timeseries")
    os.makedirs(figures_subdirectory, exist_ok=True)

    for term_list in utils.TERM_DETAILS:
        _plotter(df, term_list, figures_subdirectory)
        (
            app_logger.info(
                f"Figure saved for {term_list} in directory: {figures_subdirectory}"
            )
            if app_logger
            else print(
                f"Figure saved for {term_list} in directory: {figures_subdirectory}"
            )
        )


if __name__ == "__main__":
    results_file = "samples/sample_results.csv"
    figures_directory = "samples/Figures/"
    plot_timeseries(results_file, figures_directory)

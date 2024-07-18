#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:12:34 2022

This script reads CSV files with energy and conversion terms from the Lorenz
Energy Cycle, for each vertical level and time step, and plots them.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import os

import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import src.plots.utils as utils
from src.plots.utils import get_data_vertical_levels, read_results


def _plotter(dict_vertical, figures_subdirectory, app_logger=False):

    times = pd.to_datetime(dict_vertical["Az"].index)
    levs = dict_vertical["Az"].columns

    term_details = utils.TERM_DETAILS
    all_terms = {
        key: [term for term in details["terms"] if term in dict_vertical]
        for key, details in term_details.items()
        if key in ["energy", "conversion", "generation_dissipation"]
    }

    for term_type, terms in all_terms.items():
        if not terms:
            continue  # Skip if no terms to plot
        fig, gs = create_figure_for_term(term_type)

        for i, term in enumerate(terms):
            ax = fig.add_subplot(gs[i])
            plot_contour(times, levs, dict_vertical[term], ax, term)
            finalize_plot(ax, term)

        outfile = os.path.join(figures_subdirectory, f"hovmoller_{term_type}.png")
        plt.savefig(outfile, bbox_inches="tight")
        (
            app_logger.info(f"Figure saved to {outfile}")
            if app_logger
            else print(f"Figure saved to {outfile}")
        )


def create_figure_for_term(term_type):
    if term_type == "energy":
        fig = plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(nrows=2, ncols=2, wspace=0.5, hspace=0.3, top=0.95)
    elif term_type == "conversion":
        fig = plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(nrows=2, ncols=2, wspace=0.75, hspace=0.3, top=0.95)
    elif term_type == "generation_dissipation":
        fig = plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(nrows=1, ncols=2, right=1.2, wspace=0.4)

    return fig, gs


def plot_contour(times, levs, data, ax, term):
    imin = data.min(numeric_only=True).min()
    imax = data.max(numeric_only=True).max()
    absmax = np.amax([np.abs(imin), imax])
    cmap, norm, title = get_cmap_norm_title_for_term(term, imin, imax, absmax)

    cf = ax.contourf(
        times, levs, data[levs].transpose(), cmap=cmap, extend="neither", norm=norm
    )
    ax.contour(times, levs, data[levs].transpose(), colors="k", extend="neither")

    def fmt(x, pos):
        a, b = "{:.2e}".format(x).split("e")
        b = int(b)
        return r"${} \times 10^{{{}}}$".format(a, b)

    if (
        term in utils.TERM_DETAILS["conversion"]["terms"]
        or term in utils.TERM_DETAILS["generation_dissipation"]["terms"]
    ):
        cbar = plt.colorbar(cf, ax=ax, format=ticker.FuncFormatter(fmt))
    else:
        cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label(title, rotation=270, labelpad=20, fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    ax.invert_yaxis()


def get_cmap_norm_title_for_term(term, imin, imax, absmax):
    if term in utils.TERM_DETAILS["energy"]["terms"]:
        cmap = "cmo.amp"
        norm = colors.Normalize(vmin=imin, vmax=imax)
        title = f"{utils.TERM_DETAILS['energy']['label']} ({utils.TERM_DETAILS['energy']['unit']})"
    elif term in utils.TERM_DETAILS["conversion"]["terms"]:
        cmap = "cmo.tarn"
        norm = colors.TwoSlopeNorm(vcenter=0, vmin=-absmax, vmax=absmax)
        title = f'{utils.TERM_DETAILS["conversion"]["label"]} ({utils.TERM_DETAILS["conversion"]["unit"]})'
    elif term in utils.TERM_DETAILS["generation_dissipation"]["terms"]:
        cmap = "cmo.tarn"
        norm = colors.TwoSlopeNorm(vcenter=0, vmin=-absmax, vmax=absmax)
        title = (
            f'{utils.TERM_DETAILS["generation_dissipation"]["label"]} '
            f'({utils.TERM_DETAILS["generation_dissipation"]["unit"]})'
        )
    return cmap, norm, title


def finalize_plot(ax, term):
    locator = mdates.AutoDateLocator(minticks=5, maxticks=12)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    ax.tick_params(axis="x", labelrotation=20)
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    ax.set_title(term, fontsize=20)


def plot_hovmoller(results_file, figures_directory, app_logger=False):

    read_results(results_file, app_logger)

    results_directory = os.path.dirname(results_file)
    dict_vertical = get_data_vertical_levels(results_directory)

    figures_subdirectory = os.path.join(figures_directory, "boxplots")
    os.makedirs(figures_subdirectory, exist_ok=True)

    figures_directory = os.path.join(results_directory, "Figures")
    figures_subdirectory = os.path.join(figures_directory, "hovmollers")
    os.makedirs(figures_subdirectory, exist_ok=True)

    (
        app_logger.info("Creating hovmoller diagrams")
        if app_logger
        else print("Creating hovmoller diagrams")
    )
    _plotter(dict_vertical, figures_subdirectory, app_logger)
    (
        app_logger.info("Hovmoller diagrams created")
        if app_logger
        else print("Hovmoller diagrams created")
    )


if __name__ == "__main__":

    # Test for Reg1-Representative_fixed
    results_file = "samples/Reg1-Representative_NCEP-R2_fixed/Reg1-Representative_NCEP-R2_fixed_results.csv"
    figures_directory = "samples/Reg1-Representative_fixed/Figures/"
    plot_hovmoller(results_file, figures_directory)

    # Test for Catarina_NCEP-R2_fixed
    results_file = "samples/Catarina_NCEP-R2_fixed/Catarina_NCEP-R2_fixed_results.csv"
    figures_directory = "samples/Catarina_NCEP-R2_fixed/Figures/"
    plot_hovmoller(results_file, figures_directory)

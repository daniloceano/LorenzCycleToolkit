# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    term_output.py                                     :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2026/08/18 10:00:00 by daniloceano       #+#    #+#              #
#    Updated: 2026/08/18 10:00:00 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
This script defines the TermOutputMixin object. It holds the three policies
that every Lorenz Energy Cycle term class shares: the NaN policy, the output
unit of the term and the per-level CSV layout. Term classes in ``src/analysis``
inherit it so that those policies have a single definition.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import csv
import os

from .nan_handling import convert_units, handle_nans


def _first_row(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle))


def header_levels(path):
    """Pressure levels (Pa) listed in the header of a per-level CSV."""
    return [float(v) for v in _first_row(path)[1:]]


def align_to_header(frame, path, time_name):
    """
    Put each value of a per-level frame under the header column of its own level (E0-21a).

    Args:
        frame (pd.DataFrame): One column per pressure level (Pa), one row per time.
        path (str): The per-level CSV the frame is appended to.
        time_name (str): Name of the time coordinate, the first cell of the
            header written by the frameworks.

    Returns:
        pd.DataFrame: The frame with the columns of the file header, in its
        order; a header level absent from the profile becomes an empty cell.
        A file without that header (not created yet, or created by appending
        rows only) is not aligned: the frame is returned as it is. The
        frameworks always create the files, with their headers, first.

    Raises:
        ValueError: If the profile has a level that is not in the header, so
            that writing it would put values under the wrong levels.
    """
    if not os.path.exists(path) or _first_row(path)[:1] != [time_name]:
        return frame
    header = header_levels(path)
    unknown = [float(p) for p in frame.columns if float(p) not in header]
    if unknown:
        raise ValueError(
            f"{path}: the profile has level(s) {unknown} that are not in the file "
            f"header {header}; writing it would put values under the wrong levels"
        )
    return frame.reindex(columns=header)


class TermOutputMixin:
    """
    Mixin providing the shared output policy of the Lorenz Energy Cycle terms.

    Attributes:
        OUTPUT_UNITS (str): Units the column integral of the term is converted
            to. Subclasses override it when their term is not a flux of energy
            per unit area (for example ``J/m^2`` for the energy contents).

    Methods:
        _handle_nans: Applies the shared vertical NaN policy to an integrand.
        _convert_units: Converts a term to ``OUTPUT_UNITS``.
        _save_vertical_levels: Appends a vertical profile to its per-level CSV.

    The host class must define ``VerticalCoordIndexer``, ``TimeName``,
    ``method``, ``app_logger`` and ``results_subdirectory_vertical_levels``.
    """

    OUTPUT_UNITS = "W/m^2"

    def _handle_nans(self, function, variable_name=""):
        """
        Apply the shared vertical NaN policy to an integrand.

        Args:
            function (xr.DataArray): The integrand.
            variable_name (str): Label used in the log messages.

        Returns:
            xr.DataArray: The integrand with interior vertical gaps
            interpolated. Pressure levels are never dropped here: the vertical
            control volume is fixed once for the whole budget in BoxData.
        """
        return handle_nans(
            function,
            self.VerticalCoordIndexer,
            variable_name=variable_name,
            app_logger=self.app_logger,
        )

    def _convert_units(self, function, variable_name, target=None):
        """
        Convert a term to its output units.

        Args:
            function (xr.DataArray): The quantified term.
            variable_name (str): Label used in the error message.
            target (str, optional): Target units. Defaults to ``OUTPUT_UNITS``.

        Returns:
            xr.DataArray: The term in the target units.

        Raises:
            ValueError: If the conversion is not dimensionally possible.
        """
        return convert_units(
            function,
            target if target is not None else self.OUTPUT_UNITS,
            variable_name,
            app_logger=self.app_logger,
        )

    def _save_vertical_levels(self, function, variable_name):
        """
        Append a vertical profile to ``{variable_name}_{level}.csv``.

        Each value is written under the header column of its own level
        (``align_to_header``); the time-invariant profiles of the fixed
        framework carry their own row of levels.

        Args:
            function (xr.DataArray): The vertical profile to archive.
            variable_name (str): Term name used for the file name and column.

        Returns:
            None
        """
        path = f"{self.results_subdirectory_vertical_levels}/{variable_name}_{self.VerticalCoordIndexer}.csv"
        df = function.to_dataframe(name=variable_name)
        df.reset_index(inplace=True)

        if self.method == "fixed":
            if self.TimeName not in function.dims:
                df = df.T
            else:
                df = df.pivot(index=self.TimeName, columns=self.VerticalCoordIndexer, values=variable_name)
                df = align_to_header(df, path, self.TimeName)
        else:
            df.set_index(self.TimeName, inplace=True)
            df.index = df.index.strftime("%Y-%m-%d %H:%M:%S")
            df = df.pivot(columns=self.VerticalCoordIndexer, values=variable_name)
            df.columns.name = None
            df = align_to_header(df, path, self.TimeName)

        df.to_csv(path, mode="a", header=None)

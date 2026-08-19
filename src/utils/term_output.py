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

from .nan_handling import convert_units, handle_nans


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

        Args:
            function (xr.DataArray): The vertical profile to archive.
            variable_name (str): Term name used for the file name and column.

        Returns:
            None
        """
        df = function.to_dataframe(name=variable_name)
        df.reset_index(inplace=True)

        if self.method == "fixed":
            if self.TimeName not in function.dims:
                df = df.T
            else:
                df = df.pivot(index=self.TimeName, columns=self.VerticalCoordIndexer)
        else:
            df.set_index(self.TimeName, inplace=True)
            df.index = df.index.strftime("%Y-%m-%d %H:%M:%S")
            df = df.pivot(columns=self.VerticalCoordIndexer, values=variable_name)
            df.columns.name = None

        df.to_csv(
            f"{self.results_subdirectory_vertical_levels}/{variable_name}_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None,
        )

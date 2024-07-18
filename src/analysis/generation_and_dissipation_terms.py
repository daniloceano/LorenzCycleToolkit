# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    generation_and_dissipation_terms.py                :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2020/11/26 18:27:40 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/18 00:33:23 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
This script defines the GenerationDissipationTerms object. It uses the MetData object as an
input. The built-in functions use the input data to compute the following
partioned energy contents of the Lorenz Energy Cycle.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import logging

import numpy as np
from metpy.constants import Cp_d, g
from metpy.units import units

from ..utils.box_data import BoxData
from ..utils.calc_averages import CalcAreaAverage


class GenerationDissipationTerms:
    """
    Class to compute generation and dissipation terms of the Lorenz Energy Cycle.

    Attributes:
        method (str): The computation method used ('fixed', 'track', or 'choose').
        box_obj (BoxData): The BoxData object containing meteorological data.

    Methods:
        calc_gz: Computes Generation of Zonal Available Potential Energy (Gz).
        calc_ge: Computes Generation of Eddy Available Potential Energy (Ge).
        calc_dz: Computes Dissipation of Zonal Kinetic Energy (Dz).
        calc_de: Computes Dissipation of Eddy Kinetic Energy (De).

    Source for formulas used here:
        Brennan, F. E., & Vincent, D. G. (1980).
        Zonal and Eddy Components of the Synoptic-Scale Energy Budget
        during Intensification of Hurricane Carmen (1974),
        Monthly Weather Review, 108(7), 954-965. Retrieved Jan 25, 2022, from:
        https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_0954_zaecot_2_0_co_2.xml
    """

    def __init__(self, box_obj: BoxData, method: str, app_logger: logging.Logger):
        """Initialize the EnergyContents object with a BoxData object and a method."""
        self._initialize_attributes(box_obj, method, app_logger)

    def _initialize_attributes(self, box_obj, method, app_logger):
        """Helper method to initialize attributes from the BoxData object."""

        # Operational attributes
        self.method = method
        self.box_obj = box_obj
        self.results_subdirectory = box_obj.results_subdirectory
        self.results_subdirectory_vertical_levels = box_obj.results_subdirectory_vertical_levels
        self.app_logger = app_logger

        # Initialize spatial and temporal attributes
        self.LonIndexer = box_obj.LonIndexer
        self.LatIndexer = box_obj.LatIndexer
        self.VerticalCoordIndexer = box_obj.VerticalCoordIndexer
        self.TimeName = box_obj.TimeName

        # Initialize lengths for averaging
        self.xlength = box_obj.xlength
        self.ylength = box_obj.ylength
        self.method = method
        self.box_obj = box_obj

        # Initialize pressure data
        self.PressureData = box_obj.PressureData

        # Initialize attributes related to temperature
        self.tair = box_obj.tair
        self.tair_AE = box_obj.tair_AE
        self.tair_ZE = box_obj.tair_ZE

        # Initialize attributes related to wind components
        self.u = box_obj.u
        self.u_ZA = box_obj.u_ZA
        self.u_ZE = box_obj.u_ZE
        self.v = box_obj.v
        self.v_ZA = box_obj.v_ZA
        self.v_ZE = box_obj.v_ZE

        # Initialize attributes related to wind stress
        self.ust_ZA = box_obj.ust_ZA
        self.ust_ZE = box_obj.ust_ZE
        self.vst_ZE = box_obj.vst_ZE
        self.vst_ZA = box_obj.vst_ZA

        # Initialize attributes related to vertical velocity
        self.omega = box_obj.omega

        # Initialize attributes related to static stability
        self.sigma_AA = box_obj.sigma_AA
        self.sigma_AA = box_obj.sigma_AA

        # Initialize attributes related to the adiabatic heating term
        self.Q = box_obj.Q
        self.Q_ZA = box_obj.Q_ZA
        self.Q_AA = box_obj.Q_AA
        self.Q_ZE = box_obj.Q_ZE
        self.Q_AE = box_obj.Q_AE

    def calc_gz(self):
        """Computes Generation of Zonal Available Potential Energy (Gz)."""
        self.app_logger.debug("Computing Gz...")
        term = self.Q_AE * self.tair_AE
        function = CalcAreaAverage(term, self.ylength) / (Cp_d * self.sigma_AA)
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Gz")
        Gz = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
        )
        self._convert_units(Gz, "Gz")
        self.app_logger.debug("Ok.")
        return Gz

    def calc_ge(self):
        """Computes Generation of Eddy Available Potential Energy (Ge)."""
        self.app_logger.debug("Computing Ge...")
        term = self.Q_ZE * self.tair_ZE
        function = CalcAreaAverage(term, self.ylength, xlength=self.xlength) / (
            Cp_d * self.sigma_AA
        )
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Ge")
        Ge = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
        )
        self._convert_units(Ge, "Ge")
        self.app_logger.debug("Ok.")
        return Ge

    def calc_dz(self):
        """
        Computes Dissipation of Zonal Kinetic Energy (Dz).

        Note: Still needs to be fully implemented and tested
        """
        self.app_logger.debug("Computing Dz...")
        # Here we will use only the lowest vertical level
        term = (self.u_ZA.isel({self.VerticalCoordIndexer: 0}) * self.ust_ZA) + (
            self.v_ZA.isel({self.VerticalCoordIndexer: 0}) * self.vst_ZA
        )
        function = CalcAreaAverage(term, self.ylength) / g
        self._save_vertical_levels(function, "Dz")
        Dz = units.Pa * function
        self._convert_units(Dz, "Dz")
        self.app_logger.debug("Ok.")
        return Dz

    def calc_de(self):
        """
        Computes Dissipation of Eddy Kinetic Energy (De).

        Note: Still needs to be fully implemented and tested
        """
        self.app_logger.debug("Computing De...")
        # Here we will use only the lowest vertical level
        term = (self.u_ZE.isel({self.VerticalCoordIndexer: 0}) * self.ust_ZE) + (
            self.v_ZE.isel({self.VerticalCoordIndexer: 0}) * self.vst_ZE
        )
        function = CalcAreaAverage(term, self.ylength) / g
        self._save_vertical_levels(function, "De")
        De = units.Pa * function
        self._convert_units(De, "De")
        self.app_logger.debug("Ok.")
        return De

    def _handle_nans(self, function):
        """
        If there are any, interpolate them and drop any remaining NaN values.
        If there are still NaN values after interpolation, drop the dimensions that contain NaN values.

        Parameters:
            function (np.ndarray): The function to handle NaN values for.

        Returns:
            None
        """
        if np.isnan(function).any():
            function = (
                function.interpolate_na(dim=self.VerticalCoordIndexer)
                * function.metpy.units
            )
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)
        return function

    def _convert_units(self, function, variable_name):
        """
        Converts the units of a given function to 'J/m^2'.

        Parameters:
            function (type): The function to be converted.
            variable_name (type): The name of the variable associated with the function.

        Raises:
            ValueError: If there is a unit error in the given variable.

        Returns:
            None
        """
        try:
            function = function.metpy.convert_units("W/m^2")
        except ValueError as e:
            raise ValueError(f"Unit error in {variable_name}") from e
        return function

    def _save_vertical_levels(self, function, variable_name):
        """Save computed energy data to a CSV file."""
        df = function.to_dataframe(name=variable_name)
        df.reset_index(inplace=True)
        if self.method == "fixed":
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

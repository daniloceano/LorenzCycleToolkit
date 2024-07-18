# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    conversion_terms.py                                :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/01/31 20:15:59 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/18 00:20:21 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
This script defines the ConversionTerms object. It uses the MetData object as
an input. The built-in functions use the input data to compute the following
energy conversion terms of the Lorenz Energy Cycle.

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
from metpy.constants import Rd, Re, g

from ..utils.box_data import BoxData
from ..utils.calc_averages import CalcAreaAverage


class ConversionTerms:
    """
    Class to compute energy conversion terms of the Lorenz Energy Cycle.

    Attributes:
        method (str): The computation method used ('fixed', 'track', or 'choose').
        box_obj (BoxData): The BoxData object containing meteorological data.

    Methods:
        calc_ce: Computes the eedy energy conversion term (CE).
        calc_cz: Computes the zonal energy conversion term (CZ).
        calc_ca: Computes the available potential energy conversion term (CA).
        calc_ck: Computes the kinetic energy conversion term (CK).

    Source for formulas used here:
        Brennan, F. E., & Vincent, D. G. (1980).
        Zonal and Eddy Components of the Synoptic-Scale Energy Budget
        during Intensification of Hurricane Carmen (1974),
        Monthly Weather Review, 108(7), 954-965. Retrieved Jan 25, 2022, from:
        https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_0954_zaecot_2_0_co_2.xml
    """

    def __init__(self, box_obj: BoxData, method: str, app_logger: logging.Logger):
        """Initialize the ConversionTerms object with a BoxData object and a method."""
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
        self.tan_lats = np.tan(box_obj.tair["rlats"])
        self.VerticalCoordIndexer = box_obj.VerticalCoordIndexer
        self.TimeName = box_obj.TimeName

        # Initialize lengths for averaging
        self.xlength = box_obj.xlength
        self.ylength = box_obj.ylength

        # Initialize attributes related to temperature
        self.tair_AE = box_obj.tair_AE
        self.tair_ZE = box_obj.tair_ZE

        # Initialize attributes related to wind components
        self.u_ZA = box_obj.u_ZA
        self.u_ZE = box_obj.u_ZE
        self.v_ZA = box_obj.v_ZA
        self.v_ZE = box_obj.v_ZE

        # Initialize attirbutes related to vertical velocity
        self.omega_ZE = box_obj.omega_ZE
        self.omega_AE = box_obj.omega_AE

        # Initialize the static stability parameter
        self.sigma_AA = box_obj.sigma_AA

        # Initialize the pressure data
        self.PressureData = box_obj.PressureData

    def calc_ca(self):
        """
        Computes conversion between the two available potential energy forms (AZ and AE).

        Note: on the first term of the integral, on Brennan et al. (1980), it is missing
        a 2 in the multiplication Re * self.sigma_AA. This is confirmed by looking to the
        paper by Muench (1965).
        """
        self.app_logger.debug("Calculating CA...")

        # First term of the integral
        DelPhi_tairAE = (self.tair_AE * self.tair_AE["coslats"]).differentiate("rlats")
        term1 = (self.v_ZE * self.tair_ZE * DelPhi_tairAE) / (2 * Re * self.sigma_AA)
        term1 = CalcAreaAverage(term1, self.ylength, xlength=self.xlength)
        self._save_vertical_levels(term1, "Ca_1")

        # Second term of the integral
        DelPres_tairAE = (self.tair_AE).differentiate(
            self.VerticalCoordIndexer
        ) / self.PressureData.metpy.units
        term2 = (self.omega_ZE * self.tair_ZE) * DelPres_tairAE
        term2 = (
            CalcAreaAverage(term2, self.ylength, xlength=self.xlength) / self.sigma_AA
        )
        self._save_vertical_levels(term2, "Ca_2")

        # Process the integral and save the result
        function = -(term1 + term2)
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Ca")
        Ca = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
        )
        Ca = self._convert_units(Ca, "Ca")

        self.app_logger.debug("Done.")
        return Ca

    def calc_ce(self):
        """Computes conversion between the two eddy energy forms (AE and KE)."""
        self.app_logger.debug("Calculating CE...")

        # First term of the integral
        term1 = Rd / (self.PressureData * g)
        self._save_vertical_levels(term1, "Ce_1")

        # Second term of the integral
        omega_tair_product = self.omega_ZE * self.tair_ZE
        term2 = CalcAreaAverage(omega_tair_product, self.ylength, xlength=self.xlength)
        self._save_vertical_levels(term2, "Ce_2")

        # Process the integral and save the result
        function = -(term1 * term2)
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Ce")
        Ce = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
        )
        Ce = self._convert_units(Ce, "Ce")

        self.app_logger.debug("Done.")
        return Ce

    def calc_cz(self):
        """Computes conversion between the two zonal energy forms (ZE and KE)."""
        self.app_logger.debug("Calculating CZ...")

        # First term of the integral
        term1 = Rd / (self.PressureData * g)
        self._save_vertical_levels(term1, "Cz_1")

        # Second term of the integral
        omega_tair_product = self.omega_AE * self.tair_AE
        term2 = CalcAreaAverage(omega_tair_product, self.ylength)
        self._save_vertical_levels(term2, "Cz_2")

        # Process the integral and save the result
        function = -(term1 * term2)
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Cz")
        Cz = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
        )
        Cz = self._convert_units(Cz, "Cz")

        self.app_logger.debug("Done.")
        return Cz

    def calc_ck(self):
        """Computes conversion between the two eddy kinetic energy forms (KE and KZ)."""
        self.app_logger.debug("Calculating CK...")

        # First term of the integral
        DelPhi_uZA_cosphi = (self.u_ZA / self.u_ZA["coslats"]).differentiate("rlats")
        term1 = (self.u_ZE["coslats"] * self.u_ZE * self.v_ZE / Re) * DelPhi_uZA_cosphi
        term1 = CalcAreaAverage(term1, self.ylength, xlength=self.xlength)
        self._save_vertical_levels(term1, "Ck_1")

        # Second term of the integral
        DelPhi_vZA = (self.v_ZA).differentiate("rlats")
        term2 = ((self.v_ZE**2) / Re) * DelPhi_vZA
        term2 = CalcAreaAverage(term2, self.ylength, xlength=self.xlength)
        self._save_vertical_levels(term2, "Ck_2")

        # Third term of the integral
        term3 = (self.tan_lats * (self.u_ZE**2) * self.v_ZA) / Re
        term3 = CalcAreaAverage(term3, self.ylength, xlength=self.xlength)
        self._save_vertical_levels(term3, "Ck_3")

        # Fourth term of the integral
        DelPres_uZAp = (
            self.u_ZA.differentiate(self.VerticalCoordIndexer)
            / self.PressureData.metpy.units
        )
        term4 = self.omega_ZE * self.u_ZE * DelPres_uZAp
        term4 = CalcAreaAverage(term4, self.ylength, xlength=self.xlength)
        self._save_vertical_levels(term4, "Ck_4")

        # Fifth term of the integral
        DelPres_vZAp = (
            self.u_ZA.differentiate(self.VerticalCoordIndexer)
            / self.PressureData.metpy.units
        )
        term5 = self.omega_ZE * self.v_ZE * DelPres_vZAp
        term5 = CalcAreaAverage(term5, self.ylength, xlength=self.xlength)
        self._save_vertical_levels(term5, "Ck_5")

        # Process the integral and save the result
        function = term1 + term2 + term3 + term4 + term5
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Ck")
        Ck = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            / g
        )
        Ck = self._convert_units(Ck, "Ck")

        self.app_logger.debug("Done.")
        return Ck

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

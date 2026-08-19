# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    energy_contents.py                                 :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/01/31 20:15:59 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/18 00:19:54 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
This script defines the EnergyContent object. It uses the MetData object as an
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

from metpy.constants import g

from ..utils.box_data import BoxData
from ..utils.calc_averages import CalcAreaAverage
from ..utils.term_output import TermOutputMixin


class EnergyContents(TermOutputMixin):
    """
    Class to compute partitioned energy contents of the Lorenz Energy Cycle.

    Attributes:
        method (str): The computation method used ('fixed', 'track', or 'choose').
        box_obj (BoxData): The BoxData object containing meteorological data.

    Methods:
        calc_az: Computes Zonal Available Potential Energy (Az).
        calc_ae: Computes Eddy Available Potential Energy (Ae).
        calc_kz: Computes Zonal Kinetic Energy (Kz).
        calc_ke: Computes Eddy Kinetic Energy (Ke).

    Source for formulas used here:
        Brennan, F. E., & Vincent, D. G. (1980).
        Zonal and Eddy Components of the Synoptic-Scale Energy Budget
        during Intensification of Hurricane Carmen (1974),
        Monthly Weather Review, 108(7), 954-965. Retrieved Jan 25, 2022, from:
        https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_0954_zaecot_2_0_co_2.xml
    """

    # The energy contents are reservoirs, not fluxes: they are archived in
    # J/m^2, while every other term class keeps the mixin default of W/m^2.
    OUTPUT_UNITS = "J/m^2"

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

        # Initialize attributes related to temperature
        self.tair_AE = box_obj.tair_AE
        self.tair_ZE = box_obj.tair_ZE

        # Initialize attributes related to wind components
        self.u_ZA = box_obj.u_ZA
        self.u_ZE = box_obj.u_ZE
        self.v_ZA = box_obj.v_ZA
        self.v_ZE = box_obj.v_ZE

        # Initialize the static stability parameter
        self.sigma_AA = box_obj.sigma_AA

        # Initialize the pressure data
        self.PressureData = box_obj.PressureData

    def calc_az(self):
        """Computes Zonal Available Potential Energy (Az)."""
        self.app_logger.debug("Computing Az...")
        squared_tair = self.tair_AE**2
        function = CalcAreaAverage(squared_tair, self.ylength) / (2 * self.sigma_AA)
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Az")
        Az = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
        )
        Az = self._convert_units(Az, "Az")
        Az = Az.metpy.dequantify()
        self.app_logger.debug("Ok.")
        return Az

    def calc_ae(self):
        """Computes Eddy Available Potential Energy (Ae)."""
        self.app_logger.debug("Computing Ae...")
        squared_tair = self.tair_ZE**2
        function = CalcAreaAverage(squared_tair, self.ylength, xlength=self.xlength) / (
            2 * self.sigma_AA
        )
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Ae")
        Ae = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
        )
        Ae = self._convert_units(Ae, "Ae")
        Ae = Ae.metpy.dequantify()
        self.app_logger.debug("Ok.")
        return Ae

    def calc_kz(self):
        """Computes Zonal Kinetic Energy (Kz)."""
        self.app_logger.debug("Computing Kz...")
        function = CalcAreaAverage((self.u_ZA**2 + self.v_ZA**2), self.ylength)
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Kz")
        Kz = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            / (2 * g)
        )
        Kz = self._convert_units(Kz, "Kz")
        Kz = Kz.metpy.dequantify()
        self.app_logger.debug("Ok.")
        return Kz

    def calc_ke(self):
        """Computes Eddy Kinetic Energy (Ke)."""
        self.app_logger.debug("Computing Ke...")
        function = CalcAreaAverage(
            (self.u_ZE**2 + self.v_ZE**2), self.ylength, self.xlength
        )
        function = self._handle_nans(function)
        self._save_vertical_levels(function, "Ke")
        Ke = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            / (2 * g)
        )
        Ke = self._convert_units(Ke, "Ke")
        Ke = Ke.metpy.dequantify()
        self.app_logger.debug("Ok.")
        return Ke

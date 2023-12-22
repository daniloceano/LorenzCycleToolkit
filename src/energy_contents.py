# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    energy_contents.py                                 :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/01/31 20:15:59 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/22 13:51:19 by daniloceano      ###   ########.fr        #
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
import numpy as np
import pandas as pd
from metpy.constants import g
from metpy.units import units
from calc_averages import CalcAreaAverage
from box_data import BoxData

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class EnergyContents:
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
    
    def __init__(self, box_obj: BoxData, method: str):
        """Initialize the EnergyContents object with a BoxData object and a method."""
        self._initialize_attributes(box_obj, method)

    def _initialize_attributes(self, box_obj, method):
        """Helper method to initialize attributes from the BoxData object."""

        # Operational attributes
        self.method = method
        self.box_obj = box_obj
        self.output_dir = box_obj.output_dir
        
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
        squared_tair = self.tair_AE ** 2
        function = CalcAreaAverage(squared_tair, self.ylength) / (2 * self.sigma_AA)
        function = self._handle_nans(function)
        Az = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units
        self._convert_units(Az, 'Az')
        self._save_vertical_levels(Az, 'Az')
        Az = Az.metpy.dequantify()
        return Az
    
    def calc_ae(self):
        """Computes Eddy Available Potential Energy (Ae)."""
        squared_tair = self.tair_ZE ** 2
        function = CalcAreaAverage(squared_tair, self.ylength, xlength=self.xlength) / (2 * self.sigma_AA)
        function = self._handle_nans(function)
        Ae = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units
        self._convert_units(Ae, 'Ae')
        self._save_vertical_levels(Ae, 'Ae')
        Ae = Ae.metpy.dequantify()
        return Ae

    def calc_kz(self):
        """Computes Zonal Kinetic Energy (Kz)."""
        function = CalcAreaAverage((self.u_ZA ** 2 + self.v_ZA ** 2), self.ylength)
        function = self._handle_nans(function)
        Kz = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units / (2 * g)
        self._convert_units(Kz, 'Kz')
        self._save_vertical_levels(Kz, 'Kz')
        Kz = Kz.metpy.dequantify()
        return Kz
    
    def calc_ke(self):
        """Computes Eddy Kinetic Energy (Ke)."""
        function = CalcAreaAverage((self.u_ZE ** 2 + self.v_ZE ** 2), self.ylength, self.xlength)
        function = self._handle_nans(function)
        Ke = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units / (2 * g)
        self._convert_units(Ke, 'Ke')
        self._save_vertical_levels(Ke, 'Ke')
        Ke = Ke.metpy.dequantify()
        return Ke

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
            function = function.metpy.convert_units('J/m^2')
        except ValueError as e:
            raise ValueError(f'Unit error in {variable_name}') from e
        return function
        
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
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)
        return function
        
    def _save_vertical_levels(self, function, variable_name):
        """Save computed energy data to a CSV file."""
        if self.method == 'fixed':
            df = function.to_dataframe(name='Az').unstack()
            logging.info(f'Computed {variable_name}')
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()

        df.to_csv(f"{self.output_dir}/Az_{self.VerticalCoordIndexer}.csv",
                  mode="a", header=None)

        df.to_csv(f"{self.output_dir}/{variable_name}_{self.VerticalCoordIndexer}.csv", mode="a", header=None)
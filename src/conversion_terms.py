# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    ConversionTerms.py                                 :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/01/31 20:15:59 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/22 11:53:41 by daniloceano      ###   ########.fr        #
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
import pandas as pd
from metpy.constants import g
from metpy.constants import Rd
from metpy.constants import Re
from metpy.units import units
from box_data import BoxData
from calc_averages import CalcAreaAverage

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
    
    def __init__(self, box_obj: BoxData, method: str):
        """Initialize the ConversionTerms object with a BoxData object and a method."""
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
        # First term of the integral
        DelPhi_tairAE = (self.tair_AE * self.tair_AE["coslats"]).differentiate("rlats")
        term1 = (self.v_ZE * self.tair_ZE * DelPhi_tairAE) / (2 * Re * self.sigma_AA)
        term1 = CalcAreaAverage(term1, self.ylength, xlength=self.xlength)

        # Second term of the integral
        DelPres_tairAE = (self.tair_AE).differentiate(self.VerticalCoordIndexer) / self.PressureData.metpy.units
        term2 = (self.omega_ZE * self.tair_ZE) * DelPres_tairAE
        term2 = CalcAreaAverage(term2, self.ylength, xlength=self.xlength) / self.sigma_AA

        # Process the integral and save the result
        function = term1 + term2
        function = self._handle_nans(function)
        Ca = - function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units
        self._convert_units(Ca, 'Ca')
        self._save_vertical_levels(Ca, 'Ca')

        return Ca
        
    def calc_ce(self):
        """Computes conversion between the two eddy energy forms (AE and KE)."""
        # First term of the integral
        term1 = Rd / (self.PressureData * g)
        omega_tair_product = self.omega_ZE * self.tair_ZE

        # Second term of the integral
        term2 = - CalcAreaAverage(omega_tair_product, self.ylength, xlength=self.xlength)

        # Process the integral and save the result
        function = term1 * term2
        function = self._handle_nans(function)
        Ce = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units
        self._convert_units(Ce, 'Ce')
        self._save_vertical_levels(Ce, 'Ce')

        return Ce
    
    def calc_cz(self):
        """Computes conversion between the two zonal energy forms (ZE and KE)."""
        # First term of the integral
        term1 = Rd / (self.PressureData * g)
        omega_tair_product = self.omega_AE * self.tair_AE

        # Second term of the integral
        term2 = CalcAreaAverage(omega_tair_product, self.ylength)

        # Process the integral and save the result
        function = term1 * term2
        function = self._handle_nans(function)
        Cz = - function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units
        self._convert_units(Cz, 'Cz')
        self._save_vertical_levels(Cz, 'Cz')

        return Cz
        
    def calc_ck(self):
        """Computes conversion between the two eddy kinetic energy forms (KE and KZ)."""
        # First term of the integral
        DelPhi_uZA_cosphi = ((self.u_ZA / self.u_ZA["coslats"]) * self.u_ZA["coslats"]).differentiate("rlats")
        term1 = (self.u_ZE["coslats"] * self.u_ZE * self.v_ZE / Re) * DelPhi_uZA_cosphi
        term1 = CalcAreaAverage(term1, self.ylength, xlength=self.xlength)

        # Second term of the integral
        DelPhi_vZA = (self.v_ZA * self.v_ZA["coslats"]).differentiate("rlats")
        term2 = ((self.v_ZE ** 2) / Re) * DelPhi_vZA
        term2 = CalcAreaAverage(term2, self.ylength, xlength=self.xlength)

        # Third term of the integral
        term3 = (self.tan_lats * (self.u_ZE ** 2) * self.v_ZA) / Re
        term3 = CalcAreaAverage(term3, self.ylength, xlength=self.xlength)

        # Fourth term of the integral
        DelPres_uZAp = self.u_ZA.differentiate(self.VerticalCoordIndexer) / self.PressureData.metpy.units
        term4 = self.omega_ZE * self.u_ZE * DelPres_uZAp
        term4 = CalcAreaAverage(term4, self.ylength, xlength=self.xlength)

        # Fifth term of the integral
        DelPres_vZAp = self.u_ZA.differentiate(self.VerticalCoordIndexer) / self.PressureData.metpy.units
        term5 = self.omega_ZE * self.v_ZE * DelPres_vZAp
        term5 = CalcAreaAverage(term5, self.ylength, xlength=self.xlength)

        # Process the integral and save the result
        function = term1 + term2 + term3 + term4 + term5
        function = self._handle_nans(function)
        Ck = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units / g
        self._convert_units(Ck, 'Ck')
        self._save_vertical_levels(Ck, 'Ck')

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
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
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
            function = function.metpy.convert_units('W/m^2')
        except ValueError as e:
            raise ValueError(f'Unit error in {variable_name}') from e
    
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
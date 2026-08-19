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
from ..utils.term_output import TermOutputMixin


class ConversionTerms(TermOutputMixin):
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
        calc_c_overturning: Diagnoses domain-mean overturning pressure work.

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
        self.tair_AA = box_obj.tair_AA
        self.tair_AE = box_obj.tair_AE
        self.tair_ZE = box_obj.tair_ZE

        # Initialize attributes related to wind components
        self.u_ZA = box_obj.u_ZA
        self.u_ZE = box_obj.u_ZE
        self.v_ZA = box_obj.v_ZA
        self.v_ZE = box_obj.v_ZE

        # Initialize attirbutes related to vertical velocity
        self.omega_AA = box_obj.omega_AA
        self.omega_ZE = box_obj.omega_ZE
        self.omega_AE = box_obj.omega_AE

        # Initialize the static stability parameter
        self.sigma_AA = box_obj.sigma_AA

        # Initialize the pressure data
        self.PressureData = box_obj.PressureData

    def calc_ca(self):
        """
        Computes conversion between the two available potential energy forms (AZ and AE).

        .. math::
            C_A = -\\int_{p_t}^{p_b}\\left(
                  \\frac{1}{a\\sigma}\\overline{v'T'\\frac{\\partial T^*}{\\partial\\phi}}
                + \\frac{1}{\\sigma}\\overline{\\omega'T'\\frac{\\partial T^*}{\\partial p}}
                  \\right)dp

        Sources
        -------
        Latitude derivative: acts on ``T*`` alone.  Muench (1965),
        Norquist et al. (1977), Brennan and Vincent (1980), Michaelides (1987)
        and Michaelides et al. (1999) all print it this way.

        Leading factor: ``1/(a sigma)``, with no factor ``1/2``.  Muench (1965)
        and Michaelides (1987) print ``1/(2 a sigma)``; Norquist et al. (1977),
        who state that they follow Muench's notation, together with Brennan and
        Vincent (1980) and Michaelides et al. (1999), print it without.  The
        form used here is the one consistent with the reservoir definition
        adopted by the toolkit, ``A_E = int T'^2 / (2 sigma) dp``:
        differentiating that expression supplies a chain-rule factor 2 which
        cancels the 2 in the denominator.

        Vertical term: the reduced form ``(1/sigma) dT*/dp`` of Brennan and
        Vincent (1980), rather than the fuller
        ``p^-k d/dp (T* p^k / sigma)`` of Muench (1965) and Michaelides (1987).
        """
        self.app_logger.debug("Calculating CA...")

        # First term of the integral
        DelPhi_tairAE = self.tair_AE.differentiate("rlats")
        term1 = (self.v_ZE * self.tair_ZE * DelPhi_tairAE) / (Re * self.sigma_AA)
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

    def calc_c_overturning(self):
        r"""Diagnose domain-mean overturning pressure work.

        .. math::
            C_{\mathrm{overturning}} = -\int_{p_t}^{p_b}
                \overline{\omega}\,\alpha\,\frac{dp}{g},
            \qquad \alpha = \frac{R_d\overline{T}}{p}.

        This is an overturning-strength diagnostic, not a missing conversion
        in the Lorenz-cycle budget.  In the exact pressure-work identity it
        cancels the ``Phi_bar * omega_bar`` part of the top/bottom
        geopotential flux.  It is therefore exported for interpretation but
        correctly remains outside ``RGz`` and ``RKz``.  With the toolkit sign
        convention, mean ascent (``omega < 0``) gives positive
        ``C_overturning``.
        """
        self.app_logger.debug("Calculating C_overturning...")

        alpha = Rd * self.tair_AA / self.PressureData
        function = -(self.omega_AA * alpha) / g
        function = self._handle_nans(function, "C_overturning")
        self._save_vertical_levels(function, "C_overturning")
        result = (
            function.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
        )
        result = self._convert_units(result, "C_overturning")

        self.app_logger.debug("Done.")
        return result

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
        """Computes conversion between the two eddy kinetic energy forms (KE and KZ).

        Sources
        -------
        Muench (1965), Norquist et al. (1977, Eq. 6), Brennan and Vincent
        (1980, p. 964) and Michaelides (1987, p. 24) for the five sub-terms.
        """
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

        # Fifth term of the integral: the eddy momentum flux acts on the
        # vertical shear of the zonal-mean MERIDIONAL wind, d[v]/dp.
        # Muench (1965), Norquist et al. (1977, Eq. 6), Brennan and Vincent
        # (1980, p. 964) and Michaelides (1987, p. 24) all print d[v]/dp;
        # Michaelides et al. (1999, Eq. 12) prints d[u]/dp.
        DelPres_vZAp = (
            self.v_ZA.differentiate(self.VerticalCoordIndexer)
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

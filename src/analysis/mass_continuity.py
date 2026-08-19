# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    mass_continuity.py                                 :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2026/08/18 10:00:00 by daniloceano       #+#    #+#              #
#    Updated: 2026/08/18 10:00:00 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
This script defines the MassContinuity object. It uses the BoxData object as an
input and computes the limited-area mass-continuity residual of the Lorenz
Energy Cycle,

    M = <D_EW> + <D_NS> + d(omega_bar)/dp,

which is identically zero in the continuum but not for independently archived
reanalysis winds and omega. M therefore sets the numerical noise floor of the
budget and is reported outside the residuals.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import logging

from metpy.constants import Re, g

from ..utils.box_data import BoxData
from ..utils.calc_averages import CalcAreaAverage
from ..utils.term_output import TermOutputMixin


class MassContinuity(TermOutputMixin):
    """
    Class to compute the limited-area mass-continuity residual of the Lorenz
    Energy Cycle.

    Attributes:
        method (str): The computation method used ('fixed', 'track', or 'choose').
        box_obj (BoxData): The BoxData object containing meteorological data.

    Methods:
        calc_mass_residual_profile: Computes the vertical profile M(p), in s^-1.
        calc_mass_residual: Computes the column integral int M dp/g, in
            kg m^-2 s^-1.

    Source for formulas used here:
        docs/source/math.rst, section "Mass-continuity diagnostic".

        Brennan, F. E., & Vincent, D. G. (1980).
        Zonal and Eddy Components of the Synoptic-Scale Energy Budget
        during Intensification of Hurricane Carmen (1974),
        Monthly Weather Review, 108(7), 954-965. Retrieved Jan 25, 2022, from:
        https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_0954_zaecot_2_0_co_2.xml
    """

    OUTPUT_UNITS = "kg/m^2/s"

    def __init__(self, box_obj: BoxData, method: str, app_logger: logging.Logger):
        """Initialize the MassContinuity object with a BoxData object and a method."""
        # Operational attributes
        self.method = method
        self.box_obj = box_obj
        self.app_logger = app_logger
        self.results_subdirectory_vertical_levels = (
            box_obj.results_subdirectory_vertical_levels
        )

        # Initialize spatial and temporal attributes
        self.LonIndexer = box_obj.LonIndexer
        self.LatIndexer = box_obj.LatIndexer
        self.TimeName = box_obj.TimeName
        self.VerticalCoordIndexer = box_obj.VerticalCoordIndexer
        self.PressureData = box_obj.PressureData

        # Initialize domain limits and lengths for averaging
        self.western_limit = box_obj.western_limit
        self.eastern_limit = box_obj.eastern_limit
        self.southern_limit = box_obj.southern_limit
        self.northern_limit = box_obj.northern_limit
        self.xlength = box_obj.xlength
        self.ylength = box_obj.ylength

        # Initialize attributes related to wind and omega
        self.u = box_obj.u
        self.v_ZA = box_obj.v_ZA
        self.omega_AA = box_obj.omega_AA

    def calc_mass_residual_profile(self):
        r"""
        Computes the mass-continuity residual profile M(p), in s\ :sup:`-1`.

        The two horizontally averaged divergences are evaluated in their
        exactly telescoped boundary forms,

        ``<D_EW> = int (u_e - u_w) dphi / (a Delta(lambda) int cos(phi) dphi)``
        ``<D_NS> = ([v]cos(phi))_n - ([v]cos(phi))_s / (a int cos(phi) dphi)``,

        which remain well defined at the poles for a global domain. Both are
        normalised by the trapezoidal measure ``int cos(phi) dphi`` used by
        :func:`src.utils.calc_averages.CalcAreaAverage`, so the three summands
        of M share one area-averaging convention. Normalising the first two by
        the analytic ``sin(phi_n) - sin(phi_s)`` instead would leave an
        O(dphi^2) mismatch against ``d(omega_bar)/dp`` -- the same order of
        magnitude as the inconsistency M exists to measure.

        Returns:
            xr.DataArray: The vertical profile of M, also saved to
            ``M_{VerticalCoordIndexer}.csv``.
        """
        self.app_logger.debug("Calculating mass-continuity residual M...")

        u_faces = self.u.sel(**{self.LonIndexer: self.eastern_limit}) - self.u.sel(
            **{self.LonIndexer: self.western_limit}
        )
        self.d_ew_area = CalcAreaAverage(
            u_faces / (Re * u_faces["coslats"] * self.xlength), self.ylength
        )

        cos_lats = self.v_ZA["coslats"]
        meridional_flux = self.v_ZA * cos_lats
        self.d_ns_area = (
            meridional_flux.sel(**{self.LatIndexer: self.northern_limit})
            - meridional_flux.sel(**{self.LatIndexer: self.southern_limit})
        ) / (Re * cos_lats.integrate("rlats"))

        self.domega_dp = self.omega_AA.differentiate(self.VerticalCoordIndexer) / (
            self.PressureData.metpy.units
        )
        residual = self._handle_nans(
            self.d_ew_area + self.d_ns_area + self.domega_dp, "M"
        )
        self.residual_profile = residual
        self._save_vertical_levels(residual, "M")
        return residual

    def calc_mass_residual(self):
        r"""
        Computes the conservative column residual, in kg m\ :sup:`-2` s\ :sup:`-1`.

        The vertical derivative is integrated as the exact boundary difference
        ``omega_bottom - omega_top``. This is the finite-volume counterpart of
        ``int d(omega_bar)/dp dp`` and avoids a second differentiation-then-
        quadrature error in the scalar output.

        Returns:
            xr.DataArray: The column integral ``int M dp / g``.
        """
        self.calc_mass_residual_profile()
        horizontal_column = (
            (self.d_ew_area + self.d_ns_area).integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
        )
        vertical_column = self.omega_AA.isel(
            **{self.VerticalCoordIndexer: -1}
        ) - self.omega_AA.isel(**{self.VerticalCoordIndexer: 0})
        column = (horizontal_column + vertical_column) / g
        result = self._convert_units(column, "M")
        self.app_logger.debug("Mass-continuity residual M done.")
        return result

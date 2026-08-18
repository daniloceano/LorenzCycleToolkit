# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    boundary_terms.py                                  :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/01/31 20:15:59 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/18 00:21:16 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
This script defines the CBoundaryTerms object. It uses the MetData object as
an input. The built-in functions use the input data to compute the following
boundary terms of the Lorenz Energy Cycle.

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
from metpy.constants import Re, g

from ..utils.box_data import BoxData
from ..utils.calc_averages import CalcAreaAverage, CalcZonalAverage
from ..utils.nan_handling import convert_units, handle_nans


class BoundaryTerms:
    """
    Class to compute boundary terms of the Lorenz Energy Cycle.

    Attributes:
        method (str): The computation method used ('fixed', 'track', or 'choose').
        box_obj (BoxData): The BoxData object containing meteorological data.

    Methods:
        calc_baz(self): flux of zonal available potential energy across the boundaries (BAZ).
        calc_bae(self): flux of eddy available potential energy across the boundaries (BAE).
        calc_bkz(self): flux of zonal kinetic energy across the boundaries (BKZ).
        calc_bke(self): flux of eddy kinetic energy across the boundaries (BKE).
        calc_boz(self): appearence of zonal kinetic energywith work produced at the boundaries (BΦZ).
        calc_boe(self): appearence of eddy kinetic energy with work produced at the boundaries (BΦE).

    Sources for formulas used here:
        Michaelides, S. C. (1987).
        Limited Area Energetics of Genoa Cyclogenesis,
        Monthly Weather Review, 115(1), 13-26. From:
        https://journals.ametsoc.org/view/journals/mwre/115/1/1520-0493_1987_115_0013_laeogc_2_0_co_2.xml

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
        self.app_logger = app_logger

        # Initialize spatial and temporal attributes
        self.LonIndexer = box_obj.LonIndexer
        self.LatIndexer = box_obj.LatIndexer
        self.western_limit = box_obj.western_limit
        self.eastern_limit = box_obj.eastern_limit
        self.southern_limit = box_obj.southern_limit
        self.northern_limit = box_obj.northern_limit
        self.VerticalCoordIndexer = box_obj.VerticalCoordIndexer
        self.TimeName = box_obj.TimeName

        # Initialize lengths for averaging
        self.xlength = box_obj.xlength
        self.ylength = box_obj.ylength

        # Initialize pressure data
        self.PressureData = box_obj.PressureData

        # Initialize attributes related to temperature
        self.tair_AE = box_obj.tair_AE
        self.tair_ZE = box_obj.tair_ZE

        # Initialize attributes related to wind components
        self.u = box_obj.u
        self.u_ZA = box_obj.u_ZA
        self.u_ZE = box_obj.u_ZE
        self.v = box_obj.v
        self.v_ZA = box_obj.v_ZA
        self.v_ZE = box_obj.v_ZE

        # Initialize attributes related to vertical velocity
        self.omega = box_obj.omega
        self.omega_ZE = box_obj.omega_ZE
        self.omega_AA = box_obj.omega_AA
        self.omega_ZA = box_obj.omega_ZA
        self.omega_AE = box_obj.omega_AE

        # Initialize static stability parameter
        self.sigma_AA = box_obj.sigma_AA

        # Initialize attributes related to geopotential
        self.geopt = box_obj.geopt
        self.geopt_ZA = box_obj.geopt_ZA
        self.geopt_ZE = box_obj.geopt_ZE
        self.geopt_AE = box_obj.geopt_AE

        # Initialize operand using the notation from Michaelides (1987)
        self.c1 = -1 / (Re * self.xlength * self.ylength)
        self.c2 = -1 / (Re * self.ylength)

    def calc_baz(self):
        """
        Computes the flux of Zonal Available Potential Energy across the boundaries (BAZ).

        Note:
            On Brennan et al. (1980), the area average of term3 is missing (Muench, 1965).
        """
        self.app_logger.debug("Calculating BAZ...")

        # First term
        term1 = (
            (2 * self.tair_AE * self.tair_ZE * self.u) + (self.tair_AE**2 * self.u)
        ) / (2 * self.sigma_AA)
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(
            **{self.LonIndexer: self.western_limit}
        )
        term1 = term1.integrate("rlats")
        term1 = self._handle_nans(term1)
        term1 = (
            term1.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c1
        )

        # Second term
        term2 = self.v_ZE * self.tair_ZE
        term2 = CalcZonalAverage(term2, self.xlength) * 2 * self.tair_AE
        term2 = (term2 + ((self.tair_AE**2) * self.v_ZA)) * self.tair_AE["coslats"]
        term2 = (
            term2.sel(**{self.LatIndexer: self.northern_limit})
            - term2.sel(**{self.LatIndexer: self.southern_limit})
        ) / (2 * self.sigma_AA)
        term2 = self._handle_nans(term2)
        term2 = (
            term2.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c2
        )

        # Third term
        term3a = 2 * self.omega_ZE * self.tair_ZE
        term3a = CalcZonalAverage(term3a, self.xlength) * self.tair_AE
        term3b = self.omega_ZA * self.tair_AE**2
        term3 = term3a + term3b
        term3 = self._handle_nans(term3)
        term3 = CalcAreaAverage(term3, self.ylength) / (2 * self.sigma_AA)
        term3 = term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(
            **{self.VerticalCoordIndexer: 0}
        )

        # Combine terms and save the result
        function = term1 + term2 - term3
        function = self._handle_nans(function)
        Baz = self._convert_units(function, "BAZ")

        self.app_logger.debug("Done.")
        return Baz

    def calc_bae(self):
        """Computes the flux of Eddy Available Potential Energy across the boundaries (BAE)."""
        self.app_logger.debug("Calculating BAE...")

        # First Integral
        term1 = self.u * (self.tair_ZE**2)
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(
            **{self.LonIndexer: self.western_limit}
        )
        term1 = (term1 / (2 * self.sigma_AA)).integrate("rlats")
        term1 = self._handle_nans(term1)
        term1 = (
            term1.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c1
        )

        # Second Integral
        term2 = (
            CalcZonalAverage(self.v * self.tair_ZE**2, self.xlength)
            * self.tair_ZE["coslats"]
        )
        term2 = term2 / (2 * self.sigma_AA)
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(
            **{self.LatIndexer: self.southern_limit}
        )
        term2 = self._handle_nans(term2)
        term2 = (
            term2.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c2
        )

        # Third Term
        term3 = (self.omega * self.tair_ZE**2) / (2 * self.sigma_AA)
        term3 = CalcAreaAverage(term3, self.ylength, xlength=self.xlength)
        term3 = self._handle_nans(term3)
        term3 = term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(
            **{self.VerticalCoordIndexer: 0}
        )

        # Combine terms and save the result
        function = term1 + term2 - term3
        function = self._handle_nans(function)
        Bae = self._convert_units(function, "BAE")

        self.app_logger.debug("Done.")
        return Bae

    def calc_bkz(self):
        """Computes the flux of Zonal Kinetic Energy across the boundaries (BKZ)."""
        self.app_logger.debug("Calculating BKZ...")

        # First term
        term1 = self.u * (self.u**2 + self.v**2 - self.u_ZE**2 - self.v_ZE**2)
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(
            **{self.LonIndexer: self.western_limit}
        )
        term1 = (term1 / (2 * g)).integrate("rlats")
        term1 = self._handle_nans(term1)
        term1 = (
            term1.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c1
        )

        # Second term
        term2 = (
            (self.u**2 + self.v**2 - self.u_ZE**2 - self.v_ZE**2)
            * self.v
            * self.v["coslats"]
        )
        term2 = CalcZonalAverage(term2, self.xlength)
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(
            **{self.LatIndexer: self.southern_limit}
        )
        term2 = self._handle_nans(term2)
        term2 = (
            (term2 / (2 * g)).integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c2
        )

        # Third term
        term3 = (self.u**2 + self.v**2 - self.u_ZE**2 - self.v_ZE**2) * self.omega
        term3 = CalcAreaAverage(term3, self.ylength, xlength=self.xlength) / (2 * g)
        term3 = self._handle_nans(term3)
        term3 = term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(
            **{self.VerticalCoordIndexer: 0}
        )

        # Combine terms and save the result
        function = term1 + term2 - term3
        function = self._handle_nans(function)
        Bkz = self._convert_units(function, "BKz")

        self.app_logger.debug("Done.")
        return Bkz

    def calc_bke(self):
        """Computes the flux of Eddy Kinetic Energy across the boundaries (BKE)."""
        self.app_logger.debug("Calculating BKE...")

        # First term
        term1 = self.u * (self.u_ZE**2 + self.v_ZE**2)
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(
            **{self.LonIndexer: self.western_limit}
        )
        term1 = (term1 / (2 * g)).integrate("rlats")
        term1 = self._handle_nans(term1)
        term1 = (
            term1.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c1
        )

        # Second term
        term2 = (self.u_ZE**2 + self.v_ZE**2) * self.v * self.v["coslats"]
        term2 = CalcZonalAverage(term2, self.xlength)
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(
            **{self.LatIndexer: self.southern_limit}
        )
        term2 = self._handle_nans(term2)
        term2 = (
            (term2 / (2 * g)).integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c2
        )

        # Third term
        term3 = (self.u_ZE**2 + self.v_ZE**2) * self.omega
        term3 = CalcAreaAverage(term3, self.ylength, xlength=self.xlength) / (2 * g)
        term3 = self._handle_nans(term3)
        term3 = term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(
            **{self.VerticalCoordIndexer: 0}
        )

        # Combine terms and save the result
        function = term1 + term2 - term3
        function = self._handle_nans(function)
        Bke = self._convert_units(function, "BKE")

        self.app_logger.debug("Done.")
        return Bke

    def calc_boz(self):
        """
        Computes the rate of change of Zonal Kinetic Energy due to work done by
        pressure forces at the boundaries of the control volume (BΦZ).

        .. math::
            B\\Phi_Z = c_1\\int\\!\\!\\int
                       \\frac{[u]\\,\\Delta_{EW}\\Phi'
                             + \\Phi^*\\,\\Delta_{EW}u'}{g}\\,d\\phi\\,dp
                     + c_2\\int \\left.\\frac{[v]\\Phi^*\\cos\\phi}{g}
                       \\right|_{\\phi_s}^{\\phi_n} dp
                     - \\left.\\frac{\\overline{\\omega^*\\Phi^*}}{g}
                       \\right|_{p_t}^{p_b}

        where :math:`\\Delta_{EW}X = X|_{\\lambda_e} - X|_{\\lambda_w}`.

        Sources
        -------
        East/west walls: Brennan and Vincent (1980, Appendix) write this wall as
        the flux ``u*Phi - u'*Phi'`` evaluated from west to east.  Expanding that
        east-minus-west difference gives ``[u] * D_EW(Phi')`` plus a second piece
        carrying ``D_EW(u')``.  The form used here takes that second piece with
        the area departure ``Phi*`` in place of the zonal mean ``[Phi]``: a
        limited-area box carries a net mass flux through its walls, so with the
        full geopotential the term would change when an arbitrary constant is
        added to Phi, whereas ``Phi*`` leaves it unchanged.

        Michaelides (1987, Appendix) writes this wall as ``[v] * Phi*`` evaluated
        from west to east.  Both factors are zonal means and so do not depend on
        longitude, which makes that expression identically zero.

        North/south and top/bottom walls: Muench (1965) and Michaelides (1987),
        with the area departure ``Phi*`` used throughout.  Brennan and Vincent
        (1980) write these two walls with the full geopotential instead.
        """
        self.app_logger.debug("Calculating BΦZ...")

        # ``CalcAreaAverage`` uses the analytic spherical area in its
        # denominator but a trapezoidal numerator.  Re-centering the existing
        # BoxData anomalies by the quadrature's mean-of-one enforces the
        # defining identities mean(Phi*) = mean(omega*) = 0 to roundoff.  In
        # particular, adding any pressure-only reference Phi_ref(p) then leaves
        # BΦZ invariant at machine precision.
        geopt_star, omega_star = self._area_anomalies()

        # First term: east/west faces.  Both contributions are mean-by-eddy
        # products evaluated on the walls, and both vanish on a periodic
        # (global) domain, where the two walls coincide.
        term1 = self._east_west_pressure_work(geopt_star)

        # Second term: north/south faces, area-anomalous zonal geopotential.
        term2 = (self.v_ZA * geopt_star) * self.v_ZA["coslats"] / g
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(
            **{self.LatIndexer: self.southern_limit}
        )
        term2 = self._handle_nans(term2)
        term2 = (
            term2.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c2
        )

        # Third term: top/bottom faces, using omega* Phi*.
        term3 = CalcAreaAverage(omega_star * geopt_star, self.ylength) / g
        term3 = self._handle_nans(term3)
        term3 = term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(
            **{self.VerticalCoordIndexer: 0}
        )

        function = term1 + term2 - term3
        function = self._handle_nans(function)
        Boz = self._convert_units(function, "BOZ")

        self.app_logger.debug("Done.")
        return Boz

    def _area_anomalies(self):
        """Return quadrature-centred ``Phi*`` and ``omega*``."""
        geopt_zonal, _ = self._zonal_geopotential_decomposition()
        area_quadrature_norm = (
            geopt_zonal["coslats"].integrate("rlats") / self.ylength
        )
        # Remove a latitude-independent reference *before* the meridional
        # quadrature.  The expression is algebraically unchanged, but avoids
        # integrating an artificial O(1e6) gauge and then subtracting two large
        # nearly equal numbers on NCEP's float32 coordinates.
        geopt_relative = geopt_zonal - geopt_zonal.isel(
            **{self.LatIndexer: 0}, drop=True
        )
        geopt_star = geopt_relative - (
            CalcAreaAverage(geopt_relative, self.ylength)
            / area_quadrature_norm
        )
        omega_star = self.omega_AE - (
            CalcAreaAverage(self.omega_AE, self.ylength) / area_quadrature_norm
        )
        # At a pole-to-pole domain xarray may drop this auxiliary coordinate
        # while aligning the two separately centred arrays.  It is the shared
        # latitude weight required by CalcAreaAverage, so restore it explicitly.
        geopt_star = geopt_star.assign_coords(coslats=geopt_zonal["coslats"])
        omega_star = omega_star.assign_coords(coslats=self.omega_AE["coslats"])
        return geopt_star, omega_star

    def _zonal_geopotential_decomposition(self):
        """Return quadrature-normalised ``[Phi]`` and ``Phi'``.

        NCEP longitude coordinates are stored as float32.  Accumulating the
        trapezoidal integral over a wide box can therefore make the zonal mean
        of one differ from one by O(1e-8).  Dividing by that measured norm
        enforces the defining identities and prevents a large constant gauge
        from leaking into ``Phi'``.
        """
        zonal_norm = CalcZonalAverage(
            self.geopt.metpy.dequantify() * 0 + 1, self.xlength
        )
        geopt_zonal = self.geopt_ZA / zonal_norm
        geopt_eddy = self.geopt - geopt_zonal
        return geopt_zonal, geopt_eddy

    def _east_west_pressure_work(self, geopt_star):
        """Pressure work through the east and west faces (term I of BΦZ).

        Uses the quadrature-normalised zonal eddy geopotential, so that a
        constant added to the geopotential leaves the result unchanged.
        """
        _, geopt_eddy = self._zonal_geopotential_decomposition()
        face_flux = (self.u_ZA * geopt_eddy + geopt_star * self.u_ZE) / g
        term = face_flux.sel(**{self.LonIndexer: self.eastern_limit}) - face_flux.sel(
            **{self.LonIndexer: self.western_limit}
        )
        term = self._handle_nans(term.integrate("rlats"), "BΦZ_east_west")
        return (
            term.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c1
        )

    def calc_boe(self):
        """
        Computes the appearence of Eddy Kinetic Energy associated with work
        produced at its boundaries (BΦE).

        .. math::
            B\\Phi_E = c_1\\int\\!\\!\\int \\left.\\frac{u'\\Phi'}{g}
                       \\right|_{\\lambda_w}^{\\lambda_e} d\\phi\\,dp
                     + c_2\\int \\left.\\frac{\\overline{v'\\Phi'}\\cos\\phi}{g}
                       \\right|_{\\phi_s}^{\\phi_n} dp
                     - \\left.\\frac{\\overline{\\omega'\\Phi'}}{g}
                       \\right|_{p_t}^{p_b}

        Brennan and Vincent (1980, pp. 964-965) and Michaelides (1987, p. 25)
        agree on this expression, and so does ``docs/source/math.rst``.  Every
        factor is an eddy covariance.
        """
        self.app_logger.debug("Calculating BΦE...")

        # First term: east/west faces, zonal eddy covariance u'Phi'.
        _, geopt_eddy = self._zonal_geopotential_decomposition()
        term1 = (self.u_ZE * geopt_eddy) / g
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(
            **{self.LonIndexer: self.western_limit}
        )
        term1 = term1.integrate("rlats")
        term1 = self._handle_nans(term1)
        term1 = (
            term1.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c1
        )

        # Second term: north/south faces, zonal average of the eddy covariance.
        term2 = CalcZonalAverage(self.v_ZE * geopt_eddy, self.xlength)
        term2 = term2 * term2["coslats"] / g
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(
            **{self.LatIndexer: self.southern_limit}
        )
        term2 = self._handle_nans(term2)
        term2 = (
            term2.integrate(self.VerticalCoordIndexer)
            * self.PressureData.metpy.units
            * self.c2
        )

        # Third term: top/bottom faces.
        term3 = (
            CalcAreaAverage(
                self.omega_ZE * geopt_eddy, self.ylength, xlength=self.xlength
            )
            / g
        )
        term3 = self._handle_nans(term3)
        term3 = term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(
            **{self.VerticalCoordIndexer: 0}
        )

        function = term1 + term2 - term3
        function = self._handle_nans(function)
        Boe = self._convert_units(function, "BOE")

        self.app_logger.debug("Done.")
        return Boe

    def _handle_nans(self, function, variable_name=""):
        """Delegate to the shared NaN policy.

        Interior NaNs are interpolated along the vertical coordinate and
        reported; pressure levels are never dropped here, because the vertical
        control volume is fixed once for the whole budget in BoxData.
        """
        return handle_nans(
            function,
            self.VerticalCoordIndexer,
            variable_name=variable_name,
            app_logger=getattr(self, "app_logger", None),
        )

    def _convert_units(self, function, variable_name):
        """Delegate to the shared unit conversion.

        pint raises DimensionalityError, which subclasses TypeError, so a
        ``ValueError`` guard would not catch it.
        """
        return convert_units(
            function,
            "W/m^2",
            variable_name,
            app_logger=getattr(self, "app_logger", None),
        )

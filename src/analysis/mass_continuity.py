"""Mass-continuity diagnostics for limited-area Lorenz-cycle budgets.

The residual implemented here is the area-mean pressure-coordinate
continuity equation in pressure coordinates:

    M = <D_EW> + <D_NS> + d(omega_bar)/dp.

``M`` is identically zero in the continuum, but not for independently
archived reanalysis winds and omega.  Its pressure profile is saved with the
other vertical diagnostics.  The scalar returned to the frameworks is the
column integral ``int M dp/g``, i.e. the unresolved column mass flux.
"""

import logging

from metpy.constants import Re, g

from ..utils.box_data import BoxData
from ..utils.nan_handling import convert_units, handle_nans


class MassContinuity:
    """Compute and export the limited-area mass-continuity residual."""

    def __init__(self, box_obj: BoxData, method: str, app_logger: logging.Logger):
        self.method = method
        self.box_obj = box_obj
        self.app_logger = app_logger
        self.results_subdirectory_vertical_levels = (
            box_obj.results_subdirectory_vertical_levels
        )

        self.LonIndexer = box_obj.LonIndexer
        self.LatIndexer = box_obj.LatIndexer
        self.TimeName = box_obj.TimeName
        self.VerticalCoordIndexer = box_obj.VerticalCoordIndexer
        self.PressureData = box_obj.PressureData

        self.western_limit = box_obj.western_limit
        self.eastern_limit = box_obj.eastern_limit
        self.southern_limit = box_obj.southern_limit
        self.northern_limit = box_obj.northern_limit
        self.xlength = box_obj.xlength
        self.ylength = box_obj.ylength

        self.u = box_obj.u
        self.v_ZA = box_obj.v_ZA
        self.omega_AA = box_obj.omega_AA

    def calc_mass_residual_profile(self):
        r"""Return ``M(p)`` in s\ :sup:`-1` and save its vertical profile.

        The two horizontally averaged divergences are evaluated in their
        exactly telescoped boundary forms.  These are algebraically identical
        to area-averaging

        ``D_EW = (u_e-u_w)/(a cos(phi) Delta(lambda))`` and
        ``D_NS = d([v]cos(phi))/dphi/(a cos(phi))``,

        while remaining well defined at the poles for a global domain.
        """
        self.app_logger.debug("Calculating mass-continuity residual M...")

        u_faces = self.u.sel(**{self.LonIndexer: self.eastern_limit}) - self.u.sel(
            **{self.LonIndexer: self.western_limit}
        )
        self.d_ew_area = u_faces.integrate("rlats") / (
            Re * self.xlength * self.ylength
        )

        meridional_flux = self.v_ZA * self.v_ZA["coslats"]
        self.d_ns_area = (
            meridional_flux.sel(**{self.LatIndexer: self.northern_limit})
            - meridional_flux.sel(**{self.LatIndexer: self.southern_limit})
        ) / (Re * self.ylength)

        self.domega_dp = self.omega_AA.differentiate(self.VerticalCoordIndexer) / (
            self.PressureData.metpy.units
        )
        residual = handle_nans(
            self.d_ew_area + self.d_ns_area + self.domega_dp,
            self.VerticalCoordIndexer,
            variable_name="M",
            app_logger=self.app_logger,
        )
        self.residual_profile = residual
        self._save_vertical_levels(residual, "M")
        return residual

    def calc_mass_residual(self):
        r"""Return the conservative column residual in kg m\ :sup:`-2` s\ :sup:`-1`.

        The vertical derivative is integrated as the exact boundary
        difference ``omega_bottom - omega_top``.  This is the finite-volume
        counterpart of ``int d(omega_bar)/dp dp`` and prevents a second,
        avoidable differentiation-then-quadrature error in the scalar output.
        """
        self.calc_mass_residual_profile()
        horizontal_column = (
            (self.d_ew_area + self.d_ns_area).integrate(
                self.VerticalCoordIndexer
            )
            * self.PressureData.metpy.units
        )
        vertical_column = self.omega_AA.isel(
            **{self.VerticalCoordIndexer: -1}
        ) - self.omega_AA.isel(**{self.VerticalCoordIndexer: 0})
        column = (horizontal_column + vertical_column) / g
        result = convert_units(
            column,
            "kg/m^2/s",
            "M",
            app_logger=self.app_logger,
        )
        self.app_logger.debug("Mass-continuity residual M done.")
        return result

    def _save_vertical_levels(self, function, variable_name):
        """Append a profile using the same layout as the energy diagnostics."""
        df = function.to_dataframe(name=variable_name).reset_index()

        if self.method == "fixed":
            if self.TimeName not in function.dims:
                df = df.T
            else:
                df = df.pivot(
                    index=self.TimeName, columns=self.VerticalCoordIndexer
                )
        else:
            df.set_index(self.TimeName, inplace=True)
            df.index = df.index.strftime("%Y-%m-%d %H:%M:%S")
            df = df.pivot(
                columns=self.VerticalCoordIndexer, values=variable_name
            )
            df.columns.name = None

        df.to_csv(
            f"{self.results_subdirectory_vertical_levels}/"
            f"{variable_name}_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None,
        )

# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    box_data.py                                        :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/01/31 20:15:59 by daniloceano       #+#    #+#              #
#    Updated: 2024/07/18 00:19:08 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create an object that will store, whithin a bounding box specified by the user,
all meteorological data required by the Lorenz Energy Cycle computations.

The program uses xarray built-in function for integrations. For using it, it
would be necessary to sort the data from south to north and from west to east,
but for saving up memory, it was using a negative sign in front of the
integrations for the latitude coordinate.

Also, it was created a dimension for the radians of the latitude and longitude,
for the integration

Suffixes used:
    ZA = Zonal average (whithin the dfined box)
    AA = Area average (whithin the dfined box)
    ZE = zonal eddy (departure from zonal average)
    AE = area eddy (zonal departures from area averages)


Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com

@author: daniloceano
"""

import argparse

import numpy as np
import pandas as pd
import xarray as xr
from metpy.constants import g
from metpy.units import units

from ..utils.calc_averages import CalcAreaAverage, CalcZonalAverage
from ..utils.thermodynamics import AdiabaticHEating, StaticStability


class BoxData:
    """
    A class to store and process meteorological data within a specified bounding box
    for Lorenz Energy Cycle computations.

    Attributes:
        data (xr.Dataset): Dataset containing meteorological data.
        variable_list_df (pd.DataFrame): DataFrame mapping variables to dataset fields.
        western_limit (float): Western boundary of the bounding box.
        eastern_limit (float): Eastern boundary of the bounding box.
        southern_limit (float): Southern boundary of the bounding box.
        northern_limit (float): Northern boundary of the bounding box.
        args (argparse.Namespace): Arguments for additional processing controls.
        results_subdirectory (str): Directory path for output data.
        dTdt (xr.DataArray, optional): Dataset containing temperature gradient data.

    Methods:
        __init__: Initializes the BoxData object with the given parameters.
    """

    def __init__(
        self,
        data: xr.Dataset,
        variable_list_df: pd.DataFrame,
        western_limit: float,
        eastern_limit: float,
        southern_limit: float,
        northern_limit: float,
        args: argparse.Namespace,
        results_subdirectory: str,
        results_subdirectory_vertical_levels: str,
        dTdt: xr.DataArray = None,
    ):
        """Initialize the BoxData object with meteorological data and computation parameters."""
        self._initialize_indices(variable_list_df, data)
        self.args = args
        self.results_subdirectory = results_subdirectory
        self.results_subdirectory_vertical_levels = results_subdirectory_vertical_levels

        self.dx = float(data[self.LatIndexer][1] - data[self.LatIndexer][0])

        # Set domain limits and calculate lengths for averaging
        self._set_domain_limits(
            data, western_limit, eastern_limit, southern_limit, northern_limit
        )

        # Extract and process meteorological data within the defined box
        self._process_meteorological_data(data, variable_list_df, dTdt)

    def _initialize_indices(self, variable_list_df, data):
        """Initialize indices from the variable list DataFrame."""
        self.LonIndexer = variable_list_df.loc["Longitude"]["Variable"]
        self.LatIndexer = variable_list_df.loc["Latitude"]["Variable"]
        self.TimeName = variable_list_df.loc["Time"]["Variable"]
        self.VerticalCoordIndexer = variable_list_df.loc["Vertical Level"]["Variable"]
        self.PressureData = data[self.VerticalCoordIndexer] * units("Pa")

    def _set_domain_limits(
        self, data, western_limit, eastern_limit, southern_limit, northern_limit
    ):
        """Set domain limits for the bounding box and calculate lengths for averaging."""
        self.western_limit = self._select_nearest(data, self.LonIndexer, western_limit)
        self.eastern_limit = self._select_nearest(data, self.LonIndexer, eastern_limit)
        self.southern_limit = self._select_nearest(
            data, self.LatIndexer, southern_limit
        )
        self.northern_limit = self._select_nearest(
            data, self.LatIndexer, northern_limit
        )

        self.xlength = self.eastern_limit["rlons"] - self.western_limit["rlons"]
        self.ylength = np.sin(self.northern_limit["rlats"]) - np.sin(
            self.southern_limit["rlats"]
        )

    def _select_nearest(self, data, indexer, value):
        """Select the data point nearest to the specified value."""
        return data[indexer].sel({indexer: value}, method="nearest")

    def _process_meteorological_data(self, data, variable_list_df, dTdt):
        """Extract and process meteorological data within the defined box."""
        # Extract and process temperature data
        self._process_temperature_data(data, variable_list_df)

        # Extract and process wind component data
        self._process_wind_data(data, variable_list_df)

        # Extract and process frictional terms
        self._process_friction_terms(data, variable_list_df, self.args)

        # Extract and process vertical velocity (omega)
        self._process_omega(data, variable_list_df)

        # Extract and process geopotential
        self._process_geopotential_data(data, variable_list_df)

        # Extract and process static stability parameter
        self._process_adiabatic_heating_and_stability(dTdt)

    def _process_temperature_data(self, data, variable_list_df):
        """Extract and process temperature data."""
        self.tair = self._extract_data(data, variable_list_df, "Air Temperature", "K")
        self.tair_ZA = CalcZonalAverage(self.tair, self.xlength)
        self.tair_AA = CalcAreaAverage(self.tair_ZA, self.ylength)
        self.tair_ZE = self.tair - self.tair_ZA
        self.tair_AE = self.tair_ZA - self.tair_AA

    def _process_wind_data(self, data, variable_list_df):
        """Extract and process wind component data."""
        self.u = self._extract_data(
            data, variable_list_df, "Eastward Wind Component", "m/s"
        )
        self.u_ZA = CalcZonalAverage(self.u, self.xlength)
        self.u_AA = CalcAreaAverage(self.u_ZA, self.ylength)
        self.u_ZE = self.u - self.u_ZA
        self.u_AE = self.u_ZA - self.u_AA

        self.v = self._extract_data(
            data, variable_list_df, "Northward Wind Component", "m/s"
        )
        self.v_ZA = CalcZonalAverage(self.v, self.xlength)
        self.v_AA = CalcAreaAverage(self.v_ZA, self.ylength)
        self.v_ZE = self.v - self.v_ZA
        self.v_AE = self.v_ZA - self.v_AA

    def _process_friction_terms(self, data, variable_list_df, args):
        """Extract and process friction terms."""
        if args.residuals:
            self.ust = self.tair * np.nan
            self.vst = self.tair * np.nan

        else:
            self.ust = self._extract_data(
                data, variable_list_df, "Friction Velocity", "m/s"
            )
            self.vst = self._extract_data(
                data, variable_list_df, "Friction Velocity", "m/s"
            )

        self.ust_ZA = CalcZonalAverage(self.ust, self.xlength)
        self.ust_AA = CalcAreaAverage(self.ust_ZA, self.ylength)
        self.ust_ZE = self.ust - self.ust_ZA
        self.ust_AE = self.ust_ZA - self.ust_AA

        self.vst_ZA = CalcZonalAverage(self.vst, self.xlength)
        self.vst_AA = CalcAreaAverage(self.vst_ZA, self.ylength)
        self.vst_ZE = self.vst - self.vst_ZA
        self.vst_AE = self.vst_ZA - self.vst_AA

    def _process_omega(self, data, variable_list_df):
        """Extract and process vertical velocity."""
        self.omega = self._extract_data(
            data, variable_list_df, "Omega Velocity", "Pa/s"
        )
        self.omega_ZA = CalcZonalAverage(self.omega, self.xlength)
        self.omega_AA = CalcAreaAverage(self.omega_ZA, self.ylength)
        self.omega_ZE = self.omega - self.omega_ZA
        self.omega_AE = self.omega_ZA - self.omega_AA

    def _process_geopotential_data(self, data, variable_list_df):
        """Extract and process geopotential or geopotential height data."""
        if "Geopotential" in variable_list_df.index:
            self.geopt = self._extract_data(
                data, variable_list_df, "Geopotential", "m**2/s**2"
            )
        else:
            self.geopt = self._extract_and_convert_geopotential_height(
                data, variable_list_df
            )

        self.geopt_ZA = CalcZonalAverage(self.geopt, self.xlength)
        self.geopt_AA = CalcAreaAverage(self.geopt_ZA, self.ylength)
        self.geopt_ZE = self.geopt - self.geopt_ZA
        self.geopt_AE = self.geopt_ZA - self.geopt_AA

    def _extract_and_convert_geopotential_height(self, data, variable_list_df):
        """Extract geopotential height data and convert to geopotential."""
        geopt_height = self._extract_data(
            data,
            variable_list_df,
            "Geopotential Height",
            variable_list_df.loc["Geopotential Height"]["Units"],
        )
        return (geopt_height * g).metpy.convert_units("m**2/s**2")

    def _process_adiabatic_heating_and_stability(self, dTdt):
        """Process adiabatic heating and static stability parameter."""
        if self.args.fixed:
            self.Q = AdiabaticHEating(
                self.tair,
                self.PressureData,
                self.omega,
                self.u,
                self.v,
                self.VerticalCoordIndexer,
                self.LatIndexer,
                self.LonIndexer,
                self.TimeName,
            ).sel(
                **{
                    self.LatIndexer: slice(self.southern_limit, self.northern_limit),
                    self.LonIndexer: slice(self.western_limit, self.eastern_limit),
                }
            )
        elif self.args.track or self.args.choose:
            self.dTdt = dTdt
            self.Q = AdiabaticHEating(
                self.tair,
                self.PressureData,
                self.omega,
                self.u,
                self.v,
                self.VerticalCoordIndexer,
                self.LatIndexer,
                self.LonIndexer,
                self.TimeName,
                self.dTdt,
            ).sel(
                **{
                    self.LatIndexer: slice(self.southern_limit, self.northern_limit),
                    self.LonIndexer: slice(self.western_limit, self.eastern_limit),
                }
            )
        else:
            print("Adiabatic heating computation skipped. Check flags!")

        self.Q_ZA = CalcZonalAverage(self.Q, self.xlength)
        self.Q_AA = CalcAreaAverage(self.Q_ZA, self.ylength)
        self.Q_ZE = self.Q - self.Q_ZA
        self.Q_AE = self.Q_ZA - self.Q_AA

        self.sigma_AA = StaticStability(
            self.tair,
            self.PressureData,
            self.VerticalCoordIndexer,
            self.xlength,
            self.ylength,
        )

    def _extract_data(self, data, variable_list_df, variable_name, unit):
        """Extract data for a specific variable and convert to the specified unit."""
        var_key = variable_list_df.loc[variable_name]["Variable"]
        unit_to_convert = units(variable_list_df.loc[variable_name]["Units"])
        return (
            (data[var_key] * unit_to_convert)
            .metpy.convert_units(unit)
            .sel(
                **{
                    self.LatIndexer: slice(self.southern_limit, self.northern_limit),
                    self.LonIndexer: slice(self.western_limit, self.eastern_limit),
                }
            )
        )

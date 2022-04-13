#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 20:15:59 2022

Create an object that will store all meteorological data required by the Lorenz
Energy Cycle computations. The self functions computate the zonal and area 
averages and eddy terms for each variable, as well the static stability index.

@author: danilocoutodsouza
"""

import xarray
import calc

class MetData:
    '''
    Object containing all meteorological data required for the LEC computation
    '''
    def __init__(self,  LonIndexer: str, LatIndexer: str, TimeName: str,
                 VerticalCoordIndexer: str, TemperatureData: xarray.Dataset,
                 PressureData: xarray.Dataset,
                 UWindComponentData: xarray.Dataset,
                 VWindComponentData: xarray.Dataset,
                 OmegaData: xarray.Dataset):
        self.PressureData = PressureData
        self.LonIndexer = LonIndexer
        self.LatIndexer = LatIndexer
        self.TimeName = TimeName
        self.VerticalCoordIndexer = VerticalCoordIndexer
        
        # Temperature data values, averages and eddy terms
        self.tair = TemperatureData
        self.tair_ZA = calc.CalcZonalAverage(self.tair, self.LonIndexer)
        self.tair_AA = calc.CalcAreaAverage(self.tair, self.LonIndexer,self.LatIndexer)
        self.tair_ZE = self.tair - self.tair_ZA
        self.tair_AE = self.tair_ZA - self.tair_AA
        
        # Zonal wind component data values, averages and eddy terms
        self.u = UWindComponentData
        self.u_ZA = calc.CalcZonalAverage(self.u, self.LonIndexer)
        self.u_AA = calc.CalcAreaAverage(self.u, self.LonIndexer,self.LatIndexer)
        self.u_ZE = self.u - self.u_ZA
        self.u_AE = self.u_ZA - self.u_AA
        
        # Meridional wind component data values, averages and eddy terms
        self.v = VWindComponentData
        self.v_ZA = calc.CalcZonalAverage(self.v, self.LonIndexer)
        self.v_AA = calc.CalcAreaAverage(self.v, self.LonIndexer,self.LatIndexer)
        self.v_ZE = self.v - self.v_ZA
        self.v_AE = self.v_ZA - self.v_AA
        
        # Omega velocity (vertical velocity in pressure levels) data values,
        # averages and eddy terms
        self.omega = OmegaData
        self.omega_ZA = calc.CalcZonalAverage(self.omega, self.LonIndexer)
        self.omega_AA = calc.CalcAreaAverage(self.omega, self.LonIndexer,self.LatIndexer)
        self.omega_ZE = self.omega - self.omega_ZA
        self.omega_AE = self.omega_ZA - self.omega_AA
        
        # Static stability parameter
        self.sigma_AA = calc.StaticStability(self.tair, self.PressureData, self.VerticalCoordIndexer,
                        self.LatIndexer, self.LonIndexer)

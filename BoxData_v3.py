#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 17:05:40 2022

@author: danilocoutodsouza
"""

import xarray
import calc
import numpy as np

class BoxData:
    '''
    Object containing all meteorological data required for the LEC computation
    '''
    def __init__(self,  LonIndexer: str, LatIndexer: str, TimeName: str,
                 VerticalCoordIndexer: str, TemperatureData: xarray.Dataset,
                 PressureData: xarray.Dataset,
                 UWindComponentData: xarray.Dataset,
                 VWindComponentData: xarray.Dataset,
                 OmegaData: xarray.Dataset,
                 western_limit: float, eastern_limit: float,
                 southern_limit: float, northern_limit: float):
        self.PressureData = PressureData
        self.LonIndexer = LonIndexer
        self.LatIndexer = LatIndexer
        self.TimeName = TimeName
        self.VerticalCoordIndexer = VerticalCoordIndexer
        
        # Find data gridpoints that match the box limits
        self.BoxWest = float((TemperatureData[LonIndexer][(np.abs(TemperatureData[LonIndexer] - western_limit)).argmin()]).values)
        self.BoxEast = float((TemperatureData[LonIndexer][(np.abs(TemperatureData[LonIndexer] - eastern_limit)).argmin()]).values)
        self.BoxSouth = float((TemperatureData[LatIndexer][(np.abs(TemperatureData[LatIndexer] - southern_limit)).argmin()]).values)
        self.BoxNorth = float((TemperatureData[LatIndexer][(np.abs(TemperatureData[LatIndexer] - northern_limit)).argmin()]).values)
        
        # Suffixes used:
            # ZA = Zonal average (whithin the dfined box)
            # GZA = Global zonal average
            # AA = Area average (whithin the dfined box)
            # GA = Global Average
            # ZE = zonal eddy (departure from zonal average)
            # AE = area eddy (zonal departures from area averages)
        
        # Temperature data values, averages and eddy terms
        self.tair = TemperatureData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.tair_ZA = calc.CalcZonalAverage(self.tair, self.LonIndexer)
        self.tair_AA = calc.CalcAreaAverage(self.tair, self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.tair_GZA = calc.CalcZonalAverage(TemperatureData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth)}),self.LonIndexer)
        self.tair_GA = calc.CalcAreaAverage(TemperatureData, self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.tair_ZE = self.tair - self.tair_GZA
        self.tair_AE = self.tair_ZA - self.tair_GA
        
        # Zonal wind component data values, averages and eddy terms
        self.u = UWindComponentData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.u_ZA = calc.CalcZonalAverage(self.u, self.LonIndexer)
        self.u_AA = calc.CalcAreaAverage(self.u,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.u_GZA = calc.CalcZonalAverage(UWindComponentData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth)}),self.LonIndexer)
        self.u_GA = calc.CalcAreaAverage(UWindComponentData,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.u_ZE = self.u - self.u_GZA
        self.u_AE = self.u_ZA - self.u_GA
        
        # Meridional wind component data values, averages and eddy terms
        self.v = VWindComponentData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.v_ZA = calc.CalcZonalAverage(self.v, self.LonIndexer)
        self.v_AA = calc.CalcAreaAverage(self.v,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.v_GZA = calc.CalcZonalAverage(VWindComponentData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth)}),self.LonIndexer)
        self.v_GA = calc.CalcAreaAverage(VWindComponentData,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.v_ZE = self.v - self.v_GZA
        self.v_AE = self.v_ZA - self.v_GA
        
        # Omega velocity (vertical velocity in pressure levels) data values,
        # averages and eddy terms
        self.omega = OmegaData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.omega_ZA = calc.CalcZonalAverage(self.omega, self.LonIndexer)
        self.omega_AA = calc.CalcAreaAverage(self.omega,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.omega_GZA = calc.CalcZonalAverage(OmegaData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth)}),self.LonIndexer)
        self.omega_GA = calc.CalcAreaAverage(OmegaData,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.omega_ZE = self.omega - self.omega_GZA
        self.omega_AE = self.omega_ZA - self.omega_GA
        
        # Static stability parameter
        self.sigma_AA = calc.StaticStability(self.tair, self.PressureData, self.VerticalCoordIndexer,
                        self.LatIndexer, self.LonIndexer
                        ,self.BoxNorth, self.BoxSouth,
                        self.BoxWest, self.BoxEast)
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
        
        ## For the computations bellow, the variable is sliced for the
        # box definied for computations, whereas the zonal average is 
        # computed for the whole globe, as consistent to the Lorenz theory.
        # Likewise, the eddy terms are computed taking the 
        
        # Temperature data values, averages and eddy terms
        self.tair = TemperatureData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.tair_ZA = calc.CalcZonalAverage(self.tair,
                                self.LonIndexer).sel(**{LatIndexer: 
                                slice(self.BoxNorth, self.BoxSouth)})
        self.tair_AA = calc.CalcAreaAverage(TemperatureData,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.tair_ZE = (TemperatureData - self.tair_ZA).sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.tair_AE = (self.tair_ZA - self.tair_AA).sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth)})
        
        # Zonal wind component data values, averages and eddy terms
        self.u = UWindComponentData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.u_ZA = calc.CalcZonalAverage(UWindComponentData, self.LonIndexer)
        self.u_AA = calc.CalcAreaAverage(UWindComponentData,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.u_ZE = (UWindComponentData - self.u_ZA).sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.u_AE = (self.u_ZA - self.u_AA).sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth)})
        
        # Meridional wind component data values, averages and eddy terms
        self.v = VWindComponentData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.v_ZA = calc.CalcZonalAverage(VWindComponentData, self.LonIndexer)
        self.v_AA = calc.CalcAreaAverage(VWindComponentData,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.v_ZE = (VWindComponentData - self.v_ZA).sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.v_AE = (self.v_ZA - self.v_AA).sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth)})
        
        # Omega velocity (vertical velocity in pressure levels) data values,
        # averages and eddy terms
        self.omega = OmegaData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.omega_ZA = calc.CalcZonalAverage(OmegaData, self.LonIndexer)
        self.omega_AA = calc.CalcAreaAverage(OmegaData,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.omega_ZE = (OmegaData - self.omega_ZA).sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.omega_AE = (self.omega_ZA - self.omega_AA).sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth)})
        
        # Static stability parameter
        self.sigma_AA = calc.StaticStability(self.tair, self.PressureData, self.VerticalCoordIndexer,
                        self.LatIndexer, self.LonIndexer
                        ,self.BoxNorth, self.BoxSouth,
                        self.BoxWest, self.BoxEast)

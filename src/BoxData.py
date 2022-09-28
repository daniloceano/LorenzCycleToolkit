#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 20:15:59 2022
Create an object that will store, whithin a bounding box specified by the user,
all meteorological data required by the Lorenz Energy Cycle computations.
The recursive functions computate the zonal and area averages and eddy terms 
for each variable, as well the static stability index.
@author: danilocoutodsouza

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

"""

import xarray
import Math
from thermodynamics import StaticStability
import numpy as np
from metpy.constants import g

class BoxData:
    '''
    Object containing all meteorological data required for the LEC computation
    '''
    def __init__(self,  LonIndexer: str, LatIndexer: str,
                 TimeName: str,
                 VerticalCoordIndexer: str, 
                 TemperatureData: xarray.Dataset,
                 PressureData: xarray.Dataset,
                 UWindComponentData: xarray.Dataset,
                 VWindComponentData: xarray.Dataset,
                 AdiabaticHeatingData: xarray.Dataset,
                 ZonalWindStressData: xarray.Dataset,
                 MeridionalWindStressData: xarray.Dataset,
                 OmegaData: xarray.Dataset,
                 HgtData: xarray.Dataset,
                 western_limit: float, eastern_limit: float,
                 southern_limit: float, northern_limit: float,
                 output_dir: str):
        self.PressureData = PressureData
        self.LonIndexer = LonIndexer
        self.LatIndexer = LatIndexer
        self.TimeName = TimeName
        self.VerticalCoordIndexer = VerticalCoordIndexer
        self.output_dir = output_dir
        
        # Find data gridpoints that match the box limits
        self.BoxWest = float((TemperatureData[LonIndexer]
                              [(np.abs(TemperatureData[LonIndexer] - 
                                       western_limit)).argmin()]).values)
        self.BoxEast = float((TemperatureData[LonIndexer]
                              [(np.abs(TemperatureData[LonIndexer] - 
                                       eastern_limit)).argmin()]).values)
        self.BoxSouth = float((TemperatureData[LatIndexer]
                               [(np.abs(TemperatureData[LatIndexer] - 
                                        southern_limit)).argmin()]).values)
        self.BoxNorth = float((TemperatureData[LatIndexer]
                               [(np.abs(TemperatureData[LatIndexer] - 
                                        northern_limit)).argmin()]).values)
        
        # Temperature data values, averages and eddy terms
        self.tair = TemperatureData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.tair_ZA = Math.CalcZonalAverage(self.tair, self.LonIndexer)
        self.tair_AA = Math.CalcAreaAverage(self.tair, self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.tair_ZE = self.tair - self.tair_ZA
        self.tair_AE = self.tair_ZA - self.tair_AA
        
        # Zonal wind component data values, averages and eddy terms
        self.u = UWindComponentData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.u_ZA = Math.CalcZonalAverage(self.u, self.LonIndexer)
        self.u_AA = Math.CalcAreaAverage(self.u,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.u_ZE = self.u - self.u_ZA
        self.u_AE = self.u_ZA - self.u_AA
        
        # Meridional wind component data values, averages and eddy terms
        self.v = VWindComponentData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.v_ZA = Math.CalcZonalAverage(self.v, self.LonIndexer)
        self.v_AA = Math.CalcAreaAverage(self.v,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.v_ZE = self.v - self.v_ZA
        self.v_AE = self.v_ZA - self.v_AA
        
        # Zonal wind stress data values, averages and eddy terms
        self.ust = ZonalWindStressData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.ust_ZA = Math.CalcZonalAverage(self.ust, self.LonIndexer)
        self.ust_AA = Math.CalcAreaAverage(self.ust,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.ust_ZE = self.ust - self.ust_ZA
        self.ust_AE = self.ust_ZA - self.ust_AA
        
        # Meridional wind stress data values, averages and eddy terms
        self.vst = MeridionalWindStressData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            self.LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.vst_ZA = Math.CalcZonalAverage(self.vst, self.LonIndexer)
        self.vst_AA = Math.CalcAreaAverage(self.vst,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.vst_ZE = self.vst - self.vst_ZA
        self.vst_AE = self.vst_ZA - self.vst_AA
        
        # Omega velocity (vertical velocity in pressure levels) data values,
        # averages and eddy terms
        self.omega = OmegaData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.omega_ZA = Math.CalcZonalAverage(self.omega, self.LonIndexer)
        self.omega_AA = Math.CalcAreaAverage(self.omega,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.omega_ZE = self.omega - self.omega_ZA
        self.omega_AE = self.omega_ZA - self.omega_AA
        
        # Geopotential height data values, averages and eddy terms
        self.hgt = HgtData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.hgt_ZA = Math.CalcZonalAverage(self.hgt, self.LonIndexer)
        self.hgt_AA = Math.CalcAreaAverage(self.hgt,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.hgt_ZE = self.hgt - self.hgt_ZA
        self.hgt_AE = self.hgt_ZA - self.hgt_AA
        
        # Geopotential (g*hgt) data values, averages and eddy terms
        self.geopt = self.hgt*g
        self.geopt_ZA = Math.CalcZonalAverage(self.geopt, self.LonIndexer)
        self.geopt_AA = Math.CalcAreaAverage(self.geopt,self.LatIndexer,
                                            self.BoxSouth,self.BoxNorth,
                                            self.LonIndexer)
        self.geopt_ZE = self.geopt - self.geopt_ZA
        self.geopt_AE = self.geopt_ZA - self.geopt_AA
        
        self.Q = AdiabaticHeatingData.sel(**{LatIndexer: 
            slice(self.BoxNorth, self.BoxSouth),
            LonIndexer: slice(self.BoxWest, self.BoxEast)})
        self.Q_ZA = Math.CalcZonalAverage(self.Q,self.LonIndexer)
        self.Q_AA = Math.CalcAreaAverage(self.Q,self.LatIndexer,
                                    LonIndexer=self.LonIndexer)
        self.Q_ZE = self.Q - self.Q_ZA
        self.Q_AE = self.Q_ZA - self.Q_AA
        
        # Static stability parameter
        self.sigma_AA = StaticStability(self.tair, self.PressureData, self.VerticalCoordIndexer,
                        self.LatIndexer, self.LonIndexer
                        ,self.BoxNorth, self.BoxSouth,
                        self.BoxWest, self.BoxEast)
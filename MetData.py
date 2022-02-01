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
                 PressureData: xarray.Dataset):
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
        
        # Static stability parameter
        self.sigma_AA = calc.StaticStability(self.tair, self.PressureData, self.VerticalCoordIndexer,
                        self.LatIndexer, self.LonIndexer)
        
        
        
    # def _calc_sigma_aa(self):
    #     '''
    #     Computates the Static Stability parameter for data
    #     '''
    #     self.sigma_AA = calc.StaticStability(self.tair, self.PressureData, self.VerticalCoordIndexer,
    #                     self.LatIndexer, self.LonIndexer)
    #     return self.sigma_AA    
    # def _calc_tair_za(self):
    #     '''
    #     Computates zonal average for temperature data
    #     '''
    #     self.tair_ZA = calc.CalcZonalAverage(self.tair, self.LonIndexer)
    #     return self.tair_ZA 
    # def _calc_tair_aa(self):
    #     '''
    #     Computates area average for temperature data
    #     '''
    #     self.tair_AA = calc.CalcAreaAverage(self.tair,
    #                                         self.LonIndexer,self.LatIndexer)
    #     return self.tair_AA
    # def _calc_tair_ze(self):
    #     '''
    #     Computates departure from zonal average (zonal eddy) for temperature data
    #     '''
    #     self.tair_ZE = self.tair - self.tair_ZA
    #     return self.tair_ZE
    # def _calc_tair_ae(self):
    #     '''
    #     Computates departure from area average (area eddy) for temperature data
    #     '''
    #     self.tair_AE = self.tair_ZA - self.tair_AA
    #     return self.tair_AE
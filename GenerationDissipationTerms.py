#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:27:40 2020

@author: Danilo
"""

import numpy as np
from metpy.constants import g
from metpy.constants import Rd
from metpy.constants import Re
from metpy.units import units
from calc import (CalcAreaAverage,VerticalTrazpezoidalIntegration,Differentiate)
from BoxData import BoxData
from EnergyContents import function_to_df

class GenerationDissipationTerms:
    
    def __init__(self, box_obj: BoxData):
        self.PressureData = box_obj.PressureData
        self.LonIndexer = box_obj.LonIndexer
        self.LatIndexer = box_obj.LatIndexer
        self.TimeName = box_obj.TimeName
        self.VerticalCoordIndexer = box_obj.VerticalCoordIndexer
        self.output_dir = box_obj.output_dir
        self.tair = box_obj.tair
        self.tair_AE = box_obj.tair_AE
        self.tair_ZE = box_obj.tair_ZE
        self.u_ZA = box_obj.u_ZA
        self.u_ZE = box_obj.u_ZE
        self.v_ZA = box_obj.v_ZA
        self.v_ZE = box_obj.v_ZE
        self.sigma_AA = box_obj.sigma_AA
        self.omega_ZE = box_obj.omega_ZE
        self.omega_AE = box_obj.omega_AE
        
        self.rlats = np.deg2rad(self.tair_AE[self.LatIndexer])
        self.cos_lats = np.cos(self.rlats)
        self.tan_lats = np.tan(self.rlats)

    
    
    def AdiabaticHEating(self):
        """
        Compute the diabatic heating as a residual form the thermodynamic equation
        for all vertical levels and for the desired domain
        """
        # Temperature tendency as dT/dt
        TairTendency = self.tair.copy(deep=True).differentiate(
            self.TimeName,datetime_unit='h') 
        ## Vertical advection of temperature
        dTdp = self.tair.copy(deep=True
            ).sortby(self.VerticalCoordIndexer,ascending=True
            ).differentiate(self.VerticalCoordIndexer) / units.hPa
    
    
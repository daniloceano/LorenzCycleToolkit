#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:27:40 2020

@author: Danilo
"""

import numpy as np
from metpy.units import units
from metpy.constants import Cp_d
from calc import (CalcZonalAverage,CalcAreaAverage,VerticalTrazpezoidalIntegration,
                  Differentiate,AdiabaticHEating)
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
        self.u = box_obj.u
        self.u_ZA = box_obj.u_ZA
        self.u_ZE = box_obj.u_ZE
        self.v = box_obj.v
        self.v_ZA = box_obj.v_ZA
        self.v_ZE = box_obj.v_ZE
        self.sigma_AA = box_obj.sigma_AA
        self.omega= box_obj.omega
        self.omega_ZE = box_obj.omega_ZE
        self.omega_AE = box_obj.omega_AE
        
        self.rlons = np.deg2rad(self.tair[self.LonIndexer])
        self.cos_lons = np.cos(self.rlons)

        self.rlats = np.deg2rad(self.tair[self.LatIndexer])
        self.cos_lats = np.cos(self.rlats)
        self.tan_lats = np.tan(self.rlats)
        
        self.sigma_AA = box_obj.sigma_AA
        
        self.Q = AdiabaticHEating(self.tair,self.PressureData,self.omega,
                                    self.u,self.v,self.VerticalCoordIndexer,
                                    self.LatIndexer,self.LonIndexer,
                                    self.TimeName)
        self.Q_ZA = CalcZonalAverage(self.Q,self.LatIndexer)
        self.Q_AA = CalcAreaAverage(self.Q,self.LatIndexer,
                                    LonIndexer=self.LonIndexer)
        self.Q_ZE = self.Q - self.Q_AA
        self.Q_AE = self.Q_ZA - self.Q_AA
    
    def calc_gz(self):
        
        _ = (self.Q_AE*self.tair_AE)/(Cp_d*self.sigma_AA)
        function = CalcAreaAverage(_,self.LatIndexer,
                                    LonIndexer=self.LonIndexer)
        Gz = VerticalTrazpezoidalIntegration(function,self.PressureData,
                                             self.VerticalCoordIndexer)
        try: 
            Gz = Gz.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Gz')
            raise
        print(Gz.values*Gz.metpy.units)
        # Save Ca before vertical integration
        print('Saving Gz for each vertical level...')
        try:
            df = function_to_df(self,self.VerticalCoordIndexer,function)
            df.to_csv(self.output_dir+'/Gz_'+self.VerticalCoordIndexer+'.csv')
        except:
            raise('Could not save file with Gz for each level')
        print('Done!')
        return Gz
    
    def calc_ge(self):
        
        _ = (self.Q_ZE*self.tair_ZE)/(Cp_d*self.sigma_AA)
        function = CalcAreaAverage(_,self.LatIndexer,
                                    LonIndexer=self.LonIndexer)
        Ge = VerticalTrazpezoidalIntegration(function,self.PressureData,
                                             self.VerticalCoordIndexer)
        try: 
            Ge = Ge.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Ge')
            raise
        print(Ge.values*Ge.metpy.units)
        # Save Ca before vertical integration
        print('Saving Ge for each vertical level...')
        try:
            df = function_to_df(self,self.VerticalCoordIndexer,function)
            df.to_csv(self.output_dir+'/Ge_'+self.VerticalCoordIndexer+'.csv')
        except:
            raise('Could not save file with Ge for each level')
        print('Done!')
        return Ge
        
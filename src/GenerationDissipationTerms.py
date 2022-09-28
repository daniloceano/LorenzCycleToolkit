#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:27:40 2020

@author: Danilo
"""

from metpy.units import units
from metpy.constants import Cp_d
from metpy.constants import g
from Math import (CalcZonalAverage,CalcAreaAverage,
                  VerticalTrazpezoidalIntegration)
from BoxData import BoxData
import pandas as pd

class GenerationDissipationTerms:
    
    def __init__(self, box_obj: BoxData, method: str):
        self.method = method
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
        self.ust_ZA = box_obj.ust_ZA
        self.ust_ZE = box_obj.ust_ZE
        self.vst_ZE = box_obj.vst_ZE
        self.vst_ZA = box_obj.vst_ZA
        self.sigma_AA = box_obj.sigma_AA
        self.omega = box_obj.omega
        self.sigma_AA = box_obj.sigma_AA
        self.Q = box_obj.Q
        self.Q_ZA = CalcZonalAverage(self.Q,self.LonIndexer)
        self.Q_AA = CalcAreaAverage(self.Q,self.LatIndexer,
                                    LonIndexer=self.LonIndexer)
        self.Q_ZE = self.Q - self.Q_ZA
        self.Q_AE = self.Q_ZA - self.Q_AA
    
    def calc_gz(self):
        print('\nComputing generation of Zonal Potential Energy (Gz)...')
        _ = (self.Q_AE*self.tair_AE)/(Cp_d*self.sigma_AA)
        function = CalcAreaAverage(_,self.LatIndexer)
        Gz = VerticalTrazpezoidalIntegration(function,self.PressureData,
                                             self.VerticalCoordIndexer)
        try: 
            Gz = Gz.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Gz')
            raise
        if self.method == 'eulerian':
            print(Gz.values*Gz.metpy.units)
        # Save Gz before vertical integration
        print(Gz.values*Gz.metpy.units)
        print('Saving Gz for each vertical level...')
        if self.method == 'eulerian':
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Gz_'+self.VerticalCoordIndexer+'.csv',
                    mode="a", header=None)
        print('Done!')
        return Gz
    
    def calc_ge(self):
        print('\nComputing generation of Eddy Potential Energy (Ge)...')
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
        # Save Gz before vertical integration
        print(Ge.values*Ge.metpy.units)
        print('Saving Ge for each vertical level...')
        if self.method == 'eulerian':
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Ge_'+self.VerticalCoordIndexer+'.csv',
                    mode="a", header=None)
        print('Done!')
        return Ge
    
    def calc_dz(self):
        # Here we will use only the lowest vertical level
        _ = (self.u_ZA.isel({self.VerticalCoordIndexer:0})*self.ust_ZA) + (
            self.v_ZA.isel({self.VerticalCoordIndexer:0})*self.vst_ZA)
        Dz = units.Pa * CalcAreaAverage(_,self.LatIndexer)/g   
        return Dz
    
    def calc_de(self):
        # Here we will use only the lowest vertical level
        _ = (self.u_ZE.isel({self.VerticalCoordIndexer:0})*self.ust_ZE) + (
            self.v_ZE.isel({self.VerticalCoordIndexer:0})*self.vst_ZE)
        De = units.Pa * CalcAreaAverage(_,self.LatIndexer,
                                        LonIndexer=self.LonIndexer)/g   
        return De
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:27:40 2020

@author: Danilo
"""

from metpy.units import units
from metpy.constants import Cp_d
from metpy.constants import g
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
        self.Q_ZA = box_obj.Q_ZA
        self.Q_AA = box_obj.Q_AA
        self.Q_ZE = box_obj.Q_ZE
        self.Q_AE = box_obj.Q_AE
        self.xlength = box_obj.xlength
        self.ylength = box_obj.ylength
    
    def calc_gz(self):
        print('\nComputing generation of Zonal Potential Energy (Gz)...')
        _ = self.Q_AE*self.tair_AE
        _AA = (_*_["coslats"]).integrate("rlats")/self.ylength
        function = _AA/(Cp_d*self.sigma_AA)
        Gz = function.integrate(self.VerticalCoordIndexer
                    )* function[self.VerticalCoordIndexer].metpy.units
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
            df = function.to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.TimeName]).to_dataframe(name=time
                                                             ).transpose()
        df.to_csv(self.output_dir+'/Gz_'+self.VerticalCoordIndexer+'.csv',
                    mode="a", header=None)
        print('Done!')
        return Gz
    
    def calc_ge(self):
        print('\nComputing generation of Eddy Potential Energy (Ge)...')
        _ = self.Q_ZE*self.tair_ZE
        _ZA = _.integrate("rlons")/self.xlength
        _AA = (_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength
        function = _AA/(Cp_d*self.sigma_AA)
        Ge = function.integrate(self.VerticalCoordIndexer
                    )* function[self.VerticalCoordIndexer].metpy.units
        try: 
            Ge = Ge.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Ge')
            raise
        # Save Gz before vertical integration
        print(Ge.values*Ge.metpy.units)
        print('Saving Ge for each vertical level...')
        if self.method == 'eulerian':
            df = function.to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.TimeName]).to_dataframe(name=time
                                                             ).transpose()
        df.to_csv(self.output_dir+'/Ge_'+self.VerticalCoordIndexer+'.csv',
                    mode="a", header=None)
        print('Done!')
        return Ge
    
    def calc_dz(self):
        # Here we will use only the lowest vertical level
        _ = (self.u_ZA.isel({self.VerticalCoordIndexer:0})*self.ust_ZA) + (
            self.v_ZA.isel({self.VerticalCoordIndexer:0})*self.vst_ZA)
        _ZA = _.integrate("rlons")/self.xlength
        _AA = (_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength
        Dz = units.Pa *  _AA/g   
        return Dz
    
    def calc_de(self):
        # Here we will use only the lowest vertical level
        _ = (self.u_ZE.isel({self.VerticalCoordIndexer:0})*self.ust_ZE) + (
            self.v_ZE.isel({self.VerticalCoordIndexer:0})*self.vst_ZE)
        _ZA = _.integrate("rlons")/self.xlength
        _AA = (_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength
        De = units.Pa * _AA/g   
        return De
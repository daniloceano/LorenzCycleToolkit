#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:27:40 2020

@author: Danilo
"""
from calc_averages import CalcAreaAverage
from box_data import BoxData

from metpy.units import units
from metpy.constants import Cp_d
from metpy.constants import g

import pandas as pd
import numpy as np


class GenerationDissipationTerms:
    
    def __init__(self, box_obj: BoxData, method: str):
        self.method = method
        self.box_obj = box_obj
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
        print('Computing generation of Zonal Potential Energy (Gz)...')

        term = self.Q_AE * self.tair_AE
        function = CalcAreaAverage(term, self.ylength) / (Cp_d * self.sigma_AA)
        
        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Gz = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units

        try:
            Gz = Gz.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in Gz')
            raise ValueError('Unit error in Gz') from e

        if self.box_obj.args.verbosity:
            print(Gz.values * Gz.metpy.units)

        print('Saving Gz for each vertical level...')
        if self.method == 'fixed':
            df = function.to_dataframe(name='Ce', dim_order=[self.TimeName, self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()

        df.to_csv(self.output_dir + '/Gz_' + self.VerticalCoordIndexer + '.csv', mode="a", header=None)
        print('Done!')

        return Gz
    
    def calc_ge(self):
        print('Computing generation of Eddy Potential Energy (Ge)...')

        term = self.Q_ZE * self.tair_ZE
        function = CalcAreaAverage(term, self.ylength, xlength=self.xlength) / (Cp_d * self.sigma_AA)

        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Ge = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units

        try:
            Ge = Ge.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in Ge')
            raise ValueError('Unit error in Ge') from e

        if self.box_obj.args.verbosity:
            print(Ge.values * Ge.metpy.units)

        print('Saving Ge for each vertical level...')
        if self.method == 'fixed':
            df = function.to_dataframe(name='Ce', dim_order=[self.TimeName, self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()

        df.to_csv(self.output_dir + '/Ge_' + self.VerticalCoordIndexer + '.csv', mode="a", header=None)
        print('Done!')

        return Ge

    def calc_dz(self):
        # Here we will use only the lowest vertical level
        term = (self.u_ZA.isel({self.VerticalCoordIndexer: 0}) * self.ust_ZA) + (
            self.v_ZA.isel({self.VerticalCoordIndexer: 0}) * self.vst_ZA)
        function = CalcAreaAverage(term, self.ylength) / g
        Dz = units.Pa * function
        return Dz

    def calc_de(self):
        # Here we will use only the lowest vertical level
        term = (self.u_ZE.isel({self.VerticalCoordIndexer: 0}) * self.ust_ZE) + (
            self.v_ZE.isel({self.VerticalCoordIndexer: 0}) * self.vst_ZE)
        function = CalcAreaAverage(term, self.ylength) / g
        De = units.Pa * function
        return De
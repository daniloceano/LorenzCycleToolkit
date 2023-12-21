#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script defines the EnergyContent object. It uses the MetData object as an
input. The built-in functions use the input data to compute the following 
partioned energy contents of the Lorenz Energy Cycle:
    AZ = zonal available potential energy
    AE = eddy available potential energy
    KZ = zonal kinetic energy
    KE = eddy kinetic energy

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com


Source for formulas used here:
        Brennan, F. E., & Vincent, D. G. (1980).
        Zonal and Eddy Components of the Synoptic-Scale Energy Budget
        during Intensification of Hurricane Carmen (1974),
        Monthly Weather Review, 108(7), 954-965. Retrieved Jan 25, 2022, from:
        https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_0954_zaecot_2_0_co_2.xml
"""

from calc_averages import CalcAreaAverage
from box_data import BoxData

from metpy.constants import g
from metpy.units import units

import numpy as np
import pandas as pd
import xarray as xr

class EnergyContents:
    
    def __init__(self, box_obj: BoxData, method: str):
        self.method = method
        self.box_obj = box_obj
        self.PressureData = box_obj.PressureData
        self.LonIndexer = box_obj.LonIndexer
        self.LatIndexer = box_obj.LatIndexer
        self.TimeName = box_obj.TimeName
        self.VerticalCoordIndexer = box_obj.VerticalCoordIndexer
        self.output_dir = box_obj.output_dir
        self.tair_AE = box_obj.tair_AE
        self.tair_ZE = box_obj.tair_ZE
        self.u_ZA = box_obj.u_ZA
        self.u_ZE = box_obj.u_ZE
        self.v_ZA = box_obj.v_ZA
        self.v_ZE = box_obj.v_ZE
        self.sigma_AA = box_obj.sigma_AA
        self.xlength = box_obj.xlength
        self.ylength = box_obj.ylength
        
    def calc_az(self):
        print('Computing Zonal Available Potential Energy (Az)...')

        squared_tair = self.tair_AE ** 2
        function = CalcAreaAverage(squared_tair, self.ylength) / (2 * self.sigma_AA)

        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Az = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units

        try: 
            Az = Az.metpy.convert_units('J/m^2')
        except ValueError as e:
            raise ValueError('Unit error in Az') from e

        if self.box_obj.args.verbosity:
            print(Az.values * Az.metpy.units)

        Az = Az.metpy.dequantify()

        if self.box_obj.args.verbosity:
            print('Saving Az for each vertical level...')

        if self.method == 'fixed':
            df = function.to_dataframe(name='Az').unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()

        df.to_csv(
            f"{self.output_dir}/Az_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None
        )
        print('Done!')

        return Az
    
    def calc_ae(self):
        print('Computing Eddy Available Potential Energy (Ae)...')

        squared_tair = self.tair_ZE ** 2
        function = CalcAreaAverage(squared_tair, self.ylength, xlength=self.xlength) / (2 * self.sigma_AA)

        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Ae = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units

        try: 
            Ae = Ae.metpy.convert_units('J/m^2')
        except ValueError as e:
            print('Unit error in Ae')
            raise ValueError('Unit error in Ae') from e

        if self.box_obj.args.verbosity:
            print(Ae.values * Ae.metpy.units)

        Ae = Ae.metpy.dequantify()

        print('Saving Ae for each vertical level...')
        if self.method == 'fixed':
            df = function.to_dataframe(name='Ae').unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()
        
        df.to_csv(
            f"{self.output_dir}/Ae_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None
        )
        print('Done!')

        return Ae

    def calc_kz(self):
        print('Computing Zonal Kinetic Energy (Kz)...')

        function = CalcAreaAverage((self.u_ZA ** 2 + self.v_ZA ** 2), self.ylength)

        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Kz = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units / (2 * g)

        if self.box_obj.args.verbosity:
            print(Kz.values * Kz.metpy.units)

        print('Saving Kz for each vertical level...')

        try: 
            Kz = Kz.metpy.convert_units('J/m^2')
        except ValueError as e:
            print('Unit error in Kz')
            raise ValueError('Unit error in Kz') from e

        if self.box_obj.args.verbosity:
            print(Kz.values * Kz.metpy.units)
        
        Kz = Kz.metpy.dequantify()

        print('Saving Kz for each vertical level...')
        if self.method == 'fixed':
            df = function.to_dataframe(name='Kz').unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()
        
        df.to_csv(
            f"{self.output_dir}/Kz_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None
        )
        print('Done!')

        return Kz
    
    def calc_ke(self):
        print('Computing Eddy Kinetic Energy (Ke)...')

        function = CalcAreaAverage((self.u_ZE ** 2 + self.v_ZE ** 2), self.ylength, self.xlength)

        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Ke = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units / (2 * g)

        try: 
            Ke = Ke.metpy.convert_units('J/m^2')
        except ValueError as e:
            print('Unit error in Ke')
            raise ValueError('Unit error in Ke') from e

        if self.box_obj.args.verbosity:
            print(Ke.values * Ke.metpy.units)

        Ke = Ke.metpy.dequantify()

        print('Saving Ke for each vertical level...')
        if self.method == 'fixed':
            df = function.to_dataframe(name='Ke').unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()
        
        df.to_csv(
            f"{self.output_dir}/Ke_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None
        )
        print('Done!')

        return Ke

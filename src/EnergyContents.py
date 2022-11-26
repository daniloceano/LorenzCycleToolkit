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

from Math import (VerticalTrazpezoidalIntegration, CalcAreaAverage)
from metpy.constants import g
from BoxData import BoxData
import pandas as pd

class EnergyContents:
    
    def __init__(self, box_obj: BoxData, method: str):
        self.method = method
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
        
    def calc_az(self):
        print('\nComputing Zonal Available Potential Energy (Az)...')
        _ = CalcAreaAverage(self.tair_AE**2,self.LatIndexer)
        function = _/(2*self.sigma_AA)
        # Az = VerticalTrazpezoidalIntegration(function,self.PressureData,
        #                                      self.VerticalCoordIndexer)
        Az = -function.integrate(self.VerticalCoordIndexer
                            ) * function[self.VerticalCoordIndexer].metpy.units
        try: 
            Az = Az.metpy.convert_units('J/ m **2')
        except ValueError:
            raise ValueError('Unit error in Az')
        print(Az.values*Az.metpy.units)
        print('Saving Az for each vertical level...')
        # Save Az before vertical integration          
        if self.method == 'eulerian':
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Az_'+self.VerticalCoordIndexer+'.csv',
            mode="a", header=None)
        print('Done!')
        return Az
    
    def calc_ae(self): 
        print('\nComputing Eddy Available Potential Energy (Ae)...')
        _ = CalcAreaAverage(self.tair_ZE**2,self.LatIndexer,LonIndexer=self.LonIndexer)
        function = _/(2*self.sigma_AA)
        # Ae = VerticalTrazpezoidalIntegration(function,self.PressureData,
        #                                      self.VerticalCoordIndexer)
        Ae = -function.integrate(self.VerticalCoordIndexer
                            ) * function[self.VerticalCoordIndexer].metpy.units
        try: 
            Ae = Ae.metpy.convert_units('J/ m **2')
        except ValueError:
            print('Unit error in Ae')
            raise
        print(Ae.values*Ae.metpy.units)
        print('Saving Ae for each vertical level...')
        # Save Ae before vertical integration
        if self.method == 'eulerian':    
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Ae_'+self.VerticalCoordIndexer+'.csv',
                mode="a", header=None)
        print('Done!')
        return Ae
    
    def calc_kz(self):
        print('\nComputing Zonal Kinetic Energy (Kz)...')
        _ = (self.u_ZA**2)+(self.v_ZA**2)
        function = CalcAreaAverage(_,self.LatIndexer)
        # Kz = VerticalTrazpezoidalIntegration(function,self.PressureData,
        #                                      self.VerticalCoordIndexer)/(2*g)
        Kz = -function.integrate(self.VerticalCoordIndexer
                    ) * function[self.VerticalCoordIndexer].metpy.units/(2*g)
        print(Kz.values*Kz.metpy.units)
        print('Saving Kz for each vertical level...')
        try: 
            Kz = Kz.metpy.convert_units('J/ m **2')
        except ValueError:
            print('Unit error in Kz')
            raise
        # Save Kz before vertical integration
        if self.method == 'eulerian':
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Kz_'+self.VerticalCoordIndexer+'.csv',
                mode="a", header=None)
        print('Done!')
        return Kz
    
    def calc_ke(self):
        print('\nComputing Eddy Kinetic Energy (Ke)...')
        _ = (self.u_ZE**2)+(self.v_ZE**2)
        function = CalcAreaAverage(_,self.LatIndexer,LonIndexer=self.LonIndexer)
        # Ke = VerticalTrazpezoidalIntegration(function,self.PressureData,
        #                                      self.VerticalCoordIndexer)/(2*g)
        Ke = -function.integrate(self.VerticalCoordIndexer
                    ) * function[self.VerticalCoordIndexer].metpy.units/(2*g)
        try: 
            Ke = Ke.metpy.convert_units('J/ m **2')
        except ValueError:
            print('Unit error in Ke')
            raise
        print(Ke.values*Ke.metpy.units)
        print('Saving Ke for each vertical level...')
        # Save Ke before vertical integration
        if self.method == 'eulerian':
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Ke_'+self.VerticalCoordIndexer+'.csv',
                mode="a", header=None)
        print('Done!')
        return Ke
        
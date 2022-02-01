#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script defines the EnergyContent object. It uses the MetData object as an
input. The built-in functions use the input data to compute 
the following partioned energy contents:
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

from calc import (VerticalTrazpezoidalIntegration, CalcAreaAverage)
from metpy.constants import g
import MetData

class EnergyContents:
    def __init__(self, md_obj: MetData):
        self.PressureData = md_obj.PressureData
        self.LonIndexer = md_obj.LonIndexer
        self.LatIndexer = md_obj.LatIndexer
        self.TimeName = md_obj.TimeName
        self.VerticalCoordIndexer = md_obj.VerticalCoordIndexer
        self.tair_AE = md_obj.tair_AE
        self.tair_ZE = md_obj.tair_ZE
        self.u_ZA = md_obj.u_ZA
        self.u_ZE = md_obj.u_ZE
        self.v_ZA = md_obj.v_ZA
        self.v_ZE = md_obj.v_ZE
        self.sigma_AA = md_obj.sigma_AA
        
    def calc_az(self):
        numerator = CalcAreaAverage(self.tair_AE**2,self.LatIndexer)
        function = numerator/(2*self.sigma_AA)
        Az = VerticalTrazpezoidalIntegration(function,self.PressureData,
                                             self.VerticalCoordIndexer)
        try: 
            Az = Az.metpy.convert_units('J/ m **2')
        except ValueError:
            print('Unit error in Az')
            raise
        return Az
    
    def calc_ae(self):  
        numerator = CalcAreaAverage(self.tair_ZE**2,self.LatIndexer,self.LonIndexer)
        function = numerator/(2*self.sigma_AA)
        Az = VerticalTrazpezoidalIntegration(function,self.PressureData,
                                             self.VerticalCoordIndexer)
        try: 
            Ae = Az.metpy.convert_units('J/ m **2')
        except ValueError:
            print('Unit error in Az')
            raise
        return Ae
    
    def calc_kz(self):
        Wind_ZA = (self.u_ZA**2)+(self.v_ZA**2)
        Wind_AA = CalcAreaAverage(Wind_ZA,self.LatIndexer)
            
        Kz = VerticalTrazpezoidalIntegration(Wind_AA,self.PressureData,
                                             self.VerticalCoordIndexer)/(2*g)
        try: 
            Kz = Kz.metpy.convert_units('J/ m **2')
        except ValueError:
            print('Unit error in Kz')
            raise
        
        return Kz
    
    def calc_ke(self):
        Wind_ZE = (self.u_ZE**2)+(self.v_ZE**2)
        Wind_AA = CalcAreaAverage(Wind_ZE,self.LatIndexer,self.LonIndexer)
        Ke = VerticalTrazpezoidalIntegration(Wind_AA,self.PressureData,
                                             self.VerticalCoordIndexer)/(2*g)
        return Ke

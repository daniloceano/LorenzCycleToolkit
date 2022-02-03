#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script defines the ConversionTerms object. It uses the MetData object as
an input. The built-in functions use the input data to compute the following 
energy conversion terms of the Lorenz Energy Cycle:
    CE = Conversion between the two eddy energy terms (AE and KE)
    CZ = Conversion between the two zonal energy terms (AZ and KZ)
    CA = Conversion between the two available potential energy forms (AZ and AE)
    CK = Conversion between the two kinetic energy forms (KZ and KE)

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
import numpy as np
from metpy.constants import g
from metpy.constants import Rd
from metpy.constants import Re
from calc import (CalcAreaAverage,VerticalTrazpezoidalIntegration,Differentiate)
from MetData import MetData

class ConversionTerms:
    
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
        self.omega_ZE = md_obj.omega_ZE
        self.omega_AE = md_obj.omega_AE
        
        self.rlats = np.deg2rad(self.tair_AE[self.LatIndexer])
        self.cos_lats = np.cos(self.rlats)
        self.tan_lats = np.tan(self.rlats)
        
    def calc_ce(self):
        FirstTerm = Rd/(self.PressureData*g)
        _ = self.omega_ZE*self.tair_ZE
        SecondTerm = CalcAreaAverage(_,self.LonIndexer,self.LatIndexer)
        function = (FirstTerm*SecondTerm)
        Ce = -VerticalTrazpezoidalIntegration(function,self.PressureData,
                                             self.VerticalCoordIndexer)
        # Check if units are in accordance with expected and convert
        try: 
            Ce = Ce.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Ce')
            raise
        return Ce
    
    def calc_cz(self):
        FirstTerm = Rd/(self.PressureData*g)
        _ = self.omega_AE*self.tair_AE
        SecondTerm = CalcAreaAverage(_,self.LatIndexer)
        function = (FirstTerm*SecondTerm)
        Cz = -VerticalTrazpezoidalIntegration(function,self.PressureData,
                                             self.VerticalCoordIndexer)
        try: 
            Cz = Cz.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Cz')
            raise
        return Cz
    
    def calc_ca(self):
        # First term of the integral
        _ = (self.v_ZE*self.tair_ZE/(2*Re*self.sigma_AA)) \
            * Differentiate(self.tair_AE,self.rlats,self.LatIndexer)
        function = CalcAreaAverage(_,self.LatIndexer,self.LonIndexer)
        # Second term of the integral
        _ =  (self.omega_ZE*self.tair_ZE) \
            * Differentiate(self.tair_AE/self.sigma_AA,self.PressureData,self.VerticalCoordIndexer)
        function += CalcAreaAverage(_,self.LatIndexer,self.LonIndexer)
        # Integrate in pressure
        Ca = VerticalTrazpezoidalIntegration(function,self.PressureData,
                                             self.VerticalCoordIndexer)
        try: 
            Ca = Ca.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Ck')
            raise
        return Ca
        
    def calc_ck(self):
        # First term
        _ = ((self.cos_lats*self.u_ZE*self.v_ZE/Re)) \
            * Differentiate(self.u_ZA/self.cos_lats,self.rlats,self.LatIndexer)
        function = CalcAreaAverage(_,self.LatIndexer,self.LonIndexer)
        # Second term
        _ = ((self.v_ZE**2)/Re) \
            * Differentiate(self.v_ZA,self.rlats,self.LatIndexer)
        function += CalcAreaAverage(_,self.LatIndexer,self.LonIndexer)
        # Third term
        _ = (self.tan_lats*(self.u_ZE**2)*self.v_ZA)/Re
        function += CalcAreaAverage(_,self.LatIndexer,self.LonIndexer)
        # Fourth term
        _ = (self.omega_ZE*self.u_ZE) \
            * Differentiate(self.u_ZA,
                            self.PressureData,self.VerticalCoordIndexer)
        function += CalcAreaAverage(_,self.LatIndexer,self.LonIndexer)
        # Fifith term
        tmp1 = self.omega_ZE*self.v_ZE
        tmp2 = Differentiate(self.v_ZA,self.PressureData,self.VerticalCoordIndexer)
        function +=  CalcAreaAverage(tmp1*tmp2,self.LatIndexer,self.LonIndexer)
        # Integrate in pressure
        Ck = VerticalTrazpezoidalIntegration(function,self.PressureData,
                                             self.VerticalCoordIndexer)/g
        try: 
            Ck = Ck.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Ck')
            raise
        return Ck

# def Calc_Ce(PressureData,OmegaData,TemperatureData,LonIndexer,LatIndexer,VerticalCoordIndexer):
#     '''
#     Parameters
#     ----------
#     PressureData: xarray.Dataset
#         unit aware array for the pressure levels
#     OmegaData: xarray.Dataset
#         unit aware array for the omega velocity data
#     TemperatureData: xarray.Dataset
#         unit aware array for the temperature  
#     LonIndexer: string
#         the indexer used for the longitude variable in the xarray
#     LatIndexer: string 
#         the indexer used for the latitude variable in the xarray  
    
#     Returns
#     -------
#     Ce: xarray.DataArray
#         Unit-aware values corresponding to the conversion between eddy energy terms
#     '''        
#     FirstTerm = Rd/PressureData
#     # Zonal average terms
#     omega_ZA = CalcZonalAverage(OmegaData,LonIndexer)
#     tair_ZA = CalcZonalAverage(TemperatureData,LonIndexer)
#     # Zonal eddy terms
#     omega_ZE = OmegaData-omega_ZA
#     tair_ZE = TemperatureData-tair_ZA
#     TemperatureVerticalTransport = omega_ZE*tair_ZE
#     # Area averages
#     SecondTerm = CalcAreaAverage(TemperatureVerticalTransport,LonIndexer,LatIndexer)
#     # Integrate in pressure
#     function = FirstTerm*SecondTerm
#     Ce = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)/(-g)
#     # Check if units are in accordance with expected and convert
#     try: 
#         Ce = Ce.metpy.convert_units('W/ m **2')
#     except ValueError:
#         print('Unit error in Ce')
#         raise
#     return Ce

# def Calc_Cz(PressureData,OmegaData,TemperatureData,LonIndexer,LatIndexer,VerticalCoordIndexer):
#     '''
#     Parameters
#     ----------
#     PressureData: xarray.Dataset
#         unit aware array for the pressure levels
#     OmegaData: xarray.Dataset
#         unit aware array for the omega velocity data
#     TemperatureData: xarray.Dataset
#         unit aware array for the temperature  
#     LonIndexer: string
#         the indexer used for the longitude variable in the xarray
#     LatIndexer: string 
#         the indexer used for the latitude variable in the xarray  
    
#     Returns
#     -------
#     Cz: xarray.DataArray
#         Unit-aware values corresponding to the conversion between zonal energy terms
#     '''   
#     FirstTerm = Rd/PressureData
#     # Zonal average terms
#     omega_ZA = CalcZonalAverage(OmegaData,LonIndexer)
#     tair_ZA = CalcZonalAverage(TemperatureData,LonIndexer)
#     # Area average terms
#     omega_AA = CalcAreaAverage(omega_ZA,LatIndexer)
#     tair_AA = CalcAreaAverage(tair_ZA,LatIndexer)
#     # Area eddy terms
#     omega_AE = omega_ZA-omega_AA
#     tair_AE = tair_ZA-tair_AA
#     TemperatureVerticalTransport = omega_AE*tair_AE
#     # Area averages
#     SecondTerm = CalcAreaAverage(TemperatureVerticalTransport,LatIndexer)
#     # Integrate in pressure
#     function = FirstTerm*SecondTerm
#     Cz = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)/(-g)
#     # Check if units are in accordance with expected and convert
#     try: 
#         Cz = Cz.metpy.convert_units('W/ m **2')
#     except ValueError:
#         print('Unit error in Cz')
#         raise
#     return Cz
    
# def Calc_Ca(VWindComponentData,PressureData,OmegaData,TemperatureData,
#             LonIndexer,LatIndexer,VerticalCoordIndexer):
#     '''
    
#     This expression is a slightly modificated from the source material:
#         In the first term of the integral Brennan & Vincent write the denominator
#         as (sigma * Re), however, other sources use (2 * sigma * Re), such the
#         Muench's work (1965) and the  Michaelides' (1986). Therefore, here it
#         is adopted the latter expression.
    
#     Parameters
#     ----------
#     VWindComponentData: xarray.Dataset
#         unit aware array for the meridional component of the wind 
#     PressureData: xarray.Dataset
#         unit aware array for the pressure levels
#     OmegaData: xarray.Dataset
#         unit aware array for the omega velocity data
#     TemperatureData: xarray.Dataset
#         unit aware array for the temperature  
#     LonIndexer: string
#         the indexer used for the longitude variable in the xarray
#     LatIndexer: string 
#         the indexer used for the latitude variable in the xarray 
#     VerticalCoordIndexer: string
#         the indexer used for the xarray coordante which the integration will
#         be performed
    
#     Returns
#     -------
#     Ca: xarray.DataArray
#         Unit-aware values corresponding to the conversion between 
#         available potential energy terms
#     ''' 
#     rlats = np.deg2rad(VWindComponentData[LatIndexer])
#     # Zonal averages
#     v_ZA = CalcZonalAverage(VWindComponentData,LonIndexer)
#     tair_ZA = CalcZonalAverage(TemperatureData,LonIndexer)
#     omega_ZA = CalcZonalAverage(OmegaData,LonIndexer)
#     # Zonal eddies
#     v_ZE = VWindComponentData - v_ZA
#     tair_ZE = TemperatureData - tair_ZA
#     omega_ZE = OmegaData - omega_ZA
#     # Temperature area average
#     tair_AA = CalcAreaAverage(tair_ZA,LatIndexer)
#     # Temperature area eddy
#     tair_AE = tair_ZA-tair_AA
#     # Sigma
#     sigma_AA = StaticStability(TemperatureData,PressureData,VerticalCoordIndexer,
#                             LatIndexer,LonIndexer)
#     # First term of the integral
#     tmp1 = v_ZE*tair_ZE/(Re*2*sigma_AA)
#     tmp2 = Differentiate(tair_AE,rlats,LatIndexer)
#     FirstTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
#     # Second term of the integral
#     tmp1 =  omega_ZE*tair_ZE/(sigma_AA)
#     tmp2 = Differentiate(tair_AE,PressureData,VerticalCoordIndexer)
#     SecondTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
#     # Integrate in pressure
#     function = FirstTerm+SecondTerm
#     Ca = -VerticalTrazpezoidalIntegration(function,PressureData,
#                                           VerticalCoordIndexer)
#     # Check if units are in accordance with expected and convert
#     try: 
#         Ca = Ca.metpy.convert_units('W/ m **2')
#     except ValueError:
#         print('Unit error in Ck')
#         raise
#     return Ca

# def Calc_Ck(UWindComponentData,VWindComponentData,PressureData,OmegaData,
#             LonIndexer,LatIndexer,VerticalCoordIndexer):
#     '''
#      Parameters
#     ----------
#     VWindComponentData: xarray.Dataset
#         unit aware array for the meridional component of the wind 
#     PressureData: xarray.Dataset
#         unit aware array for the pressure levels
#     OmegaData: xarray.Dataset
#         unit aware array for the omega velocity data
#     TemperatureData: xarray.Dataset
#         unit aware array for the temperature  
#     LonIndexer: string
#         the indexer used for the longitude variable in the xarray
#     LatIndexer: string 
#         the indexer used for the latitude variable in the xarray 
#     VerticalCoordIndexer: string
#         the indexer used for the xarray coordante which the integration will
#         be performed
    
#     Returns
#     -------
#     Cz: xarray.DataArray
#         Unit-aware values corresponding to the conversion between kinetic energy terms
#     '''  
#     # Convert latitude do radians
#     rlats = np.deg2rad(UWindComponentData[LatIndexer])
#     # Zonal averages 
#     v_ZA = CalcZonalAverage(VWindComponentData,LonIndexer)
#     u_ZA = CalcZonalAverage(UWindComponentData,LonIndexer)
#     omega_ZA = CalcZonalAverage(OmegaData,LonIndexer)
#     # Zonal eddies
#     v_ZE = VWindComponentData - v_ZA
#     u_ZE = UWindComponentData - u_ZA
#     omega_ZE = OmegaData - omega_ZA
#     # First term
#     tmp1 = np.cos(rlats)*u_ZE*v_ZE/Re
#     tmp2 = Differentiate(u_ZA/np.cos(rlats),rlats,LatIndexer)
#     FirstTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
#     # Second term
#     tmp1 = (v_ZE**2)/Re
#     tmp2 = Differentiate(v_ZA,rlats,LatIndexer)
#     SecondTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
#     # Third term
#     tmp1 = (np.tan(rlats)*(u_ZE**2)*v_ZA)/Re
#     ThirdTerm = CalcAreaAverage(tmp1,LatIndexer,LonIndexer)
#     # Fourth term
#     tmp1 = omega_ZE*u_ZE
#     tmp2 = Differentiate(u_ZA,PressureData,VerticalCoordIndexer)
#     FourthTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
#     # Fifith term
#     tmp1 = omega_ZE*v_ZE
#     tmp2 = Differentiate(v_ZA,PressureData,VerticalCoordIndexer)
#     FifithTerm =  CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
#     # Integrate in pressure
#     function = (FirstTerm+SecondTerm+ThirdTerm+FourthTerm+FifithTerm)
#     Ck = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)/g
#     # Check if units are in accordance with expected and convert
#     try: 
#         Ck = Ck.metpy.convert_units('W/ m **2')
#     except ValueError:
#         print('Unit error in Ck')
#         raise
#     return Ck
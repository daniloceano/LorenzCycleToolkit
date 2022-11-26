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

Terms calculated using finite differences derivatives:
        First the original terms are copied and for the terms where
        it is needed to perform a spatial deriative, the lat/lon coordinate
        is converted from degrees to radians. It is necessary to organize the
        coordinate indexers from the smaller to the greater for the differentition
        process. For example, latitude from S to N and pressure from top to
        base. Then, it is used the Xarray's method for computating finite 
        differences derivatives using the central difference formula. Afterwards,
        the original lat/lon coords. are restaured, when needed. Lastly, the
        corresponding units are assigned, when applicable.

"""
import numpy as np
from metpy.constants import g
from metpy.constants import Rd
from metpy.constants import Re
from metpy.units import units
from Math import (CalcAreaAverage,VerticalTrazpezoidalIntegration)
from BoxData import BoxData
import pandas as pd

class ConversionTerms:
    
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
        self.omega_ZE = box_obj.omega_ZE
        self.omega_AE = box_obj.omega_AE
        
        self.rlats = np.deg2rad(box_obj.tair[self.LatIndexer])
        self.cos_lats = np.cos(self.rlats)
        self.tan_lats = np.tan(self.rlats)
        
    def calc_ce(self):
        print('\nComputing conversion between eddy energy terms (Ce)...')
        FirstTerm = Rd/(self.PressureData*g)
        _ = self.omega_ZE*self.tair_ZE
        SecondTerm = CalcAreaAverage(_,self.LatIndexer,LonIndexer=self.LonIndexer)
        function = (FirstTerm*SecondTerm)
        # Ce = VerticalTrazpezoidalIntegration(-function,self.PressureData,
        #                                      self.VerticalCoordIndexer)
        Ce = function.integrate(self.VerticalCoordIndexer
                           ) * function[self.VerticalCoordIndexer].metpy.units
        
        # Check if units are in accordance with expected and convert
        try: 
            Ce = Ce.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Ce')
            raise
        print(Ce.values*Ce.metpy.units)
        print('Saving Ce for each vertical level...')
        # Save Ce before vertical integration
        if self.method == 'eulerian':
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Ce_'+self.VerticalCoordIndexer+'.csv',
                    mode="a", header=None)
        print('Done!')
        return Ce
    
    def calc_cz(self):
        print('\nComputing conversion between zonal energy terms (Cz)...')
        FirstTerm = Rd/(self.PressureData*g)
        _ = self.omega_AE*self.tair_AE
        SecondTerm = CalcAreaAverage(_,self.LatIndexer)
        function = (FirstTerm*SecondTerm)
        # Cz = VerticalTrazpezoidalIntegration(-function,self.PressureData,
        #                                      self.VerticalCoordIndexer)
        Cz = function.integrate(self.VerticalCoordIndexer
                            ) * function[self.VerticalCoordIndexer].metpy.units
        try: 
            Cz = Cz.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Cz')
            raise
        print(Cz.values*Cz.metpy.units)
        print('Saving Cz for each vertical level...')
        # Save Cz before vertical integration
        if self.method == 'eulerian':
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Cz_'+self.VerticalCoordIndexer+'.csv',
                    mode="a", header=None)
        print('Done!')
        return Cz
    
    def calc_ca(self):
        print('\nComputing conversion between available potential energy terms (Ca)...')
        ## First term of the integral ##
        # Derivate tair_AE in respect to latitude
        DelPhi_tairAE = self.tair_AE.copy(deep=True
        ).assign_coords({self.LatIndexer:self.rlats}
                         ).sortby(self.LatIndexer,ascending=True
                        ).differentiate(self.LatIndexer).assign_coords(
                        {self.LatIndexer: self.tair_AE[self.LatIndexer]})
        _ = (self.v_ZE*self.tair_ZE/(2*Re*self.sigma_AA)) * DelPhi_tairAE
        function = CalcAreaAverage(_,self.LatIndexer,LonIndexer=self.LonIndexer)
        
        ## Second term of the integral ##
        # Derivate tair_AE in respect to pressure and divide it by sigma
        DelPres_tairAE_sigma = (self.tair_AE/self.sigma_AA).copy(deep=True
        ).sortby(self.VerticalCoordIndexer,ascending=True
        ).differentiate(self.VerticalCoordIndexer) / units.hPa
        _ =  (self.omega_ZE*self.tair_ZE) * DelPres_tairAE_sigma
        function += CalcAreaAverage(_,self.LatIndexer,LonIndexer=self.LonIndexer)
        
        ## Integrate in pressure ##
        # Ca = VerticalTrazpezoidalIntegration(-function,self.PressureData,
        #                                      self.VerticalCoordIndexer)
        Ca = function.integrate(self.VerticalCoordIndexer
                            ) * function[self.VerticalCoordIndexer].metpy.units
        try: 
            Ca = Ca.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Ca')
            raise
        print(Ca.values*Ca.metpy.units)
        print('Saving Ca for each vertical level...')
        if self.method == 'eulerian':
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Ca_'+self.VerticalCoordIndexer+'.csv',
                    mode="a", header=None)
        print('Done!')
        return Ca
        
    def calc_ck(self):
        print('\nComputing conversion between kinetic energy terms (Ck)...')
        ## First term ##
        # Divide the zonal mean of the zonal wind component (u) by the cosine
        # of the latitude (in radians) and then differentiate it in regard to
        # the latitude (also in radians)
        DelPhi_uZA_cosphi = (self.u_ZA/self.cos_lats).copy(deep=True
        ).assign_coords({self.LatIndexer:self.rlats}
                         ).sortby(self.LatIndexer,ascending=True
                        ).differentiate(self.LatIndexer).assign_coords(
                        {self.LatIndexer: self.u_ZA[self.LatIndexer]})
        _ = (self.cos_lats*self.u_ZE*self.v_ZE/Re) * DelPhi_uZA_cosphi
        function = CalcAreaAverage(_,self.LatIndexer,LonIndexer=self.LonIndexer)
        
        ## Second term ##
        # Differentiate the zonal mean of the meridional wind (v) in regard to
        # the latitude (in radians)
        DelPhi_vZA = self.v_ZA.copy(deep=True
        ).assign_coords({self.LatIndexer:self.rlats}
                          ).sortby(self.LatIndexer,ascending=True
                        ).differentiate(self.LatIndexer).assign_coords(
                        {self.LatIndexer: self.v_ZA[self.LatIndexer]})
        _ = ((self.v_ZE**2)/Re) * DelPhi_vZA
        function += CalcAreaAverage(_,self.LatIndexer,LonIndexer=self.LonIndexer)
        
        ## Third term ##
        _ = (self.tan_lats*(self.u_ZE**2)*self.v_ZA)/Re
        function += CalcAreaAverage(_,self.LatIndexer,LonIndexer=self.LonIndexer)
        
        ## Fourth term ##
        # Differentiate the zonal mean of the zonal wind (u) in regard to the
        # pressure and assign the units in hPa
        DelPres_uZAp = self.u_ZA.copy(deep=True
                         ).sortby(self.VerticalCoordIndexer,ascending=True
                        ).differentiate(self.VerticalCoordIndexer) / units.hPa
        _ = self.omega_ZE * self.u_ZE * DelPres_uZAp
        function += CalcAreaAverage(_,self.LatIndexer,LonIndexer=self.LonIndexer)
        
        ## Fifith term ##
        # Differentiate the zonal mean of the meridional wind (v) in regard to
        # the pressure and assign the units in hPa
        DelPres_vZAp = self.u_ZA.copy(deep=True
                        ).sortby(self.VerticalCoordIndexer,ascending=True
                        ).differentiate(self.VerticalCoordIndexer) / units.hPa
        _ = self.omega_ZE * self.v_ZE * DelPres_vZAp
        function +=  CalcAreaAverage(_,self.LatIndexer,LonIndexer=self.LonIndexer)
        
        ## Integrate in pressure ##
        # Ck = VerticalTrazpezoidalIntegration(function,self.PressureData,
        #                                      self.VerticalCoordIndexer)/g
        Ck = -function.integrate(self.VerticalCoordIndexer
                        ) * function[self.VerticalCoordIndexer].metpy.units/g
        try: 
            Ck = Ck.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in Ck')
            raise
        print(Ck.values*Ck.metpy.units)
        print('Saving Ck for each vertical level...')
        # Save Ck before vertical integration            
        if self.method == 'eulerian':
            df = function.drop([self.LonIndexer,self.LatIndexer]
                ).to_dataframe(name='Ce',dim_order=[
                    self.TimeName,self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop([self.LonIndexer,self.LatIndexer,self.TimeName]
                    ).to_dataframe(
                        name=time).transpose()
        df.to_csv(self.output_dir+'/Ck_'+self.VerticalCoordIndexer+'.csv',
                    mode="a", header=None)
        print('Done!')
        return Ck

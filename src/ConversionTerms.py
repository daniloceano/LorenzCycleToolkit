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
        base (this was done in the lorenz-cycle function, when opening data, 
        because other functions need this). Then, it is used the Xarray's 
        method for computating finite differences derivatives using the central
        difference formula. Lastly, the corresponding units are assigned, 
        when applicable.

"""
from calc_averages import CalcAreaAverage
import numpy as np
from metpy.constants import g
from metpy.constants import Rd
from metpy.constants import Re
from metpy.units import units
from box_data import BoxData
import pandas as pd

class ConversionTerms:
    
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
        self.omega_ZE = box_obj.omega_ZE
        self.omega_AE = box_obj.omega_AE
        self.tan_lats = np.tan(box_obj.tair["rlats"])
        self.xlength = box_obj.xlength
        self.ylength = box_obj.ylength
        
    def calc_ce(self):
        print('Computing conversion between eddy energy terms (Ce)...')

        FirstTerm = Rd / (self.PressureData * g)
        omega_tair_product = self.omega_ZE * self.tair_ZE
        SecondTerm = -CalcAreaAverage(omega_tair_product, self.ylength, xlength=self.xlength)
        function = FirstTerm * SecondTerm

        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Ce = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units

        try: 
            Ce = Ce.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in Ce')
            raise ValueError('Unit error in Ce') from e

        if self.box_obj.args.verbosity:
            print(Ce.values * Ce.metpy.units)

        print('Saving Ce for each vertical VerticalCoordIndexer...')
        if self.method == 'fixed':
            df = function.to_dataframe(name='Ce').unstack().transpose()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()

        df.to_csv(
            f"{self.output_dir}/Ce_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None
        )
        print('Done!')

        return Ce
    
    def calc_cz(self):
        print('Computing conversion between zonal energy terms (Cz)...')

        FirstTerm = Rd / (self.PressureData * g)
        omega_tair_product = self.omega_AE * self.tair_AE
        SecondTerm = CalcAreaAverage(omega_tair_product, self.ylength)
        function = FirstTerm * SecondTerm

        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Cz = -function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units

        try: 
            Cz = Cz.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in Cz')
            raise ValueError('Unit error in Cz') from e

        if self.box_obj.args.verbosity:
            print(Cz.values * Cz.metpy.units)

        print('Saving Cz for each vertical VerticalCoordIndexer...')
        if self.method == 'fixed':
            df = function.to_dataframe(name='Cz').unstack().transpose()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()

        df.to_csv(
            f"{self.output_dir}/Cz_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None
        )
        print('Done!')

        return Cz
    
    def calc_ca(self):
        print('Computing conversion between available potential energy terms (Ca)...')

        # First term of the integral
        DelPhi_tairAE = (self.tair_AE * self.tair_AE["coslats"]).differentiate("rlats")
        term1 = (self.v_ZE * self.tair_ZE * DelPhi_tairAE) / (2 * Re * self.sigma_AA)
        term1 = CalcAreaAverage(term1, self.ylength, xlength=self.xlength)

        # Second term of the integral
        DelPres_tairAE = (self.tair_AE).differentiate(self.VerticalCoordIndexer) / self.PressureData.metpy.units
        term2 = (self.omega_ZE * self.tair_ZE) * DelPres_tairAE
        term2 = CalcAreaAverage(term2, self.ylength, xlength=self.xlength) / self.sigma_AA

        function = term1 + term2

        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Ca = -function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units

        try: 
            Ca = Ca.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in Ca')
            raise ValueError('Unit error in Ca') from e

        if self.box_obj.args.verbosity:
            print(Ca.values * Ca.metpy.units)

        print('Saving Ca for each vertical VerticalCoordIndexer...')
        if self.method == 'fixed':
            df = function.to_dataframe(name='Ca', dim_order=[self.TimeName, self.VerticalCoordIndexer]).unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()

        df.to_csv(
            f"{self.output_dir}/Ca_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None
        )
        print('Done!')

        return Ca
        
    def calc_ck(self):
        print('Computing conversion between kinetic energy terms (Ck)...')

        DelPhi_uZA_cosphi = (self.u_ZA / self.u_ZA["coslats"]).differentiate("rlats")
        term1 = (self.u_ZE["coslats"] * self.u_ZE * self.v_ZE / Re) * DelPhi_uZA_cosphi
        term1 = CalcAreaAverage(term1, self.ylength, xlength=self.xlength)

        DelPhi_vZA = (self.v_ZA * self.v_ZA["coslats"]).differentiate("rlats")
        term2 = ((self.v_ZE**2) / Re) * DelPhi_vZA
        term2 = CalcAreaAverage(term2, self.ylength, xlength=self.xlength)

        term3 = (self.tan_lats * (self.u_ZE**2) * self.v_ZA) / Re
        term3 = CalcAreaAverage(term3, self.ylength, xlength=self.xlength)

        DelPres_uZAp = self.u_ZA.differentiate(self.VerticalCoordIndexer) / self.PressureData.metpy.units
        term4 = self.omega_ZE * self.u_ZE * DelPres_uZAp
        term4 = CalcAreaAverage(term4, self.ylength, xlength=self.xlength)

        DelPres_vZAp = self.u_ZA.differentiate(self.VerticalCoordIndexer) / self.PressureData.metpy.units
        term5 = self.omega_ZE * self.v_ZE * DelPres_vZAp
        term5 = CalcAreaAverage(term5, self.ylength, xlength=self.xlength)

        function = term1 + term2 + term3 + term4 + term5

        if np.isnan(function).any():
            function = function.interpolate_na(dim=self.VerticalCoordIndexer) * function.metpy.units
            if np.isnan(function).any():
                function = function.dropna(dim=self.VerticalCoordIndexer)

        Ck = function.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units / g

        try: 
            Ck = Ck.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in Ck')
            raise ValueError('Unit error in Ck') from e

        if self.box_obj.args.verbosity:
            print(Ck.values * Ck.metpy.units)

        print('Saving Ck for each vertical VerticalCoordIndexer...')
        if self.method == 'fixed':
            df = function.to_dataframe(name='Ck').unstack()
        else:
            time = pd.to_datetime(function[self.TimeName].data)
            df = function.drop(self.TimeName).to_dataframe(name=time).transpose()

        df.to_csv(
            f"{self.output_dir}/Ck_{self.VerticalCoordIndexer}.csv",
            mode="a",
            header=None
        )
        print('Done!')

        return Ck
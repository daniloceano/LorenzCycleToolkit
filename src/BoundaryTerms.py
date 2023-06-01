#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script defines the CBoundaryTerms object. It uses the MetData object as
an input. The built-in functions use the input data to compute the following 
boundary terms of the Lorenz Energy Cycle:
Computate energy fluxes across boundaries:
    BAz = zonal available potential energy
    BAe = eddy available potential energy
    BKz = zonal kinetic energy
    BKe = eddy kinetic energy
    
Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com


Source for formulas used here:
        Michaelides, S. C. (1987).
        Limited Area Energetics of Genoa Cyclogenesis,
        Monthly Weather Review, 115(1), 13-26. From:
        https://journals.ametsoc.org/view/journals/mwre/115/1/1520-0493_1987_115_0013_laeogc_2_0_co_2.xml

"""
from Math import (CalcZonalAverage, CalcAreaAverage)
import numpy as np
from metpy.constants import g
from metpy.constants import Re
from metpy.units import units
from BoxData import BoxData

def remove_nan(function, coordinate):
    function = function.interpolate_na(dim=coordinate) * function.metpy.units
    if np.isnan(function).any():
        function = function.dropna(dim=coordinate)
    return function

class BoundaryTerms:
    
    def __init__(self, box_obj: BoxData, method: str):
        self.method = method
        self.box_obj = box_obj
        self.PressureData = box_obj.PressureData
        self.LonIndexer = box_obj.LonIndexer
        self.LatIndexer = box_obj.LatIndexer
        self.TimeName = box_obj.TimeName
        self.VerticalCoordIndexer = box_obj.VerticalCoordIndexer
        self.tair_AE = box_obj.tair_AE
        self.tair_ZE = box_obj.tair_ZE
        self.u = box_obj.u
        self.u_ZA = box_obj.u_ZA
        self.u_ZE = box_obj.u_ZE
        self.v = box_obj.v
        self.v_ZA = box_obj.v_ZA
        self.v_ZE = box_obj.v_ZE
        self.omega = box_obj.omega
        self.sigma_AA = box_obj.sigma_AA
        self.omega_ZE = box_obj.omega_ZE
        self.omega_AA = box_obj.omega_AA
        self.omega_ZA = box_obj.omega_ZA
        self.omega_AE = box_obj.omega_AE
        self.geopt_ZE = box_obj.geopt_ZE
        self.geopt_AE = box_obj.geopt_AE
        self.western_limit = box_obj.western_limit
        self.eastern_limit = box_obj.eastern_limit
        self.southern_limit = box_obj.southern_limit
        self.northern_limit = box_obj.northern_limit
        self.xlength = box_obj.xlength
        self.ylength = box_obj.ylength        
        # Using the notation from Michaelides (1987)
        self.c1 = -1/(Re*self.xlength*self.ylength)
        self.c2 = -1/(Re*self.ylength)
        
    def calc_baz(self):
        print('Computing Zonal Available Potential Energy (Az) transport across boundaries (BAZ)...')

        # First Integral
        term1 = ((2 * self.tair_AE * self.tair_ZE * self.u) + (self.tair_AE ** 2 * self.u)) / (2 * self.sigma_AA)
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(**{self.LonIndexer: self.western_limit})
        term1 = term1.integrate("rlats")
        if np.isnan(term1).any():
            term1 = remove_nan(term1, self.VerticalCoordIndexer)
        term1 = term1.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c1

        # Second Integral
        term2 = self.v_ZE * self.tair_ZE
        term2 = CalcZonalAverage(term2, self.xlength) * 2 * self.tair_AE
        term2 = (term2 + (self.tair_AE ** 2 * self.v_ZA)) * self.tair_AE["coslats"]
        term2 = (term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(**{self.LatIndexer: self.southern_limit})) / (2 * self.sigma_AA)
        if np.isnan(term2).any():
            term2 = remove_nan(term2, self.VerticalCoordIndexer)
        term2 = term2.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c2

        # Third Term
        term3_1 = 2 * self.omega_ZE * self.tair_ZE
        term3_1 = CalcZonalAverage(term3_1, self.xlength) * self.tair_AE
        # Accordingly to Muench (1965), it is required to do an area average of this term and 
        # sigma is divided before subtracting the data (see next step)
        term3_2 = self.omega_ZA * self.tair_AE ** 2
        term3 = term3_1 + term3_2
        if np.isnan(term3).any():
            term3 = remove_nan(term3, self.VerticalCoordIndexer)
        term3 = CalcAreaAverage(term3, self.ylength) / (2 * self.sigma_AA)
        term3 = term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(**{self.VerticalCoordIndexer: 0})

        function = term1 + term2 - term3
        if np.isnan(function).any():
            function = remove_nan(function, self.VerticalCoordIndexer)

        try: 
            Baz = function.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in BAZ')
            raise ValueError('Unit error in BAZ') from e

        if self.box_obj.args.verbosity:
            print(Baz.values * Baz.metpy.units)

        return Baz

    def calc_bae(self):
        print('Computing Eddy Available Potential Energy (Ae) transport across boundaries (BAE)...')

        # First Integral
        term1 = self.u * self.tair_ZE ** 2
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(**{self.LonIndexer: self.western_limit})
        term1 = (term1 / (2 * self.sigma_AA)).integrate("rlats")
        if np.isnan(term1).any():
            term1 = remove_nan(term1, self.VerticalCoordIndexer)
        term1 = term1.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c1

        # Second Integral
        term2 = CalcZonalAverage(self.v * self.tair_ZE ** 2, self.xlength) * self.tair_ZE["coslats"]
        term2 = term2 / (2 * self.sigma_AA)
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(**{self.LatIndexer: self.southern_limit})
        if np.isnan(term2).any():
            term2 = remove_nan(term2, self.VerticalCoordIndexer)
        term2 = term2.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c2

        # Third Term
        term3 = (self.omega * self.tair_ZE ** 2) / (2 * self.sigma_AA)
        term3 = CalcAreaAverage(term3, self.ylength, xlength=self.xlength)
        if np.isnan(term3).any():
            term3 = remove_nan(term3, self.VerticalCoordIndexer)
        term3 = term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(**{self.VerticalCoordIndexer: 0})
        
        function = term1 + term2 - term3
        if np.isnan(term3).any():
            function = remove_nan(function, self.VerticalCoordIndexer)

        try: 
            Bae = function.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in BAE')
            raise ValueError('Unit error in BAE') from e

        if self.box_obj.args.verbosity:
            print(Bae.values * Bae.metpy.units)

        return Bae
    
    def calc_bkz(self):
        print('Computing Zonal Kinetic Energy (Kz) transport across boundaries (BKz)...')

        # First Integral
        term1 = self.u * (self.u ** 2 + self.v ** 2 - self.u_ZE ** 2 - self.v_ZE ** 2)
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(**{self.LonIndexer: self.western_limit})
        term1 = (term1 / (2 * g)).integrate("rlats")
        if np.isnan(term1).any():
            term1 = remove_nan(term1, self.VerticalCoordIndexer)
        term1 = term1.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c1

        # Second Integral
        term2 = (self.u ** 2 + self.v ** 2 - self.u_ZE ** 2 - self.v_ZE ** 2) * self.v * self.v["coslats"]
        term2 = CalcZonalAverage(term2, self.xlength)
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(**{self.LatIndexer: self.southern_limit})
        if np.isnan(term2).any():
            term2 = remove_nan(term2, self.VerticalCoordIndexer)
        term2 = (term2 / (2 * g)).integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c2

        # Third Term
        term3 = (self.u ** 2 + self.v ** 2 - self.u_ZE ** 2 - self.v_ZE ** 2) * self.omega
        term3 = CalcAreaAverage(term3, self.ylength, xlength=self.xlength) / (2 * g)
        if np.isnan(term3).any():
            term3 = remove_nan(term3, self.VerticalCoordIndexer)
        term3 = (term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(**{self.VerticalCoordIndexer: 0}))

        function = term1 + term2 - term3
        if np.isnan(function).any():
            function = remove_nan(function, self.VerticalCoordIndexer)

        try:
            Bkz = function.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in BKz')
            raise ValueError('Unit error in BKz') from e

        if self.box_obj.args.verbosity:
            print(Bkz.values * Bkz.metpy.units)

        return Bkz

    def calc_bke(self):
        print('Computing Eddy Kinetic Energy (Ke) transport across boundaries (BKe)...')

        # First Integral
        term1 = self.u * (self.u_ZE ** 2 + self.v_ZE ** 2)
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(**{self.LonIndexer: self.western_limit})
        term1 = (term1 / (2 * g)).integrate("rlats")
        if np.isnan(term1).any():
            term1 = remove_nan(term1, self.VerticalCoordIndexer)
        term1 = term1.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c1

        # Second Integral
        term2 = (self.u_ZE ** 2 + self.v_ZE ** 2) * self.v * self.v["coslats"]
        term2 = CalcZonalAverage(term2, self.xlength)
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(**{self.LatIndexer: self.southern_limit})
        if np.isnan(term2).any():
            term2 = remove_nan(term2, self.VerticalCoordIndexer)
        term2 = (term2 / (2 * g)).integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c2

        # Third Term
        term3 = (self.u_ZE ** 2 + self.v_ZE ** 2) * self.omega
        term3 = CalcAreaAverage(term3, self.ylength, xlength=self.xlength) / (2 * g)
        if np.isnan(term3).any():
            term3 = remove_nan(term3, self.VerticalCoordIndexer)
        term3 = (term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(**{self.VerticalCoordIndexer: 0}))

        function = term1 + term2 - term3
        if np.isnan(function).any():
            function = remove_nan(function, self.VerticalCoordIndexer)

        try:
            Bke = function.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in BKe')
            raise ValueError('Unit error in BKe') from e

        if self.box_obj.args.verbosity:
            print(Bke.values * Bke.metpy.units)

        return Bke
 
    def calc_boz(self):
        print('Computing Zonal Kinetic Energy (Kz) production by fluxes at the boundaries (BΦZ)...')

        # First Integral
        term1 = (self.v_ZA * self.geopt_AE) / g
        # Cannot perform eastern boundary minus western boundary
        term1 = term1.integrate("rlats")
        if np.isnan(term1).any():
            term1 = remove_nan(term1, self.VerticalCoordIndexer)
        term1 = term1.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c1

        # Second Integral
        term2 = (self.v_ZA * self.geopt_AE) * self.v_ZA["coslats"] / g
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(**{self.LatIndexer: self.southern_limit})
        if np.isnan(term2).any():
            term2 = remove_nan(term2, self.VerticalCoordIndexer)
        term2 = term2.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c2

        # Third Term
        term3 = CalcAreaAverage(self.omega_AE * self.geopt_AE, self.ylength) / g
        if np.isnan(term3).any():
            term3 = remove_nan(term3, self.VerticalCoordIndexer)
        term3 = (term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(**{self.VerticalCoordIndexer: 0}))

        function = term1 + term2 - term3
        if np.isnan(function).any():
            function = remove_nan(function, self.VerticalCoordIndexer)

        try:
            Boz = function.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in BΦZ')
            raise ValueError('Unit error in BΦZ') from e

        if self.box_obj.args.verbosity:
            print(Boz.values * Boz.metpy.units)

        return Boz

    def calc_boe(self):
        print('Computing Eddy Kinetic Energy (Ke) production by fluxes at the boundaries (BΦE)...')

        # First Integral
        term1 = (self.v_ZE * self.geopt_AE) / g
        term1 = term1.sel(**{self.LonIndexer: self.eastern_limit}) - term1.sel(**{self.LonIndexer: self.western_limit})
        term1 = term1.integrate("rlats")
        if np.isnan(term1).any():
            term1 = remove_nan(term1, self.VerticalCoordIndexer)
        term1 = term1.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c1

        # Second Integral
        term2 = (self.v_ZA * self.geopt_AE) * self.v_ZA["coslats"] / g
        if np.isnan(term2).any():
            term2 = remove_nan(term2, self.VerticalCoordIndexer)
        term2 = term2.sel(**{self.LatIndexer: self.northern_limit}) - term2.sel(**{self.LatIndexer: self.southern_limit})
        term2 = term2.integrate(self.VerticalCoordIndexer) * self.PressureData.metpy.units * self.c2

        # Third Term
        term3 = CalcAreaAverage(self.omega_AE * self.geopt_AE, self.ylength) / g
        if np.isnan(term3).any():
            term1 = remove_nan(term3, self.VerticalCoordIndexer)
        term3 = (term3.isel(**{self.VerticalCoordIndexer: -1}) - term3.isel(**{self.VerticalCoordIndexer: 0}))

        function = term1 + term2 - term3
        if np.isnan(function).any():
                    function = remove_nan(function, self.VerticalCoordIndexer)

        try:
            Boe = function.metpy.convert_units('W/m^2')
        except ValueError as e:
            print('Unit error in BΦE')
            raise ValueError('Unit error in BΦE') from e

        if self.box_obj.args.verbosity:
            print(Boe.values * Boe.metpy.units)

        return Boe
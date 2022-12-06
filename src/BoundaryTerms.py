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

import numpy as np
from metpy.constants import g
from metpy.constants import Re
from metpy.units import units
from BoxData import BoxData

class BoundaryTerms:
    
    def __init__(self, box_obj: BoxData, method: str):
        self.method = method
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
        self.c1 = -1/((Re*(
                np.deg2rad(box_obj.eastern_limit)-np.deg2rad(box_obj.western_limit))*
                (np.sin(np.deg2rad(box_obj.northern_limit))-
                np.sin(np.deg2rad(box_obj.southern_limit))))/units.radian)
        self.c2 = -1/(Re*(np.sin(np.deg2rad(box_obj.northern_limit))-
            np.sin(np.deg2rad(box_obj.southern_limit))))
        
    def calc_baz(self):
        if self.method == 'eulerian':
            print('\nComputing Zonal Available Potential Energy (Az) transport across boundaries (BAZ)...')

        ## First Integral ##
        _ = ((2*self.tair_AE*self.tair_ZE*self.u)
              + (self.tair_AE**2*self.u))/(2*self.sigma_AA)
        # Data at eastern boundary minus data at western boundary 
        _ = _.sel(**{self.LonIndexer: self.eastern_limit}) - _.sel(
            **{self.LonIndexer: self.western_limit})
        # Integrate through latitude
        _ = _.integrate("rlats")
        # Integrate through pressure levels
        function = _.integrate(self.VerticalCoordIndexer
                    ) * _[self.VerticalCoordIndexer].metpy.units * self.c1

        ## Second Integral ##
        
        _ZA = ((self.u_ZE*self.tair_ZE).integrate("rlons")/self.xlength)
        _AA = (_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength
        _ = ((2*_AA*self.tair_AE) + (self.tair_AE**2*self.v_ZA)
             * self.tair_AE["coslats"])/(2*self.sigma_AA)
        # Data at northern boundary minus data at southern boundary
        _ = _.sel(**{self.LatIndexer: self.northern_limit}) - _.sel(
            **{self.LatIndexer: self.southern_limit})
        # Integrate through pressure levels
        function += _.sortby(self.VerticalCoordIndexer,ascending=True
                             ).integrate(self.VerticalCoordIndexer
                    )* _[self.VerticalCoordIndexer].metpy.units*self.c2
        
        ## Third Term ##
        # First part
        _ZA = ((2*self.omega_ZE*self.tair_ZE).integrate("rlons")/self.xlength)
        _AA = (_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength
        _ = _AA*self.tair_AE
        # Second part
        _ZA = ((2*self.omega_AA*self.tair_ZE**2).integrate("rlons")/self.xlength)
        _AA = (_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength
        _ += _AA
        # Accordingly to Muench (1965), it is required to do an area average
        # of this term
        _AA = (_.integrate("rlats")/self.ylength)/(2*self.sigma_AA)
        # Data at the bottom minus data at top level
        function -= (_AA.isel(**{self.VerticalCoordIndexer: -1}) - 
                     _AA.isel(**{self.VerticalCoordIndexer: 0}))
        try: 
            Baz = function.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in BAZ')
            raise
        if self.method == 'eulerian':
            print(Baz.values*Baz.metpy.units)
        return Baz
    
    def calc_bae(self):
        if self.method == 'eulerian':
            print('\nComputing Eddy Available Potential Energy (Ae) transport across boundaries (BAE)...')
        ## First Integral ##
        _ = self.u*self.tair_ZE**2
         # Data at eastern boundary minus data at western boundary 
        _ = _.sel(**{self.LonIndexer: self.eastern_limit}) - _.sel(
            **{self.LonIndexer: self.western_limit}) 
        # Integrate through latitude
        _ = (_/(2*self.sigma_AA)).integrate("rlats")
        
        # Integrate through pressure levels
        function = _.integrate(self.VerticalCoordIndexer
                    )* _[self.VerticalCoordIndexer].metpy.units*self.c1
        
        ## Second Integral ##
        _ZA = ((self.v*self.tair_ZE**2).integrate("rlons")/self.xlength
               *self.tair_ZE["coslats"])/(2*self.sigma_AA)
        # Data at northern boundary minus data at southern boundary
        _ = _ZA.sel(**{self.LatIndexer: self.northern_limit}) - _ZA.sel(
            **{self.LatIndexer: self.southern_limit})
        # Integrate through pressure levels
        function += _.integrate(self.VerticalCoordIndexer
                    ) * _[self.VerticalCoordIndexer].metpy.units*self.c2
        
        ## Third Term ##
        _ZA = (self.omega*self.tair_ZE**2).integrate("rlons")/self.xlength
        _AA = (-(_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength 
               /(2*self.sigma_AA))
        _ = _AA.sortby(self.VerticalCoordIndexer)
        # Data at the bottom minus data at top level
        function -= (_.isel(**{self.VerticalCoordIndexer: -1}) - _AA.isel(
            **{self.VerticalCoordIndexer: 0}))
        try: 
            Bae = function.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in BAe')
            raise
        if self.method == 'eulerian':
            print(Bae.values*Bae.metpy.units)
        return Bae
    
    def calc_bkz(self):
        if self.method == 'eulerian':
            print('\nComputing Zonal Kinetic Energy (Kz) transport across boundaries (BKz)...')
        ## First Integral ##
        _ = self.u*(self.u**2+self.v**2-self.u_ZE**2-self.v_ZE**2)
         # Data at eastern boundary minus data at western boundary 
        _ = _.sel(**{self.LonIndexer: self.eastern_limit}) - _.sel(
            **{self.LonIndexer: self.western_limit}) 
        # Integrate through latitude
        _ = (_/(2*g)).integrate("rlats")
        # Integrate through pressure levels
        function = _.integrate(self.VerticalCoordIndexer
                    )* _[self.VerticalCoordIndexer].metpy.units*self.c1
        
        ## Second Integral ##
        _ = ((self.u**2+self.v**2-self.u_ZE**2-self.v_ZE**2)
            *self.v*self.v["coslats"]).integrate(
                "rlons")/self.xlength
        # Data at northern boundary minus data at southern boundary
        _ = _.sel(**{self.LatIndexer: self.northern_limit}) - _.sel(
            **{self.LatIndexer: self.southern_limit})
        # Integrate through pressure levels
        function += (_/(2*g)).integrate(self.VerticalCoordIndexer
                    ) * _[self.VerticalCoordIndexer].metpy.units*self.c2
        
        ## Third Term ##
        _ZA  = (((self.u**2+self.v**2-self.u_ZE**2-self.v_ZE**2)
            *self.omega).integrate("rlons")/self.xlength)
        _AA = (_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength
        _ = _AA/(2*g)
        # Data at the bottom minus data at top level
        function -= (_.isel(**{self.VerticalCoordIndexer: -1}) - _.isel(
            **{self.VerticalCoordIndexer: 0}))
        try: 
            Bkz = function.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in BKz')
            raise
        if self.method == 'eulerian':
            print(Bkz.values*Bkz.metpy.units)
        return Bkz
    
    def calc_bke(self):
        if self.method == 'eulerian':
            print('\nComputing Eddy Kinetic Energy (Ke) transport across boundaries (BKe)...')
        ## First Integral ##
        _ = self.u*(self.u_ZE**2+self.v_ZE**2)/(2*g)
         # Data at eastern boundary minus data at western boundary 
        _ = _.sel(**{self.LonIndexer: self.eastern_limit}) - _.sel(
            **{self.LonIndexer: self.western_limit}) 
        # Integrate through latitude
        _ = _.integrate("rlats")
        function = _.integrate(self.VerticalCoordIndexer
                    )* _[self.VerticalCoordIndexer].metpy.units*self.c1
        
        ## Second Integral ##
        _ = ((self.u_ZE**2+self.v_ZE**2)*self.v
               *self.v["coslats"]).integrate("rlons")/self.xlength
        # Data at northern boundary minus data at southern boundary
        _ = _.sel(**{self.LatIndexer: self.northern_limit}) - _.sel(
            **{self.LatIndexer: self.southern_limit})
        # Integrate through pressure levels
        function += (_/(2*g)).integrate(self.VerticalCoordIndexer
                    ) * _[self.VerticalCoordIndexer].metpy.units*self.c2
        
        ## Third Term ##
        _ZA = (((self.u_ZE**2+self.v_ZE**2)*self.omega
               ).integrate("rlons")/self.xlength)
        _AA = -(_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength
        _ = _AA/(2*g)
        function -= (_.isel(**{self.VerticalCoordIndexer: -1}) - _.isel(
                                   **{self.VerticalCoordIndexer: 0}))
        try: 
            Bke = function.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in BKe')
            raise
        if self.method == 'eulerian':
            print(Bke.values*Bke.metpy.units)
        return Bke 
    
    def calc_boz(self):
        """
        Used Brennan (1980) here.
        """
        if self.method == 'eulerian':
            print('\nComputing Zonal Kinetic Energy (Kz) production by fluxes at the boundaries (BΦZ)...')
        ## First Integral ##
        _ = (self.v_ZA*self.geopt_AE)/g
         # Data at eastern boundary minus data at western boundary 
        # _ = _.sel(**{self.LonIndexer: self.eastern_limit}) - _.sel(
        #     **{self.LonIndexer: self.western_limit}) 
        # Integrate through latitude
        # _ = HorizontalTrazpezoidalIntegration(_,self.LatIndexer)
        _ = _.integrate("rlats")
        # Integrate through pressure levels
        # function = VerticalTrazpezoidalIntegration(_,self.PressureData,
        #                             self.VerticalCoordIndexer)*self.c1
        function  = -_.integrate(self.VerticalCoordIndexer
                    )* _[self.VerticalCoordIndexer].metpy.units*self.c1
        
        ## Second Integral ##
        _ = (self.v_ZA*self.geopt_AE)*self.v_ZA["coslats"]/g
        # Data at northern boundary minus data at southern boundary
        _ = _.sel(**{self.LatIndexer: self.northern_limit}) - _.sel(
            **{self.LatIndexer: self.southern_limit})
        # Integrate through pressure levels
        # function += VerticalTrazpezoidalIntegration(_,self.PressureData,
        #                             self.VerticalCoordIndexer)*self.c2
        function += -_.integrate(self.VerticalCoordIndexer
                    )* _[self.VerticalCoordIndexer].metpy.units*self.c2
        
        ## Third Term ##
        # _ = CalcAreaAverage(self.omega_AE*self.geopt_AE, self.LatIndexer)/g
        _AA = (-(self.omega_AE*self.geopt_AE*self.geopt_AE["coslats"]
                 ).integrate("rlats")/self.ylength)/g  
        _ = _AA.sortby(self.VerticalCoordIndexer,ascending=False)        
        function -= _.sortby(self.VerticalCoordIndexer,ascending=False
        ).isel(**{self.VerticalCoordIndexer:-1}) - _.isel(
            **{self.VerticalCoordIndexer: 0})
        try: 
            Boz = function.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in BΦZ')
            raise
        if self.method == 'eulerian':
            print(Boz.values*Boz.metpy.units)
        return Boz 
    
    def calc_boe(self):
        if self.method == 'eulerian':
            print('\nComputing Eddy Kinetic Energy (Kz) production by fluxes at the boundaries (BΦE)...')
        ## First Integral ##
        _ = (self.u_ZE*self.geopt_ZE)/g
         # Data at eastern boundary minus data at western boundary 
        _ = _.sel(**{self.LonIndexer: self.eastern_limit}) - _.sel(
            **{self.LonIndexer: self.western_limit}) 
        # Integrate through latitude
        # _ = HorizontalTrazpezoidalIntegration(_,self.LatIndexer)
        _ = _.integrate("rlats")
        # Integrate through pressure levels
        # function = VerticalTrazpezoidalIntegration(_,self.PressureData,
        #                             self.VerticalCoordIndexer)*self.c1
        function = -_.integrate(self.VerticalCoordIndexer
                    )* _[self.VerticalCoordIndexer].metpy.units*self.c1
        
        ## Second Integral ##
        # _ = CalcZonalAverage((self.v_ZE*self.geopt_ZE),
        #                      self.LonIndexer)*self.cos_lats/g
        _ = ((((self.v_ZE*self.geopt_ZE)).integrate("rlons")/self.xlength
               )*self.v_ZE["coslats"])/g
        # Data at northern boundary minus data at southern boundary
        _ = _.sel(**{self.LatIndexer: self.northern_limit}) - _.sel(
            **{self.LatIndexer: self.southern_limit})
        # Integrate through pressure levels
        # function += VerticalTrazpezoidalIntegration(_,self.PressureData,
        #                             self.VerticalCoordIndexer)*self.c2
        function += -_.integrate(self.VerticalCoordIndexer
                    )* _[self.VerticalCoordIndexer].metpy.units*self.c2
        
        ## Third Term ##
        # _ = CalcAreaAverage(self.omega_ZE*self.geopt_ZE, self.LatIndexer,
        #                     LonIndexer=self.LonIndexer)/g
        _ZA = ((self.omega_ZE*self.geopt_ZE).integrate("rlons")/self.xlength)/g
        _AA = -(_ZA*_ZA["coslats"]).integrate("rlats")/self.ylength
        _ = _AA.sortby(self.VerticalCoordIndexer,ascending=False)  
        function -= _.isel(**{self.VerticalCoordIndexer: -1}) - _.isel(
            **{self.VerticalCoordIndexer: 0})
        # function -= _.sortby(self.VerticalCoordIndexer,ascending=False
        # ).isel(**{self.VerticalCoordIndexer: -1}) - _.isel(
        #     **{self.VerticalCoordIndexer: 0})
        try: 
            Boz = function.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in BΦZ')
            raise
        if self.method == 'eulerian':
            print(Boz.values*Boz.metpy.units)
        return Boz 
                 
        
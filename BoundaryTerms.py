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
from metpy.units import units
from calc import (CalcAreaAverage,VerticalTrazpezoidalIntegration,
                  Differentiate, HorizontalTrazpezoidalIntegration,
                  CalcZonalAverage)
from BoxData import BoxData
from EnergyContents import function_to_df

class BoundaryTerms:
    
    def __init__(self, box_obj: BoxData):
        self.PressureData = box_obj.PressureData
        self.LonIndexer = box_obj.LonIndexer
        self.LatIndexer = box_obj.LatIndexer
        self.TimeName = box_obj.TimeName
        self.VerticalCoordIndexer = box_obj.VerticalCoordIndexer
        self.output_dir = box_obj.output_dir
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
        self.omega_ZA = box_obj.omega_ZA
        self.omega_AE = box_obj.omega_AE
        
        self.rlats = np.deg2rad(self.tair_AE[self.LatIndexer])
        self.rlons = np.deg2rad(self.tair_AE[self.LonIndexer])
        self.sin_lats = np.sin(self.rlats)
        self.cos_lats = np.cos(self.rlats)
        self.tan_lats = np.tan(self.rlats)
        
        self.BoxWest = box_obj.BoxWest
        self.BoxEast = box_obj.BoxEast
        self.BoxSouth = box_obj.BoxSouth
        self.BoxNorth = box_obj.BoxNorth
        
        # Using the notation from Michaelides (1987)
        self.c1 = 1/((Re*(
                np.deg2rad(box_obj.BoxEast)-np.deg2rad(box_obj.BoxWest))*
                (np.sin(np.deg2rad(box_obj.BoxNorth))-
                np.sin(np.deg2rad(box_obj.BoxSouth))))/units.radian)
        self.c2 = 1/(Re*(np.sin(np.deg2rad(box_obj.BoxNorth))-
            np.sin(np.deg2rad(box_obj.BoxSouth))))
        
    def calc_baz(self):
        print('\nComputing Zonal Available Potential energy (Az) transport\
              across boundaries (BAZ)...')
             # needs revision
        ## First Integral ##
        _ = ((2*self.tair_AE*self.tair_ZE*self.u)
             + (self.tair_AE**2*self.u))/(2*self.sigma_AA)
        # Data at eastern boundary minus data at western boundary 
        _ = _.sel(**{self.LonIndexer: self.BoxEast}) - _.sel(
            **{self.LonIndexer: self.BoxWest})
        # Integrate through latitude
        _ = HorizontalTrazpezoidalIntegration(_,self.LatIndexer)
        # Integrate through pressure levels
        function = -VerticalTrazpezoidalIntegration(_,self.PressureData,
                                    self.VerticalCoordIndexer)*self.c1

        ## Second Integral ##
        FirstTerm = 2*self.tair_AE*CalcZonalAverage(self.u_ZE*self.tair_ZE,
                                LonIndexer=self.LonIndexer)
        SecondTerm = self.tair_AE**2*self.v_ZA
        _ = (FirstTerm + SecondTerm)*np.cos(self.rlons)/(2*self.sigma_AA)
        # Data at northern boundary minus data at southern boundary
        _ = _.sel(**{self.LatIndexer: self.BoxNorth}) - _.sel(
            **{self.LatIndexer: self.BoxSouth})
        # Integrate through pressure levels
        function += -VerticalTrazpezoidalIntegration(_,self.PressureData,
                                    self.VerticalCoordIndexer)*self.c2
        
        ## Third Term ##
        # In this case, there is an alteration from the formula originally
        # contained in Brennan, F. E., & Vincent, D. G. (1980). The original
        # formula does not indicate to do the area average, however if this
        # is not done, it is impossible to perform the required operations.
        # Consulting the paper from Muench (1965), one can see that for this
        # term this area averagin is performed, so here we just assumed that
        # the absence of the area averaging  operator is a typo.
        _ = (CalcZonalAverage(self.omega_ZE*self.tair_ZE,self.LonIndexer)*
             self.tair_AE*2) + (self.tair_AE**2*self.omega_ZA)
        _ = CalcAreaAverage(_,self.LatIndexer)/(2*self.sigma_AA)
        # Sort the data from bottom to top and then do data at the bottom
        # minus data at top level
        function += -_.sortby(self.VerticalCoordIndexer,ascending=False
        ).isel(**{self.VerticalCoordIndexer: 0}) - _.isel(
            **{self.VerticalCoordIndexer: -1})
        try: 
            Baz = function.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in BAZ')
            raise
        print(Baz.values*Baz.metpy.units)
        return Baz
    
    def calc_bae(self):
        print('\nComputing Eddy Available Potential energy (Ae) transport\
              across boundaries (BAE)...')
        ## First Integral ##
        _ = (self.u*self.tair_ZE**2)/(2*self.sigma_AA)
         # Data at eastern boundary minus data at western boundary 
        _ = _.sel(**{self.LonIndexer: self.BoxEast}) - _.sel(
            **{self.LonIndexer: self.BoxWest}) 
        # Integrate through latitude
        _ = HorizontalTrazpezoidalIntegration(_,self.LatIndexer)
        # Integrate through pressure levels
        function = -VerticalTrazpezoidalIntegration(_,self.PressureData,
                                    self.VerticalCoordIndexer)*self.c1
        
        ## Second Integral ##
        _ = CalcZonalAverage(self.v*self.tair_ZE**2,self.LonIndexer
                             )*np.cos(self.rlons)/(2*self.sigma_AA)
        # Data at northern boundary minus data at southern boundary
        _ = _.sel(**{self.LatIndexer: self.BoxNorth}) - _.sel(
            **{self.LatIndexer: self.BoxSouth})
        # Integrate through pressure levels
        function += -VerticalTrazpezoidalIntegration(_,self.PressureData,
                                    self.VerticalCoordIndexer)*self.c2
        
        ## Third Term ##
        _ = CalcAreaAverage(self.omega*self.tair_ZE**2,self.LatIndexer,
                            LonIndexer=self.LonIndexer)/(2*self.sigma_AA)
        function += -_.sortby(self.VerticalCoordIndexer,ascending=False
        ).isel(**{self.VerticalCoordIndexer: 0}) - _.isel(
            **{self.VerticalCoordIndexer: -1})
        try: 
            Bae = function.metpy.convert_units('W/ m **2')
        except ValueError:
            print('Unit error in BAZ')
            raise
        print(Bae.values*Bae.metpy.units)
        return Bae
                 
        
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 20:15:59 2022

Create an object that will store, whithin a bounding box specified by the user,
all meteorological data required by the Lorenz Energy Cycle computations.

The program uses xarray built-in function for integrations. For using it, it 
would be necessary to sort the data from south to north and from west to east,
but for saving up memory, it was using a negative sign in front of the 
integrations for the latitude coordinate.

Also, it was created a dimension for the radians of the latitude and longitude,
for the integration

Suffixes used:
    ZA = Zonal average (whithin the dfined box)
    AA = Area average (whithin the dfined box)
    ZE = zonal eddy (departure from zonal average)
    AE = area eddy (zonal departures from area averages)


Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
    
@author: daniloceano
"""

import xarray
from Math import (CalcZonalAverage, CalcAreaAverage)
import argparse
from thermodynamics import StaticStability, AdiabaticHEating
from metpy.constants import g
from metpy.units import units
import numpy as np
import pandas as pd

class BoxData:
    '''
    Object containing all meteorological data required for the LEC computation
    '''
    def __init__(self, data: xarray.Dataset, dfVars: pd.DataFrame,
                 western_limit: float, eastern_limit: float,
                 southern_limit: float, northern_limit: float,
                 args: argparse.Namespace, output_dir: str,
                 Q=None):
        self.LonIndexer = dfVars.loc['Longitude']['Variable']
        self.LatIndexer = dfVars.loc['Latitude']['Variable']
        self.TimeName = dfVars.loc['Time']['Variable']
        self.VerticalCoordIndexer = dfVars.loc['Vertical Level']['Variable']
        self.PressureData = data[self.VerticalCoordIndexer]*units(
             data[self.VerticalCoordIndexer].units).to('Pa')
        self.args = args
        self.output_dir = output_dir
        
        self.dx = float(data[self.LatIndexer][1]-data[self.LatIndexer][0])
        
        # Domain limits 
        self.western_limit = data[self.LonIndexer].sel(
            {self.LonIndexer:western_limit}, method='nearest')
        self.eastern_limit = data[self.LonIndexer].sel(
            {self.LonIndexer:eastern_limit}, method='nearest')
        self.southern_limit = data[self.LatIndexer].sel(
            {self.LatIndexer:southern_limit}, method='nearest')
        self.northern_limit = data[self.LatIndexer].sel(
            {self.LatIndexer:northern_limit}, method='nearest')
        
        # Set length for doing averages
        self.xlength = self.eastern_limit['rlons']-self.western_limit['rlons']
        self.ylength = np.sin(self.northern_limit['rlats']
                              ) - np.sin(self.southern_limit['rlats'])
        
        # Temperature data values, averages and eddy terms
        self.tair = (data[dfVars.loc['Air Temperature']['Variable']] \
              * units(dfVars.loc['Air Temperature']['Units']).to('K')).sel(
                  **{self.LatIndexer:slice(self.southern_limit, self.northern_limit),
                  self.LonIndexer: slice(self.western_limit, self.eastern_limit)})
        self.tair_ZA = CalcZonalAverage(self.tair,self.xlength)
        self.tair_AA = CalcAreaAverage(self.tair_ZA,self.ylength)
        self.tair_ZE = self.tair - self.tair_ZA
        self.tair_AE = self.tair_ZA - self.tair_AA 
        
        # Zonal wind component data values, averages and eddy terms
        self.u = (data[dfVars.loc['Eastward Wind Component']['Variable']] \
             * units(dfVars.loc['Eastward Wind Component']['Units']).to('m/s')
             ).sel(**{self.LatIndexer:slice(self.southern_limit, self.northern_limit),
                 self.LonIndexer: slice(self.western_limit, self.eastern_limit)})
        self.u_ZA = CalcZonalAverage(self.u,self.xlength)
        self.u_AA = CalcAreaAverage(self.u_ZA,self.ylength)
        self.u_ZE = self.u - self.u_ZA
        self.u_AE = self.u_ZA - self.u_AA
        
        # Meridional wind component data values, averages and eddy terms
        self.v = (data[dfVars.loc['Northward Wind Component']['Variable']] \
             * units(dfVars.loc['Northward Wind Component']['Units']).to('m/s')
             ).sel(**{self.LatIndexer:slice(self.southern_limit, self.northern_limit),
                 self.LonIndexer: slice(self.western_limit, self.eastern_limit)})
        self.v_ZA = CalcZonalAverage(self.v,self.xlength)
        self.v_AA = CalcAreaAverage(self.v_ZA,self.ylength)
        self.v_ZE = self.v - self.v_ZA
        self.v_AE = self.v_ZA - self.v_AA
        
        # Zonal and Meridional wind stress data values, averages and eddy terms
        if args.residuals:
            self.ust = self.tair*np.nan
            self.vst = self.tair*np.nan
        else:
            self.ust = (data[dfVars.loc['Zonal Wind Stress']['Variable']] \
                 * units(dfVars.loc['Zonal Wind Stress']['Units'])
                 ).sel(**{self.LatIndexer:slice(self.southern_limit, self.northern_limit),
                     self.LonIndexer: slice(self.western_limit, self.eastern_limit)})
            self.vst = (data[dfVars.loc['Meridional Wind Stress']['Variable']] \
                 * units(dfVars.loc['Meridional Wind Stress']['Units'])
                 ).sel(**{self.LatIndexer:slice(self.southern_limit, self.northern_limit),
                     self.LonIndexer: slice(self.western_limit, self.eastern_limit)})
                          
        self.ust_ZA = CalcZonalAverage(self.ust,self.xlength)
        self.ust_AA = CalcAreaAverage(self.ust_ZA,self.ylength)
        self.vst_ZA = CalcZonalAverage(self.vst,self.xlength)
        self.vst_AA = CalcAreaAverage(self.vst_ZA,self.ylength)
        self.ust_ZE = self.ust - self.ust_ZA
        self.ust_AE = self.ust_ZA - self.ust_AA
        self.vst_ZE = self.vst - self.vst_ZA
        self.vst_AE = self.vst_ZA - self.vst_AA
        
        # Omega velocity (vertical velocity in pressure levels) data values,
        # averages and eddy terms
        self.omega = (data[dfVars.loc['Omega Velocity']['Variable']] \
             * units(dfVars.loc['Omega Velocity']['Units']).to('Pa/s')
             ).sel(**{self.LatIndexer:slice(self.southern_limit, self.northern_limit),
                 self.LonIndexer: slice(self.western_limit, self.eastern_limit)})
        self.omega_ZA = CalcZonalAverage(self.omega,self.xlength)
        self.omega_AA = CalcAreaAverage(self.omega_ZA,self.ylength)
        self.omega_ZE = self.omega - self.omega_ZA
        self.omega_AE = self.omega_ZA - self.omega_AA
        
        # Geopotential (g*hgt) height data values data values, a
        # verages and eddy terms
        if args.geopotential:
            self.geopt = (data[dfVars.loc['Geopotential']['Variable']] \
                 * units(dfVars.loc['Geopotential']['Units']).to('m**2/s**2')
                 ).sel(**{self.LatIndexer:slice(self.southern_limit, self.northern_limit),
                     self.LonIndexer: slice(self.western_limit, self.eastern_limit)})
        else:
            self.geopt = (data[dfVars.loc['Geopotential Height']['Variable']]*g\
             * units(dfVars.loc['Geopotential Height']['Units'])
             ).sel(**{self.LatIndexer:slice(self.southern_limit,
                                            self.northern_limit),
             self.LonIndexer: slice(self.western_limit, self.eastern_limit)}
                      ).metpy.convert_units('m**2/s**2')        
        self.geopt_ZA = CalcZonalAverage(self.geopt,self.xlength)
        self.geopt_AA = CalcAreaAverage(self.geopt_ZA,self.ylength)
        self.geopt_ZE = self.geopt - self.geopt_ZA
        self.geopt_AE = self.geopt_ZA - self.geopt_AA
        
        # Adiaatic heating
        if args.fixed:
            self.Q = AdiabaticHEating(self.tair,self.tair[self.VerticalCoordIndexer],
                self.omega, self.u,self.v,self.VerticalCoordIndexer,
                self.LatIndexer,self.LonIndexer,self.TimeName).sel(
                    **{self.LatIndexer:slice(self.southern_limit, self.northern_limit),
                self.LonIndexer: slice(self.western_limit, self.eastern_limit)})
        elif args.track or args.choose:
            self.Q = Q
        else:
            print("could not compute Q. Check flags!")
        
        self.Q_ZA = CalcZonalAverage(self.Q,self.xlength)
        self.Q_AA = CalcAreaAverage(self.Q_ZA,self.ylength)
        self.Q_ZE = self.Q - self.Q_ZA
        self.Q_AE = self.Q_ZA - self.Q_AA
        
        # Static stability parameter
        self.sigma_AA = StaticStability(self.tair, self.PressureData, self.VerticalCoordIndexer,
                        self.xlength, self.ylength)
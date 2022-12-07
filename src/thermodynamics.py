#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 28 11:33:22 2022

Script for thermodynamics calculations necessary for computate
Lorenz Energy Cycle such static stability parameter, and all temrs of the
thermodynamic equation.

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com
"""

import numpy as np
from metpy.units import units
from metpy.constants import Rd
from metpy.constants import Cp_d
from metpy.constants import g
from metpy.constants import Re
from metpy.calc import potential_temperature


# def StaticStability(TemperatureData,PressureData,VerticalCoordIndexer,
#                     LatIndexer,LonIndexer,BoxNorth,BoxSouth,BoxWest, BoxEast):
def StaticStability(TemperatureData,PressureData,VerticalCoordIndexer,
                    xlength,ylength):    
    
    """
    Compute the static stability parameter sigma for all vertical levels
    and for the desired domain
    
    Source:
        Michaelides, S. C. (1987). 
        Limited Area Energetics of Genoa Cyclogenesis,
        Monthly Weather Review, 115(1), 13-26. Retrieved Jan 24, 2022, from:
        https://journals.ametsoc.org/view/journals/mwre/115/1/1520-0493_1987_115_0013_laeogc_2_0_co_2.xml
    
    Parameters
    ----------
    temp: xarray.Dataset
        temperature data in Kelvin
    min_lon, max_lon, minlat, min_lon: float
        minimum and maximum longitude/latitude to be calculated the zonal average

    Returns
    -------
    sigma: xarray.Dataset
        Dataset containing sigma values for all pressure levels and for the
        box specyfied by min_lon, max_lon, min_lat and max_lat    
    
    """
    FirstTerm = g*TemperatureData/Cp_d
    SecondTerm = (PressureData*g/Rd)
    ThirdTerm = TemperatureData.differentiate(VerticalCoordIndexer
                            )/TemperatureData[VerticalCoordIndexer].metpy.units
    function = (FirstTerm-(SecondTerm*ThirdTerm))
    sigma_ZA = function.integrate("rlons")/xlength
    sigma_AA = (sigma_ZA*sigma_ZA["coslats"]).integrate("rlats")/ylength
    return sigma_AA

def TairTendency(TemperatureData,TimeName,LatIndexer,LonIndexer,
                 VerticalCoordIndexer):
    delTdelt = TemperatureData.differentiate(
            TimeName,datetime_unit='s') / units('seconds')

    return delTdelt

def HorizontalTemperatureAdvection(TemperatureData,LatIndexer,LonIndexer,
                                   UWindComponentData,VWindComponentData,
                                   VerticalCoordIndexer):
    lons,lats = TemperatureData[LonIndexer],TemperatureData[LatIndexer]
    cos_lats = np.cos(np.deg2rad(lats))
    ## Horizontal advection of temperature ##
    # Differentiate temperature in respect to longitude and latitude
    dTdlambda = TemperatureData.copy(deep=True).differentiate(LonIndexer)
    dTdphi = TemperatureData.copy(deep=True).differentiate(LatIndexer)
    # Get the values for dx and dy in meters
    dx = np.deg2rad(lons.differentiate(LonIndexer))*cos_lats*Re
    dy = np.deg2rad(lats.differentiate(LatIndexer))*Re
    AdvHT = -1* ((UWindComponentData*dTdlambda/dx)+(
                                VWindComponentData*dTdphi/dy))

    
    return AdvHT

def AdiabaticHEating(TemperatureData,PressureData, OmegaData,
                     UWindComponentData,VWindComponentData,
                     VerticalCoordIndexer,LatIndexer,LonIndexer,TimeName):
        """
        Compute the diabatic heating as a residual form the thermodynamic 
        equation for all vertical levels and for the desired domain
        """
        # Temperature tendency as dT/dt
        delTdelt = TairTendency(TemperatureData,TimeName,LatIndexer,LonIndexer,
                 VerticalCoordIndexer)
        
        # Horizonal temperature advection
        AdvHT = HorizontalTemperatureAdvection(TemperatureData,LatIndexer,
                 LonIndexer,UWindComponentData,VWindComponentData,
                                   VerticalCoordIndexer)
        
        # Static stability parameter (here we need a slight modified version
        # from calc.py so the units can match)
        theta = potential_temperature(PressureData,TemperatureData)
        sigma = -(TemperatureData/theta) *theta.differentiate(
            VerticalCoordIndexer)/units.hPa
        
        ## Diabatic reating as a residual ##
        Q = delTdelt-AdvHT-(sigma*OmegaData)
        J = Q*Cp_d
        return J
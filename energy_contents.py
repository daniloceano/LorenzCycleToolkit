#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computate partioned energy contents:
    AZ = zonal available potential energy
    AE = eddy available potential energy
    KZ = zonal kinetic energy
    KE = eddy kinetic energy


Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com
"""

from calc import (VerticalTrazpezoidalIntegration,
                       CalcZonalAverage, CalcAreaAverage,
                       StaticStability)
from metpy.constants import g


def Calc_Az(TemperatureData,PressureData,LonIndexer,LatIndexer,VerticalCoordIndexer):
    '''
    Parameters
    ----------
    TemperatureData: xarray.Dataset
        unit aware array for the temperature  
    PressureData: xarray.Dataset
        unit aware array for the pressure levels
    LonIndexer: string
        the indexer used for the longitude variable in the xarray
    LatIndexer: string 
        the indexer used for the latitude variable in the xarray         
    VerticalCoordIndexer: string
        the indexer used for the xarray coordante which the integration will
        be performed
    
    Returns
    -------
    Az: xarray.DataArray
        Unit-aware values corresponding of the Zonal Available Potential Energy
    '''        

    ## Temperature averages
    tair_ZA = CalcZonalAverage(TemperatureData,LonIndexer) # Zonal Average
    tair_AA = CalcAreaAverage(TemperatureData,LonIndexer,LatIndexer) # Area Average
    
    ## Temperature area eddy
    tair_AE = tair_ZA-tair_AA # Area Eddy
    
    sigma_AA = StaticStability(TemperatureData,PressureData,VerticalCoordIndexer,
                            LatIndexer,LonIndexer)
    
    area_ave = CalcAreaAverage(tair_AE**2,LatIndexer)
    function = area_ave/(2*sigma_AA)
    Az = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)
    
    try: 
        Az = Az.metpy.convert_units('J/ m **2')
    except ValueError:
        print('Unit error in Az')
        raise
    
    return Az
    

def Calc_Ae(TemperatureData,PressureData,LonIndexer,LatIndexer,VerticalCoordIndexer):
    '''
    Parameters
    ----------
    TemperatureData: xarray.Dataset
        unit aware array for the temperature  
    PressureData: xarray.Dataset
        unit aware array for the pressure levels
    LonIndexer: string
        the indexer used for the longitude variable in the xarray
    LatIndexer: string 
        the indexer used for the latitude variable in the xarray         
    VerticalCoordIndexer: string
        the indexer used for the xarray coordante which the integration will
        be performed
    
    Returns
    -------
    Az: xarray.DataArray
        Unit-aware values corresponding of the Eddy Available Potential Energy
    '''   

    ## Temperature zonal averages
    tair_ZA = CalcZonalAverage(TemperatureData,LonIndexer) # Zonal Average
    
    ## Temperature zonal eddy 
    tair_ZE = TemperatureData-tair_ZA # Zonal Eddy
    
    sigma_AA = StaticStability(TemperatureData,PressureData,VerticalCoordIndexer,
                            LatIndexer,LonIndexer)
    
    area_ave = CalcAreaAverage(tair_ZE**2,LatIndexer,LonIndexer)
    function = area_ave/(2*sigma_AA)
    Az = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)
    
    try: 
        Ae = Az.metpy.convert_units('J/ m **2')
    except ValueError:
        print('Unit error in Az')
        raise
    
    return Ae

def Calc_Kz(UWindComponentData,VWindComponentData,PressureData,LonIndexer,LatIndexer,VerticalCoordIndexer):
    '''
    Parameters
    ----------
    UWindComponentData: xarray.Dataset
        unit aware array for the meridional component of the wind  
    VWindComponentData: xarray.Dataset
        unit aware array for the zonal component of the wind 
    PressureData: xarray.Dataset
        unit aware array for the pressure levels
    LonIndexer: string
        the indexer used for the longitude variable in the xarray
    LatIndexer: string 
        the indexer used for the latitude variable in the xarray         
    VerticalCoordIndexer: string
        the indexer used for the xarray coordante which the integration will
        be performed
    
    Returns
    -------
    Kz: xarray.DataArray
        Unit-aware values corresponding of the Zonal Kinectic Energy
    '''
    
    # Wind Component Averages
    u_ZA = CalcZonalAverage(UWindComponentData,LonIndexer)
    v_ZA = CalcZonalAverage(VWindComponentData,LonIndexer)

    Wind_ZA = (u_ZA**2)+(v_ZA**2)
    Wind_AA = CalcAreaAverage(Wind_ZA,LatIndexer)
    
    function = Wind_AA/(2*g)
    
    Kz = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)
    
    
    ### TEST NUMPY TRAPZ INTEGRATION!!! ####
    
    try: 
        Kz = Kz.metpy.convert_units('J/ m **2')
    except ValueError:
        print('Unit error in Kz')
        raise
        
    # Kz = -np.trapz(function,PressureData)
    
    return Kz

def Calc_Ke(UWindComponentData,VWindComponentData,PressureData,LonIndexer,LatIndexer,VerticalCoordIndexer): 
    '''
    Parameters
    ----------
    UWindComponentData: xarray.Dataset
        unit aware array for the meridional component of the wind  
    VWindComponentData: xarray.Dataset
        unit aware array for the zonal component of the wind 
    PressureData: xarray.Dataset
        unit aware array for the pressure levels
    LonIndexer: string
        the indexer used for the longitude variable in the xarray
    LatIndexer: string 
        the indexer used for the latitude variable in the xarray         
    VerticalCoordIndexer: string
        the indexer used for the xarray coordante which the integration will
        be performed
    
    Returns
    -------
    Ke: xarray.DataArray
        Unit-aware values corresponding of the Eddy Kinectic Energy
    '''
    
    # Wind Component Averages
    u_ZA = CalcZonalAverage(UWindComponentData,LonIndexer)
    v_ZA = CalcZonalAverage(VWindComponentData,LonIndexer)
    
    # Wind Eddy Terms
    u_ZE = UWindComponentData-u_ZA
    v_ZE = VWindComponentData-v_ZA

    Wind_ZE = (u_ZE**2)+(v_ZE**2)
    Wind_AA = CalcAreaAverage(Wind_ZE,LatIndexer,LonIndexer)
    
    function = Wind_AA/(2*g)
    
    Ke = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)
    
    return Ke

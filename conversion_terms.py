#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for computating the conversion terms as part of the Lorenz Energy Cycle

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
from calc import (CalcZonalAverage, CalcAreaAverage,
                 VerticalTrazpezoidalIntegration, Differentiate,
                 StaticStability)


def Calc_Ce(PressureData,OmegaData,TemperatureData,LonIndexer,LatIndexer,VerticalCoordIndexer):
    '''
    Parameters
    ----------
    PressureData: xarray.Dataset
        unit aware array for the pressure levels
    OmegaData: xarray.Dataset
        unit aware array for the omega velocity data
    TemperatureData: xarray.Dataset
        unit aware array for the temperature  
    LonIndexer: string
        the indexer used for the longitude variable in the xarray
    LatIndexer: string 
        the indexer used for the latitude variable in the xarray  
    
    Returns
    -------
    Ce: xarray.DataArray
        Unit-aware values corresponding to the conversion between eddy energy terms
    '''        
    FirstTerm = Rd/PressureData
    # Zonal average terms
    omega_ZA = CalcZonalAverage(OmegaData,LonIndexer)
    tair_ZA = CalcZonalAverage(TemperatureData,LonIndexer)
    # Zonal eddy terms
    omega_ZE = OmegaData-omega_ZA
    tair_ZE = TemperatureData-tair_ZA
    TemperatureVerticalTransport = omega_ZE*tair_ZE
    # Area averages
    SecondTerm = CalcAreaAverage(TemperatureVerticalTransport,LonIndexer,LatIndexer)
    # Integrate in pressure
    function = FirstTerm*SecondTerm
    Ce = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)/(-g)
    # Check if units are in accordance with expected and convert
    try: 
        Ce = Ce.metpy.convert_units('W/ m **2')
    except ValueError:
        print('Unit error in Ce')
        raise
    return Ce

def Calc_Cz(PressureData,OmegaData,TemperatureData,LonIndexer,LatIndexer,VerticalCoordIndexer):
    '''
    Parameters
    ----------
    PressureData: xarray.Dataset
        unit aware array for the pressure levels
    OmegaData: xarray.Dataset
        unit aware array for the omega velocity data
    TemperatureData: xarray.Dataset
        unit aware array for the temperature  
    LonIndexer: string
        the indexer used for the longitude variable in the xarray
    LatIndexer: string 
        the indexer used for the latitude variable in the xarray  
    
    Returns
    -------
    Cz: xarray.DataArray
        Unit-aware values corresponding to the conversion between zonal energy terms
    '''   
    FirstTerm = Rd/PressureData
    # Zonal average terms
    omega_ZA = CalcZonalAverage(OmegaData,LonIndexer)
    tair_ZA = CalcZonalAverage(TemperatureData,LonIndexer)
    # Area average terms
    omega_AA = CalcAreaAverage(omega_ZA,LatIndexer)
    tair_AA = CalcAreaAverage(tair_ZA,LatIndexer)
    # Area eddy terms
    omega_AE = omega_ZA-omega_AA
    tair_AE = tair_ZA-tair_AA
    TemperatureVerticalTransport = omega_AE*tair_AE
    # Area averages
    SecondTerm = CalcAreaAverage(TemperatureVerticalTransport,LatIndexer)
    # Integrate in pressure
    function = FirstTerm*SecondTerm
    Cz = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)/(-g)
    # Check if units are in accordance with expected and convert
    try: 
        Cz = Cz.metpy.convert_units('W/ m **2')
    except ValueError:
        print('Unit error in Cz')
        raise
    return Cz
    
def Calc_Ca(VWindComponentData,PressureData,OmegaData,TemperatureData,
            LonIndexer,LatIndexer,VerticalCoordIndexer):
    '''
    
    This expression is a slightly modificated from the source material:
        In the first term of the integral Brennan & Vincent write the denominator
        as (sigma * Re), however, other sources use (2 * sigma * Re), such the
        Muench's work (1965) and the  Michaelides' (1986). Therefore, here it
        is adopted the latter expression.
    
    Parameters
    ----------
    VWindComponentData: xarray.Dataset
        unit aware array for the meridional component of the wind 
    PressureData: xarray.Dataset
        unit aware array for the pressure levels
    OmegaData: xarray.Dataset
        unit aware array for the omega velocity data
    TemperatureData: xarray.Dataset
        unit aware array for the temperature  
    LonIndexer: string
        the indexer used for the longitude variable in the xarray
    LatIndexer: string 
        the indexer used for the latitude variable in the xarray 
    VerticalCoordIndexer: string
        the indexer used for the xarray coordante which the integration will
        be performed
    
    Returns
    -------
    Ca: xarray.DataArray
        Unit-aware values corresponding to the conversion between 
        available potential energy terms
    ''' 
    rlats = np.deg2rad(VWindComponentData[LatIndexer])
    # Zonal averages
    v_ZA = CalcZonalAverage(VWindComponentData,LonIndexer)
    tair_ZA = CalcZonalAverage(TemperatureData,LonIndexer)
    omega_ZA = CalcZonalAverage(OmegaData,LonIndexer)
    # Zonal eddies
    v_ZE = VWindComponentData - v_ZA
    tair_ZE = TemperatureData - tair_ZA
    omega_ZE = OmegaData - omega_ZA
    # Temperature area average
    tair_AA = CalcAreaAverage(tair_ZA,LatIndexer)
    # Temperature area eddy
    tair_AE = tair_ZA-tair_AA
    # Sigma
    sigma_AA = StaticStability(TemperatureData,PressureData,VerticalCoordIndexer,
                            LatIndexer,LonIndexer)
    # First term of the integral
    tmp1 = v_ZE*tair_ZE/(Re*2*sigma_AA)
    tmp2 = Differentiate(tair_AE,rlats,LatIndexer)
    FirstTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
    # Second term of the integral
    tmp1 =  omega_ZE*tair_ZE/(sigma_AA)
    tmp2 = Differentiate(tair_AE,PressureData,VerticalCoordIndexer)
    SecondTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
    # Integrate in pressure
    function = FirstTerm+SecondTerm
    Ca = -VerticalTrazpezoidalIntegration(function,PressureData,
                                          VerticalCoordIndexer)*9.8
    # Check if units are in accordance with expected and convert
    try: 
        Ca = Ca.metpy.convert_units('W/ m **2')
    except ValueError:
        print('Unit error in Ck')
        raise
    return Ca

def Calc_Ck(UWindComponentData,VWindComponentData,PressureData,OmegaData,
            LonIndexer,LatIndexer,VerticalCoordIndexer):
    '''
     Parameters
    ----------
    VWindComponentData: xarray.Dataset
        unit aware array for the meridional component of the wind 
    PressureData: xarray.Dataset
        unit aware array for the pressure levels
    OmegaData: xarray.Dataset
        unit aware array for the omega velocity data
    TemperatureData: xarray.Dataset
        unit aware array for the temperature  
    LonIndexer: string
        the indexer used for the longitude variable in the xarray
    LatIndexer: string 
        the indexer used for the latitude variable in the xarray 
    VerticalCoordIndexer: string
        the indexer used for the xarray coordante which the integration will
        be performed
    
    Returns
    -------
    Cz: xarray.DataArray
        Unit-aware values corresponding to the conversion between kinetic energy terms
    '''  
    # Convert latitude do radians
    rlats = np.deg2rad(UWindComponentData[LatIndexer])
    # Zonal averages 
    v_ZA = CalcZonalAverage(VWindComponentData,LonIndexer)
    u_ZA = CalcZonalAverage(UWindComponentData,LonIndexer)
    omega_ZA = CalcZonalAverage(OmegaData,LonIndexer)
    # Zonal eddies
    v_ZE = VWindComponentData - v_ZA
    u_ZE = UWindComponentData - u_ZA
    omega_ZE = OmegaData - omega_ZA
    # First term
    tmp1 = np.cos(rlats)*u_ZE*v_ZE/Re
    tmp2 = Differentiate(u_ZA/np.cos(rlats),rlats,LatIndexer)
    FirstTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
    # Second term
    tmp1 = (v_ZE**2)/Re
    tmp2 = Differentiate(v_ZA,rlats,LatIndexer)
    SecondTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
    # Third term
    tmp1 = (np.tan(rlats)*(u_ZE**2)*v_ZA)/Re
    ThirdTerm = CalcAreaAverage(tmp1,LatIndexer,LonIndexer)
    # Fourth term
    tmp1 = omega_ZE*u_ZE
    tmp2 = Differentiate(u_ZA,PressureData,VerticalCoordIndexer)
    FourthTerm = CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
    # Fifith term
    tmp1 = omega_ZE*v_ZE
    tmp2 = Differentiate(v_ZA,PressureData,VerticalCoordIndexer)
    FifithTerm =  CalcAreaAverage(tmp1*tmp2,LatIndexer,LonIndexer)
    # Integrate in pressure
    function = (FirstTerm+SecondTerm+ThirdTerm+FourthTerm+FifithTerm)
    Ck = VerticalTrazpezoidalIntegration(function,PressureData,VerticalCoordIndexer)/g
    # Check if units are in accordance with expected and convert
    try: 
        Ck = Ck.metpy.convert_units('W/ m **2')
    except ValueError:
        print('Unit error in Ck')
        raise
    return Ck
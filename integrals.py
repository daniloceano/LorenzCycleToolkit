#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculating zonal average and area average of some variable
Used for computate Lorenz Energy Cycle

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com

"""

# from derivatives import calc_delvar_delp
import numpy as np
# import math

def TrazpezoidalIntegration(VariableData,dimension):
    """
    Integrate data using the trapezoidal rule 
    
    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated
    dimension: string
        the indexer used for the xarray coordante which the integration will
        be performed
   
    Returns
    -------
    trapz: xarray.Dataset
        Arrays with the integrated data
    """
    DimensionArray = VariableData[dimension].values # array containing either lats or lons 
    trapz = 0
    delta = np.deg2rad(np.abs((DimensionArray[-1]-DimensionArray[0])/(len(DimensionArray)-1)))
    for i in DimensionArray:
        if i == DimensionArray[0] or i == DimensionArray[-1]:
            trapz += VariableData.sel({dimension:i})
        else:
            trapz += 2*VariableData.sel({dimension:i})
    return trapz*(delta/2)

def CalcZonalAverage(VariableData,LonName):
    """
    Computates zonal averages of a variable for all z levels and times 
    
    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated
    LonName: string
        the indexer used for the longitude variable in the xarray       
   
    Returns
    -------
    zonal_ave: xarray.Dataset
        Arrays of zonal avreages for all longitudes from the passed Dataset
    """
    # Get latitude and logitude data
    lons = VariableData[LonName]
    rlons = np.deg2rad(lons)
    # Integrate through longitude
    trapz = TrazpezoidalIntegration(VariableData,LonName)
    # Take the zonal average
    zonal_ave = trapz/(rlons[-1]-rlons[0]).values
    return zonal_ave

def CalcAreaAverage(VariableData,LatName,LonName=None):
    """
    Computates the Area Average of a function.
    
    The default is to computate the zonal average and then a meridional average.
    If the input data is already some sort of zonal quantity (average or not),
    simply set LonName to None
    
    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated
    LatName: string 
        the indexer used for the latitude variable in the xarray 
    LonName: string (optional)
        the indexer used for the longitude variable in the xarray 
   
    Returns
    -------
    zonal_ave: xarray.Dataset
        Arrays of area avreages for all latitudes and longitudes from
        the passed Dataset
    """
    # Get latitude and logitude data
    lats = VariableData[LatName]
    rlats = np.deg2rad(lats)
    if LonName:
        # If LonName is provided, get zonal ave
        zonal_ave = CalcZonalAverage(VariableData,LonName)
    else:
        zonal_ave = VariableData
    # Take the area avearge
    trapz = TrazpezoidalIntegration(zonal_ave*np.cos(rlats),LatName)
    area_ave = trapz/(np.sin(rlats[0])-np.sin(rlats[-1])).values
    return area_ave

# def calc_static_stability(temp,min_lon,max_lon,min_lat,max_lat):
#     """
#     Computates the static stability parameter sigma for all vertical levels
#     and for the desired domain
    
#     Parameters
#     ----------
#     temp: xarray.Dataset
#         temperature data in Kelvin
#     min_lon, max_lon, minlat, min_lon: float
#         minimum and maximum longitude/latitude to be calculated the zonal average

#     Returns
#     -------
#     sigma: xarray.Dataset
#         Dataset containing sigma values for all pressure levels and for the
#         box specyfied by min_lon, max_lon, min_lat and max_lat    
    
#     """
#     R = 287     # ideal gas cte.
#     g = 9.81    # gravity
#     cp = 1004.0 # specific heat of air at cte. pressure
    
#     # breaks the equation into two terms
#     # 1) g * T/cp
#     tmp1 = temp*0
#     for lev in temp.level.values:
#         tmp1.loc[dict(level=lev)] = g * temp.sel(level=lev)/cp

#     # 2)  (g*p/R) * (delT/delP)
#     tmp2 = calc_delvar_delp(temp)
#     tmp3 = temp*0
#     for lev in temp.level.values:
#         tmp3.loc[dict(level=lev)] = (g * lev/R) * tmp2.sel(level=lev)
    
#     # first term - second term
#     var = tmp1 - tmp3
    
#     sigma = calc_area_ave(var,min_lon,max_lon,min_lat,max_lat)
    
#     return sigma
    
#def calc_static_stability(temp,min_lon,max_lon,min_lat,max_lat):
#    """
#    Computates the static stability parameter sigma for all vertical levels
#    and for the desired domain
#    
#    Parameters
#    ----------
#    temp: xarray.Dataset
#        temperature data in Kelvin
#    min_lon, max_lon, minlat, min_lon: float
#        minimum and maximum longitude/latitude to be calculated the zonal average
#
#    Returns
#    -------
#    sigma: xarray.Dataset
#        Dataset containing sigma values for all pressure levels and for the
#        box specyfied by min_lon, max_lon, min_lat and max_lat    
#    
#    """
#    R = 287     # ideal gas cte.
#    g = 9.81    # gravity
#    cp = 1004.0 # specific heat of air at cte. pressure
#    
#    za_temp = calc_zonal_ave(temp,min_lon,max_lon,min_lat,max_lat)
#    
#    ave_temp = calc_area_ave(za_temp,min_lon,max_lon,min_lat,max_lat)
#    
#    # breaks the equation into two terms
#    # 1) R * T/ cp * p
#    tmp1 = ave_temp*0
#    for p in temp.level.values:
#        tmp1.loc[dict(level=p)] = (R * ave_temp.sel(level=p))/ (cp * p)
#
#    # 2)  delT / delP
#    tmp2 = calc_delvar_delp(ave_temp)
#    
#    #g * ( first term - second term )
#    sigma = g*(tmp1 - tmp2)
#    
#    
#    return sigma
    

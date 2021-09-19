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

from derivatives import calc_delvar_delp
import numpy as np
import math

def calc_zonal_ave(var,min_lon,max_lon,min_lat,max_lat):
    """
    Makes zonal averages of a variable for all z levels and times 
    
    Parameters
    ----------
    var: xarray.Dataset
        the variable to be integrated
    lat: float
        the latitude to be calculated the zonal average
    lons: array
        the corresponding array of longitudes for the variable    
    min_lon, max_lon: float
        minimum and maximum longitude to be calculated the zonal average        
   
    Returns
    -------
    zonal_ave: xarray.Dataset
        Arrays of zonal avreages for all latitudes in the passed varieble
    """
        
    # trimm data for desired lons
    var = var.sel(lon=slice(min_lon,max_lon),lat=slice(max_lat,min_lat))
    
    # # integrate using trapezoidal rule
    # lons, trapz =var.lon.values, 0
    # dx = (lons[-1]-lons[0])/(len(lons))
    # for lon in lons:
    #     if lon == lons[0] or lon == lons[-1]:
    #         trapz += var.sel(lon=lon)
    #     else:
    #         trapz += 2*var.sel(lon=lon)
    # trapz = trapz*(dx/2)
    # # take the zonal average
    # zonal_ave = trapz/(lons[0]-lons[1])
    
    # integrate using numpy trapezoidal function
    lons = var.lon.values
    trapz = np.trapz(var,lons)
    zonal_ave = trapz/(lons[0]-lons[1])
    
    # # makes zonal average 
    # zonal_ave = var*0
    # for lon in var.lon.values:
    #     zonal_ave.loc[dict(lon=lon)] =  var.mean('lon')
        
    return zonal_ave

def calc_area_ave(var,min_lon,max_lon,min_lat,max_lat):
    """
    Makes an area average of a variable for all z levels and times
    
    Parameters
    ----------
    var: xarray.Dataset
        the variable to be integrated
    min_lon, max_lon, minlat, min_lon: float
        minimum and maximum longitude/latitude to be calculated the zonal average

    Returns
    -------
    area_ave: xarray.Dataset
        Dataset of the passed varieble avreages trough the area delimited by
        min_lon, max_lon, min_lat and max_lat
        Note that this function returns the CZ value in all grid points    
    """
    
    # # integrate using trapezoidal rule
    # trapz = 0
    # lons = var.lon
    # dx = (lons[-1]-lons[0])/(len(lons))
    # for lon in lons:
    #     if lon == lons[0] or lon == lons[-1]:
    #         trapz += (var.sel(lon=lon)*np.sin(var.sel(lon=lon)))
    #     else:
    #         trapz += 2*(var.sel(lon=lon)*np.sin(var.sel(lon=lon)))
    # trapz = trapz*(dx/2)
    
    # trimm data for desired lats and lons
    var = var.sel(lon=slice(min_lon,max_lon),lat=slice(max_lat,min_lat))
    # calc zonal average
    za = calc_zonal_ave(var,min_lon,max_lon,min_lat,max_lat)
    lats = var.lat.values
    trapz = np.trapz(za*np.cos(za))
    area_ave = trapz/(lats[0]-lats[1])
    
    # # trimm data for desired lats and lons
    # x = var.sel(lon=slice(min_lon,max_lon),lat=slice(max_lat,min_lat))
    # # calc zonal average
    # za = calc_zonal_ave(var,min_lon,max_lon,min_lat,max_lat)
    # # makes area average
    # area_ave = x*0
    # for lat in x.lat.values:
    #     area_ave.loc[dict(lat=lat)] = \
    #             za.mean('lat')
    
    return area_ave


def interga_P(var):
    """
    Intergate the provided variable trough its pressure levels
    
    Parameters
    ----------
    var: xarray.Dataset
        the variable to be integrated
    """
    
    # get pressure levels
    levs = var.level.values
    # sum trough pressure levels
    # and 
    sum_levs = 0
    for lev in range((len(levs)-1),0,-1):
        # multiply by 100 to get pressur elevels in pascals
        dp = (levs[lev-1] - levs[lev])
        #
        sum_levs = sum_levs + (
                    ((var.sel(level=levs[lev-1]) + \
                    var.sel(level=levs[lev]))/2) * dp)
        
    return  sum_levs

def interga_phi(var,min_lat,max_lat):
    """
    Intergate the provided variable trough latitude
    
    Parameters
    ----------
    var: xarray.Dataset
        the variable to be integrated
    """
    
    # get pressure levels
    lats = var.lat.sel(lat=slice(max_lat,min_lat)).values
    # sum trough pressure levels
    # and 
    sum_lats = 0
    for lat in range((len(lats)-1),0,-1):
        # multiply by 100 to get pressur elevels in pascals
        dy = (lats[lat-1] - lats[lat])
        #
        sum_lats = sum_lats + (
                    ((var.sel(lat=lats[lat-1]) + \
                    var.sel(lat=lats[lat]))/2) * dy)
        
    return  sum_lats



def calc_static_stability(temp,min_lon,max_lon,min_lat,max_lat):
    """
    Computates the static stability parameter sigma for all vertical levels
    and for the desired domain
    
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
    R = 287     # ideal gas cte.
    g = 9.81    # gravity
    cp = 1004.0 # specific heat of air at cte. pressure
    
    # breaks the equation into two terms
    # 1) g * T/cp
    tmp1 = temp*0
    for lev in temp.level.values:
        tmp1.loc[dict(level=lev)] = g * temp.sel(level=lev)/cp

    # 2)  (g*p/R) * (delT/delP)
    tmp2 = calc_delvar_delp(temp)
    tmp3 = temp*0
    for lev in temp.level.values:
        tmp3.loc[dict(level=lev)] = (g * lev/R) * tmp2.sel(level=lev)
    
    # first term - second term
    var = tmp1 - tmp3
    
    sigma = calc_area_ave(var,min_lon,max_lon,min_lat,max_lat)
    
    return sigma
    
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
    

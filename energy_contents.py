#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computate partioned energy contents:
    AZ = zonal available potential energy
    AE = eddy available potential energy
    KZ = zonal kinetic energy
    KE = eddy kinetic energy

Parameters
----------
temp, uwnd, vwnd: xarray.Dataset
    array containing temperature/wind data
min_lon, max_lon, minlat, min_lon: float
    lons and lats for the box bounding the computation domain  

Returns
-------
AZ/AE/KZ/KE: xarray.Dataset
    Values corresponding of the energy content term
    Note that this function returns the term value for all grid points

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com
"""

g = 9.81 # gravity acceleration 

from integrals import (calc_area_ave, calc_zonal_ave,
                       calc_static_stability, interga_P)

def calc_AZ(temp,min_lon,max_lon,min_lat,max_lat):
    
    za_temp = calc_zonal_ave(temp,min_lon,max_lon,min_lat,max_lat)
    aa_temp = calc_area_ave(temp,min_lon,max_lon,min_lat,max_lat)
    da_temp  = za_temp - aa_temp
    
    sigma = calc_static_stability(temp,min_lon,max_lon,min_lat,max_lat)

    tmp = calc_area_ave((da_temp**2)/(2*sigma),min_lon,max_lon,min_lat,max_lat)
    
    AZ = interga_P(tmp)
    
    return AZ.mean('lon').mean('lat')
    
def calc_AE(temp,min_lon,max_lon,min_lat,max_lat):
    
    za_temp = calc_zonal_ave(temp,min_lon,max_lon,min_lat,max_lat)
    dz_temp  = temp - za_temp
    
    sigma = calc_static_stability(temp,min_lon,max_lon,min_lat,max_lat)

    tmp = calc_area_ave((dz_temp**2)/(2*sigma),min_lon,max_lon,min_lat,max_lat)
    
    AE = interga_P(tmp)
    
    return AE.mean('lon').mean('lat')

def calc_KZ(uwnd,vwnd,min_lon,max_lon,min_lat,max_lat):
    
    za_uwnd = calc_zonal_ave(uwnd,min_lon,max_lon,min_lat,max_lat)
    za_vwnd = calc_zonal_ave(vwnd,min_lon,max_lon,min_lat,max_lat)
        
    tmp = calc_area_ave((za_uwnd**2) + (za_vwnd**2),min_lon,max_lon,min_lat,max_lat)
    
    KZ = interga_P(tmp/(2*g))
    
    return KZ.mean('lon').mean('lat')

def calc_KE(uwnd,vwnd,min_lon,max_lon,min_lat,max_lat):
    
    za_uwnd = calc_zonal_ave(uwnd,min_lon,max_lon,min_lat,max_lat)
    za_vwnd = calc_zonal_ave(vwnd,min_lon,max_lon,min_lat,max_lat)
    
    dz_uwnd = uwnd - za_uwnd
    dz_vwnd = vwnd - za_vwnd
        
    tmp = calc_area_ave((dz_uwnd**2) + (dz_vwnd**2),min_lon,max_lon,min_lat,max_lat)
    
    KE = interga_P(tmp/(2*g))
    
    return KE.mean('lon').mean('lat')

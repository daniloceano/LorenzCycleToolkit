#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computate energy fluxes across boundaries:
    BAz = zonal available potential energy
    BAe = eddy available potential energy
    BKz = zonal kinetic energy
    BKe = eddy kinetic energy

Parameters
----------
temp, uwnd, vwnd, omg: xarray.Dataset
    array containing temperature/wind/omega data
min_lon, max_lon, minlat, min_lon: float
    lons and lats for the box bounding the computation domain  

Returns
-------
BAz/BAe/BKz/BKe: xarray.Dataset
    Values corresponding of the energy content term
    Note that this function returns the term value for all grid points

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com
"""

import numpy as np
from integrals import (calc_area_ave, calc_zonal_ave,
                       calc_static_stability, interga_P, interga_phi)

# constants
ra = 6378000    # Earth's radius
g = 9.81        # gravity acceleration 


def calc_BAz(temp,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat):

    print('')    
    print('------------------------------------------')
    print('Computating BAz....')
    
    # zonal_aves
    za_vwnd = calc_zonal_ave(vwnd,min_lon,max_lon,min_lat,max_lat)
    za_omg = calc_zonal_ave(omg,min_lon,max_lon,min_lat,max_lat)
    za_temp = calc_zonal_ave(temp,min_lon,max_lon,min_lat,max_lat)

    # area_ave
    aa_temp = calc_area_ave(temp,min_lon,max_lon,min_lat,max_lat)

    # dev from zonal_ave
    dz_vwnd = vwnd - za_vwnd
    dz_omg = omg - za_omg
    dz_temp = temp - za_temp
    
    # dev from area_ave
    da_temp = za_temp - aa_temp
    
    # sigma
    sigma = calc_static_stability(temp,min_lon,max_lon,min_lat,max_lat)
        
    # 1) First integral
    
    tmp = (2*da_temp*dz_temp*uwnd) + ((da_temp**2)*uwnd)
    
    lamb_w = tmp.sel(lon=min_lon)/2*sigma.sel(lon=min_lon)
    lamb_e = tmp.sel(lon=max_lon)/2*sigma.sel(lon=max_lon)
    
    tmp = (lamb_e - lamb_w)
    
    tmp = interga_phi(tmp,min_lat,max_lat)
    
    tmp = interga_P(tmp)
    
    tmp2 = (ra*(max_lon-min_lon)*\
                      (np.sin(max_lat)-np.sin(min_lat)))
    
    first_int = tmp/tmp2 
    
    # 2) Second integral
    
    tmp = calc_zonal_ave(dz_vwnd*dz_temp,min_lon,max_lon,min_lat,max_lat)

    tmp = (2*tmp*da_temp)+((da_temp**2)*za_vwnd)
    
    phi_n = ((tmp.sel(lat=max_lat)*np.cos(max_lat))/(2*sigma)).sel(lat=max_lat)
    phi_s = ((tmp.sel(lat=min_lat)*np.cos(min_lat))/(2*sigma)).sel(lat=min_lat)
    
    tmp = phi_n - phi_s
    
    tmp = interga_P(tmp)

    tmp2 = (ra*(np.sin(max_lat)-np.sin(min_lat))) 

    second_int = (tmp/tmp2).mean('lon')
    
    # 3) Third term
    
    tmp = calc_zonal_ave(dz_omg*dz_temp,min_lon,max_lon,min_lat,max_lat)
    
    tmp = (2*tmp*da_temp)+((da_temp**2)*za_omg)
    
    tmp = calc_area_ave(tmp,min_lon,max_lon,min_lat,max_lat)/(2*sigma)
    
    levs = tmp.level.values
    tmp_b = tmp.sel(level=np.amax(levs))
    tmp_t = tmp.sel(level=np.amin(levs))
    
    tmp = tmp_b-tmp_t
    
    tmp = tmp.mean('lon').mean('lat')
    
    BAz = -first_int-second_int-tmp
    
    print('Done!')
    
    return BAz

def calc_BAe(temp,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat):
    
    print('')    
    print('------------------------------------------')
    print('Computating BAe....')
    
    sigma = calc_static_stability(temp,min_lon,max_lon,min_lat,max_lat)
    za_temp = calc_zonal_ave(temp,min_lon,max_lon,min_lat,max_lat)
    dz_temp = temp - za_temp

    # 1) first integral
    
    tmp = uwnd*(dz_temp**2)
    
    lamb_w = tmp.sel(lon=min_lon)/2*sigma.sel(lon=min_lon)
    lamb_e = tmp.sel(lon=max_lon)/2*sigma.sel(lon=max_lon)
    
    tmp = (lamb_e-lamb_w)
    
    tmp = interga_phi(tmp,min_lat,max_lat)
    
    tmp = interga_P(tmp)
    
    tmp2 = (ra*(max_lon-min_lon)*\
                      (np.sin(max_lat)-np.sin(min_lat)))
    
    first_int = tmp/tmp2 

    # 2) second ingtegral

    tmp =  calc_zonal_ave(vwnd*(dz_temp**2),min_lon,max_lon,min_lat,max_lat)      
    
    phi_n = (tmp.sel(lat=max_lat)*np.cos(max_lat))/(2*sigma).sel(lat=max_lat)
    phi_s = (tmp.sel(lat=min_lat)*np.cos(min_lat))/(2*sigma).sel(lat=min_lat)
    
    tmp = (phi_n - phi_s)
    
    tmp = interga_P(tmp)

    tmp2 = (ra*(np.sin(max_lat)-np.sin(min_lat))) 

    second_int = (tmp/tmp2).mean('lon')
    
    # third term
    
    tmp = calc_area_ave(omg*(dz_temp**2),min_lon,max_lon,min_lat,max_lat)/(2*sigma)
    
    levs = tmp.level.values
    tmp_b = tmp.sel(level=np.amax(levs))/(2*sigma).sel(level=np.amax(levs))
    tmp_t = tmp.sel(level=np.amin(levs))/(2*sigma).sel(level=np.amin(levs))
    
    tmp = (tmp_b-tmp_t).mean('lon').mean('lat')
    
    BAe = -first_int-second_int-tmp
    
    print('Done!')
    
    return BAe

def calc_BKz(temp,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat):

    print('')    
    print('------------------------------------------')
    print('Computating BKz....')    
    
    za_uwnd = calc_zonal_ave(uwnd,min_lon,max_lon,min_lat,max_lat)
    za_vwnd = calc_zonal_ave(vwnd,min_lon,max_lon,min_lat,max_lat)

    dz_uwnd = uwnd - za_uwnd
    dz_vwnd = vwnd - za_vwnd
    
    # 1) first integral
    
    tmp = (uwnd**2 + vwnd**2 - dz_uwnd**2 - dz_vwnd**2)*vwnd
    
    tmp = (calc_zonal_ave(tmp,min_lon,max_lon,min_lat,max_lat))/(2*g)
    
    lamb_e = tmp.sel(lon=max_lon)
    lamb_w = tmp.sel(lon=min_lon)
    
    tmp = (lamb_e-lamb_w)/(2*g)
    
    tmp = interga_phi(tmp,min_lat,max_lat)
    
    tmp = interga_P(tmp)
    
    tmp2 = (ra*(max_lon-min_lon)*\
                      (np.sin(max_lat)-np.sin(min_lat)))
    
    first_int = tmp/tmp2 
    
    # 2) second integral
    
    tmp = ((uwnd**2 + vwnd**2 - dz_uwnd**2 - dz_vwnd**2)*vwnd)/(2*g)
    
    tmp = (calc_zonal_ave(tmp,min_lon,max_lon,min_lat,max_lat))/(2*g)
    
    phi_n = tmp.sel(lat=max_lat)*np.cos(max_lat)
    phi_s = tmp.sel(lat=min_lat)*np.cos(min_lat)
    
    tmp = (phi_n-phi_s)
    
    tmp = interga_P(tmp)
    
    tmp2 = (ra*(np.sin(max_lat)-np.sin(min_lat)))  
    
    second_int = (tmp/tmp2).mean('lon')
    
    # 3) Third term
    
    tmp = ((uwnd**2 + vwnd**2 - dz_uwnd**2 - dz_vwnd**2)*omg)/(2*g)
        
    tmp = calc_area_ave(tmp,min_lon,max_lon,min_lat,max_lat)   
    
    levs = tmp.level.values
    tmp_b = tmp.sel(level=np.amax(levs))
    tmp_t = tmp.sel(level=np.amin(levs))
    
    tmp = (tmp_b-tmp_t).mean('lon').mean('lat')
    
    BKz = -first_int-second_int-tmp 
    
    print('Done!')
    
    return BKz
    
    
def calc_BKe(temp,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat):

    print('')    
    print('------------------------------------------')
    print('Computating BKe....')
    
    za_uwnd = calc_zonal_ave(uwnd,min_lon,max_lon,min_lat,max_lat)
    za_vwnd = calc_zonal_ave(vwnd,min_lon,max_lon,min_lat,max_lat)

    dz_uwnd = uwnd - za_uwnd
    dz_vwnd = vwnd - za_vwnd

    # 1) First integral

    tmp = (dz_uwnd**2 + dz_vwnd**2)*uwnd
    
    tmp = (calc_zonal_ave(tmp,min_lon,max_lon,min_lat,max_lat))/(2*g)
    
    lamb_e = tmp.sel(lon=max_lon)
    lamb_w = tmp.sel(lon=min_lon)
    
    tmp = (lamb_e-lamb_w)/(2*g)
    
    tmp = interga_phi(tmp,min_lat,max_lat)
    
    tmp = interga_P(tmp)    
    
    tmp2 = (ra*(max_lon-min_lon)*\
                      (np.sin(max_lat)-np.sin(min_lat)))
    
    first_int = tmp/tmp2 

    #  2) Second integral

    tmp = ((dz_uwnd**2 + dz_vwnd**2)*vwnd)
    
    tmp = calc_zonal_ave(tmp,min_lon,max_lon,min_lat,max_lat)    
    
    phi_n = tmp.sel(lat=max_lat) * np.cos(max_lat)
    phi_s = tmp.sel(lat=min_lat) * np.cos(min_lat)
    
    tmp = (phi_n-phi_s)/(2*g)
    
    tmp = interga_P(tmp)    
        
    tmp2 = ra*(np.sin(max_lat)-np.sin(min_lat))
    
    second_int = (tmp/tmp2).mean('lon')
    
    # 3) third term
    
    tmp = ((dz_uwnd**2 + dz_vwnd**2)*omg)/(2*g)
    
    tmp = calc_area_ave(tmp,min_lon,max_lon,min_lat,max_lat)/(2*g)    
    
    levs = tmp.level.values
    tmp_b = tmp.sel(level=np.amax(levs))
    tmp_t = tmp.sel(level=np.amin(levs))
    
    tmp = (tmp_b-tmp_t).mean('lon').mean('lat')
    
    BKe = -first_int-second_int-tmp 
    
    print('Done!')

    return BKe    
    
def calc_BOz(hgt,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat):

    print('')    
    print('------------------------------------------')
    print('Computating BOz....')
    
    za_hgt = calc_zonal_ave(hgt,min_lon,max_lon,min_lat,max_lat)
    za_vwnd = calc_zonal_ave(vwnd,min_lon,max_lon,min_lat,max_lat)
    za_uwnd = calc_zonal_ave(uwnd,min_lon,max_lon,min_lat,max_lat)
    za_omg = calc_zonal_ave(omg,min_lon,max_lon,min_lat,max_lat)

    dz_hgt = hgt - za_hgt
    dz_uwnd = uwnd - za_uwnd
    
    aa_hgt = calc_area_ave(hgt,min_lon,max_lon,min_lat,max_lat)
    aa_omg = calc_area_ave(omg,min_lon,max_lon,min_lat,max_lat)

    da_hgt = za_hgt - aa_hgt
    da_omg = za_omg - aa_omg
    
    # 1) First integral
    
    tmp = (uwnd*hgt) - (dz_uwnd-dz_hgt)
    
    lamb_e = tmp.sel(lon=max_lon)
    lamb_w = tmp.sel(lon=min_lon)
    
    tmp = (lamb_e - lamb_w)/g
    
    tmp = interga_phi(tmp,min_lat,max_lat)
    
    tmp = interga_P(tmp)    
    
    tmp2 = (ra*(max_lon-min_lon)*\
                      (np.sin(max_lat)-np.sin(min_lat)))
    
    first_int = tmp/tmp2 

    # 2) Second integral

    tmp = za_vwnd * za_hgt
    
    phi_n = tmp.sel(lat=max_lat) * np.cos(max_lat)
    phi_s = tmp.sel(lat=min_lat) * np.cos(min_lat)
    
    tmp = (phi_n - phi_s)/g
    
    tmp = interga_P(tmp)    
        
    tmp2 = ra*(np.sin(max_lat)-np.sin(min_lat))
    
    second_int = (tmp/tmp2).mean('lon')     

    # 3) Third term
    
    tmp = ((aa_hgt*aa_omg) + \
           (calc_area_ave(da_hgt+da_omg,min_lon,max_lon,min_lat,max_lat)))/g

    levs = tmp.level.values
    tmp_b = tmp.sel(level=np.amax(levs))
    tmp_t = tmp.sel(level=np.amin(levs))
    
    tmp = (tmp_b-tmp_t).mean('lon').mean('lat')  
    
    BOz = -first_int-second_int-tmp 

    print('Done!')

    return BOz

def calc_BOe(hgt,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat):

    print('')    
    print('------------------------------------------')
    print('Computating BOe....')    

    za_hgt = calc_zonal_ave(hgt,min_lon,max_lon,min_lat,max_lat)
    za_vwnd = calc_zonal_ave(vwnd,min_lon,max_lon,min_lat,max_lat)
    za_uwnd = calc_zonal_ave(uwnd,min_lon,max_lon,min_lat,max_lat)
    za_omg = calc_zonal_ave(omg,min_lon,max_lon,min_lat,max_lat)

    dz_hgt = hgt - za_hgt
    dz_uwnd = uwnd - za_uwnd
    dz_vwnd = vwnd - za_vwnd
    dz_omg = omg - za_omg

    # 1) First integral
    
    tmp = (dz_uwnd*dz_hgt)
    
    lamb_e = tmp.sel(lon=max_lon)
    lamb_w = tmp.sel(lon=min_lon)
    
    tmp = (lamb_e - lamb_w)/g
    
    tmp = interga_phi(tmp,min_lat,max_lat)
    
    tmp = interga_P(tmp)    
    
    tmp2 = (ra*(max_lon-min_lon)*\
                      (np.sin(max_lat)-np.sin(min_lat)))
    
    first_int = tmp/tmp2 

    # 2) Second integral

    tmp = calc_zonal_ave(dz_hgt*dz_vwnd,min_lon,max_lon,min_lat,max_lat)
    
    phi_n = tmp.sel(lat=max_lat) * np.cos(max_lat)
    phi_s = tmp.sel(lat=min_lat) * np.cos(min_lat)
    
    tmp = (phi_n - phi_s)/g
    
    tmp = interga_P(tmp)    
        
    tmp2 = ra*(np.sin(max_lat)-np.sin(min_lat))
    
    second_int = (tmp/tmp2).mean('lon')    

    # 3) Third term

    tmp = calc_area_ave(dz_omg*dz_hgt,min_lon,max_lon,min_lat,max_lat)  
    
    levs = tmp.level.values
    tmp_b = tmp.sel(level=np.amax(levs))
    tmp_t = tmp.sel(level=np.amin(levs))
    
    tmp = (tmp_b-tmp_t).mean('lon').mean('lat') 
    
    BOe = -first_int-second_int-tmp     
    
    print('Done!')

    return BOe
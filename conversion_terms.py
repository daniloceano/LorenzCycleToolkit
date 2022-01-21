#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for conversion terms used for computate the Lorenz Energy Cycle

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com

"""
import numpy as np
from integrals import (TrazpezoidalIntegration, CalcZonalAverage,
                       CalcAreaAverage)
import metpy.calc as mpcalc
from derivatives import (calc_delvar_delp, calc_delvar_delphi)

# from integrals import (calc_area_ave, calc_zonal_ave,
#                        calc_static_stability, interga_P)


# constants
R = 287         # ideal gas cte.
g = 9.81        # gravity acceleration 
ra = 6378000    # Earth's radius
cp = 1004.0     # specific heat of air at cte. pressure



def calc_cz(omg,temp,min_lon,max_lon,min_lat,max_lat):
    """
    Computate conversion term of potencial to kinetic energy (zonal, both)
    
    Parameters
    ----------
    omg, temp: xarray.Dataset
        array containing omega e temperature data
    min_lon, max_lon, minlat, min_lon: float
        lons and lats for the box bounding the computation domain  
    
    Returns
    -------
    cz: xarray.Dataset
        Values corresponding of the conversion term
        Note that this function returns the CZ value in all grid points
    """
    
    print('')    
    print('------------------------------------------')
    print('Computating CZ....')
    
    # 1) calculate zonal devation from the area average
    # dev_area = zonal_ave - area_ave
    za_omg = calc_zonal_ave(omg,min_lon,max_lon,min_lat,max_lat)
    aa_omg = calc_area_ave(omg,min_lon,max_lon,min_lat,max_lat)
    da_omg = za_omg - aa_omg
    
    # 2) calculate devation from the area average - for temperature
    za_temp = calc_zonal_ave(temp,min_lon,max_lon,min_lat,max_lat)
    aa_temp = calc_area_ave(temp,min_lon,max_lon,min_lat,max_lat)
    da_temp = za_temp - aa_temp

    # 3) calcuate the area average of 1) * 2)    
    tmp = calc_area_ave(da_omg*da_temp,min_lon,max_lon,min_lat,max_lat)

    # 4) divide by pressure
    p = (tmp*0)+tmp.level
    tmp = tmp * R/p
       
    # 5) Integrate through all pressure levels
    cz = interga_P(tmp)
        
    # divide by -g
    cz = cz/-g
    
    print('Done!')
    
    return cz.mean('lon').mean('lat')

def calc_ca(vwnd,temp,omg,min_lon,max_lon,min_lat,max_lat):
    """
    Computate conversion term of potencial zonal to eddy kinetic energy
    
    Parameters
    ----------
    omg, temp: xarray.Dataset
        array containing omega e temperature data
    min_lon, max_lon, minlat, min_lon: float
        lons and lats for the box bounding the computation domain  
    
    Returns
    -------
    ce: xarray.Dataset
        Values corresponding of the conversion term
        Note that this function returns the CE value in all grid points
    """
    
    print('')    
    print('------------------------------------------')
    print('Computating CA....')
   
    #-----------------------------------
    # 1) First term divided in two steps   
    
    # 1.1) (v)*(T)/2*sigma*ra
    # breakdown of the first step:
    
    # 1.1.1) (v) = v - [v]
    dz_vwnd = vwnd - calc_zonal_ave(vwnd,
                                           min_lon,max_lon,min_lat,max_lat)
    # 1.1.2) (T) = T - [T]
    dz_temp = temp - calc_zonal_ave(temp,
                                           min_lon,max_lon,min_lat,max_lat)
    # 1.1.3)  sigma
    sigma = calc_static_stability(temp,min_lon,max_lon,min_lat,max_lat)
    
    
    # 1.2) del T*/del phi
    
    # 1.2.1) dev_area_T = zonal_T - area_T
    za_temp = calc_zonal_ave(temp,min_lon,max_lon,min_lat,max_lat)
    aa_temp = calc_area_ave(temp,min_lon,max_lon,min_lat,max_lat)
    da_temp  = za_temp - aa_temp
    # 1.2.2) derivate
    delphi_temp = calc_delvar_delphi(da_temp) 
    
    # 1.3) area average of those terms  multiplied
    tmp = calc_area_ave(dz_vwnd*dz_temp*delphi_temp,min_lon,max_lon,min_lat,max_lat)
    
    # get first term of the integral    
    first_int = tmp/(ra*sigma)
    
    #-----------------------------------
    # 2) Second term 
    
    # 2.1) (omega)*(T)
    # (omega) = omega - [omega] ((T) was previously calculated)
    dz_omg = omg - calc_zonal_ave(omg,
                                           min_lon,max_lon,min_lat,max_lat)
    # 2.2) del dev_areaT/ del_p
    delp_datemp = calc_delvar_delp(da_temp)
    
    # 2.3) area_ave of 2.1) * 2.2)
    tmp = calc_area_ave(dz_omg*dz_temp*dz_temp*delp_datemp,min_lon,max_lon,min_lat,max_lat)    
    
    # 2.3) get second term of the integral    
    
    second_int = tmp/sigma
    
    #-----------------------------------
    # 3) Sum and integrate for all levs
    
    ca = -interga_P(first_int+second_int)
    
    print('Done!')
    
    return ca.mean('lon').mean('lat')

def calc_ce(omg,temp,min_lon,max_lon,min_lat,max_lat):
    """
    Computate conversion term of potencial to kinetic energy (eddy, both)
    
    Parameters
    ----------
    omg, temp: xarray.Dataset
        array containing omega e temperature data
    min_lon, max_lon, minlat, min_lon: float
        lons and lats for the box bounding the computation domain  
    
    Returns
    -------
    ce: xarray.Dataset
        Values corresponding of the conversion term
        Note that this function returns the CE value in all grid points
    """
    
    print('')    
    print('------------------------------------------')
    print('Computating CE....')
    
    # 1) calculate devation from the area average
    # dev_area = zonalave - area_ave
    za_omg = calc_zonal_ave(omg,min_lon,max_lon,min_lat,max_lat)
    dz_omg = omg - za_omg 
    
    # 2) calculate devation from the area average - for temperature
    za_temp = calc_zonal_ave(temp,min_lon,max_lon,min_lat,max_lat)
    dz_temp = temp - za_temp
    
    # 3) calculate area_ave of 1) * 2)
    tmp = calc_area_ave(dz_omg*dz_temp,min_lon,max_lon,min_lat,max_lat)

    # 4) divide by pressure
    p = (tmp*0)+tmp.level
    tmp = tmp * R/p
        
    # 5) Integrate through all pressure levels
    ce = interga_P(tmp)
        
    # divide by -g
    ce = ce/-g
    
    print('Done!')
    
    return ce.mean('lon').mean('lat')

def calc_ck(uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat):
    """
    Computate conversion term of eddy to zonal kinetic energy
    
    Parameters
    ----------
    omg, temp: xarray.Dataset
        array containing omega e temperature data
    min_lon, max_lon, minlat, min_lon: float
        lons and lats for the box bounding the computation domain  
    
    Returns
    -------
    ce: xarray.Dataset
        Values corresponding of the conversion term
        Note that this function returns the CK value in all grid points
    """
    
    print('')    
    print('------------------------------------------')
    print('Computating CK....')
    
    # Computate averages that will be used
    za_uwnd = calc_zonal_ave(uwnd,min_lon,max_lon,min_lat,max_lat)
    za_vwnd = calc_zonal_ave(vwnd,min_lon,max_lon,min_lat,max_lat)
    za_omg = calc_zonal_ave(omg,min_lon,max_lon,min_lat,max_lat)
    # deviations from zonal averages
    dz_uwnd = uwnd - za_uwnd
    dz_vwnd = vwnd - za_vwnd
    dz_omg = omg - za_omg
   
   
    #-----------------------------------
    # 1) First intergal divided in two steps   
    
    # 1.1) (v)*(u)* cos(phi)/ra
    tmp1 = dz_uwnd*dz_vwnd*np.cos(dz_uwnd.lat)/ra
    
    # 1.2) del[u]/cos(phi)/del phi
    tmp2 = calc_delvar_delphi(za_uwnd/np.cos(za_uwnd.lat))
    
    # area avearge
    var = calc_area_ave(tmp1*tmp2,min_lon,max_lon,min_lat,max_lat)
    
    first_int = interga_P(var)
    
    #-----------------------------------
    # 2) Second intergal divided in two steps
    
    # 2.1) (v)**2/ra
    tmp1 = (dz_uwnd**2)/ra
    
    # 2.2) del [v]/ delphi
    tmp2 = calc_delvar_delphi(za_vwnd)

    # area avearge
    var = calc_area_ave(tmp1*tmp2,min_lon,max_lon,min_lat,max_lat)

    second_int = interga_P(var)
    
    #-----------------------------------
    # 3) Third intergal 

    # 1.1) [v]*(u)**2 * tan(phi)/ra
    tmp1 = za_uwnd * (dz_uwnd**2) * np.tan(dz_uwnd.lat)/ra

    # area avearge    
    var = calc_area_ave(tmp1,min_lon,max_lon,min_lat,max_lat)
    
    third_int = interga_P(var)
   
    #-----------------------------------
    # 4) Fourth intergal

    # 2.1) (omega)*(u)
    tmp1 = dz_omg*dz_uwnd
    
    # 2.2) del [u]/ delp
    tmp2 = calc_delvar_delp(za_uwnd)

    # area avearge
    var = calc_area_ave(tmp1*tmp2,min_lon,max_lon,min_lat,max_lat)

    fourth_int = interga_P(var)     
    
    #-----------------------------------
    # 5) Fifth intergal

    # 2.1) (omega)*(v)
    tmp1 = dz_omg*dz_uwnd
    
    # 2.2) del [u]/ delp
    tmp2 = calc_delvar_delp(za_vwnd)

    # zonal avearge
    var = calc_area_ave(tmp1*tmp2,min_lon,max_lon,min_lat,max_lat)

    fifth_int = interga_P(var)  
    
    # sum all integrals

    ck = (first_int + second_int + third_int + fourth_int + fifth_int)/g
       
    return ck.mean('lon').mean('lat')
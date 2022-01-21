#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:27:40 2020

@author: Danilo
"""

from integrals import (calc_area_ave, calc_zonal_ave,
                       calc_static_stability, interga_P)
import numpy as np

# constants
R = 287         # ideal gas cte.
g = 9.81        # gravity acceleration 
rho = 1.25      # air density
Cd = 2.6e-8     # drag coeficient

def calc_GE(omg,temp,min_lon,max_lon,min_lat,max_lat):
    
    # 1) zonal_dev of omega and T
    za_temp = calc_zonal_ave(temp,min_lon,max_lon,min_lat,max_lat)
    za_omg = calc_zonal_ave(omg,min_lon,max_lon,min_lat,max_lat)
    zd_temp = temp - za_temp
    zd_omg = omg - za_omg
    
    # 2) area_ave of  dev_temp * dev_omg
    tmp = calc_area_ave(zd_temp*zd_omg,min_lon,max_lon,min_lat,max_lat)
        
    # 4) divide by pressure
    p = (tmp*0)+tmp.level
    tmp = tmp * R/p
    
    # 5) Integrate through all pressure levels
    GE = interga_P(tmp)
    
    # divide by -g
    GE = GE/-g
    
    return GE.mean('lon').mean('lat')

def calc_drag(uwnd,vwnd):
    
    W  = np.sqrt(uwnd**2 + vwnd**2)
    taux = Cd*rho*W*uwnd
    tauy = Cd*rho*W*vwnd
    
    return taux, tauy
    

def calc_CE(uwnd,vwnd,temp,min_lon,max_lon,min_lat,max_lat):
      
    za_vwnd = calc_zonal_ave(vwnd,min_lon,max_lon,min_lat,max_lat)
    za_uwnd = calc_zonal_ave(uwnd,min_lon,max_lon,min_lat,max_lat)

    dz_uwnd = uwnd - za_uwnd
    dz_vwnd = vwnd - za_vwnd
    
    fric = calc_drag(uwnd,vwnd)[0]   
    fric_lamb, fric_phi = fric[0], fric[1]
    
    za_frclamb = calc_zonal_ave(fric_lamb,min_lon,max_lon,min_lat,max_lat)
    za_frcphi = calc_zonal_ave(fric_phi,min_lon,max_lon,min_lat,max_lat)

    dz_friclamb = fric_lamb -  za_frclamb  
    dz_fricphi = fric_phi -  za_frcphi      
    
    tmp = (dz_uwnd*dz_friclamb) + (dz_vwnd*dz_fricphi)
    
    tmp = calc_area_ave(tmp,min_lon,max_lon,min_lat,max_lat)
    
    tmp = interga_P(tmp)
    
    DE = -tmp/g
    
    return DE
    
    
    
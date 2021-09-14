#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for computating Lorenz Energy Cycle

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com

"""

import aux
from conversion_terms import (calc_cz, calc_ca, calc_ce, calc_ck)
from energy_contents import (calc_AZ, calc_AE, calc_KZ, calc_KE)
from boundary_terms import (calc_BAz, calc_BAe, calc_BKz, calc_BKe,
                            calc_BOz, calc_BOe)
from derivatives import calc_delvar_delt

# ---------------------------------
# Get data and set variables

## define spatial domain ##
min_lon, max_lon, min_lat, max_lat = -60, -30, -42.5, -17.5
aux.plot_area(min_lon, max_lon, min_lat, max_lat)

## open file and get vars ##
data = aux.get_data('ncepR2.nc')
uwnd, vwnd, omg, temp, rhum, hgt = data[0],data[1],data[2], data[3], data[4], data[5]



# -----------------------------------------
# Computate conversion terms and plot values

Cz = calc_cz(omg,temp,min_lon,max_lon,min_lat,max_lat)
Ca = calc_ca(vwnd,temp,omg,min_lon,max_lon,min_lat,max_lat)
Ce = calc_ce(omg,temp,min_lon,max_lon,min_lat,max_lat)
Ck = calc_ck(uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat)

aux.plot_timeseries(Cz,Ca,Ce,Ck,1)


# ----------------------------------------------
# Computate energy content terms and plot values

Az = calc_AZ(temp,min_lon,max_lon,min_lat,max_lat)/10e2
Ae = calc_AE(temp,min_lon,max_lon,min_lat,max_lat)/10e2
Kz = calc_KZ(uwnd,vwnd,min_lon,max_lon,min_lat,max_lat)/10e2
Ke = calc_KE(uwnd,vwnd,min_lon,max_lon,min_lat,max_lat)/10e2

aux.plot_timeseries(Az,Ae,Kz,Ke,2)

# ----------------------------------------------
# Computate boundary fluxes and plot values

BAz = calc_BAz(temp,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat)
BAe = calc_BAe(temp,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat)
BKz = calc_BKz(temp,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat)
BKe = calc_BKe(temp,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat)
BOz = calc_BOz(hgt,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat)
BOe = calc_BOe(hgt,uwnd,vwnd,omg,min_lon,max_lon,min_lat,max_lat)

aux.plot_timeseries(BAz,BAe,BKz,BKe,3)


# ----------------------------------------------
# Computate energy terms and plot values
del_Az = calc_delvar_delt(Az) 
del_Ae = calc_delvar_delt(Ae)
del_Kz = calc_delvar_delt(Kz)
del_Ke = calc_delvar_delt(Ke)


# ----------------------------------------------
# Computate generation/dissipation terms and plot values

Gz = del_Az + Ca + Cz - BAz
Ge = del_Ae - Ca + Ce - BAe
Dz = del_Kz - Ck - Cz - BKz - BOz
De = del_Ke - Ck - Ce - BKe - BOe


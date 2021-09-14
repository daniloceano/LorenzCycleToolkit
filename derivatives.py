#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for computating vertical integrations and partial derivatives

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com

"""

import numpy as np


def calc_delvar_delp(var):
    """
    Computates the partial derivative of some variable by each pressure level as:
        del var / del p
    using a second-order finite difference scheme for intermediate levels
    and a first-order forward/backward scheme for bottom/top levels
    
    Parameters
    ----------
    var: xarray.Dataset
        the variable to be integrated    
    """
    # get pressure levels
    levs = var.level.values
    # make array for derivatives
    del_var = var*0
    for lev in range(len(levs)):
        if lev == 0:
            dp = (levs[lev+1] - levs[lev])
            del_var.loc[dict(level=levs[lev])] = \
            (var.sel(level=levs[lev+1]) - var.sel(level=levs[lev]))/dp
        elif lev == (len(levs)-1):
            dp = (levs[lev] - levs[lev-1])
            del_var.loc[dict(level=levs[lev])] = \
            (var.sel(level=levs[lev]) - var.sel(level=levs[lev-1]))/dp
        else:
            dp = (levs[lev+1] - levs[lev-1])
            del_var.loc[dict(level=levs[lev])] = \
            (var.sel(level=levs[lev+1]) - var.sel(level=levs[lev-1]))/dp

    return del_var

def calc_delvar_delphi(var):
    """
    Computates the partial derivative of some variable by each latitude as:
        del var / del phi
    using a second-order finite difference scheme for intermediate levels
    and a first-order forward/backward scheme for left/right boundaries
    
    Parameters
    ----------
    var: xarray.Dataset
        the variable to be integrated    
    """
    # get latitudes
    lats = var.lat.values
    # make array for derivatives
    del_var = var*0
    for lat in range(len(lats)):
        if lat == 0:
            dy = (lats[lat+1] - lats[lat])
            del_var.loc[dict(lat=lats[lat])] = \
            (var.sel(lat=lats[lat+1]) - var.sel(lat=lats[lat]))/dy
        elif lat == (len(lats)-1):
            dy = (lats[lat] - lats[lat-1])
            del_var.loc[dict(lat=lats[lat])] = \
            (var.sel(lat=lats[lat]) - var.sel(lat=lats[lat-1]))/dy
        else:
            dy = (lats[lat+1] - lats[lat-1])
            del_var.loc[dict(lat=lats[lat])] = \
            (var.sel(lat=lats[lat+1]) - var.sel(lat=lats[lat-1]))/dy

    return del_var

def calc_delvar_delt(var):
    """
    Computates the partial derivative of some variable by each time as:
        del var / del t
    using a second-order finite difference scheme for intermediate times
    and a first-order forward/backward scheme for initial/final time steps
    
    Parameters
    ----------
    var: xarray.Dataset
        the variable to be integrated    
    """
    # get latitudes
    times = var.time.values
    # make array for derivatives
    del_var = var*0
    for t in range(len(times)):
        if t == 0:
            dt = (times[t+1] - times[t])/np.timedelta64(1, 's')
            del_var.loc[dict(time=times[t])] = \
            (var.sel(time=times[t+1]) - var.sel(time=times[t]))/dt
        elif t == (len(times)-1):
            dt = (times[t] - times[t-1])/np.timedelta64(1, 's')
            del_var.loc[dict(time=times[t])] = \
            (var.sel(time=times[t]) - var.sel(time=times[t-1]))/dt
        else:
            dt = (times[t+1] - times[t-1])/np.timedelta64(1, 's')
            del_var.loc[dict(time=times[t])] = \
            (var.sel(time=times[t+1]) - var.sel(time=times[t-1]))/dt

    return del_var

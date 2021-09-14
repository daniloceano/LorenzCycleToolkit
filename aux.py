#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Axuliary functions 

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com

"""

import cartopy.crs as ccrs
import xarray as xr
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon
import numpy as np
import matplotlib.dates as mdates


# ---------------------------

def get_data(filename):
    """
    Get data and variables from a NCEP Reanalysis 2 netCDF file
    """    
    
    file = convert_lon(xr.open_dataset('ncepR2.nc'))
    
    uwnd = file.uwnd
    vwnd = file.vwnd 
    omg = file.omega
    temp = file.air
    rhum = file.rhum
    hgt = file.hgt
    
    return uwnd, vwnd, omg, temp, rhum, hgt

# ---------------------------
# auxilliary functions for plot maps

def convert_lon(df):
    
    """
    Convert longitudes from 0:360 range to -180:180
    """
    
    df.coords['lon'] = (df.coords['lon'] + 180) % 360 - 180
    df = df.sortby(df.lon)
    
    return df

def plot_area(min_lon, max_lon, min_lat, max_lat):
    
    """
    Plot the area limited by the lons and lats values,
    that will be used for the computations
    """
    
    datacrs = ccrs.PlateCarree() # projection
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes([0, 0, 0.95, 0.95], projection=datacrs,
                  frameon=True)
    ax.set_extent([-90, 0, 0, -60], crs=datacrs)
    ax.coastlines(zorder = 1)
    ax.stock_img()
    # plot selected domain
    # create a sample polygon, `pgon`
    pgon = Polygon(((min_lon, min_lat),
            (min_lon, max_lat),
            (max_lon, max_lat),
            (max_lon, min_lat),
            (min_lon, min_lat)))
    ax.add_geometries([pgon], crs=datacrs, 
                      facecolor='red', edgecolor='k', linewidth = 3,
                      alpha=0.5, zorder = 3)
    ax.gridlines(draw_labels=True,zorder=2)    

    plt.title('Box defined for compuations \n', fontsize = 22)

# -----------------------------
# plotting time series of terms

def align_yaxis(ax1, v1, ax2, v2):
    
    """
    adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1
    """
    
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)
    
def align_yaxis_np(ax1, ax2):
    
    """
    Align zeros of the two axes, zooming them out by same ratio
    """
    
    axes = np.array([ax1, ax2])
    extrema = np.array([ax.get_ylim() for ax in axes])
    tops = extrema[:,1] / (extrema[:,1] - extrema[:,0])
    # Ensure that plots (intervals) are ordered bottom to top:
    if tops[0] > tops[1]:
        axes, extrema, tops = [a[::-1] for a in (axes, extrema, tops)]

    # How much would the plot overflow if we kept current zoom levels?
    tot_span = tops[1] + 1 - tops[0]

    extrema[0,1] = extrema[0,0] + tot_span * (extrema[0,1] - extrema[0,0])
    extrema[1,0] = extrema[1,1] + tot_span * (extrema[1,0] - extrema[1,1])
    [axes[i].set_ylim(*extrema[i]) for i in range(2)]    


def plot_timeseries(var1,var2,var3,var4,opt):
    
    """
    Plot time series of variables 1, 2, 3 and 4
    
    Parameters
    ----------
    var1,var2,var3,var4: xarray.DataArray
        One dimensional arrays containig the time series of some variable
    opt: integer (1 or 2)
        1 for plotting the conversion terms
        2 for the energy contents
        3 for boundaries
        4 for energy derivatives
    """  
    
    time = var1.time
    
    # create figure
    fig = plt.figure(figsize=(8,5))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    
    
    # plot 0 line
    ax2.plot(time,np.zeros(len(var1)),
             linewidth=2,color='grey', alpha = 0.8)
    
    if opt == 1:
    
        # plot CZ in a separated axis and other terms in the same
        lns1 = ax1.plot(time,var1,c='#473BF0', linewidth=5, label='CZ') 
        lns2 = ax2.plot(time,var2,c='#384A0F', linewidth=5, label='CA')
        lns3 = ax1.plot(time,var3,c='#C9B857', linewidth=5, label='CE')
        lns4 = ax2.plot(time,var4,c='#A53860', linewidth=5, label='CK')

        # set y-axis labels
        ax1.set_ylabel('Energy flux for CZ and CE (w $m^{-2}$ )', fontsize = 18)                
        ax2.set_ylabel('CA and CK (w $m^{-2}$ )',
                       rotation=270, fontsize = 18, labelpad=25)
        
        plt.title('Conversion terms', fontsize=22)

    if opt == 2:

        # plot AZ and AE in one axis and KZ and KE in another        
        lns1 = ax1.plot(time,var1,c='#473BF0', linewidth=5, label='AZ')
        lns2 = ax1.plot(time,var2,c='#384A0F', linewidth=5, label='AE')
        lns3 = ax1.plot(time,var3,c='#C9B857', linewidth=5, label='KZ')
        lns4 = ax1.plot(time,var4,c='#A53860', linewidth=5, label='KE')

        # set y-axis labels
        ax1.set_ylabel(' (J $m^{-2}$ )', fontsize = 18)                
#        ax2.set_ylabel('Eddy (J $m^{-2}$ )',
#                       rotation=270, fontsize = 18, labelpad=25)     

        plt.title('Energy contents', fontsize=22)
        
    if opt == 3:

        # plot AZ and AE in one axis and KZ and KE in another        
        lns1 = ax1.plot(time,var1,c='#473BF0', linewidth=5, label='BAz')
        lns2 = ax1.plot(time,var2,c='#384A0F', linewidth=5, label='BAe')
        lns3 = ax2.plot(time,var3,c='#C9B857', linewidth=5, label='BKz')
        lns4 = ax2.plot(time,var4,c='#A53860', linewidth=5, label='BKe')

        # set y-axis labels
        ax1.set_ylabel('Potential (w $m^{-2}$ )', fontsize = 18)                
        ax2.set_ylabel('Kinetic (w $m^{-2}$ )',
                       rotation=270, fontsize = 18, labelpad=25)     

        plt.title('Boundary terms', fontsize=22)        

       
    # plot params     
#    plt.xlim(0,len(var1)-1)
    plt.xlim(np.amin(time),np.amax(time))    
    ax1.grid(axis='x')
    plt.xlabel('Time step', fontsize=18)
    plt.rc('xtick',labelsize=16)
    plt.rc('ytick',labelsize=16)
    plt.xticks(time,label=var1.time)    
    myFmt = mdates.DateFormatter('%d/%m')
    ax1.xaxis.set_major_formatter(myFmt)
    fig.autofmt_xdate()
    start, end = ax1.get_xlim()
    ax1.xaxis.set_ticks(np.arange(start, end, 1))    

    
    # legend
    lns = lns1+lns2+lns3+lns4
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs)
    # align 0 from both axis
    align_yaxis_np(ax1, ax2)
    
# ---------------------------
    

cols = ['#1724B2', '#473BF0', '#648DE5', '#9EB7E5', '#E8E5DA',
        '#384A0F', '#384A0F', '#527A0D', '#54733B', '#4B5E44',
        '#786604', '#C9B857', '#BDA419', '#D6BB1E', '#E5F77D',
        '#450920', '#A53860', '#DA627D'] 
        
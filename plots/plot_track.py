#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 19:48:17 2022

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
    
"""


import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon
import pandas as pd
import sys
import cmocean
import matplotlib.colors as colors
import numpy as np
from bisect import bisect_left


def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before

# Plot the area limited by the lons and lats values that will be used
# for the computations
def main():

    track = pd.read_csv('./track',parse_dates=[0],delimiter=';',index_col='time')
    outfile = sys.argv[1]
    ResultsSubDirectory = '/'.join(outfile.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/'
        
    min_lon,max_lon,min_lat,max_lat=min(track['Lon']),max(track['Lon']),\
        min(track['Lat']),max(track['Lat'])
    
    # DataFrame with results
    df = pd.read_csv(outfile)
    df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')
    
    plt.close('all')
    datacrs = ccrs.PlateCarree() # projection
    # If the figure lenght > height the fig size and legend position have
    # ot be adjusted otherwise they won't appear
    lenx, leny = max_lon - min_lon, max_lat - min_lat
    if lenx <= leny:
        fig = plt.figure(figsize=(12, 9))
        labelx = 0.73
    else:
        fig = plt.figure(figsize=(14, 9))
        labelx = 0.67
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection=datacrs,
                  frameon=True)
    ax.set_extent([min_lon-10, max_lon+10, max_lat+10, min_lat-10], crs=datacrs)
    ax.coastlines(zorder = 1)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("lightblue"))
    gl = ax.gridlines(draw_labels=True,zorder=2,linestyle='dashed',alpha=0.8,
                 color='#383838')
    gl.xlabel_style = {'size': 14, 'color': '#383838'}
    gl.ylabel_style = {'size': 14, 'color': '#383838'}
    gl.bottom_labels = None
    gl.right_labels = None
    
    for i in range(len(df['Datetime'])):
        
        # Model timestep
        itime = str(df['Datetime'].iloc[i])
        
        itime_track = track.index[track.index==itime]
        lon = track.loc[itime_track]['Lon'].values
        lat = track.loc[itime_track]['Lat'].values
        
        min_lon, max_lon = lon-7.5,lon+7.5
        min_lat, max_lat = lat-7.5,lat+7.5
        
        c = '#383838'
        
        # plot selected domain
        if i == 0:
            pgon = Polygon(((min_lon, min_lat),
                    (min_lon, max_lat),
                    (max_lon, max_lat),
                    (max_lon, min_lat),
                    (min_lon, min_lat)))
            ax.add_geometries([pgon], crs=datacrs, 
                              facecolor='none', edgecolor=c,
                              linewidth = 5,
                              alpha=0.75, zorder=1+i)
        
    
    # use potential energy for dot color
    Ae = df['Ae']
    # dots size as kinetic energy
    msizes = [200,400,600,800,1000]
    
    intervals = [3e5,4e5,5e5,6e5]
    data = df['Ke']
    title = 'Eddy Kinect\n    Energy\n(Ke - '+r' $J\,m^{-2})$'
    sizes = []
    for val in data:
        if val <= intervals[0]:
            sizes.append(msizes[0])
        elif val > intervals[0] and val <= intervals[1]:
            sizes.append(msizes[1])
        elif val > intervals[1] and val <= intervals[2]:
            sizes.append(msizes[2])
        elif val > intervals[2] and val <= intervals[3]:
            sizes.append(msizes[3])
        else:
            sizes.append(msizes[4])
    df['sizes'] = sizes
    # Plot legend
    labels = ['< '+str(intervals[0]),
              '< '+str(intervals[1]),
              '< '+str(intervals[2]),
              '< '+str(intervals[3]),
              '> '+str(intervals[3])]
    l1 = plt.scatter([],[],c='k', s=msizes[0],label=labels[0])
    l2 = plt.scatter([],[], c='k', s=msizes[1],label=labels[1])
    l3 = plt.scatter([],[],c='k', s=msizes[2],label=labels[2])
    l4 = plt.scatter([],[],c='k', s=msizes[3],label=labels[3])
    l5 = plt.scatter([],[],c='k', s=msizes[4],label=labels[4])
    leg = plt.legend([l1, l2, l3, l4, l5], labels, ncol=1, frameon=False,
                     fontsize=10, handlelength = 0.3, handleheight = 4,
                     borderpad = 1.5, scatteryoffsets = [0.1], framealpha = 1,
                handletextpad=1.5, title=title,
                scatterpoints = 1, loc = 1,
                bbox_to_anchor=(labelx, -0.57, 0.5, 1),labelcolor = '#383838')
    leg._legend_box.align = "center"
    plt.setp(leg.get_title(), color='#383838')
    plt.setp(leg.get_title(),fontsize=12)
    for i in range(len(leg.legendHandles)):
        leg.legendHandles[i].set_color('#383838')
        leg.legendHandles[i].set_edgecolor('gray')
    
    plt.plot(track['Lon'], track['Lat'],c='#383838')
    norm = colors.TwoSlopeNorm(vmin=min(Ae), vcenter=np.mean(Ae), vmax=max(Ae))
    dots = plt.scatter(track['Lon'].loc[df['Datetime']],
                       track['Lat'].loc[df['Datetime']],
                       c=Ae,cmap=cmocean.cm.matter,s=df['sizes'],zorder=100, 
                       edgecolors='grey', norm=norm)
    
    # colorbar
    cax = fig.add_axes([ax.get_position().x1+0.02,
                    ax.get_position().y0+0.38,0.02,ax.get_position().height/1.74])
    cbar = plt.colorbar(dots, extend='both',cax=cax)
    cbar.ax.set_ylabel('Eddy Potential Energy (Ae - '+r' $J\,m^{-2})$',
                   rotation=270,fontsize=12,verticalalignment='bottom',
                   c='#383838',labelpad=15)
    
    # mark beginning and end of system
    start = [track.loc[df['Datetime'][0]]['Lon'],
             track.loc[df['Datetime'][0]]['Lat']]
    end =  [track.loc[df['Datetime'].iloc[-1]]['Lon'],
             track.loc[df['Datetime'].iloc[-1]]['Lat']]
    ax.text(*start,'A',zorder=101,fontsize=24,horizontalalignment='center',
            verticalalignment='center')
    ax.text(*end,'Z',zorder=101,fontsize=24,horizontalalignment='center',
            verticalalignment='center')
    
    # plt.title('System track and boxes defined for compuations \n', fontsize = 22)
    plt.savefig(FigsDir+'track_boxes.png',bbox_inches='tight')
    print('\nCreated figure with track and boxes defined for computations')
    
if __name__ == "__main__":
    main()
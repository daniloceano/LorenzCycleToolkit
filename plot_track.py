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
from LorenzPhaseSpace import MarkerSizeKe
import cmocean
import matplotlib.colors as colors
import numpy as np


# Plot the area limited by the lons and lats values that will be used
# for the computations
def main():

    track = pd.read_csv('./track',parse_dates=[0],delimiter=';',
                        names=['time','Lat','Lon'],index_col='time')
    outfile = sys.argv[1]
    ResultsSubDirectory = '/'.join(outfile.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/'
        
    min_lon,max_lon,min_lat,max_lat=min(track['Lon']),max(track['Lon']),\
        min(track['Lat']),max(track['Lat'])
    plt.close('all')
    datacrs = ccrs.PlateCarree() # projection
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection=datacrs,
                  frameon=True)
    ax.set_extent([min_lon-20, max_lon+20, max_lat+20, min_lat-20], crs=datacrs)
    ax.coastlines(zorder = 1)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("lightblue"))
    # ax.stock_img()
    gl = ax.gridlines(draw_labels=True,zorder=2,linestyle='dashed',alpha=0.8,
                 color='#383838')    
    gl.bottom_labels = None
    gl.right_labels = None
    
    for i in range(len(track)):
    
        itime = str(track.iloc[i].name)
        lon,lat = track.loc[itime]['Lon'], track.loc[itime]['Lat']
        
        min_lon, max_lon = lon-7.5,lon+7.5
        min_lat, max_lat = lat-7.5,lat+7.5
        
        if i == 0 or i == len(track)-1:
            # plot selected domain
            # create a sample polygon, `pgon`
            pgon = Polygon(((min_lon, min_lat),
                    (min_lon, max_lat),
                    (max_lon, max_lat),
                    (max_lon, min_lat),
                    (min_lon, min_lat)))
            ax.add_geometries([pgon], crs=datacrs, 
                              facecolor='none', edgecolor='#383838',
                              linewidth = 5,
                              alpha=0.75, zorder=1+i)
        
    df = pd.read_csv(outfile)
    # use potential energy for dot color
    Ae = df['Ae']
    # dots size as kinetic energy
    s = MarkerSizeKe(df,1)['sizes']
    plt.plot(track['Lon'], track['Lat'],c='#383838')
    norm = colors.TwoSlopeNorm(vmin=min(Ae), vcenter=np.mean(Ae), vmax=max(Ae))
    dots = plt.scatter(track['Lon'], track['Lat'],c=Ae,cmap=cmocean.cm.matter,
                       s=s,zorder=100, edgecolors='grey', norm=norm)
    # colorbar
    cax = fig.add_axes([ax.get_position().x1+0.02,
                    ax.get_position().y0+0.38,0.02,ax.get_position().height/1.74])
    cbar = plt.colorbar(dots, extend='both',cax=cax)
    cbar.ax.set_ylabel('Eddy Potential Energy (Ae - '+r' $J\,m^{-2})$',
                   rotation=270,fontsize=12,verticalalignment='bottom',
                   c='#383838',labelpad=15)
    
    # mark beginning and end of system
    start = [track.iloc[0]['Lon'],track.iloc[0]['Lat']]
    end = [track.iloc[-1]['Lon'],track.iloc[-1]['Lat']]
    ax.text(*start,'A',zorder=101,fontsize=24,horizontalalignment='center',
            verticalalignment='center')
    ax.text(*end,'Z',zorder=101,fontsize=24,horizontalalignment='center',
            verticalalignment='center')
    
    # plt.title('System track and boxes defined for compuations \n', fontsize = 22)
    plt.savefig(FigsDir+'track_boxes.png')
    print('\nCreated figure with track and boxes defined for computations')
    
if __name__ == "__main__":
    main()
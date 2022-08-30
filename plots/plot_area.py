#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:58:51 2022

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
    
"""


import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon
import sys
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature, COASTLINE
from cartopy.feature import BORDERS

def map_features(ax):
    ax.add_feature(COASTLINE,edgecolor='#283618',linewidth=1)
    ax.add_feature(BORDERS,edgecolor='#283618',linewidth=1)
    return ax

def Brazil_states(ax):    
    
    _ = ax.add_feature(cfeature.NaturalEarthFeature('physical',
                        'land', '50m', edgecolor='face', facecolor='#a4ab98'))
    
    states = NaturalEarthFeature(category='cultural', scale='50m', 
                                 facecolor='none',
                                  name='admin_1_states_provinces_lines')
    _ = ax.add_feature(states, edgecolor='#283618',linewidth=1)
    
    cities = NaturalEarthFeature(category='cultural', scale='50m',
                                 facecolor='none',
                                  name='populated_places')
    _ = ax.add_feature(cities, edgecolor='#283618',linewidth=1)
    
    
    
    
# Plot the area limited by the lons and lats values that will be used
# for the computations
def main():

    min_lon, max_lon,min_lat, max_lat = [float(i) for i in sys.argv[1:-1]]
    outdir = sys.argv[-1]
    
    plt.close('all')
    datacrs = ccrs.PlateCarree() # projection
    fig = plt.figure(figsize=(8, 8.5))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.83], projection=datacrs,
                  frameon=True)
    ax.set_extent([min_lon-20, max_lon+20, max_lat+20, min_lat-20], crs=datacrs)
    Brazil_states(ax)
    map_features(ax)
    
    # plot selected domain
    # create a sample polygon, `pgon`
    pgon = Polygon(((min_lon, min_lat),
            (min_lon, max_lat),
            (max_lon, max_lat),
            (max_lon, min_lat),
            (min_lon, min_lat)))
    ax.add_geometries([pgon], crs=datacrs, 
                      facecolor='None', edgecolor='#BF3D3B', linewidth = 3,
                      alpha=1, zorder = 3)
    gl = ax.gridlines(draw_labels=True,zorder=2)    
    gl.xlabel_style = {'size': 16}
    gl.ylabel_style = {'size': 16}

    plt.title('Box defined for compuations \n', fontsize = 22)
    plt.savefig(outdir+'Figures/box.png')
    print('\nCreated figure with box defined for computations')
    
if __name__ == "__main__":
    main()
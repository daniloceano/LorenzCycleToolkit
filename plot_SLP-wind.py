#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 17:12:04 2022

@author: daniloceano
"""

import xarray as xr
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import numpy as np 
import matplotlib.gridspec as gridspec
from matplotlib import colors
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature, COASTLINE
from cartopy.feature import BORDERS
from celluloid import Camera
import pandas as pd
import sys
from main import check_create_folder

def convert_lon(df,LonIndexer):
    df.coords[LonIndexer] = (df.coords[LonIndexer] + 180) % 360 - 180
    df = df.sortby(df[LonIndexer])
    return df

def grid_labels_params(ax):
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5,linestyle='--')
    gl.top_labels = False
    gl.left_labels = False
    gl.xlabel_style = {'size': 14, 'color': 'gray'}
    gl.ylabel_style = {'size': 14, 'color': 'gray'}
    ax.spines['geo'].set_edgecolor('gray')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax

def map_features(ax):
    ax.add_feature(COASTLINE)
    ax.add_feature(BORDERS, edgecolor='white')
    return ax

def Brazil_states(ax):    
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                  name='admin_1_states_provinces_lines')
    _ = ax.add_feature(states, edgecolor='white')
    
    cities = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                  name='populated_places')
    _ = ax.add_feature(cities)

def main():
    path = '/'.join(file.split('/')[:-1])
    infile_prefix = ''.join(file.split('/')[-1].split('.')[:-1])
    gif_folder = path+'/Surface_'+infile_prefix
    check_create_folder(gif_folder)
    
    try:
        dfVars = pd.read_csv('./fvars',sep= ';',index_col=0,header=0)
    except:
        raise SystemExit('Could not open fvars file. Check if file exists.\
 Please check path and file format')    
 
    LonIndexer,LatIndexer,TimeIndexer,LevelIndexer = \
      dfVars.loc['Longitude']['Variable'],dfVars.loc['Latitude']['Variable'],\
      dfVars.loc['Time']['Variable'],dfVars.loc['Vertical Level']['Variable']
    
    try:
        data = convert_lon(xr.open_dataset(file),'lon_2')
    except:
        raise SystemExit('Could not open specified file.\
 Please check path and file format')
 
    data = data.sel(lat_2=slice(0,-60),lon_2=slice(-100,10))
    slp,u,v = data[dfVars.loc['Sea Level Pressure']['Variable']],\
              data[dfVars.loc['Eastward Wind Component']['Variable']],\
              data[dfVars.loc['Northward Wind Component']['Variable']]
    sfc_lev = float(max(data[LevelIndexer]).values)
    
    lon,lat = data[LonIndexer], data[LatIndexer]
    proj = ccrs.PlateCarree() 
    
    cmap = cmo.balance
    levels = np.arange(1000,1036,2)
    norm = colors.TwoSlopeNorm(vmin=1010, vcenter=1014, vmax=1030)
    lims = [-60, -20, -35, -10]
    

    for i in range(len(data[TimeIndexer])):
        plt.close('all')
        fig = plt.figure(constrained_layout=False,figsize=(10,10))
        gs = gridspec.GridSpec(1, 1, hspace=0.1, wspace=0,
                                   left=0, right=0.95)
        ax = fig.add_subplot(gs[0], projection=proj)
        ax.set_extent(lims) 
        islp = slp.isel({TimeIndexer:i})/100
        iu, iv = u.isel({TimeIndexer:i}).sel({LevelIndexer:sfc_lev}),\
            v.isel({TimeIndexer:i}).sel({LevelIndexer:sfc_lev})
        
        cf1 = ax.contourf(lon, lat, islp, levels=levels, cmap=cmap,
                          norm=norm) 
        ax.contour(lon, lat, islp, cf1.levels,colors='grey', linewidths=1) 
        qv = ax.quiver(lon[::20],lat[::20],iu[::20,::20],
                       iv[::20,::20],color= 'k')
        ax.quiverkey(qv,-0.3, 1.07, 10, r'$10 \frac{m}{s}$', labelpos = 'E',
                           coordinates='axes', labelsep = 0.05,
                           fontproperties={'size': 14, 'weight': 'bold'})
        grid_labels_params(ax)
        map_features(ax)
        Brazil_states(ax)
        timestr = pd.to_datetime(str(data[TimeIndexer][i].values))
        date = timestr.strftime('%Y-%m-%dT%H%MZ')
        ax.text(0.05,1.01,date, transform=ax.transAxes, fontsize=16)
        ofile = gif_folder+'/'+date+'.png'
        plt.savefig(ofile)
        print(ofile+' created')

    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Error! No file specified.')
        quit()
    file = sys.argv[1]
    print('making surface charts for file: '+file)
    main()
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 15:03:35 2022

@author: danilocoutodsouza
"""

import xarray as xr
import sys
import os

if __name__ == "__main__":
    
    file3D,file2D = sys.argv[1],sys.argv[2]
    output = sys.argv[3]
    
    data3D = xr.open_dataset(file3D)
    data2D = xr.open_dataset(file2D)
        
    
    for coord in data3D.coords:
        if data3D.coords[coord].long_name in ['latitude','lat']:
            lat1 = coord
        elif data3D.coords[coord].long_name in ['longitude','lon']:
            lon1 = coord
        if 'level' in  data3D.coords[coord].long_name:
            lev1 = coord
        if 'time' in  data3D.coords[coord].long_name:
            time1 = coord  
    
    for coord in data2D.coords:
        if data2D.coords[coord].long_name in ['latitude','lat']:
            lat2 = coord
        elif data2D.coords[coord].long_name in ['longitude','lon']:
            lon2 = coord
        if 'time' in  data2D.coords[coord].long_name:
            time2 = coord  
        if 'level' in  data2D.coords[coord].long_name:
            lev2 = coord  
        else:
            lev2 = False
         
    # TO DO: If the data already contains vertical levels do something else
    if lev2 != False:
        pass
    # If data is 2D
    else:
        data2D = data2D.interp({lat2:data3D[lat1].values,
                lon2:data3D[lon1].values})
        # In some cases, the 2D data comes with extra time steps, so let's get
        # rid of them
        try:
            data2D  = data2D.sel({time2:slice(data3D[time1][0],data3D[time1][-1])})
        except:
            pass
        
    # Rename data2D coordinates so they match data3D
    data2D = data2D.rename({lon2:lon1}).rename({lat2:lat1})
    
    full_data = xr.merge([data3D,data2D])
    
    # Save
    DataDir = os.path.dirname(os.path.realpath(file3D))
    full_data.to_netcdf(path=DataDir+'/'+output+'.nc')
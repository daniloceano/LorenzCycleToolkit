#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculating zonal average and area average of some variable
Used for computate Lorenz Energy Cycle

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com

"""

import numpy as np
from metpy.units import units


def HorizontalTrazpezoidalIntegration(VariableData,dimension):
    """
    Integrate data using the trapezoidal rule 
    
    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated
    dimension: string
        the indexer used for the xarray coordante which the integration will
        be performed
   
    Returns
    -------
    trapz: xarray.Dataset
        Arrays with the integrated data
    """
    DimensionArray = VariableData[dimension] # array containing either lats or lons 
    trapz = 0
    delta = np.deg2rad(np.abs((DimensionArray[-1]-DimensionArray[0])/(len(DimensionArray)-1)))
    delta = units('radians')*delta
    for i in DimensionArray:
        if i == DimensionArray[0] or i == DimensionArray[-1]:
            trapz += VariableData.sel({dimension:i})
        else:
            trapz += 2*(VariableData.sel({dimension:i}))
    return trapz*(delta/2)

def VerticalTrazpezoidalIntegration(VariableData,VerticalAxis,VerticalCoordIndexer):
    """
    Integrate data using the trapezoidal rule 
    
    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated
    VerticalCoordIndexer: string
        the indexer used for the xarray coordante which the integration will
        be performed
   
    Returns
    -------
    trapz: xarray.Dataset
        Arrays with the integrated data
    """
    # # Convert do Pa if pressure is in hPa
    # DimensionArray = VariableData[VerticalCoordIndexer]
    # if DimensionArray.metpy.units == units('hPa'):
    #     DimensionArray = DimensionArray.metpy.convert_units('Pa')
    #     DimensionArray = DimensionArray*units('Pa')
    # Sort from top to bottom
    if VerticalAxis.metpy.units == units('Pa'):
        VerticalAxis = VerticalAxis.sortby(VerticalCoordIndexer, ascending=True)
    elif VerticalAxis.metpy.units == units('m'):
        VerticalAxis = VerticalAxis.sortby(VerticalCoordIndexer, ascending=False)
    # Integrate
    trapz = 0
    for i in range(1,len(VerticalAxis)):
        dz = VerticalAxis[i]-VerticalAxis[i-1]
        indexer = VerticalAxis[i].metpy.convert_units('hPa')
        fnext = VariableData.sel({VerticalCoordIndexer:indexer})
        indexer = VerticalAxis[i-1].metpy.convert_units('hPa')
        fprev = VariableData.sel({VerticalCoordIndexer:indexer})
        trapz += (fprev+fnext)*(dz/2)
        
    ## TO DO: WORK WITH HEIGHT AS VERTICAL COORDINATE

    return trapz


def CalcZonalAverage(VariableData,LonIndexer):
    """
    Computates zonal averages of a variable for all z levels and times 
    
    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated
    LonIndexer: string
        the indexer used for the longitude variable in the xarray       
   
    Returns
    -------
    zonal_ave: xarray.Dataset
        Arrays of zonal avreages for all longitudes from the passed Dataset
    """
    # Get latitude and logitude data
    lons = VariableData[LonIndexer]
    rlons = units('radians')*np.deg2rad(lons)
    # Integrate through longitude
    trapz = HorizontalTrazpezoidalIntegration(VariableData,LonIndexer)
    # Take the zonal average
    zonal_ave = trapz/(rlons[-1]-rlons[0])
    return zonal_ave

def CalcAreaAverage(VariableData,LatIndexer,LonIndexer=None):
    """
    Computates the Area Average of a function.
    
    The default is to computate the zonal average and then a meridional average.
    If the input data is already some sort of zonal quantity (average or not),
    simply set LonIndexer to None
    
    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated
    LatIndexer: string 
        the indexer used for the latitude variable in the xarray 
    LonIndexer: string (optional)
        the indexer used for the longitude variable in the xarray 
   
    Returns
    -------
    zonal_ave: xarray.Dataset
        Arrays of area avreages for all latitudes and longitudes from
        the passed Dataset
    """

    # Get latitude and logitude data
    lats = VariableData[LatIndexer]
    rlats = np.deg2rad(lats)
    rlats = units('radians')*rlats
    if LonIndexer:
        # If LonIndexer is provided, get zonal ave
        zonal_ave = CalcZonalAverage(VariableData,LonIndexer)
    else:
        zonal_ave = VariableData
    # Take the area avearge
    cosines = np.cos(rlats)
    trapz = HorizontalTrazpezoidalIntegration(zonal_ave*cosines,LatIndexer)
    length = ((np.sin(rlats[0])-np.sin(rlats[-1])))
    area_ave = trapz/length
    area_ave  = area_ave/units('radians')
    return area_ave

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 10:05:52 2022

Main script for the Lorenc Enery Cycle Program.
It will take user inputs to open the data and specify the bounding box for the
calculations and there are some functions to certify that the user input is
valid.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com



"""
from EnergyContents import EnergyContents
from ConversionTerms import ConversionTerms
from BoundaryTerms import BoundaryTerms
from GenerationDissipationTerms import GenerationDissipationTerms
from BoxData import BoxData
from metpy.units import units
import pandas as pd
import xarray as xr
import dataclasses
import sys
from typing import List
import os
import numpy as np

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon


# Arguments passed by user
file3D = sys.argv[1]
file2D = sys.argv[2]
varlist = sys.argv[3]
min_lon = float(sys.argv[4])
max_lon = float(sys.argv[5])
min_lat = float(sys.argv[6])
max_lat = float(sys.argv[7])
output = sys.argv[8]

print(sys.argv)

# This will be printed as an error message
USAGE = f"Usage: python {sys.argv[0]} [--help] | file file2D fvar\
 min_lon, max_lon, min_lat, max_lat output\
     \n\
     \n    Arguments:\
     \n    file3D: file containing 3D fileds\
     \n    file2D: file containing 2D fields (wind stress)\
     \n    fvar: csv file containing variable and dimension indexers\
     \n    min_lon: westernmost limit of the box\
     \n    max_lon: easternmost limit of the box \
     \n    min_lat: southernmost limit of the box \
     \n    max_lat: northernmost limit of the box\
     \n    output: name that will be used as prefix for saving results"

# Object with the inputs given by the user
@dataclasses.dataclass
class Arguments:
    file3D: str
    file2D: str
    fvar: str
    min_lon: float
    max_lon: float
    min_lat: float
    max_lat: float
    output: str

# This fucntion will print the inputs type and what was expected
# It calls the 'validate' function that do the actual variable check
def check_type(obj):
    for field in dataclasses.fields(obj):
        value = getattr(obj, field.name)
        print(
            f"Value: {value}, "
            f"Expected type {field.type} for {field.name}, "
            f"got {type(value)}"
        )
        # checks if provided type matches the expected type
        if type(value) != field.type:
            print("Type Error")
        else:
            print("Type Ok")

# Checks if the inputs matches the expected type and if fails, prints the 
# usage message to the user
def validate(args: List[str]):
    if len(args) > 2:
        for i in range(3,len(args)-1):
                args[i] = float(args[i])
    # Attempts to construct the object with the user inputs
    try:
        arguments = Arguments(*args)
    except TypeError:
        raise SystemExit(USAGE)
    check_type(arguments)

# Checks if the user actually inputed arguments in the first place
def tests_args() -> None:
    args = sys.argv[1:]
    # if no inputs were provided, exits and calls the usage message
    if not args:
        raise SystemExit(USAGE)
    if args[0] == "--help":
        print(USAGE)
    # if the inputs were provided, see if they have the expected type
    else:
        validate(args)
        
def check_create_folder(DirName):
    if not os.path.exists(DirName):
                os.makedirs(DirName)
                print(DirName+' created')
    else:
        print(DirName+' directory exists')
 
        
# Convert longitudes from 0:360 range to -180:180
def convert_lon(df,LonIndexer):
    df.coords[LonIndexer] = (df.coords[LonIndexer] + 180) % 360 - 180
    df = df.sortby(df[LonIndexer])
    return df

# Plot the area limited by the lons and lats values that will be used
# for the computations
def plot_area(min_lon, max_lon,min_lat, max_lat, outdir) :

    plt.close('all')
    datacrs = ccrs.PlateCarree() # projection
    fig = plt.figure(figsize=(8, 8.5))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection=datacrs,
                  frameon=True)
    ax.set_extent([min_lon-20, max_lon+20, max_lat+20, min_lat-20], crs=datacrs)
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
    plt.savefig(outdir+'/Figures/box.png')
    print('\nCreated figure with box defined for computations')

# Function for opening the data
def get_data(file3D, file2D, varlist, min_lon, max_lon, min_lat, max_lat):   
    print('Variables specified by the user in: '+varlist)
    print('Attempting to read '+varlist+' file...')
    try:
        dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)
    except:
        raise SystemExit('ERROR!!!!!')
    # Get data indexers
    LonIndexer = dfVars.loc['Longitude']['Variable']
    LatIndexer = dfVars.loc['Latitude']['Variable']
    TimeIndexer = dfVars.loc['Time']['Variable']
    LevelIndexer = dfVars.loc['Vertical Level']['Variable']
    print('Ok!')
    # Check if the merged file already exists
    DataDir = os.path.dirname(os.path.realpath(file3D))
    merged_file = DataDir+'/'+output+'.nc'
    if not os.path.exists(DataDir+'/'+output+'.nc'):
        try:
            print('Attempting to merge '+file3D+' and '+file2D+' into '
                  + merged_file)
            os.system("python merge_2d_into_3d.py "
                       +file3D+' '+file2D+' '+output)
        except:
            raise('Could not create '+output+' file')
    else:
        print(merged_file+' already exists!')
    print('Attempting to open '+merged_file)
    try:
        full_data = convert_lon(xr.open_dataset(merged_file),LonIndexer)
    except:
        raise SystemExit('ERROR!!!!!\n COULD NOT OPEN OR MERGE DATA')
    print('Ok!')
    # Sort data coordinates - data from distinc sources might have different
    # arrangements, which could affect the results from the integrations
    full_data = full_data.sortby(LonIndexer).sortby(LevelIndexer,ascending=False).sortby(LatIndexer,ascending=False)
    # Fill missing values with 0
    full_data = full_data.fillna(0)
    try:
        full_data = full_data.where(full_data.apply(np.isfinite)).fillna(0.0)
    except:
        full_data = full_data.fillna(0)
    # load data into memory (code optmization)
    data = full_data.load()
    # Stores data as separated variables
    tair = data[dfVars.loc['Air Temperature']['Variable']] \
        * units(dfVars.loc['Air Temperature']['Units']).to('K')
    hgt = data[dfVars.loc['Geopotential Height']['Variable']]\
        *units(dfVars.loc['Geopotential Height']['Units']).to('gpm')
    rhum = data[dfVars.loc['Relative Humidity']['Variable']]*\
        units(dfVars.loc['Relative Humidity']['Units']).to('percent')
    omega = data[dfVars.loc['Omega Velocity']['Variable']]*\
        units(dfVars.loc['Omega Velocity']['Units']).to('Pa/s')
    u = data[dfVars.loc['Eastward Wind Component']['Variable']]*\
        units(dfVars.loc['Eastward Wind Component']['Units'])
    v = data[dfVars.loc['Northward Wind Component']['Variable']]*\
        units(dfVars.loc['Northward Wind Component']['Units'])
    u_stress = data[dfVars.loc['Zonal Wind Stress']['Variable']]*\
        units(dfVars.loc['Zonal Wind Stress']['Units'])
    v_stress = data[dfVars.loc['Meridional Wind Stress']['Variable']]*\
        units(dfVars.loc['Meridional Wind Stress']['Units'])
    # slp = data[dfVars.loc['Sea Level Pressure']['Variable']]
    # Print variables for the user
    print('List of variables found:')
    print(dfVars)
    return LonIndexer, LatIndexer, TimeIndexer, LevelIndexer, tair, hgt,\
        rhum, omega, u, v, u_stress, v_stress

# Compute the budget equation for the energy terms (Az, Ae, Kz and Ke) using
# finite differences method (used for estimating generation, disspation and
# boundary work terms as residuals)
def calc_budget_diff(df,time):
    # get time delta in seconds
    dt = float((time[1]-time[0]).astype('timedelta64[h]'
                                    ) / np.timedelta64(1, 's'))
    # Estimate budget values for all energy terms
    for term in ['Az','Ae','Kz','Ke']:
        name = '∂'+term+'/∂t (finite diff.)'
        print('\nEstimating '+name)
        # forward finite difference for the first value
        forward = (df[term].iloc[1]-df[term].iloc[0])/dt
        # central finite differentes for the second and the penultimate values
        central_second = (df[term].iloc[2]-df[term].iloc[0])/dt
        central_penultimate = (df[term].iloc[-1]-df[term].iloc[-3])/dt
        # fourth order for the third to antepenultimate value
        fourthorder1 = (4/3)*(
            df[term].iloc[3:-1].values-df[term].iloc[1:-3].values)/(2*dt)
        fourthorder2 = (1/3)*(
            df[term].iloc[4:].values-df[term].iloc[1:-3].values)/(4*dt)
        fourthorder = fourthorder1-fourthorder2
        # backward finite difference for the last value
        backward = (df[term].iloc[-1]-df[term].iloc[-2])/dt
        # put all values together
        df[name] = [forward,central_second,*fourthorder,
                    central_penultimate,backward]
        print(df[name].values*units('W/ m **2'))
    return df

# Compute the residuals RGz, RKz, RGe and RKe using the budget terms estimated
# via finite differences
def calc_residuals(df):
    print('\nResiduals ('+str((1*units('W/ m **2')).units)+'):')
    df['RGz'] = df['∂Az/∂t (finite diff.)'] + df['Cz'] + df['Ca'] - df['BAz']
    df['RGe'] = df['∂Ae/∂t (finite diff.)'] - df['Ca'] + df['Ce'] - df['BAe']
    df['RKz'] = df['∂Kz/∂t (finite diff.)'] - df['Cz'] - df['Ck'] - df['BKz']
    df['RKe'] = df['∂Ke/∂t (finite diff.)'] - df['Ce'] + df['Ck'] - df['BKe']
    print(df[['RGz', 'RKz', 'RGe', 'RKe']])
    return df

# The main function. It will open the data, read the variables and calls the
# functions for making the calculations 
def main():
    
    # 1) Checks if the inputs are correct.
    tests_args()
    print('')
    
    # 2) Open the data
    data = get_data(file3D, file2D, varlist, min_lon, max_lon, min_lat, max_lat)   
    # Data indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = data[0],data[1],data[2], data[3]
    # Data variables
    tair = data[4]#*units(data[4].units).to('K')
    hgt = data[5]#*units(data[5].units).to('gpm')
    rhum = data[6]#*units(data[6].units)
    omega = data[7]#*units(data[7].units).to('Pa/s')
    u = data[8]#*units(data[8].units).to('m/s')
    v = data[9]#*units(data[9].units).to('m/s')
    u_stress = data[10]
    v_stress = data[11]
    pres = tair[VerticalCoordIndexer]*units(tair[VerticalCoordIndexer].units).to('Pa')
    #
    print('\n Parameters spcified for the bounding box:')
    print('min_lon, max_lon, min_lat, max_lat: '+str([min_lon,
                                                     max_lon, min_lat, max_lat]))
    
    # 3) Create folder to save results
    # Convert box limits to strings for apprending to file name
    lims = ''
    for i in [min_lon, max_lon]:
        i = str(int(i))
        if i[0] == '-':
            j = i.replace('-','')+'W'
        else:
            j = i+'E'
        lims += j
    for i in [min_lat, max_lat]:
        i = str(int(i))
        if i[0] == '-':
            j = i.replace('-','')+'S'
        else:
            j = i+'N'
        lims += j
    # Directory where results will be stored
    ResultsMainDirectory = '../LEC_Results'
    # Append data limits to outfile name
    outfile_name = output+'_'+lims
    # Each dataset of results have its own directory, allowing to store results
    # from more than one experiment at each time
    ResultsSubDirectory = ResultsMainDirectory+'/'+outfile_name
    # Check if the LEC_Figures directory exists. If not, creates it
    check_create_folder(ResultsMainDirectory)
    # Check if a directory for current data exists. If not, creates it
    check_create_folder(ResultsSubDirectory)
                
    # 4) 
    print('Computing zonal and area averages and eddy terms for each variable')
    print('and the static stability parameter...')
    try:
        box_obj = BoxData(TemperatureData=tair, PressureData=pres,
                     UWindComponentData=u, VWindComponentData=v,
                     ZonalWindStressData=u_stress,
                     MeridionalWindStressData=v_stress,
                     OmegaData=omega, HgtData=hgt,
                     LonIndexer=LonIndexer, LatIndexer=LatIndexer,
                     TimeName=TimeName,
                     VerticalCoordIndexer=VerticalCoordIndexer,
                     western_limit=min_lon, eastern_limit=max_lon,
                     southern_limit=min_lat, northern_limit=max_lat,
                     output_dir=ResultsSubDirectory)
    except:
        raise SystemExit('ERROR!!!!!')
    print('Ok!')
    
    # 5) 
    print('\n------------------------------------------------------------------------')
    print('Computing zonal and eddy kinectic and available potential energy terms')
    try:
        ec_obj = EnergyContents(box_obj)
        EnergyList = [ec_obj.calc_az(), ec_obj.calc_ae(),
                      ec_obj.calc_kz(),ec_obj.calc_ke()]
    except:
        raise SystemExit('ERROR!!!!!')
    print('Ok!')
    
    # 6)
    print('\n------------------------------------------------------------------------')
    print('Computing the conversion terms between energy contents') 
    try:
        ct_obj = ConversionTerms(box_obj)
        ConversionList = [ct_obj.calc_cz(),ct_obj.calc_ca(),
                          ct_obj.calc_ck(),ct_obj.calc_ce()]
    except:
        raise SystemExit('ERROR!!!!!')
    print('Ok!')
    
    # 7)
    print('\n------------------------------------------------------------------------')
    print('Computing the boundary terms') 
    try:
        bt_obj = BoundaryTerms(box_obj)
        BoundaryList = [bt_obj.calc_baz(),bt_obj.calc_bae(),
                        bt_obj.calc_bkz(),bt_obj.calc_bke(),
                        bt_obj.calc_boz(),bt_obj.calc_boe()]
    except:
        raise SystemExit('ERROR!!!!!')
    print('Ok!')
    
    # 8)
    print('\n------------------------------------------------------------------------')
    print('Computing generation and disspiation terms') 
    try:
        gdt_obj = GenerationDissipationTerms(box_obj)
        GenDissList = [gdt_obj.calc_gz(),gdt_obj.calc_ge(),
                       gdt_obj.calc_dz(),gdt_obj.calc_de()]
    except:
        raise SystemExit('ERROR!!!!!')
    print('Ok!')
    
    # 9)
    print('\nOrganising results in a Pandas DataFrame')
    # First, extract dates to construct a dataframe
    dates = tair[TimeName].values
    days = dates.astype('datetime64[D]')
    hours = pd.to_datetime(dates).hour
    df = pd.DataFrame(data=[*days],columns=['Date'])
    df['Hour'] = hours
    # Then adds the data to the DataFrame
    for i,j,k,l in zip(range(4),['Az','Ae','Kz','Ke'],
                     ['Cz','Ca','Ck','Ce'],
                     ['Gz','Ge','Dz','De']):
        df[j] = EnergyList[i]
        df[k] = ConversionList[i]
        df[l] = GenDissList[i]
    for i,l in zip(range(6),['BAz','BAe','BKz','BKe','BΦZ','BΦE']):
        df[l] = BoundaryList[i]       
    # 10) 
    print('\n------------------------------------------------------------------------')
    print('Estimating budget terms (∂X/∂t) using finite differences ')
    df = calc_budget_diff(df,dates) 
    print('Ok!')
    
    # 11) 
    print('\n------------------------------------------------------------------------')
    print('Computing residuals RGz, RKz, RGe and RKe')
    df = calc_residuals(df)
    print('Ok!')
    
    # 12) save file
    print('\nCreating a csv to store results...')
    outfile = ResultsSubDirectory+'/'+outfile_name+'.csv'
    df.to_csv(outfile)
    print(outfile+' created')
    print('All done!')
    
    # 13) Make figures
    FigsDirectory = ResultsSubDirectory+'/Figures'
    check_create_folder(FigsDirectory)
    os.system("python plot_timeseries.py "+outfile)
    os.system("python plot_vertical.py "+ResultsSubDirectory)
    os.system("python plot_boxplot.py "+ResultsSubDirectory)
    plot_area(min_lon, max_lon,min_lat,max_lat, ResultsSubDirectory)

if __name__ == "__main__":
    main()
    
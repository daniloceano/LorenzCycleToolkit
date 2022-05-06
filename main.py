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
from BoxData import BoxData
from metpy.units import units
import pandas as pd
import xarray as xr
import dataclasses
import sys
from typing import List
import os
import numpy as np


# This will be printed as an error message
USAGE = f"Usage: python {sys.argv[0]} [--help] | file fvar min_lon, max_lon, min_lat, max_lat]"

# Object with the inputs given by the user
@dataclasses.dataclass
class Arguments:
    file: str
    fvar: str
    min_lon: float
    max_lon: float
    min_lat: float
    max_lat: float

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
        for i in range(2,len(args)):
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
        
# Convert longitudes from 0:360 range to -180:180
def convert_lon(df,LonIndexer):
    df.coords[LonIndexer] = (df.coords[LonIndexer] + 180) % 360 - 180
    df = df.sortby(df[LonIndexer])
    return df

# Function for opening the data
def get_data(file,varlist,min_lon, max_lon, min_lat, max_lat):   
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
    # Open actual data
    print('Attempting to open '+file+' file...')
    try:
        full_data = convert_lon(xr.open_dataset(file),LonIndexer)
    except:
        raise SystemExit('ERROR!!!!!')
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
    # slp = data[dfVars.loc['Sea Level Pressure']['Variable']]
    # Print variables for the user
    print('List of variables found:')
    print(dfVars)
    return LonIndexer, LatIndexer, TimeIndexer, LevelIndexer, tair, hgt, rhum, omega, u, v#, slp


# The main function. It will open the data, read the variables and calls the
# functions for making the calculations 
def main():
    ##################
    # 1) Checks if the inputs are correct.
    tests_args()
    print('')
    ################## 
    # 2) Open the data
    data = get_data(*sys.argv[1:])   
    # Data indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = data[0],data[1],data[2], data[3]
    # Data variables
    tair = data[4]#*units(data[4].units).to('K')
    hgt = data[5]#*units(data[5].units).to('gpm')
    rhum = data[6]#*units(data[6].units)
    omega = data[7]#*units(data[7].units).to('Pa/s')
    u = data[8]#*units(data[8].units).to('m/s')
    v = data[9]#*units(data[9].units).to('m/s')
    # slp = data[10]*units(data[10].units).to('Pa')
    pres = tair[VerticalCoordIndexer]*units(tair[VerticalCoordIndexer].units).to('Pa')
    #
    print('\n Parameters spcified for the bounding box:')
    print('min_lon, max_lon, min_lat, max_lat: '+str(sys.argv[3:]))
    ##################
    # 3) Create folder to save results
    # Convert box limits to strings for apprending to file name
    lims = ''
    for i in sys.argv[3:5]:
        if i[0] == '-':
            j = i.replace('-','')+'W'
        else:
            j = i+'E'
        lims += j
    for i in sys.argv[5:]:
        if i[0] == '-':
            j = i.replace('-','')+'S'
        else:
            j = i+'N'
        lims += j
    # Directory where results will be stored
    ResultsDirectory = 'LEC_Results'
    infile_name = sys.argv[1]
    # Append data limits to outfile name
    outfile_name = infile_name.split('.nc')[0]+'_'+lims
    # Each dataset of results have its own directory, allowing to store results
    # from more than one experiment at each time
    DataDirectory = ResultsDirectory+'/'+outfile_name
    # Check if the LEC_Figures directory exists. If not, creates it
    if not os.path.exists(ResultsDirectory):
                os.makedirs(ResultsDirectory)
                print(ResultsDirectory+' created')
    else:
        print(ResultsDirectory+' directory exists')
    # Check if a directory for current data exists. If not, creates it
    if not os.path.exists(DataDirectory):
                os.makedirs(DataDirectory)
                print(DataDirectory+' created')
    ##################        
    # 4) 
    print('Computing zonal and area averages and eddy terms for each variable')
    print('and the static stability parameter...')
    try:
        box_obj = BoxData(TemperatureData=tair, PressureData=pres,
                     UWindComponentData=u, VWindComponentData=v,
                     OmegaData=omega,
                     LonIndexer=LonIndexer, LatIndexer=LatIndexer, TimeName=TimeName,
                     VerticalCoordIndexer=VerticalCoordIndexer,
                     western_limit=float(sys.argv[3]), eastern_limit=float(sys.argv[4]),
                     southern_limit=float(sys.argv[5]), northern_limit=float(sys.argv[6]),
                     output_dir=DataDirectory)
    except:
        raise SystemExit('ERROR!!!!!')
    print('Ok!')
    ##################
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
    ##################
    # 6)
    print('\n------------------------------------------------------------------------')
    print('Computing the conversion terms between energy contents') 
    try:
        ct_obj = ConversionTerms(box_obj)
        ConversionList = [ct_obj.calc_cz(),ct_obj.calc_ca(),ct_obj.calc_ck(),ct_obj.calc_ce()]
    except:
        raise SystemExit('ERROR!!!!!')
    print('Ok!')
    ##################
    # 7)
    # First, extract dates to construct the dataframe
    print('\nCreating a csv to store results...')
    dates = tair[TimeName].values
    days = dates.astype('datetime64[D]')
    hours = pd.to_datetime(dates).hour
    df = pd.DataFrame(data=[*days],columns=['Date'])
    df['Hour'] = hours
    # Then adds the data to the DataFrame
    for i,j,k in zip(range(4),['Az','Ae','Kz','Ke'],['Cz','Ca','Ck','Ce']):
        df[j] = EnergyList[i]
        df[k] = ConversionList[i]
    # Lastly, save file
    outfile = DataDirectory+'/'+outfile_name+'.csv'
    df.to_csv(outfile)
    print(outfile+' created')
    print('All done!')

if __name__ == "__main__":
    main()
    
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
from compute_terms import calc_budget_diff,calc_residuals
from metpy.units import units
import pandas as pd
import xarray as xr
import os
import numpy as np
import argparse


        
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


# Function for opening the data
def get_data(infile, varlist, min_lon, max_lon, min_lat, max_lat):   
    print('Variables specified by the user in: '+varlist)
    print('Attempting to read '+varlist+' file...')
    try:
        dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)
    except:
        raise SystemExit("Error: verify that there is a 'fvar' text file\
 located on your working directory")
    # Print variables for the user
    print('List of variables found:')
    print(dfVars)
    # Get data indexers
    LonIndexer = dfVars.loc['Longitude']['Variable']
    LatIndexer = dfVars.loc['Latitude']['Variable']
    TimeIndexer = dfVars.loc['Time']['Variable']
    LevelIndexer = dfVars.loc['Vertical Level']['Variable']
    print('Ok!')
    print('Opening inpyt data...')
    try:
        full_data = convert_lon(xr.open_dataset(infile),LonIndexer)
    except:
        raise SystemExit('ERROR!!!!!\n COULD NOT OPEN INPUT DATA')
    print('Ok!')
    # Sort data coordinates - data from distinc sources might have different
    # arrangements, which could affect the results from the integrations
    full_data = full_data.sortby(LonIndexer).sortby(LevelIndexer,
                            ascending=False).sortby(LatIndexer,ascending=False)
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
    omega = data[dfVars.loc['Omega Velocity']['Variable']]*\
        units(dfVars.loc['Omega Velocity']['Units']).to('Pa/s')
    u = data[dfVars.loc['Eastward Wind Component']['Variable']]*\
        units(dfVars.loc['Eastward Wind Component']['Units'])
    v = data[dfVars.loc['Northward Wind Component']['Variable']]*\
        units(dfVars.loc['Northward Wind Component']['Units'])
    if not args.residuals:
        u_stress = data[dfVars.loc['Zonal Wind Stress']['Variable']]*\
            units(dfVars.loc['Zonal Wind Stress']['Units'])
        v_stress = data[dfVars.loc['Meridional Wind Stress']['Variable']]*\
            units(dfVars.loc['Meridional Wind Stress']['Units'])
        return LonIndexer, LatIndexer, TimeIndexer, LevelIndexer, tair, hgt,\
            omega, u, v, u_stress, v_stress
    else:
        return LonIndexer, LatIndexer, TimeIndexer, LevelIndexer, tair, hgt,\
            omega, u, v



# The main function. It will open the data, read the variables and calls the
# functions for making the calculations 
def main():
    
    print('')
    
    # 2) Open the data
    data = get_data(args.infile, varlist, min_lon, max_lon, min_lat, max_lat)   
    # Data indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = data[0],data[1],data[2], data[3]
    # Data variables
    tair = data[4]
    hgt = data[5]
    omega = data[6]
    u = data[7]
    v = data[8]
    if not args.residuals:
        u_stress = data[9]
        v_stress = data[10]
    else:
        u_stress = v*np.nan
        v_stress = v*np.nan
    pres = tair[VerticalCoordIndexer]*units(
        tair[VerticalCoordIndexer].units).to('Pa')
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
    outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_'+lims
    # Each dataset of results have its own directory, allowing to store results
    # from more than one experiment at each time
    ResultsSubDirectory = ResultsMainDirectory+'/'+outfile_name+'/'
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
        raise SystemExit('Error on creating the box for the computations')
    print('Ok!')
    
    # 5) 
    print('\n------------------------------------------------------------------------')
    print('Computing zonal and eddy kinectic and available potential energy terms')
    try:
        ec_obj = EnergyContents(box_obj)
        EnergyList = [ec_obj.calc_az(), ec_obj.calc_ae(),
                      ec_obj.calc_kz(),ec_obj.calc_ke()]
    except:
        raise SystemExit('Error on computing Energy Contents')
    print('Ok!')
    
    # 6)
    print('\n------------------------------------------------------------------------')
    print('Computing the conversion terms between energy contents') 
    try:
        ct_obj = ConversionTerms(box_obj)
        ConversionList = [ct_obj.calc_cz(),ct_obj.calc_ca(),
                          ct_obj.calc_ck(),ct_obj.calc_ce()]
    except:
        raise SystemExit('Error on computing Conversion Terms')
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
        raise SystemExit('Error on computing Boundary Terms')
    print('Ok!')
    
    # 8)
    print('\n------------------------------------------------------------------------')
    print('Computing generation and disspiation terms') 
    try:
        gdt_obj = GenerationDissipationTerms(box_obj)
        if args.residuals:
            GenDissList = [gdt_obj.calc_gz(),gdt_obj.calc_ge()]
        else:
             GenDissList = [gdt_obj.calc_gz(),gdt_obj.calc_ge(),
                            gdt_obj.calc_dz(),gdt_obj.calc_de()]
    except:
        raise SystemExit('Error on computing generation/Dissipation Terms')
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
    for i,j in zip(range(4),['Az','Ae','Kz','Ke']):
        df[j] = EnergyList[i]
    for i,k in zip(range(4),['Cz','Ca','Ck','Ce']):
        df[k] = ConversionList[i]
    if args.residuals:
        for i,l in zip(range(4),['BAz','BAe','BKz','BKe']):
            df[l] = BoundaryList[i]
        for i,m in zip(range(4),['Gz','Ge']):
            df[m] = GenDissList[i]  
    else:
        for i,l in  zip(range(4),['Gz','Ge','Dz','De']):
            df[l] = GenDissList[i]
        for i,m in zip(range(6),['BAz','BAe','BKz','BKe','BΦZ','BΦE']):
            df[m] = BoundaryList[i]  
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
    if args.residuals:
        flag = ' -r'
    else:
        flag = ' '
    os.system("python plot_timeseries.py "+outfile+flag)
    os.system("python plot_vertical.py "+ResultsSubDirectory)
    os.system("python plot_boxplot.py "+ResultsSubDirectory+flag)
    os.system("python draw_cycle.py "+outfile)
    os.system("python LorenzPhaseSpace.py "+outfile)
    cmd = "python plot_area.py {0} {1} {2} {3} {4}".format(min_lon, max_lon,min_lat,max_lat, ResultsSubDirectory)
    os.system(cmd)

if __name__ == "__main__":
    
    # Box limits used for compuations
    dfbox = pd.read_csv('./box_limits',header=None,delimiter=';',index_col=0)
    
    # Arguments passed by user
    # infile = sys.argv[1]
    varlist = './fvars'
    min_lon = float(dfbox.loc['min_lon'].values)
    max_lon = float(dfbox.loc['max_lon'].values)
    min_lat = float(dfbox.loc['min_lat'].values)
    max_lat = float(dfbox.loc['max_lat'].values)
    
    # Warns user if limits are wrong
    if min_lon > max_lon:
        raise ValueError('Error in box_limits: min_lon > max_lon')
        quit()
    if min_lat > max_lat:
        raise ValueError('Error in box_limits: min_lat > max_lat')
        quit()
    
    parser = argparse.ArgumentParser(description = "\
Lorenz Energy Cycle program. \n \
It takes an input NetCDF file and crop the data for the bounding box specified\
 at the 'box_lims' file. An auxilliary 'fvars' file is also needed, where it is\
 specified the names used for each variable. It computates the energy, \
 conversion, boundary, generation and dissipation terms. The results are stored\
 as csv files in the 'LEC_results' directory on ../ and it also creates figures\
 for visualising the results.")
    parser.add_argument("infile", help = "Input .nc file with temperature,\
 geopotential and meridional, zonal and vertical components of the wind,\
 in pressure levels")
    parser.add_argument("-r", "--residuals", default = False, action='store_true',
    help = "Flag for computing the Dissipation and Generation terms as\
  residuals (needs to provide the 3D friction terms in the infile)")
    args = parser.parse_args()
    infile = args.infile
        
    main()
    
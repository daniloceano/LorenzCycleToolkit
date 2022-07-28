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
from thermodynamics import AdiabaticHEating
from compute_terms import calc_budget_diff,calc_residuals
from metpy.units import units
import pandas as pd
import xarray as xr
import os
import numpy as np
import argparse
from metpy.constants import g

import traceback
import logging

import time



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
def get_data(infile, varlist):   
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
    LonIndexer,LatIndexer,TimeIndexer,LevelIndexer = \
      dfVars.loc['Longitude']['Variable'],dfVars.loc['Latitude']['Variable'],\
      dfVars.loc['Time']['Variable'],dfVars.loc['Vertical Level']['Variable']
    print('Ok!')
    print('Opening inpyt data...')
    try:
        full_data = convert_lon(xr.open_dataset(infile),LonIndexer)
    except:
        raise SystemExit('ERROR!!!!!\n Could not open data. Check if path is\
 correct, fvars file and file format (should be .nc)')
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
    # Stores data as separated variables and give them correct units
    tair = data[dfVars.loc['Air Temperature']['Variable']] \
        * units(dfVars.loc['Air Temperature']['Units']).to('K')
    omega = data[dfVars.loc['Omega Velocity']['Variable']]*\
        units(dfVars.loc['Omega Velocity']['Units']).to('Pa/s')
    u = data[dfVars.loc['Eastward Wind Component']['Variable']]*\
        units(dfVars.loc['Eastward Wind Component']['Units'])
    v = data[dfVars.loc['Northward Wind Component']['Variable']]*\
        units(dfVars.loc['Northward Wind Component']['Units'])
    if args.geopotential:
         hgt = (data[dfVars.loc['Geopotential']['Variable']] \
        * units(dfVars.loc['Geopotential']['Units'])/g).metpy.convert_units('gpm')
    else:
        hgt = data[dfVars.loc['Geopotential Height']['Variable']]\
            *units(dfVars.loc['Geopotential Height']['Units']).to('gpm')
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

#---------------------------------------------------------------------------
# Computes the Lorenz Energy Cycle using an eulerian framework.
# It requires the box_lims file with the limits for the box used for the
# computations, which is fixed in time. IN this framework we can analyse how
# the eddies contribute for the local energy cycle.
def LEC_eulerian():
    print('Computing energetics using eulerian framework')
    # Box limits used for compuations
    dfbox = pd.read_csv('./box_limits',header=None,delimiter=';',index_col=0)
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
    # 2) Open the data
    data = get_data(infile, varlist)   
    # Data indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = \
        data[0], data[1], data[2], data[3]
    # Data variables
    tair, hgt, omega, u, v = data[4], data[5], data[6], data[7], data[8]
    pres = tair[VerticalCoordIndexer]*units(
        tair[VerticalCoordIndexer].units).to('Pa')
    Q = AdiabaticHEating(tair,pres,omega,u,v,VerticalCoordIndexer,
                                LatIndexer,LonIndexer,TimeName)
    if not args.residuals:
        u_stress = data[9]
        v_stress = data[10]
    else:
        u_stress = v*np.nan
        v_stress = v*np.nan
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
    # Create csv files for storing vertical results.
    # When those terms are computed data is simply appended to csv
    for term in ['Az','Ae','Kz','Ke','Cz','Ca','Ck','Ce','Ge','Gz']:
        tmp = pd.DataFrame(columns=[TimeName,
        *[float(i)/100 for i in pres.values]])
        tmp.to_csv(ResultsSubDirectory+term+'_'+VerticalCoordIndexer+'.csv',
                   index=None)
    # 4) 
    print('Computing zonal and area averages and eddy terms for each variable')
    print('and the static stability parameter...')
    try:
        box_obj = BoxData(TemperatureData=tair, PressureData=pres,
        UWindComponentData=u, VWindComponentData=v,
        ZonalWindStressData=u_stress,MeridionalWindStressData=v_stress,
        OmegaData=omega, HgtData=hgt,LonIndexer=LonIndexer, 
        AdiabaticHeatingData=Q,
        LatIndexer=LatIndexer, VerticalCoordIndexer=VerticalCoordIndexer,
        TimeName=TimeName,western_limit=min_lon, eastern_limit=max_lon,
        southern_limit=min_lat, northern_limit=max_lat,
        output_dir=ResultsSubDirectory)
    except:
        raise SystemExit('Error on creating the box for the computations')
    print('Ok!')
    # 5) 
    print('\n------------------------------------------------------------------------')
    print('Computing zonal and eddy kinectic and available potential energy terms')
    try:
        ec_obj = EnergyContents(box_obj,method='eulerian')
        EnergyList = [ec_obj.calc_az(), ec_obj.calc_ae(),
                      ec_obj.calc_kz(),ec_obj.calc_ke()]
    except:
        raise SystemExit('Error on computing Energy Contents')
    print('Ok!')
    # 6)
    print('\n------------------------------------------------------------------------')
    print('Computing the conversion terms between energy contents') 
    try:
        ct_obj = ConversionTerms(box_obj,method='eulerian')
        ConversionList = [ct_obj.calc_cz(),ct_obj.calc_ca(),
                          ct_obj.calc_ck(),ct_obj.calc_ce()]
    except:
        raise SystemExit('Error on computing Conversion Terms')
    print('Ok!')
    # 7)
    print('\n------------------------------------------------------------------------')
    print('Computing the boundary terms') 
    try:
        bt_obj = BoundaryTerms(box_obj,method='eulerian')
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
        gdt_obj = GenerationDissipationTerms(box_obj,method='eulerian')
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
    os.system("python plots/plot_timeseries.py "+outfile+flag)
    os.system("python plots/plot_vertical.py "+ResultsSubDirectory)
    os.system("python plots/plot_boxplot.py "+ResultsSubDirectory+flag)
    os.system("python plots/plot_LEC.py "+outfile)
    os.system("python plots/plot_LPS.py "+outfile)
    cmd = "python plots/plot_area.py {0} {1} {2} {3} {4}".format(min_lon, max_lon,min_lat,max_lat, ResultsSubDirectory)
    os.system(cmd)
    
#---------------------------------------------------------------------------
# Computes the Lorenz Energy Cycle using an eulerian framework.
# It requires the box_lims file with the limits for the box used for the
# computations, which is fixed in time. IN this framework we can analyse how
# the eddies contribute for the local energy cycle.
def LEC_lagrangian():
    print('Computing energetics using eulerian framework')
    # 1) Get data
    data = get_data(infile, varlist)   
    # Data indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = \
        data[0], data[1], data[2], data[3]
    # Data variables
    tair, hgt, omega, u, v = data[4], data[5], data[6], data[7], data[8]
    pres = tair[VerticalCoordIndexer]*units(
        tair[VerticalCoordIndexer].units).to('Pa')
    Q = AdiabaticHEating(tair,pres,omega,u,v,VerticalCoordIndexer,
                                LatIndexer,LonIndexer,TimeName)
    if not args.residuals:
        u_stress = data[9]
        v_stress = data[10]
    else:
        u_stress = v*np.nan
        v_stress = v*np.nan    
    # Directory where results will be stored
    ResultsMainDirectory = '../LEC_Results'
    # Append data limits to outfile name
    outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_lagranigan'
    # Each dataset of results have its own directory, allowing to store results
    # from more than one experiment at each time
    ResultsSubDirectory = ResultsMainDirectory+'/'+outfile_name+'/'
    # Check if the LEC_Figures directory exists. If not, creates it
    check_create_folder(ResultsMainDirectory)
    # Check if a directory for current data exists. If not, creates it
    check_create_folder(ResultsSubDirectory)
    # Create csv files for storing vertical results.
    # When those terms are computed data is simply appended to csv
    for term in ['Az','Ae','Kz','Ke','Cz','Ca','Ck','Ce','Ge','Gz']:
        tmp = pd.DataFrame(columns=[TimeName,
        *[float(i)/100 for i in pres.values]])
        tmp.to_csv(ResultsSubDirectory+term+'_'+VerticalCoordIndexer+'.csv',
                   index=None)  
    # Track file
    trackfile = './track'
    track = pd.read_csv(trackfile,parse_dates=[0],delimiter=';',index_col='time')
    # Create dict for store results
    TermsDict = {}
    energy = ['Az','Ae','Kz','Ke']
    conversion = ['Cz','Ca','Ck','Ce']
    boundary = ['BAz','BAe','BKz','BKe','BΦZ','BΦE']
    gendiss = ['Gz','Ge']
    if not args.residuals:  
        gendiss = ['Gz','Ge','Dz','De']
    for term in [*energy,*conversion,*boundary,*gendiss]:
        TermsDict[term] = []
    
    # Loop for each time step:
    times = tair[TimeName]
    for t in times:
        itair, ihgt, iomega, iu, iv, iu_stress, iv_stress, iQ  = \
            tair.sel({TimeName:t}), \
                hgt.sel({TimeName:t}), omega.sel({TimeName:t}),\
                    u.sel({TimeName:t}), v.sel({TimeName:t}),\
                        u_stress.sel({TimeName:t}),\
                            v_stress.sel({TimeName:t}),\
                                Q.sel({TimeName:t})              
        
        # Get current time and box limits
        itime = str(t.values)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        min_lon, max_lon = track.loc[itime]['Lon']-7.5,track.loc[itime]['Lon']+7.5
        min_lat, max_lat = track.loc[itime]['Lat']-7.5,track.loc[itime]['Lat']+7.5
        print('\nComputing terms for '+datestr+'...')
        print('Box limits (lon/lat): '+str(max_lon)+'/'+str(max_lat),
              ' '+str(min_lon)+'/'+str(min_lat))
    
        # Create box object
        try:
            box_obj = BoxData(TemperatureData=itair, PressureData=pres,
            UWindComponentData=iu, VWindComponentData=iv,
            ZonalWindStressData=iu_stress,MeridionalWindStressData=iv_stress,
            OmegaData=iomega, HgtData=ihgt,LonIndexer=LonIndexer, 
            AdiabaticHeatingData=iQ,
            LatIndexer=LatIndexer, VerticalCoordIndexer=VerticalCoordIndexer,
            TimeName=TimeName,western_limit=min_lon, eastern_limit=max_lon,
            southern_limit=min_lat, northern_limit=max_lat,
            output_dir=ResultsSubDirectory)
        except:
            raise SystemExit('Error on creating the box for the computations')
        # Compute energy terms
        try:
            ec_obj = EnergyContents(box_obj,method='lagrangian')
            TermsDict['Az'].append(ec_obj.calc_az())
            TermsDict['Ae'].append(ec_obj.calc_ae())
            TermsDict['Kz'].append(ec_obj.calc_kz())
            TermsDict['Ke'].append(ec_obj.calc_ke())
        except Exception as e:
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Energy Contents')
        # Compute conversion terms
        try:
            ct_obj = ConversionTerms(box_obj,method='lagrangian')
            TermsDict['Ca'].append(ct_obj.calc_ca())
            TermsDict['Ce'].append(ct_obj.calc_ce())
            TermsDict['Ck'].append(ct_obj.calc_ck())
            TermsDict['Cz'].append(ct_obj.calc_cz())
        except Exception as e:
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Conversion Terms')
        # Compute boundary terms
        try:
            bt_obj = BoundaryTerms(box_obj,method='lagrangian')
            TermsDict['BAe'].append(bt_obj.calc_bae().values)
            TermsDict['BAz'].append(bt_obj.calc_baz())
            TermsDict['BKe'].append(bt_obj.calc_bke())
            TermsDict['BKz'].append(bt_obj.calc_bkz())
            TermsDict['BΦE'].append(bt_obj.calc_boe())
            TermsDict['BΦZ'].append(bt_obj.calc_boz())
        except Exception as e:
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Boundary Terms')
        # Compute generation/dissipation terms
        try:
            gdt_obj = GenerationDissipationTerms(box_obj,method='lagrangian')
            TermsDict['Ge'].append(gdt_obj.calc_ge())
            TermsDict['Gz'].append(gdt_obj.calc_gz())
            if not args.residuals:
                TermsDict['De'].append(gdt_obj.calc_de())
                TermsDict['Dz'].append(gdt_obj.calc_dz())
        except Exception as e:
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Generation Terms')
    print('\nOrganising results in a Pandas DataFrame')
    df = pd.DataFrame.from_dict(TermsDict, orient ='columns',dtype = float)
    days = times.values.astype('datetime64[D]')
    hours = pd.to_datetime(times.values).hour
    df['Date'],df['Hour'] = days, hours
    
    # Print results for user
    for term in TermsDict.keys():
        print('\n'+term)
        print(df[term].tolist())
    
    print('\n------------------------------------------------------------------------')
    print('Estimating budget terms (∂X/∂t) using finite differences ')
    df = calc_budget_diff(df,times) 
    print('Ok!')
    
    print('\n------------------------------------------------------------------------')
    print('Computing residuals RGz, RKz, RGe and RKe')
    df = calc_residuals(df)
    print('Ok!')
    
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
    os.system("python plots/plot_timeseries.py "+outfile+flag)
    os.system("python plots/plot_vertical.py "+ResultsSubDirectory)
    os.system("python plots/plot_boxplot.py "+ResultsSubDirectory+flag)
    os.system("python plots/plot_LEC.py "+outfile)
    os.system("python plots/plot_LPS.py "+outfile)
    os.system("python plots/plot_track.py "+outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "\
Lorenz Energy Cycle (LEC) program. \n \
The program can compute the LEC using two distinct frameworks:\
    1) Lagragian framework. A box is definid in the box_lims' file and then the \
       energetics are computed for a fixed domain.\
    2) Eulerian framework. The domain is not fixed and follows the system using \
       the track file.\
 Both frameworks can be applied at the same time, given the required files are\
 provided. An auxilliary 'fvars' file is also needed for both frameworks. \
 It contains thespecified names used for each variable. The program computates \
 the energy,conversion, boundary, generation anddissipation terms. The results \
 are stored as csv files in the 'LEC_results' directory on ../ and it also \
 creates figures for visualising the results. \
 The use of -r flag is required while the computation of friction parameters \
 is not implemented")
    parser.add_argument("infile", help = "input .nc file with temperature,\
 geopotential and meridional, zonal and vertical components of the wind,\
 in pressure levels")
    parser.add_argument("-r", "--residuals", default = False, action='store_true',
    help = "compute the Dissipation and Generation terms as\
 residuals (needs to provide the 3D friction terms in the infile)")
    parser.add_argument("-g", "--geopotential", default = False,
    action='store_true', help = "use the geopotential data instead of\
 geopotential height. The file fvars must be adjusted for doing so.")
    parser.add_argument("-e", "--eulerian", default = False,
    action='store_true', help = "compute the energetics for a fixed domain\
 specified by the box_lims file.")
    parser.add_argument("-l", "--lagrangian", default = False,
    action='store_true', help = "compute the energetics for a fixed domain\
 specified by the box_lims file.")
    args = parser.parse_args()
    infile  = args.infile
    # box_limits = args.box_limits
    varlist = './fvars'
    # Run the program
    start_time = time.time()
    if args.eulerian:
        LEC_eulerian()
    print("--- %s seconds running eulerian framework ---" % (time.time() - start_time))
    start_time = time.time()
    if args.lagrangian:
        LEC_lagrangian()
    print("--- %s seconds for running lagrangian framework ---" % (time.time() - start_time))
    
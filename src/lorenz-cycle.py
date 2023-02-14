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
from BudgetResidual import calc_budget_diff,calc_residuals
from thermodynamics import AdiabaticHEating

from metpy.units import units
from metpy.calc import vorticity
from metpy.calc import wind_speed

from select_area import slice_domain
from select_area import draw_box_map

from scipy.signal import savgol_filter    
import pandas as pd
import xarray as xr
import os
import numpy as np
import argparse
import sys
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
    print('Opening input data...')
    try:
        full_data = convert_lon(xr.open_dataset(infile),LonIndexer)
    except:
        raise SystemExit('ERROR!!!!!\n Could not open data. Check if path is\
 correct, fvars file and file format (should be .nc)')
    print('Ok!')
    # load data into memory (code optmization)
    data = full_data.load()
    # Assign lat and lon as radians, for calculations
    data = data.assign_coords({"rlats": np.deg2rad(data[LatIndexer])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data[LatIndexer]))})
    data = data.assign_coords({"rlons": np.deg2rad(data[LonIndexer])})
    # Sort data coordinates as data from distinc sources might have different
    # arrangements, which could affect the results from the integrations
    data = data.sortby(LonIndexer).sortby(LevelIndexer,
                ascending=True).sortby(LatIndexer,ascending=True)
                        
    # Fill missing values with 0
    data = data.fillna(0)
    try:
        data = data.where(data.apply(np.isfinite)).fillna(0.0)
    except:
        data = data.fillna(0)
    
    return data

#---------------------------------------------------------------------------
# Computes the Lorenz Energy Cycle using an fixed framework.
# It requires the box_lims file with the limits for the box used for the
# computations, which is fixed in time. IN this framework we can analyse how
# the eddies contribute for the local energy cycle.
def LEC_fixed(data):
    print('Computing energetics using fixed framework')
    # Box limits used for compuations
    dfbox = pd.read_csv('../inputs/box_limits',header=None,delimiter=';',index_col=0)
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
    dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)
    LonIndexer,LatIndexer,TimeName,VerticalCoordIndexer = \
      dfVars.loc['Longitude']['Variable'],dfVars.loc['Latitude']['Variable'],\
      dfVars.loc['Time']['Variable'],dfVars.loc['Vertical Level']['Variable']
    pres = data[VerticalCoordIndexer]*units(
         data[VerticalCoordIndexer].units).to('Pa')
    
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
        box_obj = BoxData(data=data, dfVars=dfVars, args=args,
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
        ec_obj = EnergyContents(box_obj,method='fixed')
        EnergyList = [ec_obj.calc_az(), ec_obj.calc_ae(),
                      ec_obj.calc_kz(),ec_obj.calc_ke()]
    except:
        raise SystemExit('Error on computing Energy Contents')
    print('Ok!')
    # 6)
    print('\n------------------------------------------------------------------------')
    print('Computing the conversion terms between energy contents') 
    try:
        ct_obj = ConversionTerms(box_obj,method='fixed')
        ConversionList = [ct_obj.calc_cz(),ct_obj.calc_ca(),
                          ct_obj.calc_ck(),ct_obj.calc_ce()]
    except:
        raise SystemExit('Error on computing Conversion Terms')
    print('Ok!')
    # 7)
    print('\n------------------------------------------------------------------------')
    print('Computing the boundary terms') 
    try:
        bt_obj = BoundaryTerms(box_obj,method='fixed')
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
        gdt_obj = GenerationDissipationTerms(box_obj,method='fixed')
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
    dates = data[TimeName].values
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
    os.system("python ../plots/plot_timeseries.py "+outfile+flag)
    os.system("python ../plots/plot_vertical.py "+ResultsSubDirectory)
    os.system("python ../plots/plot_boxplot.py "+ResultsSubDirectory+flag)
    os.system("python ../plots/plot_LEC.py "+outfile)
    os.system("python ../plots/plot_LPS.py "+outfile)
    cmd = "python ../plots/plot_area.py {0} {1} {2} {3} {4}".format(min_lon, max_lon,min_lat,max_lat, ResultsSubDirectory)
    os.system(cmd)
    
#---------------------------------------------------------------------------
# Computes the Lorenz Energy Cycle using an fixed framework.
# It requires the box_lims file with the limits for the box used for the
# computations, which is fixed in time. IN this framework we can analyse how
# the eddies contribute for the local energy cycle.
def LEC_moving(data):
    print('Computing energetics using moving framework')
    # 2) Open the data
    # data = get_data(infile, varlist)  
    # Indexers
    dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)
    LonIndexer,LatIndexer,TimeName,VerticalCoordIndexer = \
      dfVars.loc['Longitude']['Variable'],dfVars.loc['Latitude']['Variable'],\
      dfVars.loc['Time']['Variable'],dfVars.loc['Vertical Level']['Variable']
    pres = data[VerticalCoordIndexer]*units(
         data[VerticalCoordIndexer].units).to('Pa')
    # Get variables for computing Adiabatic Heating
    tair =  (data[dfVars.loc['Air Temperature']['Variable']] * 
             units(dfVars.loc['Air Temperature']['Units']).to('K'))
    omega = (data[dfVars.loc['Omega Velocity']['Variable']] *
             units(dfVars.loc['Omega Velocity']['Units']).to('Pa/s'))
    u = (data[dfVars.loc['Eastward Wind Component']['Variable']] *
         units(dfVars.loc['Eastward Wind Component']['Units']).to('m/s'))
    v = (data[dfVars.loc['Northward Wind Component']['Variable']] *
         units(dfVars.loc['Northward Wind Component']['Units']).to('m/s'))
    hgt = (data[dfVars.loc['Geopotential Height']['Variable']] *
         units(dfVars.loc['Geopotential Height']['Units']).to('gpm'))
    Q = AdiabaticHEating(tair, pres, omega ,u, v,
            VerticalCoordIndexer,LatIndexer,LonIndexer,TimeName)
    
    # Create csv files for storing vertical results.
    # When those terms are computed data is simply appended to csv
    for term in ['Az','Ae','Kz','Ke','Cz','Ca','Ck','Ce','Ge','Gz']:
        tmp = pd.DataFrame(columns=[TimeName,
        *[float(i)/100 for i in pres.values]])
        tmp.to_csv(ResultsSubDirectory+term+'_'+VerticalCoordIndexer+'.csv',
                   index=None)  
    # Track file
    if args.track:
        trackfile = '../inputs/track'
        track = pd.read_csv(trackfile,parse_dates=[0],
                            delimiter=';',index_col='time')

    # Dictionary for saving system position and attributes
    position = {}
    results_keys = ['time', 'central_lat', 'central_lon', 'length', 'width',
            'min_zeta_850','min_hgt_850','max_wind_850']
    for key in results_keys:
        position[key] =  []
    
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
    
    times = pd.to_datetime(data[TimeName].values)
    # Slice the time array so the first and the last timestep will be the same
    # as in the track file
    if args.track:
        times = times[(times>=track.index[0]) & (times<=track.index[-1])]
        if len(times) == 0:
            print("Mismatch between trackfile and data! Check that and try again!")
            sys.exit(1)
        
    # Loop for each time step:
    for t in times:
        
        idata = data.sel({TimeName:t})
        iQ = Q.sel({TimeName:t})
        # Get current time and box limits
        itime = str(t)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d-%H%M')
        
        # Open those variables for saving 
        iu_850 = u.sel({TimeName:t}).sel({VerticalCoordIndexer:850})
        iv_850 = v.sel({TimeName:t}).sel({VerticalCoordIndexer:850})
        ight_850 = hgt.sel({TimeName:t}).sel({VerticalCoordIndexer:850})
        zeta = vorticity(iu_850, iv_850).metpy.dequantify()
        # Apply filter when using high resolution gridded data
        dx = float(iv_850[LonIndexer][1]-iv_850[LonIndexer][0])
        if dx < 1:
            zeta = zeta.to_dataset(name='vorticity'
                ).apply(savgol_filter,window_length=31, polyorder=2).vorticity
            
        lat, lon = iu_850[LatIndexer], iu_850[LonIndexer]
        
        if args.track:
            # Get current time and box limits
            if 'width'in track.columns:
                width, length = track.loc[itime]['width'],track.loc[itime]['length']
            else:
                width, length = 15, 15
                
            # Track timestep closest to the model timestep, just in case
            # the track file has a poorer temporal resolution
            track_itime = track.index[track.index.get_loc(
                t, method='nearest')].strftime('%Y-%m-%d %HZ')
            min_lon = track.loc[track_itime]['Lon']-(width/2)
            max_lon = track.loc[track_itime]['Lon']+(width/2)
            min_lat = track.loc[track_itime]['Lat']-(length/2)
            max_lat = track.loc[track_itime]['Lat']+(length/2)
            limits = {'min_lon':min_lon,'max_lon':max_lon,
                      'min_lat':min_lat,'max_lat':max_lat}
        
        elif args.choose:
            # Draw maps and ask user to specify corners for specifying the box
            limits = draw_box_map(iu_850, iv_850, zeta, ight_850,
                                  lat, lon, itime)
            min_lon, max_lon = limits['min_lon'],  limits['max_lon']
            min_lat, max_lat = limits['min_lat'],  limits['max_lat']
        
        # Store system position and attributes
        central_lat = (limits['max_lat'] + limits['min_lat'])/2
        central_lon = (limits['max_lon'] + limits['min_lon'])/2
        length = limits['max_lat'] - limits['min_lat']
        width = limits['max_lon'] - limits['min_lon']
        min_zeta = float(zeta.min())
        min_hgt = float(ight_850.min())
        max_wind = float(wind_speed(iu_850, iv_850).max())
        
        values = [datestr, central_lat, central_lon, length, width,
                  min_zeta, min_hgt, max_wind]
        for key,val in zip(results_keys,values):
            position[key].append(val)
        
        print('\nTime: ',datestr)
        print('Box min_lon, max_lon: '+str(min_lon)+'/'+str(max_lon))
        print('Box min_lat, max_lat: '+str(min_lat)+'/'+str(max_lat))
        print('Box size (longitude): '+str(width))
        print('Box size (latitude): '+str(length))
        print('Minimum vorticity at 850 hPa:',min_zeta)
        print('Minimum geopotential height at 850 hPa:',min_hgt)
        print('Maximum wind speed at 850 hPa:',max_wind)
    
        # Create box object
        try:
            box_obj = BoxData(data=idata, dfVars=dfVars, args=args,
            western_limit=min_lon, eastern_limit=max_lon,
            southern_limit=min_lat, northern_limit=max_lat,
            output_dir=ResultsSubDirectory, Q=iQ)
        except:
            raise SystemExit('Error ondatestr creating the box for the computations')
        # Compute energy terms
        try:
            ec_obj = EnergyContents(box_obj,method='moving')
            TermsDict['Az'].append(ec_obj.calc_az())
            TermsDict['Ae'].append(ec_obj.calc_ae())
            TermsDict['Kz'].append(ec_obj.calc_kz())
            TermsDict['Ke'].append(ec_obj.calc_ke())
        except Exception as e:
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Energy Contents')
        # Compute conversion terms
        try:
            ct_obj = ConversionTerms(box_obj,method='moving')
            TermsDict['Ca'].append(ct_obj.calc_ca())
            TermsDict['Ce'].append(ct_obj.calc_ce())
            TermsDict['Ck'].append(ct_obj.calc_ck())
            TermsDict['Cz'].append(ct_obj.calc_cz())
        except Exception as e:
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Conversion Terms')
        # Compute boundary terms
        try:
            bt_obj = BoundaryTerms(box_obj,method='moving')
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
            gdt_obj = GenerationDissipationTerms(box_obj,method='moving')
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
    
    # Save system position as a csv file for replicability
    track = pd.DataFrame.from_dict(position)
    track = track.rename(columns={'central_lat':'Lat','central_lon':'Lon'})
    track.to_csv(ResultsSubDirectory+outfile_name+'_track',
                 index=False, sep=";")
    
    # 13) Make figures
    FigsDirectory = ResultsSubDirectory+'/Figures'
    check_create_folder(FigsDirectory)
    
    
    if args.residuals:
        flag = ' -r'
    else:
        flag = ' '
    os.system("python ../plots/plot_timeseries.py "+outfile+flag)
    os.system("python ../plots/plot_vertical.py "+ResultsSubDirectory)
    os.system("python ../plots/plot_boxplot.py "+ResultsSubDirectory+flag)
    os.system("python ../plots/plot_LEC.py "+outfile)
    os.system("python ../plots/plot_LPS.py "+outfile)
    os.system("python ../plots/plot_track.py "+outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "\
Lorenz Energy Cycle (LEC) program. \n \
The program can compute the LEC using two distinct frameworks:\
    1) moving framework. A box is definid in the box_lims' file and then the \
       energetics are computed for a fixed domain.\
    2) fixed framework. The domain is not fixed and follows the system using \
       the track file.\
 Both frameworks can be applied at the same time, given the required files are\
 provided. An auxilliary 'fvars' file is also needed for both frameworks: it\
 contains the specified names used for each variable. The program computes \
 the energy, conversion, boundary, generation and dissipation terms. The results \
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
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fixed", default = False,
    action='store_true', help = "compute the energetics for a Fixed domain\
 specified by the 'inputs/box_lims' file.")
    group.add_argument("-t", "--track", default = False,
    action='store_true', help = "define the box using a track file specified \
by the 'inputs/track' file. The track indicate the central point of the system\
and a arbitraty box of 15°x15° is constructed.")
    group.add_argument("-c", "--choose", default = False,
    action='store_true', help = "For each time step, the user can choose the\
domain by clicking on the screen.")
    parser.add_argument("-o", "--outname", default = False, type=str,
    help = "choose a name for saving results (default is\
 the same as infile)")
 
    args = parser.parse_args()
    
    infile  = args.infile
    varlist = '../inputs/fvars'
    
    # Open data
    data = get_data(infile, varlist) 
    
    # Slice data so the code runs faster
    data, method = slice_domain(data, args, varlist)
    
    # Directory where results will be stored
    ResultsMainDirectory = '../LEC_Results'
    # Append data limits to outfile name
    if args.outname:
        outfile_name = args.outname
    else:
        outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_'+method
    # Each dataset of results have its own directory, allowing to store results
    # from more than one experiment at each time
    ResultsSubDirectory = ResultsMainDirectory+'/'+outfile_name+'/'
    # Check if the LEC_Figures directory exists. If not, creates it
    check_create_folder(ResultsMainDirectory)
    # Check if a directory for current data exists. If not, creates it
    check_create_folder(ResultsSubDirectory)
    
    # Run the program
    start_time = time.time()
    if args.fixed:
        LEC_fixed(data)
        print("--- %s seconds running fixed framework ---" % (time.time() - start_time))
    start_time = time.time()
    if args.track or args.choose:
        LEC_moving(data)
        print("--- %s seconds for running moving framework ---" % (time.time() - start_time))
    

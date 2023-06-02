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
from metpy.constants import g

from select_area import slice_domain
from select_area import draw_box_map
from select_area import plot_domain_attributes

from determine_periods import determine_periods

from scipy.signal import savgol_filter 
from numpy.linalg import LinAlgError

import pandas as pd
import xarray as xr
import os
import numpy as np
import argparse
import dask
import time
import logging

logging.getLogger('dask').setLevel(logging.ERROR)
dask.config.set(scheduler="processes")

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

def get_data(infile: str, varlist: str) -> xr.Dataset:
    """
    Opens a NetCDF file and extracts the variables specified in a CSV file.

    Args:
        infile (str): The file path of the NetCDF file to be opened.
        varlist (str): The file path of the CSV file containing the list of variables to extract.

    Returns:
        xr.Dataset: A dataset containing the extracted variables.

    Raises:
        SystemExit: If the specified CSV file or NetCDF file is not found or if an error occurs while opening the NetCDF file.
    """
    print(f"Variables specified by the user in: {varlist}")
    print(f"Attempting to read {varlist} file...")
    try:
        dfVars = pd.read_csv(varlist, sep=';', index_col=0, header=0)
    except FileNotFoundError:
        raise SystemExit("Error: The 'fvar' text file could not be found.")
    except pd.errors.EmptyDataError:
        raise SystemExit("Error: The 'fvar' text file is empty.")

    print("List of variables found:")
    print(dfVars)

    LonIndexer = dfVars.loc["Longitude"]["Variable"]
    LatIndexer = dfVars.loc["Latitude"]["Variable"]
    LevelIndexer = dfVars.loc["Vertical Level"]["Variable"]

    print("Opening input data...")
    try:
        with dask.config.set(array={'slicing': {'split_large_chunks': True}}):
            data = convert_lon(
                xr.open_dataset(infile),
                dfVars.loc['Longitude']['Variable']
            )
    except FileNotFoundError:
        raise SystemExit("ERROR: Could not open file. Check if path, fvars file, and file format (.nc) are correct.")
    except:
        raise

    print("Assigning geospatial coordinates in radians...")
    data = data.assign_coords({"rlats": np.deg2rad(data[LatIndexer])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data[LatIndexer]))})
    data = data.assign_coords({"rlons": np.deg2rad(data[LonIndexer])})

    levels_Pa = (data[LevelIndexer] * units(str(data[LevelIndexer].units))).metpy.convert_units("Pa")
    data = data.assign_coords({LevelIndexer: levels_Pa})
    
    data = data.sortby(LonIndexer).sortby(LevelIndexer, ascending=True).sortby(LatIndexer, ascending=True)

    print("Data opened successfully.")
    return data

def find_extremum_coordinates(data, lat, lon, variable):
    """
    Finds the indices of the extremum values for a given variable.

    Args:
    data: An xarray DataArray containing the data to compute the energy cycle.
    lat: An xarray DataArray containing the latitudes of the data.
    lon: An xarray DataArray containing the longitudes of the data.
    variable: A string containing the name of the variable to find the indices for.

    Returns:
    A tuple containing the indices of the extremum values for the specified variable.
    """
    
    lat_values = lat.values
    lon_values = lon.values

    if variable == 'min_zeta':
        index = np.unravel_index(data.argmin(), data.shape)
    elif variable == 'min_hgt':
        index = np.unravel_index(data.argmin(), data.shape)
    elif variable == 'max_wind':
        index = np.unravel_index(data.argmax(), data.shape)
    else:
        raise ValueError("Invalid variable specified.")

    lat_index = index[0]
    lon_index = index[1]

    lat_value = lat_values[lat_index]
    lon_value = lon_values[lon_index]

    return lat_value, lon_value

def LEC_fixed(data):
    """
    Computes the Lorenz Energy Cycle using a fixed framework.
    
    Args:
    data: A Xarray Dataset containing the data to compute the energy cycle.
    
    Returns:
    None
    """
    
    print('Computing energetics using fixed framework')
    
    print('Opening data into memory..')
    data = data.compute()
    print('done.')
    
    # Read the box limits from box_limits file
    dfbox = pd.read_csv('../inputs/box_limits',header=None,delimiter=';',index_col=0)
    min_lon = float(dfbox.loc['min_lon'].values)
    max_lon = float(dfbox.loc['max_lon'].values)
    min_lat = float(dfbox.loc['min_lat'].values)
    max_lat = float(dfbox.loc['max_lat'].values)
    
    # Raise a ValueError if the box limits are invalid
    if min_lon > max_lon:
        raise ValueError('Error in box_limits: min_lon > max_lon')
        quit()
    if min_lat > max_lat:
        raise ValueError('Error in box_limits: min_lat > max_lat')
        quit()
        
    dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)

    LonIndexer,LatIndexer,TimeName,VerticalCoordIndexer = (
      dfVars.loc['Longitude']['Variable'],
      dfVars.loc['Latitude']['Variable'],
      dfVars.loc['Time']['Variable'],
      dfVars.loc['Vertical Level']['Variable'])
    
    pres = data[VerticalCoordIndexer]
    pres = pres * pres.metpy.units
    
    print('\n Parameters spcified for the bounding box:')
    print('min_lon, max_lon, min_lat, max_lat: '+str([min_lon,
                                                     max_lon, min_lat, max_lat]))
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
    df = calc_budget_diff(df,dates, args) 
    print('Ok!')
    # 11) 
    print('\n------------------------------------------------------------------------')
    print('Computing residuals RGz, RKz, RGe and RKe')
    df = calc_residuals(df, args)
    print('Ok!')
    # 12) save file
    print('\nCreating a csv to store results...')
    outfile = ResultsSubDirectory+'/'+outfile_name+'.csv'
    df.to_csv(outfile)
    print(outfile+' created') 
    print('All done!')
    # 13) Make figures
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
    

def LEC_moving(data, dfVars, dTdt, ResultsSubDirectory, FigsDirectory):
    """
    Computes the Lorenz Energy Cycle using a moving framework.
    
    Args:
    data: A Xarray Dataset containing the data to compute the energy cycle.
    
    Returns:
    None
    """
    print('Computing energetics using moving framework')

    # Indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = (
      dfVars.loc['Longitude']['Variable'],
      dfVars.loc['Latitude']['Variable'],
      dfVars.loc['Time']['Variable'],
      dfVars.loc['Vertical Level']['Variable']
      )
    
    PressureData = data[VerticalCoordIndexer]
    
    # Create csv files for storing vertical results
    for term in ['Az', 'Ae', 'Kz', 'Ke', 'Cz', 'Ca', 'Ck', 'Ce', 'Ge', 'Gz']:
        columns = [TimeName] + [float(i) for i in PressureData.values]
        df = pd.DataFrame(columns=columns)
        file_name = term + '_' + VerticalCoordIndexer + '.csv'
        file_path = os.path.join(ResultsSubDirectory, file_name)
        df.to_csv(file_path, index=None) 
        print(file_path+' created (but still empty)')

    # Track file
    if args.track:
        trackfile = '../inputs/track'
        track = pd.read_csv(trackfile,parse_dates=[0],
                            delimiter=';',index_col='time')

    # Dictionary for saving system position and attributes
    results_keys = ['datestr', 'central_lat', 'central_lon', 'length', 'width',
                'min_zeta_850', 'min_hgt_850', 'max_wind_850']
    out_track = pd.DataFrame(columns=results_keys)

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
    
    # Slice the time array so the first and the last timestep will be the same
    # as in the track file
    times = pd.to_datetime(data[TimeName].values)
    if args.track:
        times = times[(times>=track.index[0]) & (times<=track.index[-1])]
        if len(times) == 0:
            raise ValueError("Mismatch between trackfile and data! Check that and try again!")
        
    # times = times[:10]
    for t in times:

        idata = data.sel({TimeName: t})
        idTdt = dTdt.sel({TimeName: t})

        # Sometimes there are errors on indexing the time coordinate on the NetCDF file
        #  and I can't do anything about it. So I will have to discard one of the timesteps.
        # Any complaints can be sent directly to the original owner of the .nc file.
        if idata[TimeName].shape != ():
            idata = data.isel({TimeName: 1})
            idTdt = dTdt.isel({TimeName: 1})

        iu = (idata[dfVars.loc['Eastward Wind Component']['Variable']].compute() *
            units(dfVars.loc['Eastward Wind Component']['Units']).to('m/s'))
        iv = (idata[dfVars.loc['Northward Wind Component']['Variable']].compute() *
            units(dfVars.loc['Northward Wind Component']['Units']).to('m/s'))
        if args.geopotential:
            ihgt = ((idata[dfVars.loc['Geopotential']['Variable']].compute() *
            units(dfVars.loc['Geopotential']['Units']))/g).metpy.convert_units('gpm')
        else:
            ihgt = (idata[dfVars.loc['Geopotential Height']['Variable']].compute() *
            units(dfVars.loc['Geopotential Height']['Units']).to('gpm'))

        iu_850, iv_850, ight_850 = (
            iu.sel({VerticalCoordIndexer: 85000}),
            iv.sel({VerticalCoordIndexer: 85000}),
            ihgt.sel({VerticalCoordIndexer: 85000})
        )

        iwspd_850, izeta_850 = wind_speed(iu_850, iv_850), vorticity(iu_850, iv_850).metpy.dequantify()

        lat, lon = idata[LatIndexer], idata[LonIndexer]

        # Get current time and box limits
        itime = str(t)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d-%H%M')

        # Apply filter when using high resolution gridded data
        dx = float(idata[LonIndexer][1] - idata[LonIndexer][0])
        if dx < 1:
            try:
                izeta_850 = izeta_850.to_dataset(name='vorticity').apply(savgol_filter, window_length=31, polyorder=2).vorticity
            except LinAlgError:
                izeta_850 = izeta_850.fillna(0).to_dataset(name='vorticity').apply(savgol_filter, window_length=31, polyorder=2).vorticity
            except Exception as e:
                raise e

        if args.track:
            # Get current time and box limits
            if 'width'in track.columns:
                width, length = track.loc[itime]['width'],track.loc[itime]['length']
            else:
                width, length = 15, 15
                
            # Track timestep closest to the model timestep, just in case
            # the track file has a poorer temporal resolution
            closest_index = np.searchsorted(track.index, t)
            if closest_index > 0 and (closest_index == len(track.index
                                ) or abs(t - track.index[closest_index-1]
                                ) < abs(t - track.index[closest_index])):
                closest_index -= 1
            track_itime = track.index[closest_index]
                # Store system position and attributes
            central_lon, central_lat = track.loc[track_itime]['Lon'], track.loc[track_itime]['Lat']
            min_lon = central_lon-(width/2)
            max_lon = central_lon+(width/2)
            min_lat = central_lat-(length/2)
            max_lat = central_lat+(length/2)
            limits = {'min_lon':min_lon,'max_lon':max_lon,
                        'min_lat':min_lat,'max_lat':max_lat}
            
            # Slice data for defined box
            izeta_850_slice = izeta_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            ight_850_slice = ight_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            iwspd_850_slice = iwspd_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})

                # Check if 'min_zeta_850', 'min_hgt_850' and 'max_wind_850' columns exists in the track file.
                # If they exist, then retrieve and convert the value from the track file.
                # If they do not exist, calculate them.
            try:
                min_zeta = float(track.loc[track_itime]['min_zeta_850'])
            except KeyError:
                if args.zeta:
                    min_zeta_unformatted = izeta_850.sel(latitude=central_lat, longitude=central_lon, method='nearest')
                else:
                    min_zeta_unformatted = izeta_850_slice.min()
                min_zeta = float(np.nanmin(min_zeta_unformatted))
            try:
                min_hgt = float(track.loc[track_itime]['min_hgt_850'])
            except KeyError:
                min_hgt = float(ight_850_slice.min())
            try:
                max_wind = float(track.loc[track_itime]['max_wind_850'])
            except KeyError:
                max_wind = float(iwspd_850_slice.max())

        elif args.choose:

            # Draw maps and ask user to specify corners for specifying the box
            limits = draw_box_map(iu_850, iv_850, izeta_850, ight_850,
                                    lat, lon, itime)
            
                # Store system position and attributes
            min_lon, max_lon = limits['min_lon'],  limits['max_lon']
            min_lat, max_lat = limits['min_lat'],  limits['max_lat']
            width, length = limits['max_lon'] - limits['min_lon'], limits['max_lat'] - limits['min_lat']
            central_lat = (limits['max_lat'] + limits['min_lat'])/2
            central_lon = (limits['max_lon'] + limits['min_lon'])/2

            # Slice data for defined box and find extremes
            izeta_850_slice = izeta_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            ight_850_slice = ight_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            iwspd_850_slice = iwspd_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            min_zeta = float(izeta_850_slice.min())
            min_hgt = float(ight_850_slice.min())
            max_wind = float(iwspd_850_slice.max())

        # Find position of the extremes
        lat_slice, lon_slice = izeta_850_slice[LatIndexer], izeta_850_slice[LonIndexer]
        min_zeta_lat, min_zeta_lon = find_extremum_coordinates(izeta_850_slice, lat_slice, lon_slice, 'min_zeta')
        min_hgt_lat, min_hgt_lon = find_extremum_coordinates(ight_850_slice, lat_slice, lon_slice, 'min_hgt')
        max_wind_lat, max_wind_lon = find_extremum_coordinates(iwspd_850_slice, lat_slice, lon_slice, 'max_wind')

        # Store the results in a dictionary for plotting purposes
        data850 = {
            'min_zeta': {
                'latitude': min_zeta_lat,
                'longitude': min_zeta_lon,
                'data': izeta_850
            },
            'min_hgt': {
                'latitude': min_hgt_lat,
                'longitude': min_hgt_lon,
                'data': ight_850
            },
            'max_wind': {
                'latitude': max_wind_lat,
                'longitude': max_wind_lon,
                'data': iwspd_850
            },
            'lat': lat,
            'lon': lon,
        }

        position = {
        'datestr': datestr,
        'central_lat': central_lat,
        'central_lon': central_lon,
        'length': length,
        'width': width,
        'min_zeta_850': min_zeta,
        'min_hgt_850': min_hgt,
        'max_wind_850': max_wind
        }

        out_track = out_track.append(position, ignore_index=True)

        plot_domain_attributes(data850, position, FigsDirectory)

        print(f'\nTime: {datestr}')
        print(f'central lat/lon: {central_lon}, {central_lat}')
        print(f'Box min_lon, max_lon: {min_lon}/{max_lon}')
        print(f'Box min_lat, max_lat: {min_lat}/{max_lat}')
        print(f'Box size (longitude): {width}')
        print(f'Box size (latitude): {length}')
        print(f'Minimum vorticity at 850 hPa: {min_zeta}')
        print(f'Minimum geopotential height at 850 hPa: {min_hgt}')
        print(f'Maximum wind speed at 850 hPa: {max_wind}')
    
        # Create box object
        try:
            box_obj = BoxData(
                data=idata.compute(),
                dfVars=dfVars,
                args=args,
                western_limit=min_lon,
                eastern_limit=max_lon,
                southern_limit=min_lat,
                northern_limit=max_lat,
                output_dir=ResultsSubDirectory,
                dTdt=idTdt
            )
        except Exception as e:
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error creating the box for computations')
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
            TermsDict['BAe'].append(bt_obj.calc_bae())
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
    if args.verbosity:
        for term in TermsDict.keys():
            print('\n'+term)
            print(df[term].tolist())
    
    print('\n------------------------------------------------------------------------')
    print('Estimating budget terms (∂X/∂t) using finite differences ')
    df = calc_budget_diff(df,times, args) 
    print('Ok!')
    
    print('\n------------------------------------------------------------------------')
    print('Computing residuals RGz, RKz, RGe and RKe')
    df = calc_residuals(df, args)
    print('Ok!')
    
    print('\nCreating a csv to store results...')
    outfile = ResultsSubDirectory+'/'+outfile_name+'.csv'
    df.to_csv(outfile)
    print(outfile+' created') 
    print('All done!')
    
    # Save system position as a csv file for replicability
    # df = pd.DataFrame.from_dict(position, orient='index').T
    out_track = out_track.rename(columns={'datestr':'time','central_lat':'Lat','central_lon':'Lon'})
    output_trackfile =  ResultsSubDirectory+outfile_name+'_track'
    out_track.to_csv(output_trackfile, index=False, sep=";")
    
    # 13) Make figures directory
    FigsDirectory = ResultsSubDirectory+'/Figures'
    check_create_folder(FigsDirectory)
    
    # Determine periods
    determine_periods(output_trackfile, ResultsSubDirectory)
    
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
    
    parser = argparse.ArgumentParser(description="Lorenz Energy Cycle (LEC) program.")
    parser.add_argument("infile", help="Input .nc file with temperature, geopotential, and wind components.")
    parser.add_argument("-r", "--residuals", action='store_true', help="Compute the Dissipation and Generation terms as residuals.")
    parser.add_argument("-g", "--geopotential", action='store_true', help="Use geopotential data instead of geopotential height.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fixed", action='store_true', help="Compute the energetics for a fixed domain specified by the 'box_lims' file.")
    group.add_argument("-t", "--track", action='store_true', help="Define the box using a track file specified by the 'track' file.")
    group.add_argument("-c", "--choose", action='store_true', help="Choose the domain for each time step by clicking on the screen.")
    parser.add_argument("-z", "--zeta", action='store_true', help="Use this flag if the track file was created using vorticity.")
    parser.add_argument("-o", "--outname", type=str, help="Choose a name for saving results.")
    parser.add_argument("-v", "--verbosity", action='store_true', help="Increase output verbosity.")

    # args = parser.parse_args()

    # Debuggings:
    # args = parser.parse_args(['../samples/Reg1-Representative_NCEP-R2.nc', '-r', '-t'])
    args = parser.parse_args(['/p1-nemo/danilocs/mpas/MPAS-BR/post_proc/py/interpolations/Catarina-2403-2903_MPAS.nc',
     '-t', '-g', '-r'])

    infile  = args.infile
    varlist = '../inputs/fvars'
    
    # Open data
    data = get_data(infile, varlist) 
    
    # Slice data so the code runs faster
    data, method = slice_domain(data, args, varlist)
    print(f'Data sliced: \n{data.coords}')

    # Compute terms with temporal dependencies
    print('Computing terms with temporal dependencies..')
    dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)
    dTdt =  data[dfVars.loc['Air Temperature']['Variable']].differentiate(
                dfVars.loc['Time']['Variable'],datetime_unit='s') * units('K/s')
    print('ok')
    
    # Append data limits to outfile name
    if args.outname:
        outfile_name = args.outname
    else:
        outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_'+method

    # Directory where results will be stored
    ResultsMainDirectory = '../LEC_Results'
    # Each dataset of results have its own directory, allowing to store results
    # from more than one experiment at each time
    ResultsSubDirectory = ResultsMainDirectory+'/'+outfile_name+'/'
    FigsDirectory = ResultsSubDirectory+'/Figures'
    check_create_folder(ResultsMainDirectory)
    check_create_folder(ResultsSubDirectory)
    check_create_folder(FigsDirectory)
    
    # Run the program
    start_time = time.time()
    if args.fixed:
        LEC_fixed(data)
        print("--- %s seconds running fixed framework ---" % (time.time() - start_time))
    start_time = time.time()
    if args.track or args.choose:
        LEC_moving(data, dfVars, dTdt, ResultsSubDirectory, FigsDirectory)
        print("--- %s seconds for running moving framework ---" % (time.time() - start_time))
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 17:17:29 2023

@author: daniloceano
"""
import glob
import matplotlib 

import xarray as xr
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import cmocean.cm as cmo

from metpy.calc import vorticity
from metpy.units import units

from scipy.signal import savgol_filter 
from scipy.ndimage.measurements import label, find_objects

from sklearn.preprocessing import normalize   
from sklearn import preprocessing

import cartopy.crs as ccrs
import cartopy

def convert_lon(xr,LonIndexer):
    xr.coords[LonIndexer] = (xr.coords[LonIndexer] + 180) % 360 - 180
    xr = xr.sortby(xr[LonIndexer])
    return xr 

def filter_var(variable):
    window_lenght = round(len(variable)/2)
    if (window_lenght % 2) == 0:
        window_lenght += 1
    return savgol_filter(variable, window_lenght, 3, mode="nearest")

def normalise_var(variable):
    return (variable - variable.min()) / (variable.max() - variable.min()) 

def array_vorticity(df):
    
    da = df.to_xarray()
    
    # Filter vorticity twice
    zeta_fil = xr.DataArray(filter_var(da.zeta), coords={'time':df.index})
    da = da.assign(variables={'zeta_fil':zeta_fil})
    zeta_fil2 = xr.DataArray(filter_var(zeta_fil),
                                 coords={'time':df.index})
    da = da.assign(variables={'zeta_fil2':zeta_fil2})
    
    # derivatives of the double-filtered vorticity
    da = da.assign(variables={'dzfil2_dt':
                da.zeta_fil2.differentiate('time',datetime_unit='h')})
        
    da = da.assign(variables={'dzfil2_dt2':
                da.dzfil2_dt.differentiate('time',datetime_unit='h')})
        
    da = da.assign(variables={'dzfil2_dt3':
                da.dzfil2_dt2.differentiate('time',datetime_unit='h')})
        
    # filter derivatives
    da = da.assign(variables={'dz_dt_fil2':
        xr.DataArray(filter_var(da.dzfil2_dt), coords={'time':df.index})})
    da = da.assign(variables={'dz_dt2_fil2':
        xr.DataArray(filter_var(da.dzfil2_dt2), coords={'time':df.index})})
    da = da.assign(variables={'dz_dt3_fil2':
        xr.DataArray(filter_var(da.dzfil2_dt3), coords={'time':df.index})})

    return da

def plot_track(da, fname):
    plt.close('all')
    datacrs = ccrs.PlateCarree()
    fig = plt.figure(figsize=(15, 12))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[:, :], projection=datacrs,frameon=True)
    ax.set_extent([-80, 40, 0, -70], crs=datacrs) 
    ax.coastlines(zorder = 1)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("lightblue"))
    gl = ax.gridlines(draw_labels=True,zorder=2,linestyle='dashed',alpha=0.8,
                 color='#383838')
    gl.xlabel_style = {'size': 14, 'color': '#383838'}
    gl.ylabel_style = {'size': 14, 'color': '#383838'}
    gl.bottom_labels = None
    gl.right_labels = None
    ax.plot(da.lon,da.lat,'-',c='k')
    scatter = ax.scatter(da.lon,da.lat,zorder=100,cmap=cmo.deep_r,c=da.zeta)
    plt.colorbar(scatter, pad=0.07, orientation='vertical',fraction=0.026,
                      label=' 850 hPa vorticity (ζ)')
    outname = '../vorticity_analysis/tracks/'+fname+'_track.png'
    plt.savefig(outname,dpi=500)
    print(outname,'saved')

def plot_timeseries(ax, x, *args, **kwargs):
    colors = ['#1d3557', '#d62828', '#606c38', '#f77f00']
    ls = [2.5, 2, 2, 2]
    for i in range(len(args)):
        ax.plot(x,args[i],c=colors[i],linewidth=ls[i],
                label=kwargs['labels'][i])
    plt.legend()
    plt.grid(linewidth=0.5, alpha=0.5)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()
    
def get_peaks_valleys(function):
    maxs = function.where(function > function.quantile(0.75))
    mins = function.where(function < function.quantile(0.25))
    # get indices of continuous data
    labels_peak, num_features_peak = label(~np.isnan(maxs))
    labels_valley, num_features_valley = label(~np.isnan(mins))
    slices_peak = find_objects(labels_peak)
    slices_valley = find_objects(labels_valley)
    # Get continuous series for dz_dt3 values higher than the 0.85 quantile
    # and lower than the 0.15 quantile
    continuous_max_dz = [maxs[sl] for sl in slices_peak if sl is not None]
    continuous_min_dz = [mins[sl] for sl in slices_valley if sl is not None]
    # Now, get maximum and mininum values (peaks and valleys)
    peaks, valleys = [], []
    times_peaks, times_valleys = [], []
    for data_peak in continuous_max_dz:
        peaks.append(data_peak.max())
        times_peaks.append(data_peak.idxmax())
    for data_valley in continuous_min_dz:
        valleys.append(data_valley.min())
        times_valleys.append(data_valley.idxmin())
    peaks = [float(p.values) for p in peaks]
    times_peaks = [t.values for t in times_peaks]
    df_peaks = pd.Series(peaks,index=times_peaks)
    valleys = [float(v.values) for v in valleys]
    times_valleys = [t.values for t in times_valleys]
    df_valleys = pd.Series(valleys,index=times_valleys)
    return df_peaks, df_valleys

def filter_peaks_valleys(series):
    # Remove peaks and valleys that occured very close to each other
    min_time_diff = pd.Timedelta(days=1)
    valid_indices = [series.index[0]]
    for i in range(1, len(series)):
        if series.index[i] - series.index[i-1] > min_time_diff:
            valid_indices.append(series.index[i])
    return series.loc[valid_indices]

def mature_stage(dz3):
    # Mature stage will be defined as all points between a consecutive
    # maximum and minimum of the third vorticity!
    dz3_peaks, dz3_valleys = get_peaks_valleys(dz3)
    dz3_peaks_filtered = filter_peaks_valleys(dz3_peaks)
    dz3_valleys_filtered = filter_peaks_valleys(dz3_valleys)
    
    peaks = dz3_peaks_filtered.replace(dz3_peaks_filtered.values,
                                       'peak')
    valleys = dz3_valleys_filtered.replace(dz3_valleys_filtered.values,
                                           'valley')
    dz3_peaks_valleys = pd.concat([peaks,valleys]).sort_index()
    
    dt = dz3.time[1] - dz3.time[0]
    dt = pd.Timedelta(dt.values)
        
    mature = []
    for i in range(len(dz3_peaks_valleys[:-1])):
        if (dz3_peaks_valleys[i] == 'peak') and \
            (dz3_peaks_valleys[i+1] == 'valley') :
            mature.append(pd.date_range(dz3_peaks_valleys.index[i],
                        dz3_peaks_valleys.index[i+1],
                        freq=f'{int(dt.total_seconds() / 3600)} H'))
    return mature


def intensification_incipient_stages(mature, dz2):
    
    dz2_peaks, dz2_valleys = get_peaks_valleys(dz2)
    dt = dz2.time[1] - dz2.time[0]
    dt = pd.Timedelta(dt.values)
    
    # intensification stage is defined as the period between a local minimum
    # of dz_dt2 and the beginning of the mature stage~
    intensification = []
    for series, i in zip(mature, range(len(mature))):
        mature_start = series[0] 
        if i == 0:
            # If thre is a minimum of the second derivative before the start of 
            # mature phase, the intensification starts there. Else, there is an
            # incipient stage
            if dz2_valleys.index[i] < mature_start:
                intensification_start = dz2_valleys.index[i]
                incipient = pd.date_range(dz2.time[0].values,intensification_start-dt,
                                freq=f'{int(dt.total_seconds() / 3600)} H')
            else:
                intensification_start = (dz2.time[0]).values
                incipient = []
        else:
            if dz2_valleys.index[0] < mature[0][0]:
                intensification_start = dz2_valleys.index[i]
            else:
                intensification_start = dz2_valleys.index[i-1]
        intensification_end = mature_start-dt
        intensification.append(pd.date_range(intensification_start,
                        intensification_end, 
                        freq=f'{int(dt.total_seconds() / 3600)} H'))   
    return intensification, incipient
    
def decaying_stage(mature, z, dz2):
    
    dz2_peaks, dz2_valleys = get_peaks_valleys(dz2)
    dt = z.time[1] - z.time[0]
    dt = pd.Timedelta(dt.values)
    
    decaying = []
    
    # If there's only one mature stage, all time steps after it are considered
    # as decaying stage
    if len(mature) == 1:
        decaying_start = mature[0][-1] + dt
        decaying_end = z.time[-1].values
        decaying.append(pd.date_range(decaying_start, decaying_end, 
                        freq=f'{int(dt.total_seconds() / 3600)} H'))
    
    # Else, decaying is defined as the period between the end of the mature
    # stage and the minimum of the second  derivative of the vorticity 
    else:
        for series, i in zip(mature, range(len(mature))):
            mature_end = series[-1]
            sliced_z = z.sel(time=slice(mature_end, z.time.max()))
            decaying_end = dz2_valleys[
                dz2_valleys.index >= sliced_z.time.min().values].index.min()
            # If there isn't a peak of dz_dt2 after the end of the mature
            # period, the decaying period goest until the end of the system
            if pd.isnull(decaying_end):
                decaying_end = z.time[-1].values            
            decaying_start = mature_end+dt
            decaying.append(pd.date_range(decaying_start, decaying_end, 
                            freq=f'{int(dt.total_seconds() / 3600)} H'))
            
    return decaying

def get_phases(da, outfile_name):
    
    z = da.zeta_fil2
    dz2 = da.dz_dt2_fil2
    dz3 = da.dz_dt3_fil2
    
    mature = mature_stage(dz3)

    intensification, incipient = intensification_incipient_stages(mature, dz2)
    
    decaying = decaying_stage(mature, z, dz2)
    
    # Add six hours in the end and beggining of each period, so there is an
    # overlap between each period, acting as a confidence interval
    six_hours = pd.Timedelta('6H')
    dt = z.time[1] - z.time[0]
    dt = pd.Timedelta(dt.values)
    
    phases = {}
    for phase, key in zip([[incipient], intensification, mature, decaying],
                    ['incipient', 'intensification', 'mature', 'decay']):
        tmp = []
        for period in phase:
            if len(period) > 0:
                tmp.append(pd.date_range(
                    period[0]-six_hours, period[-1]+six_hours,
                              freq=f'{int(dt.total_seconds() / 3600)} H'))
        phases[key] = tmp
        
    # Extract the first and last elements from each list
    df_dict = {}
    for k, v in phases.items():
        if len(v) == 0:
            df_dict[k] = []
        else:
            df_dict[k] = [v[0][0], v[-1][-1]]
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame.from_dict(df_dict, orient='index',
                                columns=['start', 'end'])
    
    df.to_csv(outfile_name+'.csv')
    print(outfile_name+'.csv saved')
     
    return phases
    
def plot_periods(da, periods, fname):
    
    z = da.zeta
    
    plt.close('all')
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111,frameon=True)
    colors = ['k', '#134074', '#d62828', '#f7b538', '#5b8e7d',]
    
    # Plot periods
    five_hours = pd.Timedelta('5H')
    six_hours = pd.Timedelta('6H')
    for phase, c in zip(periods.keys(),
                          ['#65a1e6','#d62828','#f7b538','#9aa981']):
        for series in periods[phase]:
            if len(series) > 0:
                ax.fill_betweenx((z.min(),z.max()+1e-5),
                            series[0]+six_hours,series[-1]-five_hours,
                            facecolor=c, alpha=0.4)
    
    ax.plot(da.zeta.time, da.zeta,c='#6b675b',
            linewidth=4,label=r'$ζ$', alpha=0.8)
    ax.plot(da.zeta_fil2.time, da.zeta_fil2,c=colors[0],
            linewidth=6,label=r'$ζ_f$')
    
    y = np.arange(z.min(),z.max()+5e-6,1e-5)
    plt.ylim(y[0],y[-1]+5e-6)
    plt.xlim(da.zeta_fil.time[0].values, da.zeta_fil.time[-1].values-five_hours)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.grid(linewidth=0.5, alpha=0.5)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()
    
    plt.tight_layout()
    
    outname = fname+'.png'
    plt.savefig(outname,dpi=500)
    print(outname,'saved')
    


def get_periods(track_file, varlist, outfile_name):
    
    dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)
    
    track = pd.read_csv(track_file, parse_dates=[0],delimiter=';',index_col=[0])
    
    df_zeta = pd.DataFrame(track['min_hgt_850'].rename('zeta'))
        
    da = array_vorticity(df_zeta)
    
    periods  = get_phases(da, outfile_name)
        
    plot_periods(da, periods, outfile_name)


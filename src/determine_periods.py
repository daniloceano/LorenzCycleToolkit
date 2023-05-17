#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 17:17:29 2023

@author: daniloceano
"""
import os
import re

import xarray as xr
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import cmocean.cm as cmo

from scipy.signal import savgol_filter 
from scipy.ndimage.measurements import label, find_objects


def check_create_folder(DirName):
    if not os.path.exists(DirName):
                os.makedirs(DirName)
                print(DirName+' created')
    else:
        print(DirName+' directory exists')

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
    """
    Calculate vorticity of an input dataframe

    Args:
    df: pandas DataFrame

    Returns:
    xarray DataArray
    """
    # Convert dataframe to xarray
    da = df.to_xarray()

    # Filter vorticity twice
    zeta_filt = xr.DataArray(filter_var(da.zeta), coords={'time':df.index})
    da = da.assign(variables={'zeta_filt': zeta_filt})
    zeta_filt2 = xr.DataArray(filter_var(zeta_filt), coords={'time':df.index})
    da = da.assign(variables={'zeta_filt2': zeta_filt2})

    # Calculate derivatives of the double-filtered vorticity
    dzfilt2_dt = da.zeta_filt2.differentiate('time', datetime_unit='h')
    dzfilt2_dt2 = dzfilt2_dt.differentiate('time', datetime_unit='h')
    dzfilt2_dt3 = dzfilt2_dt2.differentiate('time', datetime_unit='h')

    # Filter derivatives
    dz_dt_filt2 = xr.DataArray(filter_var(dzfilt2_dt), coords={'time':df.index})
    dz_dt2_filt2 = xr.DataArray(filter_var(dzfilt2_dt2), coords={'time':df.index})
    dz_dt3_filt2 = xr.DataArray(filter_var(dzfilt2_dt3), coords={'time':df.index})

    # Assign variables to xarray
    da = da.assign(variables={'dzfilt2_dt': dzfilt2_dt,
                              'dzfilt2_dt2': dzfilt2_dt2,
                              'dzfilt2_dt3': dzfilt2_dt3,
                              'dz_dt_filt2': dz_dt_filt2,
                              'dz_dt2_filt2': dz_dt2_filt2,
                              'dz_dt3_filt2': dz_dt3_filt2})

    return da
    
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
    filtered_series = series.copy()

    # Find the indices of peaks or valleys
    indices = filtered_series.index

    # Iterate through the indices
    for i in range(len(indices) - 1):
        # Check if the time difference between consecutive indices is less than 24 hours
        if indices[i + 1] - indices[i] < pd.Timedelta(hours=24):
            # Remove the peak or valley that occurs later
            if filtered_series[indices[i + 1]] > filtered_series[indices[i]]:
                filtered_series = filtered_series.drop(indices[i + 1])
            else:
                filtered_series = filtered_series.drop(indices[i])

    return filtered_series

def plot_didactic(vorticity: xr.Dataset, periods, outfile_name_didactic: str) -> None:
    """
    Plots a series of subplots illustrating the different stages of a cyclone development.

    Args:
        vorticity: A xarray Dataset containing vorticity data and its derivatives
        outfile_name_didactic: The name (full path) of the file where the plot will be saved
    """

    def plot_subplot(ax: plt.Axes, z: pd.Series, z_filt: pd.Series, dz: pd.Series,
                      dz2: pd.Series, dz3: pd.Series) -> None:
        ax.plot(vorticity.time, dz3, c=colors_phases["incipient"],
                 linewidth=0.75, label=r"$\frac{∂^{3}ζ}{∂t^{3}}$")
        ax.plot(vorticity.time, dz2, c=colors_phases["intensification"],
                 linewidth=0.75, label=r"$\frac{∂^{2}ζ}{∂t^{2}}$")
        ax.plot(vorticity.time, dz, c=colors_phases["mature"], 
                linewidth=0.75, label=r"$\frac{∂ζ}{∂t}$")
        ax.plot(vorticity.time, z, c="gray", linewidth=0.75, label=r"$ζ_{filt}$")
        ax.plot(vorticity.time, z_filt, c="k", linewidth=2, label="ζ")
        ax.set_title("Plot ζ and its derivatives")
        ax.legend(loc="upper center", bbox_to_anchor=(1.67, 1.45), ncol=5, fontsize=16)

    def plot_peaks_valleys(ax, peaks, valleys, peaks_color, valleys_color):
        # Plot peaks
        if peaks is not None:
            if isinstance(peaks, list):
                s = 0
                for i, peak in enumerate(peaks):
                    peaks_time, peaks_z = find_matching_z_values(peak, z)
                    peak_color = peaks_color[i] if i < len(peaks_color) else None
                    peak_size = sizes[s] if s < len(sizes) else 150
                    ax.scatter(peaks_time, peaks_z, facecolor=peak_color, edgecolor='k', s=peak_size)
                    s += 1
            else:
                peaks_time, peaks_z = find_matching_z_values(peaks, z)
                ax.scatter(peaks_time, peaks_z, facecolor=peaks_color, edgecolor='k', s=150)

        # Plot valleys
        if valleys is not None:
            if isinstance(valleys, list):
                s = 0
                for i, valley in enumerate(valleys):
                    valleys_time, valleys_z = find_matching_z_values(valley, z)
                    valley_color = valleys_color[i] if i < len(valleys_color) else None
                    valley_size = sizes[s] if s < len(sizes) else 150
                    ax.scatter(valleys_time, valleys_z, facecolor=valley_color, s=valley_size)
                    s += 1
            else:
                valleys_time, valleys_z = find_matching_z_values(valleys, z)
                ax.scatter(valleys_time, valleys_z, facecolor=valleys_color, s=150)


    def plot_stage(ax, stage_name=None):
        # Plot vorticity data
        ax.plot(vorticity.time, z, c="k", linewidth=2)

        if stage_name is not None:
            # Compute stage periods
            stage_periods = periods[periods.index.str.contains(stage_name)]

            # Plot stage data
            for _, row in stage_periods.iterrows():
                stage_start = row['start']
                stage_end = row['end']
                stage_data = z.sel(time=slice(stage_start, stage_end))

                ax.fill_betweenx((z.min(), z.max()), stage_data.time.min(), stage_data.time.max(),
                                facecolor=colors_phases[stage_name], alpha=0.6)

        # Set plot title
        if stage_name is not None:
            ax.set_title(stage_name.capitalize() + ' stage')


    def find_matching_z_values(dz2_valleys, z):
        matching_times = []
        matching_z_values = []

        for time in dz2_valleys.index:
            if time in z.time:
                matching_times.append(time)
                matching_z_values.append(z.sel(time=time))

        return matching_times, matching_z_values

    colors_phases = {'incipient': '#65a1e6', 'intensification': '#f7b538',
          'mature': '#d62828', 'decay': '#9aa981'}
    
    colors_derivatives = {'dz':'#d62828', 'dz2':'#f7b538', 'dz3':'#65a1e6'}

    z = vorticity.zeta_filt2
    dz = vorticity.dz_dt_filt2 * 50
    dz2 = vorticity.dz_dt2_filt2 * 500
    dz3 = vorticity.dz_dt3_filt2 * 5000

    dz_peaks, dz_valleys = get_peaks_valleys(dz)
    dz2_peaks, dz2_valleys = get_peaks_valleys(dz2)
    dz3_peaks, dz3_valleys = get_peaks_valleys(dz3)

    plt.close("all")
    fig = plt.figure(figsize=(15, 18))
    gs = gridspec.GridSpec(3, 3)

    # 1) Plot vorticity and its derivatives
    ax_list = []
    for i in range(3):
        for j in range(3):
            ax_list.append(fig.add_subplot(gs[i, j], frameon=True))
    plot_subplot(ax_list[0], vorticity.zeta, z, dz, dz2, dz3)

    # 2) Plot dz peaks and valleys
    ax2 = fig.add_subplot(gs[0, 1], frameon=True)
    col = list(colors_derivatives.keys()) 
    sizes = [350, 250, 150]
    s = 0
    for function, color_key in zip([dz, dz2, dz3], col):
        peaks, valleys = get_peaks_valleys(function)
        ax2.scatter(peaks.index, peaks, facecolor=colors_derivatives[color_key],
                     edgecolor='k', linewidth=2, s=sizes[s])
        ax2.scatter(valleys.index, valleys, facecolor=colors_derivatives[color_key],
                     linewidth=2, s=sizes[s])
        ax2.plot(function.time, function, c=colors_derivatives[color_key],
                  linewidth=1, linestyle='-')
        s += 1
    ax2.plot(vorticity.time, z, c="k", linewidth=2)
    ax2.set_title('Get peaks and valleys in derivatives')
    
    # 3) Filter peaks and valeys closer than 1 day
    ax3 = fig.add_subplot(gs[0, 2], frameon=True)
    dz_peaks_filtered, dz_valleys_filtered = filter_peaks_valleys(dz_peaks), filter_peaks_valleys(dz_valleys)
    dz2_peaks_filtered, dz2_valleys_filtered = filter_peaks_valleys(dz2_peaks), filter_peaks_valleys(dz2_valleys)
    dz3_peaks_filtered, dz3_valleys_filtered = filter_peaks_valleys(dz3_peaks), filter_peaks_valleys(dz3_valleys)
    plot_peaks_valleys(ax3, [dz_peaks_filtered, dz2_peaks_filtered, dz3_peaks_filtered],
                       [dz_valleys_filtered, dz2_valleys_filtered, dz3_valleys_filtered],
                       [colors_derivatives["dz"], colors_derivatives["dz2"],colors_derivatives["dz3"]],
                       [colors_derivatives["dz"], colors_derivatives["dz2"],colors_derivatives["dz3"]])
    ax3.plot(vorticity.time, z, c="k", linewidth=2)
    ax3.set_title('Filter peaks closer than 1 day\n and transpose to ζ series')

    # 4) Plot mature stages
    ax4 = fig.add_subplot(gs[1, 0],frameon=True)  
    plot_stage(ax4, "mature")
    plot_peaks_valleys(ax4, dz_peaks, dz_valleys,
                       colors_phases["mature"], colors_phases["mature"])

    # Plot intensification stages
    ax5 = fig.add_subplot(gs[1, 1],frameon=True)
    plot_stage(ax5, "intensification")
    plot_peaks_valleys(ax5, peaks=[], valleys=[dz_valleys, dz2_valleys], peaks_color=[],
                         valleys_color=[colors_derivatives["dz"], colors_derivatives["dz2"]])
    
    # Plot decaying stages
    ax6 = fig.add_subplot(gs[1, 2],frameon=True)
    plot_stage(ax6, "decay")    
    plot_peaks_valleys(ax6, dz_peaks, dz2_valleys,
                        colors_derivatives["dz"], colors_derivatives["dz2"])

    # Plot incipient stage, if there's any
    ax7 = fig.add_subplot(gs[2, 0], frameon=True)
    if 'incipient' in periods.index:
        plot_stage(ax7, "incipient")
        plot_peaks_valleys(ax7, dz2_peaks, dz2_valleys,
                        colors_derivatives["dz2"], colors_derivatives["dz2"])
    else:
        ax7.plot(vorticity.time, z, c="k", linewidth=2)
        ax7.set_title("Incipient stage not found")
    
    # Plot everything together
    ax8 = fig.add_subplot(gs[2, 1:3],frameon=True)
    
    for _, row in periods.iterrows():
        period_start = row['start']
        period_end = row['end']
        period_data = z.sel(time=slice(period_start, period_end))

        # Extract the base phase name using regular expression
        phase_name = re.search(r'^(\w+)', row.name).group(1)

        ax8.fill_betweenx((z.min(), z.max()), period_data.time.min(), period_data.time.max(),
                        facecolor=colors_phases[phase_name], alpha=0.6)

    ax8.plot(vorticity.time, z, c="k", linewidth=2)
    
    ax8.set_title('Add overlaps for six hours between periods as confidence interval')
    
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()

    outfile_name_didactic+='.png'
    plt.savefig(outfile_name_didactic,dpi=500)
    print(outfile_name_didactic,'saved')

def mature_stage(dz):

    # Mature stage will be defined as all points between a consecutive
    # minimum and maximum of the derivative of the vorticity.

    dz_peaks, dz_valleys = get_peaks_valleys(dz)
    dz_peaks_filtered = filter_peaks_valleys(dz_peaks)
    dz_valleys_filtered = filter_peaks_valleys(dz_valleys)
    
    peaks = dz_peaks_filtered.replace(dz_peaks_filtered.values,
                                       'peak')
    valleys = dz_valleys_filtered.replace(dz_valleys_filtered.values,
                                           'valley')
    dz_peaks_valleys = pd.concat([peaks,valleys]).sort_index()
    
    dt = dz.time[1] - dz.time[0]
    dt = pd.Timedelta(dt.values)
        
    mature = []
    for i in range(len(dz_peaks_valleys[:-1])):
        if (dz_peaks_valleys[i] == 'valley') and \
            (dz_peaks_valleys[i+1] == 'peak') :
            mature.append(pd.date_range(dz_peaks_valleys.index[i],
                        dz_peaks_valleys.index[i+1],
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
                # If there's only one valley for dz2, there's no intensification
                if len(dz2_valleys) <= 1:
                    intensification = []
                    return intensification, incipient
                else:
                    intensification_start = dz2_valleys.index[i-1]
        intensification_end = mature_start-dt
        intensification.append(pd.date_range(intensification_start,
                        intensification_end, 
                        freq=f'{int(dt.total_seconds() / 3600)} H'))   
    
    return intensification, incipient

def decaying_stage(dz, dz2):
    # Decaying stage will be defined as all points between a consecutive peak of dz and a valley of dz2.

    dz_peaks, _ = get_peaks_valleys(dz)
    _, dz2_valleys = get_peaks_valleys(dz2)
    dz_peaks_filtered = filter_peaks_valleys(dz_peaks)
    dz2_valleys_filtered = filter_peaks_valleys(dz2_valleys)

    peaks = dz_peaks_filtered.replace(dz_peaks_filtered.values, 'peak')
    valleys = dz2_valleys_filtered.replace(dz2_valleys_filtered.values, 'valley')
    dz_peaks_valleys = pd.concat([peaks, valleys]).sort_index()

    dt = dz.time[1] - dz.time[0]
    dt = pd.Timedelta(dt.values)

    decaying = []
    i = 0
    while i < len(dz_peaks_valleys) - 1:
        if dz_peaks_valleys[i] == 'peak':
            j = i + 1
            while j < len(dz_peaks_valleys) and dz_peaks_valleys[j] == 'valley':
                decaying.append(pd.date_range(dz_peaks_valleys.index[i], dz_peaks_valleys.index[j], freq=f'{int(dt.total_seconds() / 3600)} H'))
                j += 1
            i = j - 1
        i += 1

    return decaying


def get_formatted_phases(phases):
    new_phases = {}
   
    for phase_key in phases:
        phase_list = phases[phase_key]
       
        if len(phase_list) == 1:
            new_phases[phase_key] = phase_list[0]
        else:
            for idx, phase in enumerate(phase_list):
                new_key = f"{phase_key} {idx + 1}" if idx > 0 else phase_key
                new_phases[new_key] = phase

    return new_phases

def get_phases(vorticity, output_directory):
    
    z = vorticity.zeta_filt2
    dz = vorticity.dz_dt_filt2
    dz2 = vorticity.dz_dt2_filt2
    dz3 = vorticity.dz_dt3_filt2
    
    mature = mature_stage(dz)

    intensification, incipient = intensification_incipient_stages(mature, dz2)
    
    decaying = decaying_stage(dz, dz2)
    
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
    
    phases = get_formatted_phases(phases)
    
    # Extract the first and last elements from each list
    df_dict = {}
    for k, v in phases.items():
        if len(v) == 0:
            df_dict[k] = []
        else:
            df_dict[k] = [v[0], v[-1]]
            
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame.from_dict(df_dict, orient='index',
                                columns=['start', 'end'])
    
    df.to_csv(output_directory+'periods.csv')
    print(output_directory+'periods.csv created')
     
    return df

def plot_periods(da, periods, outfile_name):
    
    z = da.zeta
    
    plt.close('all')
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111,frameon=True)
    colors = {'incipient': '#134074', 'intensification': '#f7b538',
          'mature': '#d62828', 'decay': '#9aa981'}
    
    six_hours = pd.Timedelta('5h 30 min')
    
    periods = periods.dropna()
    # Iterate over the periods and plot each phase
    for i in range(len(periods)):
        phase = periods.iloc[i].name
        # Remove any variation from the phase name
        phase_name = phase.split()[0]
        # Get the color for the phase from the color dictionary
        color = colors.get(phase_name, '#000000')
        ax.fill_betweenx((z.min(),z.max()+1e-5),
                        periods.loc[phase].start+six_hours,
                        periods.loc[phase].end-six_hours,
                        label=phase, facecolor=color, alpha=0.4)
        
    ax.plot(da.zeta.time, da.zeta,c='#6b675b',
            linewidth=4,label=r'$ζ$', alpha=0.8)
    ax.plot(da.zeta_filt2.time, da.zeta_filt2,c='k',
            linewidth=6,label=r'$ζ_f$')
    
    y = np.arange(z.min(),z.max()+5e-6,1e-5)
    plt.ylim(y[0],y[-1]+5e-6)
    plt.xlim(da.zeta_filt.time[0].values, da.zeta_filt.time[-1].values)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.grid(linewidth=0.5, alpha=0.5)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()
    
    plt.tight_layout()    

    outname = outfile_name+'.png'
    plt.savefig(outname,dpi=500)
    print(outname,'saved')

import pandas as pd

def get_periods(track_file, output_directory):
    # Set the output file names
    periods_outfile_path = output_directory + 'periods'
    periods_didatic_outfile_path = output_directory + 'periods_didatic'

    # Read the track file and extract the vorticity data
    track = pd.read_csv(track_file, parse_dates=[0], delimiter=';', index_col=[0])
    zeta_df = pd.DataFrame(track['min_zeta_850'].rename('zeta'))        
    vorticity = array_vorticity(zeta_df)

    # Determine the periods
    periods = get_phases(vorticity, output_directory)

    # Create plots
    plot_periods(vorticity, periods, periods_outfile_path)
    plot_didactic(vorticity, periods, periods_didatic_outfile_path)
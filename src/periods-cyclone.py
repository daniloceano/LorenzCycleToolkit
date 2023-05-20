#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    periods-cyclone.py                                 :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: danilocs <danilocs@student.42.fr>          +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/05/19 19:06:47 by danilocs          #+#    #+#              #
#    Updated: 2023/05/19 19:06:47 by danilocs         ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import re

import xarray as xr
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import cmocean.cm as cmo

from scipy.signal import argrelextrema
from scipy.signal import savgol_filter 


def check_create_folder(DirName):
    if not os.path.exists(DirName):
                os.makedirs(DirName)
                print(DirName+' created')
    else:
        print(DirName+' directory exists')

def filter_var(variable):
    window_lenght = round(len(variable)/2)
    if (window_lenght % 2) == 0:
        window_lenght += 1
    return savgol_filter(variable, window_lenght, 3, mode="nearest")

def array_vorticity(df):
    """
    Calculate derivatives of the vorticity and filter the resulting series

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
    dz_dt3_filt2 = xr.DataArray(filter_var(dz_dt3_filt2), coords={'time':df.index})

    # Assign variables to xarray
    da = da.assign(variables={'dzfilt2_dt': dzfilt2_dt,
                              'dzfilt2_dt2': dzfilt2_dt2,
                              'dzfilt2_dt3': dzfilt2_dt3,
                              'dz_dt_filt2': dz_dt_filt2,
                              'dz_dt2_filt2': dz_dt2_filt2,
                              'dz_dt3_filt2': dz_dt3_filt2})

    return da
    

def find_peaks_valleys(series):
    """
    Find peaks and valleys in a pandas series

    Args:
    series: pandas Series

    Returns:
    result: pandas Series with nans and "peaks" or "valley" in their respective positions
    """
    # Extract the values of the series
    data = series.values

    # Find peaks and valleys
    peaks = argrelextrema(data, np.greater_equal)[0]
    valleys = argrelextrema(data, np.less_equal)[0]

    # Create a series of NaNs
    result = pd.Series(index=series.index, dtype=object)
    result[:] = np.nan

    # Label the peaks and valleys
    result.iloc[peaks] = 'peak'
    result.iloc[valleys] = 'valley'

    # Don't assume that first and last data are either a peak or a valley
    result.iloc[0] = np.nan
    result.iloc[-1] = np.nan

    return result

def filter_peaks_valleys(result):

    filtered_result = result.copy()
    peaks = result[result == 'peak'].index
    valleys = result[result == 'valley'].index

    for peak_idx in peaks:
        next_valley = valleys[valleys > peak_idx].min()
        if next_valley is not pd.NaT and (next_valley - peak_idx) <= pd.Timedelta(days=1):
            filtered_result[next_valley] = np.nan

    for valley_idx in valleys:
        previous_peak = peaks[peaks < valley_idx].max()
        if previous_peak is not pd.NaT and (valley_idx - previous_peak) <= pd.Timedelta(days=1):
            filtered_result[valley_idx] = np.nan

    return filtered_result

def find_mature_stage(df):

    # Find valleys and peaks indices
    valleys = df[df['dz_peaks_valleys'] == 'valley'].index
    peaks = df[df['dz_peaks_valleys'] == 'peak'].index

    # Fill periods between valleys and peaks with "mature"
    for valley, peak in zip(valleys, peaks):
        df.loc[valley:peak, 'periods'] = 'mature'

    return df

def find_intensification_period(df):
    # Find dz and dz2 valleys indices
    dz_valleys = df[df['dz_peaks_valleys'] == 'valley'].index
    dz2_valleys = df[df['dz2_peaks_valleys'] == 'valley'].index

    # Find the corresponding dz2 valley for each dz valley
    dz2_valley_indices = []
    for dz_valley in dz_valleys:
        idx = np.argmin(np.abs(dz2_valleys - dz_valley))
        dz2_valley_indices.append(idx)

    # Fill periods between dz2 and dz valleys with "intensification"
    for dz2_valley_idx, dz_valley in zip(dz2_valley_indices, dz_valleys):
        dz2_valley = dz2_valleys[dz2_valley_idx]
        df.loc[dz2_valley:dz_valley, 'periods'] = 'intensification'

    return df

def find_decay_period(df):
    # Find dz peaks and dz2 valleys indices
    dz_peaks = df[df['dz_peaks_valleys'] == 'peak'].index
    dz2_valleys = df[df['dz2_peaks_valleys'] == 'valley'].index

    # Find the corresponding dz2 valley for each dz peak
    dz2_valley_indices = []
    for dz_peak in dz_peaks:
        idx = np.argmin(np.abs(dz2_valleys - dz_peak))
        dz2_valley_indices.append(idx)

    # Fill periods between dz peak and dz2 valley with "decay"
    for dz2_valley_idx, dz_peak in zip(dz2_valley_indices, dz_peaks):
        dz2_valley = dz2_valleys[dz2_valley_idx]
        df.loc[dz_peak:dz2_valley, 'periods'] = 'decay'

    return df

def plot_phase(df, phase):
    # Create a copy of the DataFrame
    df_copy = df.copy()

    colors_phases = {'incipient': '#65a1e6', 'intensification': '#f7b538',
          'mature': '#d62828', 'decay': '#9aa981'}

    # Find the start and end indices of "mature" periods
    mature_starts = df_copy[(df_copy['periods'] == phase) &
                             (df_copy['periods'].shift(1) != phase)].index
    mature_ends = df_copy[(df_copy['periods'] == phase) &
                           (df_copy['periods'].shift(-1) != phase)].index

    # Iterate over the "mature" periods and fill the area
    for start, end in zip(mature_starts, mature_ends):
        plt.fill_between(df_copy.index, df_copy['z'], where=(df_copy.index >= start) &
                          (df_copy.index <= end), alpha=0.7, color=colors_phases[phase])

    # Plot the "z" series
    plt.plot(df_copy.index, df_copy['z'], c='k')

    # Set labels and title
    plt.xlabel('Time')
    plt.ylabel('z')
    plt.title(f'{phase} phase')

    # Show the plot
    plt.show()

def plot_series_with_peaks_valleys(df):
    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the "z" series
    ax.plot(df.index, df['z'], label='z', color='k')

    # Define the series and colors for plotting
    series_names = ['dz', 'dz2', 'dz3']
    series_colors = ['#d62828', '#f7b538', '#65a1e6']
    marker_sizes = [80, 60, 40]
    peaks_valleys_columns = ['dz_peaks_valleys', 'dz2_peaks_valleys', 'dz3_peaks_valleys']
    scaling_factors = [100, 1000, 10000]

    # Plot the series and their peaks/valleys
    for series_name, series_color, peaks_valleys_col, marker_size, scaling_factor in zip(series_names, series_colors,
                                                                                         peaks_valleys_columns,
                                                                                         marker_sizes,
                                                                                         scaling_factors):
        ax.plot(df.index, df[series_name] * scaling_factor, color=series_color, label=series_name)
        ax.scatter(df.index[df[peaks_valleys_col] == 'peak'], df[series_name][df[peaks_valleys_col] == 'peak'] * scaling_factor,
                   color=series_color, marker='o', s=marker_size)
        ax.scatter(df.index[df[peaks_valleys_col] == 'valley'], df[series_name][df[peaks_valleys_col] == 'valley'] * scaling_factor,
                   color=series_color, marker='o', facecolors='none', s=marker_size)

    # Set labels and title
    ax.set_xlabel('Time')
    ax.set_ylabel('Series')
    ax.set_title('Series with Peaks and Valleys')

    # Move the legend outside the figure and to the top
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4)

    # Show the plot
    plt.show()

def plot_peaks_valleys_series(series, *peaks_valleys_series_list):
    # Plot the series
    plt.figure(figsize=(10, 6))
    plt.plot(series, color='k')

    # Plot peaks and valleys
    colors = ['#d62828', '#f7b538', '#65a1e6']  # List of colors for differentiating multiple series
    markers = ['o', 'o']  # List of markers for peaks and valleys
    labels = ['Peaks', 'Valleys']  # List of labels for the legend

    for i, peaks_valleys_series in enumerate(peaks_valleys_series_list):
        mask_notna = peaks_valleys_series.notna()
        mask_peaks = peaks_valleys_series == 'peak'

        # Calculate decreasing marker size
        marker_size = 300 - (i * 70)

        # Plot peaks
        plt.scatter(series.index[mask_notna & mask_peaks], series[mask_notna & mask_peaks],
                    marker=markers[0], color=colors[i], s=marker_size)

        # Plot valleys
        plt.scatter(series.index[mask_notna & ~mask_peaks], series[mask_notna & ~mask_peaks],
                    marker=markers[1], facecolors='none', edgecolors=colors[i], s=marker_size)

    plt.title('Series with Peaks and Valleys')
    plt.show()

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

    df = z.to_dataframe().rename(columns={'zeta_filt2':'z'})
    df['dz'] = dz.to_dataframe()
    df['dz2'] = dz2.to_dataframe()
    df['dz3'] = dz3.to_dataframe()

    df['dz_peaks_valleys'] = find_peaks_valleys(df['dz'])
    df['dz2_peaks_valleys'] = find_peaks_valleys(df['dz2'])
    df['dz3_peaks_valleys'] = find_peaks_valleys(df['dz3'])

    # Initialize periods column as NaN
    df['periods'] = np.nan

    # First step: identify peaks and valleys of vorticity
    # plot_series_with_peaks_valleys(df)

    # Second step: look for patterns
    # plot_peaks_valleys_series(df['z'], df['dz_peaks_valleys'], df['dz2_peaks_valleys'], df['dz3_peaks_valleys'])
    
    # Mature phase: between consecutive valley and peak of dz
    df = find_mature_stage(df)
    # plot_phase(df, "mature")

    # Intensification phase: between consecutive valleys of dz2 and dz
    df = find_intensification_period(df)
    # plot_phase(df, "intensification")

    # Decay phase: between consecutive peak of dz and valley of dz2
    df = find_decay_period(df)
    plot_phase(df, "decay")

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
    # plot_periods(vorticity, periods, periods_outfile_path)
    # plot_didactic(vorticity, periods, periods_didatic_outfile_path)

# Testing #
if __name__ == "__main__":

    track_file = '../inputs/track-test-periods'
    output_directory = './'
    get_periods(track_file, output_directory)


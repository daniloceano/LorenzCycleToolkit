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
#    Updated: 2023/05/22 13:26:04 by danilocs         ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import glob

import xarray as xr
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import matplotlib.ticker as ticker


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

    # # Don't assume that first and last data are either a peak or a valley
    # result.iloc[0] = np.nan
    # result.iloc[-1] = np.nan

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

import pandas as pd
import numpy as np

def find_mature_stage(df):

    valleys = df[df['dz_peaks_valleys'] == 'valley'].index
    peaks = df[df['dz_peaks_valleys'] == 'peak'].index

    # Iterate over valleys and find corresponding peaks
    for i in range(len(valleys)):
        valley = valleys[i]

        # Find the peaks that occur after the current valley
        corresponding_peaks = peaks[peaks > valley]

        if len(corresponding_peaks) > 0:
            # Find the first peak that occurs after the current valley
            corresponding_peak = corresponding_peaks[0]

            # Mature stage will only be defined if the valley starts at negative values
            if df.loc[valley, 'dz'] < 0:

                duration = corresponding_peak - valley

                # Mature stage needs to be at least 1 day long
                if duration >= pd.Timedelta(hours=24):

                    df.loc[valley:corresponding_peak, 'periods'] = 'mature'

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

        if (df.loc[dz2_valley, 'dz2'] < 0) and (df.loc[dz_valley, 'dz'] < 0):

            df.loc[dz2_valley:dz_valley, 'periods'] = 'intensification'

    return df

def find_decay_period(df):
    # Find dz peaks and dz2 valleys indices
    dz_peaks = df[df['dz_peaks_valleys'] == 'peak'].index
    dz2_valleys = df[df['dz2_peaks_valleys'] == 'valley'].index

    # Iterate over dz peaks
    for dz_peak in dz_peaks:
        # Find dz2 valleys with index greater than dz_peak
        valid_dz2_valleys = dz2_valleys[dz2_valleys > dz_peak]

        if len(valid_dz2_valleys) > 0:
            # Find dz2 valleys that occurs after dz_peaks
            for valid_dz2_valley in valid_dz2_valleys:

                # Check if there are no dz peaks between dz_peak and dz2_valley
                # dz_peak needs to be positive (in the mean) to be considered decay period
                if (not any(df['dz_peaks_valleys'].loc[dz_peak:valid_dz2_valley][1:] == 'peak')
                    ) and (df.loc[dz_peak:valid_dz2_valley, 'dz'].mean() > 0):
                    # Fill the period between dz_peak and valid_dz2_valley with 'decay'
                    df.loc[dz_peak:valid_dz2_valley, 'periods'] = 'decay'
            
    return df

def find_incipient_period(df):
    df['periods'].fillna('incipient', inplace=True)
    return df

import pandas as pd

def clean_periods(df):

    df['time'] = df.index

    # Calculate the duration of each period
    df['period_duration'] = df.groupby((df['periods'] != df['periods'].shift()).cumsum())['time'].transform(lambda x: x.max() - x.min())

    # Filter out periods with duration <= 6 hours
    df = df[df['period_duration'] > pd.Timedelta(hours=6)].copy()

    # Remove the extra columns
    df.drop('period_duration', axis=1, inplace=True)
    df.drop('time', axis=1, inplace=True)

    return df

def periods_to_dict(df):
    periods_dict = {}

    # Find the start and end indices of each period
    period_starts = df[df['periods'] != df['periods'].shift()].index
    period_ends = df[df['periods'] != df['periods'].shift(-1)].index

    # Iterate over the periods and create keys in the dictionary
    for i in range(len(period_starts)):
        period_name = df.loc[period_starts[i], 'periods']
        start = period_starts[i]
        end = period_ends[i]

        # Check if the period name already exists in the dictionary
        if period_name in periods_dict:
            # Append a suffix to the period name
            suffix = len(periods_dict[period_name]) + 1
            new_period_name = f"{period_name} {suffix}"
            periods_dict[new_period_name] = (start, end)
        else:
            periods_dict[period_name] = (start, end)

    return periods_dict


def plot_phase(df, phase, ax=None, show_title=True):
    # Create a copy of the DataFrame
    df_copy = df.copy()

    colors_phases = {'incipient': '#65a1e6', 'intensification': '#f7b538',
                     'mature': '#d62828', 'decay': '#9aa981'}

    # Find the start and end indices of the period
    phase_starts = df_copy[(df_copy['periods'] == phase) &
                            (df_copy['periods'].shift(1) != phase)].index
    phase_ends = df_copy[(df_copy['periods'] == phase) &
                          (df_copy['periods'].shift(-1) != phase)].index

    # Use the provided axes or create new ones
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    ax.axhline(0, c='gray', linewidth=0.5)

    # Iterate over the periods and fill the area
    for start, end in zip(phase_starts, phase_ends):
        ax.fill_between(df_copy.index, df_copy['z'], where=(df_copy.index >= start) &
                        (df_copy.index <= end), alpha=0.7, color=colors_phases[phase])

    ax.plot(df_copy.index, df_copy['z'], c='k')

    if show_title:
        title = ax.set_title(f'{phase} phase')
        title.set_position([0.55, 1.05])  # Adjust the title position as needed


    if ax is None:
        plt.show()

def plot_specific_peaks_valleys(df, key1, key2, ax):
    # Define the series and colors for plotting
    series_colors = {'dz':'#d62828', 'dz2':'#f7b538', 'dz3':'#65a1e6'}
    marker_sizes = {'dz': 250, 'dz2': 200, 'dz3': 150}

    for key in [key1, key2]:

        key_name =  key.split('_')[0]
        peak_or_valley = key.split('_')[1][:-1]

        peaks_valleys_series =  df[f"{key_name}_peaks_valleys"]

        color  = series_colors[key_name]
        marker_size = marker_sizes[key_name]

        # Plot the specified peaks or valleys
        if peak_or_valley == 'peak':
            ax.scatter(df.index[peaks_valleys_series == peak_or_valley],
                       df[peaks_valleys_series == 'peak'].z, color=color, marker='o', s=marker_size)
        elif peak_or_valley == 'valley':
            ax.scatter(df.index[peaks_valleys_series == peak_or_valley],
                       df[peaks_valleys_series == 'valley'].z, color=color, marker='o', s=marker_size,
                         facecolors='none', linewidth=2)

def plot_vorticity(ax, vorticity):

    ax.axhline(0, c='gray', linewidth=0.5)
    
    ax.plot(vorticity.time, vorticity.zeta, c='gray', linewidth=0.75, label='ζ')

    ax.plot(vorticity.time, vorticity.zeta_filt2, c='k', linewidth=2, label=r"$ζ_{filt}$")

    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))

    ax.legend(loc='best')

    ax.set_title("filter ζ")

def plot_series_with_peaks_valleys(df, ax):

    ax.axhline(0, c='gray', linewidth=0.5)

    ax.plot(df.index, df['z'], label='ζ', color='k')

    # Define the series and colors for plotting
    series_names = ['dz', 'dz2', 'dz3']
    labels = [r"$\frac{∂ζ_{filt}}{∂t} \times 10^{3}$",
               r"$\frac{∂^{2}ζ_{filt}}{∂t^{2}} \times 10^{4}$",
                 r"$\frac{∂^{3}ζ_{filt}}{∂t^{3}} \times 10^{5}$"]
    series_colors = ['#d62828', '#f7b538', '#65a1e6']
    marker_sizes = [80, 60, 40]
    peaks_valleys_columns = ['dz_peaks_valleys', 'dz2_peaks_valleys', 'dz3_peaks_valleys']
    scaling_factors = [100, 1000, 10000]

    # Plot the series and their peaks/valleys
    for series_name, label, series_color, peaks_valleys_col, marker_size, scaling_factor in zip(series_names,
                                                                                         labels,
                                                                                         series_colors,
                                                                                         peaks_valleys_columns,
                                                                                         marker_sizes,
                                                                                         scaling_factors):
        ax.plot(df.index, df[series_name] * scaling_factor, color=series_color, label=label)
        ax.scatter(df.index[df[peaks_valleys_col] == 'peak'], df[series_name][df[peaks_valleys_col] == 'peak'] * scaling_factor,
                   color=series_color, marker='o', s=marker_size)
        ax.scatter(df.index[df[peaks_valleys_col] == 'valley'], df[series_name][df[peaks_valleys_col] == 'valley'] * scaling_factor,
                   color=series_color, marker='o', facecolors='none', s=marker_size, linewidth=2)

    ax.set_title('derivate ζ')

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.5), ncol=4)

def plot_peaks_valleys_series(series, ax, *peaks_valleys_series_list):

    ax.axhline(0, c='gray', linewidth=0.5)
    
    # Plot the series
    ax.plot(series, color='k')

    # Plot peaks and valleys
    colors = ['#d62828', '#f7b538', '#65a1e6']  # List of colors for differentiating multiple series
    markers = ['o', 'o']  # List of markers for peaks and valleys

    for i, peaks_valleys_series in enumerate(peaks_valleys_series_list):
        mask_notna = peaks_valleys_series.notna()
        mask_peaks = peaks_valleys_series == 'peak'

        # Calculate decreasing marker size
        marker_size = 250 - (i * 50)

        # Plot peaks
        ax.scatter(series.index[mask_notna & mask_peaks], series[mask_notna & mask_peaks],
                    marker=markers[0], color=colors[i], s=marker_size)

        # Plot valleys
        ax.scatter(series.index[mask_notna & ~mask_peaks], series[mask_notna & ~mask_peaks],
                    marker=markers[1], facecolors='none', edgecolors=colors[i], s=marker_size, linewidth=2)

    ax.set_title('peaks and valleys')
    ax.title.set_position([0.55, 1.05])

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

def get_phases(vorticity, output_directory, i):
    
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

    # First step: filter vorticity data
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(331)
    plot_vorticity(ax1, vorticity)

    # Second step: identify peaks and valleys of vorticity
    ax2 = fig.add_subplot(332)
    plot_series_with_peaks_valleys(df, ax2)

    # Third step: look for patterns
    ax3 = fig.add_subplot(333)
    plot_peaks_valleys_series(
        df['z'], ax3,
        df['dz_peaks_valleys'],
        df['dz2_peaks_valleys'], 
        df['dz3_peaks_valleys']
        )
    
    # Mature phase: between consecutive valley and peak of dz
    df = find_mature_stage(df)
    ax4 = fig.add_subplot(334)
    plot_phase(df, "mature", ax4)
    plot_specific_peaks_valleys(df, "dz_valleys", "dz_peaks", ax4)

    # Intensification phase: between consecutive valleys of dz2 and dz
    df = find_intensification_period(df)
    ax5 = fig.add_subplot(335)
    plot_phase(df, "intensification", ax5)
    plot_specific_peaks_valleys(df, "dz_valleys", "dz2_valleys", ax5)

    # Decay phase: between consecutive peak of dz and valley of dz2
    df = find_decay_period(df)
    ax6 = fig.add_subplot(336)
    plot_phase(df, "decay", ax6)
    plot_specific_peaks_valleys(df, "dz_peaks", "dz2_valleys", ax6)

    # Put everything together
    ax7 = fig.add_subplot(337)
    plot_phase(df, "intensification", ax7, show_title=False)
    plot_phase(df, "mature", ax7, show_title=False)
    plot_phase(df, "decay", ax7, show_title=False)
    ax7.set_title("combine everythng")
    ax7.title.set_position([0.55, 1.05])

    # Incipient phase
    df = find_incipient_period(df)
    ax8 = fig.add_subplot(338)
    plot_phase(df, "intensification", ax8, show_title=False)
    plot_phase(df, "mature", ax8, show_title=False)
    plot_phase(df, "decay", ax8, show_title=False)
    plot_phase(df, "incipient", ax8, show_title=False)
    ax8.set_title("add incipient phase")
    ax8.title.set_position([0.55, 1.05])

    # Post processing of periods
    df = clean_periods(df)
    periods_dict = periods_to_dict(df)

    # Set y-axis labels in scientific notation (power notation) and change date format to "%m%d"
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
        date_format = mdates.DateFormatter("%m-%d")
        ax.xaxis.set_major_formatter(date_format)
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    plt.subplots_adjust(hspace=0.6)
    
    plt.savefig(f'{output_directory}/periods-test_{i}.png', dpi=500)
    
    # df.to_csv(output_directory+'periods.csv')
    # print(output_directory+'periods.csv created')
     
    return df


def get_periods(track_file, output_directory, i):
    # Set the output file names
    periods_outfile_path = output_directory + 'periods'
    periods_didatic_outfile_path = output_directory + 'periods_didatic'

    # Read the track file and extract the vorticity data
    track = pd.read_csv(track_file, parse_dates=[0], delimiter=';', index_col=[0])
    zeta_df = pd.DataFrame(track['min_zeta_850'].rename('zeta'))        
    vorticity = array_vorticity(zeta_df)

    # Determine the periods
    periods = get_phases(vorticity, output_directory, i)

    # Create plots
    # plot_periods(vorticity, periods, periods_outfile_path)
    # plot_didactic(vorticity, periods, periods_didatic_outfile_path)

# Testing #
if __name__ == "__main__":

    # track_file = '../inputs/track-test-periods'
    # output_directory = './'
    # get_periods(track_file, output_directory)

    files = glob.glob('../../SWSA-cyclones_energetic-analysis/LEC_results/*/*track')

    for track_file in files:

        filename = os.path.basename(track_file)
        file_id = filename.split('_')[0]
        
        # if file_id == '19840620':

        output_directory = '../../SWSA-cyclones_energetic-analysis/figures/periods_test'
        get_periods(track_file, output_directory, file_id)
        print(file_id)


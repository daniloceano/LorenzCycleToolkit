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
#    Updated: 2023/05/29 13:29:04 by danilocs         ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import csv

import xarray as xr
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
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

    # Filter derivatives
    dz_dt_filt2 = xr.DataArray(filter_var(dzfilt2_dt), coords={'time':df.index})
    dz_dt2_filt2 = xr.DataArray(filter_var(dzfilt2_dt2), coords={'time':df.index})

    # Assign variables to xarray
    da = da.assign(variables={'dzfilt2_dt': dzfilt2_dt,
                              'dzfilt2_dt2': dzfilt2_dt2,
                              'dz_dt_filt2': dz_dt_filt2,
                              'dz_dt2_filt2': dz_dt2_filt2})

    return da
    

def find_peaks_valleys(series):
    """
    Find peaks, valleys, and zero locations in a pandas series

    Args:
    series: pandas Series

    Returns:
    result: pandas Series with nans, "peak", "valley", and 0 in their respective positions
    """
    # Extract the values of the series
    data = series.values

    # Find peaks, valleys, and zero locations
    peaks = argrelextrema(data, np.greater_equal)[0]
    valleys = argrelextrema(data, np.less_equal)[0]
    zeros = np.where(data == 0)[0]

    # Create a series of NaNs
    result = pd.Series(index=series.index, dtype=object)
    result[:] = np.nan

    # Label the peaks, valleys, and zero locations
    result.iloc[peaks] = 'peak'
    result.iloc[valleys] = 'valley'
    result.iloc[zeros] = 0

    return result

import pandas as pd
import numpy as np

def find_mature_stage(df):
    dz_peaks = df[df['dz_peaks_valleys'] == 'peak'].index
    dz_valleys = df[df['dz_peaks_valleys'] == 'valley'].index
    z_valleys = df[df['z_peaks_valleys'] == 'valley'].index

    # Iterate over z valleys
    for z_valley in z_valleys:
        # Find the previous and next dz valleys relative to the current z valley
        previous_dz_valley = dz_valleys[dz_valleys < z_valley]

        # Check if there is a previous dz valley
        if len(previous_dz_valley) == 0:
            continue

        previous_dz_valley = previous_dz_valley[-1]

        # Find the next dz peak after the z valley
        next_dz_peak = dz_peaks[dz_peaks > z_valley]

        # Check if there is a next dz peak
        if len(next_dz_peak) == 0:
            continue

        next_dz_peak = next_dz_peak[0]

        # Calculate the distances between z valley and the previous/next dz valleys
        distance_to_previous_dz_valley = z_valley - previous_dz_valley
        distance_to_next_dz_peak = next_dz_peak - z_valley

        # Calculate the 3/4 distances from z valley to the previous/next dz valleys
        three_fourth_previous = z_valley - (1 / 4) * distance_to_previous_dz_valley
        three_fourth_next = z_valley + (1 / 4) * distance_to_next_dz_peak

        # Fill the period between three_fourth_previous and three_fourth_next with 'mature'
        df.loc[three_fourth_previous:three_fourth_next, 'periods'] = 'mature'

    return df

def find_intensification_period(df):
    # Find z and dz2 valleys indices
    z_valleys = df[df['z_peaks_valleys'] == 'valley'].index
    dz2_valleys = df[df['dz2_peaks_valleys'] == 'valley'].index
    dz_peaks = df[df['dz_peaks_valleys'] == 'peak'].index
    dz_valleys = df[df['dz_peaks_valleys'] == 'valley'].index

    # Remove z_valley if it is the last value of the series,
    #  here we assume that the cyclone cannot end intensifying
    z_valleys = z_valleys[z_valleys != df.index[-1]]

    # Filter dz_peaks that are at least 12 hours apart from dz_valleys
    # This is to avoid false interruption of intensification periods
    filtered_dz_peaks = dz_peaks.copy()
    for dz_peak in dz_peaks:
        time_diffs = abs(dz_peak - dz_valleys)
        closest_dz_valley_idx = time_diffs.argmin()
        closest_dz_valley = dz_valleys[closest_dz_valley_idx]
        if abs(dz_peak - closest_dz_valley) < pd.Timedelta(hours=12):
            filtered_dz_peaks = filtered_dz_peaks.drop(dz_peak)

    # Find the periods between dz2 and z valleys, excluding positive dz peaks,
    # negative dz valleys, and peaks more than 6 hours apart
    for dz2_valley in dz2_valleys:
        next_z_valley = z_valleys[z_valleys > dz2_valley].min()
        valid_dz_peaks = [peak for peak in filtered_dz_peaks if dz2_valley < peak < next_z_valley]
        if not valid_dz_peaks or (df.loc[valid_dz_peaks, 'dz'] > 0).all():
            if next_z_valley is not pd.NaT:
                intensification_end = next_z_valley
                df.loc[dz2_valley:intensification_end, 'periods'] = 'intensification'

    return df

def find_decay_period(df):
    # Find z and dz valleys and dz2 valleys indices
    z_valleys = df[df['z_peaks_valleys'] == 'valley'].index
    z_peaks = df[df['z_peaks_valleys'] == 'peak'].index
    dz2_valleys = df[df['dz2_peaks_valleys'] == 'valley'].index
    dz_valleys = df[df['dz_peaks_valleys'] == 'valley'].index

    # remove z_peaks that match the last index of the series
    z_peaks = z_peaks[:-1] if z_peaks[-1] == df.index[-1] else z_peaks

    for z_valley in z_valleys:
        next_z_peak = z_peaks[z_peaks > z_valley].min()
        if next_z_peak is not pd.NaT:
            valid_dz2_valleys = dz2_valleys[(dz2_valleys > z_valley) & (dz2_valleys < next_z_peak)]
        else:
            valid_dz2_valleys = dz2_valleys[dz2_valleys > z_valley]

        for valid_dz2_valley in valid_dz2_valleys:
                df.loc[z_valley:valid_dz2_valley, 'periods'].fillna('decay', inplace=True)

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

        # Add 6 hours to the beginning and to the end of the period as confidence intervals,
        # except for the beginning adn the end of the series
        if start != df.index[0]:
            start -= pd.Timedelta(hours=6)
        if end != df.index[0]:
            end += pd.Timedelta(hours=6)

        # Check if the period name already exists in the dictionary
        if period_name in periods_dict.keys():
            # Append a suffix to the period name
            suffix = len(periods_dict[period_name]) + 1 if len(periods_dict[period_name]) > 2 else 2
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

def plot_specific_peaks_valleys(df, ax, key1, key2, key3=None):
    # Define the series and colors for plotting
    series_colors = {'z': 'k', 'dz':'#d62828', 'dz2':'#f7b538'}
    marker_sizes = {'z': 190, 'dz': 120, 'dz2': 50}

    for key in [key1, key2, key3]:

        if key != None:
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

    # Define the series and colors for plotting
    series_names = ['z', 'dz', 'dz2']
    labels = ['ζ',
              r"$\frac{∂ζ_{filt}}{∂t} \times 10^{3}$",
               r"$\frac{∂^{2}ζ_{filt}}{∂t^{2}} \times 20^{4}$"]
    series_colors = ['k', '#d62828', '#f7b538']
    marker_sizes = [190, 120, 50]
    peaks_valleys_columns = ['z_peaks_valleys', 'dz_peaks_valleys', 'dz2_peaks_valleys']
    scaling_factors = [1, 100, 2000]

    # Plot the series and their peaks/valleys
    for series_name, label, series_color, peaks_valleys_col, marker_size, scaling_factor in zip(series_names,
                                                                                         labels,
                                                                                         series_colors,
                                                                                         peaks_valleys_columns,
                                                                                         marker_sizes,
                                                                                         scaling_factors):
        ax.plot(df.index, df[series_name] * scaling_factor, color=series_color, label=label)
        ax.scatter(df.index[df[peaks_valleys_col] == 'peak'],
                    df[series_name][df[peaks_valleys_col] == 'peak'] * scaling_factor,
                   color=series_color, marker='o', s=marker_size)
        ax.scatter(df.index[df[peaks_valleys_col] == 'valley'],
                    df[series_name][df[peaks_valleys_col] == 'valley'] * scaling_factor,
                   color=series_color, marker='o', facecolor='None', linewidth=2, s=marker_size)

    ax.set_title('derivate ζ')

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.5), ncol=4)

def plot_peaks_valleys_series(series, ax, *peaks_valleys_series_list):

    # Plot the series
    ax.plot(series, color='k')

    # Plot peaks and valleys
    colors = ['k', '#d62828', '#f7b538','#9aa981'] 
    marker_size = [190, 120, 50]
    zorder = [99, 100, 101]
    for i, peaks_valleys_series in enumerate(peaks_valleys_series_list):
        mask_notna = peaks_valleys_series.notna()
        mask_peaks = peaks_valleys_series == 'peak'

        # Plot peaks
        ax.scatter(series.index[mask_notna & mask_peaks], series[mask_notna & mask_peaks],
                    marker='o', color=colors[i], s=marker_size[i], zorder=zorder[i])

        # Plot valleys
        ax.scatter(series.index[mask_notna & ~mask_peaks], series[mask_notna & ~mask_peaks],
                    marker='o', edgecolors=colors[i], facecolors='None', s=marker_size[i],
                      linewidth=2,  zorder=zorder[i])
        
    # Plot vertical lines at valleys of z, valleys of dz, and peaks of dz
    valleys_z = peaks_valleys_series_list[0][peaks_valleys_series_list[0] == 'valley']
    valleys_dz = peaks_valleys_series_list[1][peaks_valleys_series_list[1] == 'valley']
    peaks_dz = peaks_valleys_series_list[1][peaks_valleys_series_list[1] == 'peak']

    for x in valleys_z.index:
        ax.axvline(x=x, color=colors[1], linestyle='-', linewidth=1.5, zorder=10)
    for x in valleys_dz.index:
        ax.axvline(x=x, color=colors[2], linestyle='-', linewidth=1.5, zorder=11)
    for x in peaks_dz.index:
        ax.axvline(x=x, color=colors[3], linestyle='-', linewidth=1.5, zorder=12)

    ax.set_title('find center of stages')
    ax.title.set_position([0.55, 1.05])


def plot_all_periods(phases_dict, df, ax=None, vorticity=None, periods_outfile_path=None):
    colors_phases = {'incipient': '#65a1e6', 'intensification': '#f7b538', 'mature': '#d62828', 'decay': '#9aa981'}

    # Create a new figure if ax is not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 5))

    if vorticity is not None:
        ax.plot(vorticity.time, vorticity, linewidth=0.75, color='gray', label=r'ζ')

    # Plot the vorticity data
    ax.plot(df.index, df['z'], c='k', label=r'$ζ_{f}$')

    legend_labels = set()  # To store unique legend labels

    # Shade the areas between the beginning and end of each period
    for phase, (start, end) in phases_dict.items():
        # Extract the base phase name (without suffix)
        base_phase = phase.split()[0]

        # Access the color based on the base phase name
        color = colors_phases[base_phase]

        # Fill between the start and end indices with the corresponding color
        ax.fill_between(df.index, df['z'], where=(df.index >= start) & (df.index <= end),
                        alpha=0.7, color=color, label=base_phase)

        # Add the base phase name to the legend labels set
        legend_labels.add(base_phase)

    # Add legend labels for Vorticity and ζ
    legend_labels.add('ζ')
    legend_labels.add(r'$ζ_{f}$')

    # Set the title
    ax.set_title('Vorticity Data with Periods')

    if periods_outfile_path is not None:
        # Remove duplicate labels from the legend
        handles, labels = ax.get_legend_handles_labels()
        unique_labels = []
        for label in labels:
            if label not in unique_labels and label in legend_labels:
                unique_labels.append(label)

        ax.legend(handles, unique_labels, loc='upper right', bbox_to_anchor=(1.5, 1))

        ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
        date_format = mdates.DateFormatter("%Y-%m-%d")
        ax.xaxis.set_major_formatter(date_format)
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        plt.tight_layout()

        fname = f"{periods_outfile_path}.png"
        plt.savefig(fname, dpi=500)
        print(f"{fname} created.")

def get_periods(vorticity):

    z = vorticity.zeta_filt2
    dz = vorticity.dz_dt_filt2
    dz2 = vorticity.dz_dt2_filt2

    df = z.to_dataframe().rename(columns={'zeta_filt2':'z'})
    df['dz'] = dz.to_dataframe()
    df['dz2'] = dz2.to_dataframe()

    df['z_peaks_valleys'] = find_peaks_valleys(df['z'])
    df['dz_peaks_valleys'] = find_peaks_valleys(df['dz'])
    df['dz2_peaks_valleys'] = find_peaks_valleys(df['dz2'])

    df['periods'] = np.nan

    # Intensification phase: between consecutive valleys of dz2 and zeta.
    df = find_intensification_period(df)

    # Decay phase: between consecutive valleys of z and dz2.
    # The mean dz value between consecutive peak of dz and valley of dz2 must be negtive
    df = find_decay_period(df)

    # Calculate the periods# Mature stage: between consecutive valleys and peaks of dz.
    # Valleys of dz must be negative and the phase needs to be at least 1 day long.
    df = find_mature_stage(df)

    # Incipient phase: all times not classified previously
    df = find_incipient_period(df)

    # Remove noisy periods that are less than 6 hours long
    df = clean_periods(df)

    # Pass the periods to a dictionary with each period's name as key
    #  and their corresponding start and end times as values.
    # Also, add extra 6 hours to the start and end of the periods as "confidence intervals"
    periods_dict = periods_to_dict(df)

    return periods_dict, df

def plot_didactic(periods_dict, df, vorticity, output_directory):
    
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
        df['z_peaks_valleys'],
        df['dz_peaks_valleys'],
        df['dz2_peaks_valleys'], 
        )
    
    # Intensification phase
    ax4 = fig.add_subplot(334)
    plot_phase(df, "intensification", ax4)
    plot_specific_peaks_valleys(df, ax4, "z_valleys", "dz_valleys", "dz2_valleys")

    # Decay phase
    ax5 = fig.add_subplot(335)
    plot_phase(df, "decay", ax5)
    plot_specific_peaks_valleys(df, ax5, "z_valleys", "dz_peaks", "dz2_valleys")

    # Mature phase
    ax6 = fig.add_subplot(336)
    plot_phase(df, "mature", ax6)
    plot_specific_peaks_valleys(df, ax6, "dz_valleys", "dz_peaks", )

    # Put everything together
    ax7 = fig.add_subplot(337)
    plot_phase(df, "intensification", ax7, show_title=False)
    plot_phase(df, "mature", ax7, show_title=False)
    plot_phase(df, "decay", ax7, show_title=False)
    ax7.set_title("combine everythng")
    ax7.title.set_position([0.55, 1.05])

    # Incipient phase
    ax8 = fig.add_subplot(338)
    plot_phase(df, "intensification", ax8, show_title=False)
    plot_phase(df, "mature", ax8, show_title=False)
    plot_phase(df, "decay", ax8, show_title=False)
    plot_phase(df, "incipient", ax8, show_title=False)
    ax8.set_title("add incipient phase")
    ax8.title.set_position([0.55, 1.05])

    # Post processing of periods
    ax9 = fig.add_subplot(339)
    plot_all_periods(periods_dict, df, ax9)
    ax9.set_title("add buffer times (6h)")
    ax9.title.set_position([0.55, 1.05])

    # Set y-axis labels in scientific notation (power notation) and change date format to "%m%d"
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
        date_format = mdates.DateFormatter("%m-%d")
        ax.xaxis.set_major_formatter(date_format)
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    plt.subplots_adjust(hspace=0.6)
    
    outfile = f'{output_directory}.png'
    plt.savefig(outfile, dpi=500)
    print(f"{outfile} created.")

def export_periods_to_csv(phases_dict, periods_outfile_path):

    filepath = f"{periods_outfile_path}.csv"

    # Extract phase names, start dates, and end dates from the periods dictionary
    data = [(phase, start, end) for phase, (start, end) in phases_dict.items()]
    
    # Write the data to a CSV file
    with open(filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['', 'start', 'end'])  # Write the header
        writer.writerows(data)  # Write the data rows

    print(f"{filepath} written.")

def determine_periods(track_file, output_directory):
    # Set the output file names
    periods_outfile_path = output_directory + 'periods'
    periods_didatic_outfile_path = output_directory + 'periods_didatic'

    # Read the track file and extract the vorticity data
    track = pd.read_csv(track_file, parse_dates=[0], delimiter=';', index_col=[0])
    zeta_df = pd.DataFrame(track['min_zeta_850'].rename('zeta'))        
    vorticity = array_vorticity(zeta_df)

    # Determine the periods
    periods_dict, df = get_periods(vorticity)

    # Create plots
    plot_all_periods(periods_dict, df, ax=None, vorticity=vorticity.zeta, periods_outfile_path=periods_outfile_path)
    plot_didactic(periods_dict, df, vorticity, periods_didatic_outfile_path)
    export_periods_to_csv(periods_dict, periods_outfile_path)
    

# Testing #
if __name__ == "__main__":

    track_file = '../inputs/track-test-periods'
    output_directory = './'
    determine_periods(track_file, output_directory)

# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    determine_periods.py                               :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo  <danilo.oceano@gmail.com>          +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/05/19 19:06:47 by danilocs          #+#    #+#              #
#    Updated: 2023/08/24 10:57:17 by Danilo           ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
Version: 1.2.0

This script processes vorticity data, identifies different phases of the cyclone     
and plots the identified periods on periods.png and periods_didatic.png   

The input track file should have the following columns:

    time;Lat;Lon;length;width;min_zeta_850;min_hgt_850;max_wind_850

The current version only works for negative vorciticity values
"""

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
from scipy.signal import convolve

from sklearn.preprocessing import MinMaxScaler

def check_create_folder(DirName, verbosity=False):
    if not os.path.exists(DirName):
                os.makedirs(DirName)
                print(DirName+' created')
    else:
        if verbosity:
            print(DirName+' directory exists')

def pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1.0, n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2.0 * np.pi * cutoff * k) / (np.pi * k)
    w[n - 1 : 0 : -1] = firstfactor * sigma
    w[n + 1 : -1] = firstfactor * sigma
    return w[1:-1]

def lanczos_filter(variable, window_length_lanczo, frequency):
    weights = pass_weights(window_length_lanczo, 1.0 / frequency)
    filtered_variable = convolve(variable, weights, mode="same")
    return filtered_variable

def pass_weights_bandpass(window, cutoff_low, cutoff_high):
    """Calculate weights for a bandpass Lanczos filter.

    Args:
    window: int
        The length of the filter window.

    cutoff_low: float
        The low cutoff frequency in inverse time steps.

    cutoff_high: float
        The high cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * (cutoff_high - cutoff_low)
    k = np.arange(1.0, n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = (
        np.sin(2.0 * np.pi * cutoff_high * k) / (np.pi * k)
        - np.sin(2.0 * np.pi * cutoff_low * k) / (np.pi * k)
    )
    w[n - 1 : 0 : -1] = firstfactor * sigma
    w[n + 1 : -1] = firstfactor * sigma
    return w[1:-1]

def lanczos_bandpass_filter(variable, window_length_lanczo, cutoff_low, cutoff_high):
    weights = pass_weights_bandpass(window_length_lanczo, cutoff_low, cutoff_high)
    filtered_variable = convolve(variable, weights, mode="same")
    return filtered_variable

def array_vorticity(zeta_df):
    """
    Calculate derivatives of the vorticity and filter the resulting series

    Args:
    df: pandas DataFrame

    Returns:
    xarray DataArray
    """

    # Parameters for filtering and smoothing processes
    frequency = 24.0
    savgol_polynomial = 3
    window_length_lanczo = len(zeta_df) // 2 
    window_length_savgol = len(zeta_df) | 1
    if pd.Timedelta(zeta_df.index[-1] - zeta_df.index[0]) > pd.Timedelta('8D'):
        window_length_savgol_2nd = window_length_savgol // 2 | 1
    else:
        window_length_savgol_2nd = window_length_savgol // 4 | 1
    if window_length_savgol_2nd < savgol_polynomial:
        window_length_savgol_2nd = 3
    cutoff_low = 1.0 / (7 * 24.0)
    cutoff_high = 1.0 / 48.0  # 24 hours
    
    # Convert dataframe to xarray
    da = zeta_df.to_xarray()

    # Apply Lanczos filter to vorticity 
    zeta_filtred = lanczos_bandpass_filter(da['zeta'].copy(), window_length_lanczo, cutoff_low, cutoff_high)
    zeta_filtred_low_pass = lanczos_filter(da.zeta.copy(), window_length_lanczo, frequency)
    zeta_filtred = xr.DataArray(zeta_filtred, coords={'time':zeta_df.index})
    da = da.assign(variables={'zeta_filt': zeta_filtred})

    num_samples = len(zeta_filtred)
    num_copy_samples = int(0.05 * num_samples)
    zeta_filtred.data[:num_copy_samples] = zeta_filtred_low_pass.data[:num_copy_samples]
    zeta_filtred.data[-num_copy_samples:] = zeta_filtred_low_pass.data[-num_copy_samples:]

    if pd.Timedelta(zeta_df.index[-1] - zeta_df.index[0]) > pd.Timedelta('8D'):
        window_length_savgol = window_length_savgol // 2 | 1
        
    # Smooth filtered vorticity with Savgol filter
    zeta_smoothed = xr.DataArray(
        savgol_filter(zeta_filtred, window_length_savgol//2|1, savgol_polynomial, mode="nearest"),
        coords={'time':zeta_df.index})

    zeta_smoothed2 = xr.DataArray(
            savgol_filter(zeta_smoothed, window_length_savgol_2nd, savgol_polynomial, mode="nearest"),
            coords={'time':zeta_df.index})
    
    da = da.assign(variables={'zeta_smoothed': zeta_smoothed})
    da = da.assign(variables={'zeta_smoothed2': zeta_smoothed2})

    # Calculate vorticity derivatives
    dz_dt = da.zeta.differentiate('time', datetime_unit='h')
    dz_dt2 = dz_dt.differentiate('time', datetime_unit='h')
    
    # Calculate the smoothed vorticity derivatives 
    dzfilt_dt = zeta_smoothed2.differentiate('time', datetime_unit='h')
    dzfilt_dt2 = dzfilt_dt.differentiate('time', datetime_unit='h')

    # Filter derivatives
    dz_dt_filt = xr.DataArray(
        savgol_filter(dzfilt_dt, window_length_savgol, savgol_polynomial, mode="nearest"),
        coords={'time':zeta_df.index})
    dz_dt2_filt = xr.DataArray(
        savgol_filter(dzfilt_dt2, window_length_savgol, savgol_polynomial, mode="nearest"),
        coords={'time':zeta_df.index})
    
    dz_dt_smoothed2 = xr.DataArray(
        savgol_filter(dz_dt_filt, window_length_savgol, savgol_polynomial, mode="nearest"),
        coords={'time':zeta_df.index})
    dz_dt2_smoothed2 = xr.DataArray(
        savgol_filter(dz_dt2_filt, window_length_savgol, savgol_polynomial, mode="nearest"),
        coords={'time':zeta_df.index})

    # Assign variables to xarray
    da = da.assign(variables={'dz_dt': dz_dt,
                              'dz_dt2': dz_dt2,
                              'dz_dt_filt': dz_dt_filt,
                              'dz_dt2_filt': dz_dt2_filt,
                              'dz_dt_smoothed2': dz_dt_smoothed2,
                              'dz_dt2_smoothed2': dz_dt2_smoothed2})

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

def find_mature_stage(df):
    dz_peaks = df[df['dz_peaks_valleys'] == 'peak'].index
    dz_valleys = df[df['dz_peaks_valleys'] == 'valley'].index
    z_valleys = df[df['z_peaks_valleys'] == 'valley'].index
    z_peaks = df[df['z_peaks_valleys'] == 'peak'].index

    series_length = df.index[-1] - df.index[0]
    dt = df.index[1] - df.index[0]

    # Iterate over z valleys
    for z_valley in z_valleys:
        # Find the previous and next dz valleys relative to the current z valley
        next_z_peak = z_peaks[z_peaks > z_valley]
        previous_z_peak =  z_peaks[z_peaks < z_valley]

        # Check if there is a previous or next z_peak
        if len(previous_z_peak) == 0 or len(next_z_peak) == 0:
            continue

        previous_z_peak = previous_z_peak[-1]
        next_z_peak = next_z_peak[0]

        # Calculate the distances between z valley and the previous/next dz valleys
        distance_to_previous_z_peak = z_valley - previous_z_peak
        distance_to_next_z_peak = next_z_peak - z_valley

        mature_distance_previous = 0.125 * distance_to_previous_z_peak
        mature_distance_next = 0.125 * distance_to_next_z_peak

        mature_start = z_valley - mature_distance_previous
        mature_end = z_valley + mature_distance_next

        # Mature stage needs to be at least 3% of total length
        mature_indexes = df.loc[mature_start:mature_end].index
        if mature_indexes[-1] - mature_indexes[0] > 0.03 * series_length:
            # Fill the period between mature_start and mature_end with 'mature'
            df.loc[mature_start:mature_end, 'periods'] = 'mature'

    # Check if all mature stages are preceeded by an intensification
    mature_periods = df[df['periods'] == 'mature'].index
    if len(mature_periods) > 0:
        blocks = np.split(mature_periods, np.where(np.diff(mature_periods) != dt)[0] + 1)
        for block in blocks:
            block_start, block_end = block[0], block[-1]
            if df.loc[block_start - dt, 'periods'] != 'intensification':
                df.loc[block_start:block_end, 'periods'] = np.nan

    return df


def find_intensification_period(df):
    # Find z peaks and valleys
    z_peaks = df[df['z_peaks_valleys'] == 'peak'].index
    z_valleys = df[df['z_peaks_valleys'] == 'valley'].index

    length = df.index[-1] - df.index[0]
    dt = df.index[1] - df.index[0]

    # Find intensification periods between z peaks and valleys
    for z_peak in z_peaks:
        next_z_valley = z_valleys[z_valleys > z_peak].min()
        if next_z_valley is not pd.NaT:
            intensification_start = z_peak
            intensification_end = next_z_valley

            # Intensification needs to be at least 7.5% of the total series length
            if intensification_end-intensification_start > length*0.12:
                df.loc[intensification_start:intensification_end, 'periods'] = 'intensification'
    
    # Check if there are multiple blocks of consecutive intensification periods
    intensefication_periods = df[df['periods'] == 'intensification'].index
    blocks = np.split(intensefication_periods, np.where(np.diff(intensefication_periods) != dt)[0] + 1)

    for i in range(len(blocks) - 1):
        block_end = blocks[i][-1]
        next_block_start = blocks[i+1][0]
        gap = next_block_start - block_end

        # If the gap between blocks is smaller than 7.5%, fill with intensification
        if gap < length*0.075:
            df.loc[block_end:next_block_start, 'periods'] = 'intensification'

    return df

def find_decay_period(df):
    # Find z peaks and valleys
    z_peaks = df[df['z_peaks_valleys'] == 'peak'].index
    z_valleys = df[df['z_peaks_valleys'] == 'valley'].index

    length = df.index[-1] - df.index[0]
    dt = df.index[1] - df.index[0]

    # Find decay periods between z valleys and peaks
    for z_valley in z_valleys:
        next_z_peak = z_peaks[z_peaks > z_valley].min()
        if next_z_peak is not pd.NaT:
            decay_start = z_valley
            decay_end = next_z_peak
        else:
            decay_start = z_valley
            decay_end = df.index[-1]  # Last index of the DataFrame

        # Decay needs to be at least 12% of the total series length
        if decay_end - decay_start > length*0.12:
            df.loc[decay_start:decay_end, 'periods'] = 'decay'

    # Check if there are multiple blocks of consecutive decay periods
    decay_periods = df[df['periods'] == 'decay'].index
    blocks = np.split(decay_periods, np.where(np.diff(decay_periods) != dt)[0] + 1)

    for i in range(len(blocks) - 1):
        block_end = blocks[i][-1]
        next_block_start = blocks[i+1][0]
        gap = next_block_start - block_end

        # If the gap between blocks is smaller than 7.5%, fill with decay
        if gap < length*0.075:
            df.loc[block_end:next_block_start, 'periods'] = 'decay'

    return df

def find_residual_period(df):
    unique_phases = [item for item in df['periods'].unique() if pd.notnull(item)]
    num_unique_phases = len(unique_phases)

    if num_unique_phases == 1:
        phase_to_fill = unique_phases[0]

        last_phase_index = df[df['periods'] == phase_to_fill].index[-1]
        dt = df.index[1] - df.index[0]
        df.loc[last_phase_index + dt:, 'periods'].fillna('residual', inplace=True)
    else:
        mature_periods = df[df['periods'] == 'mature'].index
        decay_periods = df[df['periods'] == 'decay'].index
        intensification_periods = df[df['periods'] == 'intensification'].index

        # Find residual periods where there is no decay stage after the mature stage
        for mature_period in mature_periods:
            if len(unique_phases) > 2:
                next_decay_period = decay_periods[decay_periods > mature_period].min()
                if next_decay_period is pd.NaT:
                    df.loc[mature_period:, 'periods'] = 'residual'
                    
        # Update mature periods
        mature_periods = df[df['periods'] == 'mature'].index

        # Fills with residual period intensification stage if there isn't a mature stage after it
        # but only if there's more than two periods
        if len(unique_phases) > 2:
            for intensification_period in intensification_periods:
                next_mature_period = mature_periods[mature_periods > intensification_period].min()
                if next_mature_period is pd.NaT:
                    df.loc[intensification_period:, 'periods'] = 'residual'

        # Fill NaNs after decay with residual if there is a decay, else, fill the NaNs after mature
        if 'decay' in unique_phases:
            last_decay_index = df[df['periods'] == 'decay'].index[-1]
        elif 'mature' in unique_phases:
            last_decay_index = df[df['periods'] == 'mature'].index[-1]
        dt = df.index[1] - df.index[0]
        df.loc[last_decay_index + dt:, 'periods'].fillna('residual', inplace=True)

    return df

def find_incipient_period(df):

    periods = df['periods']
    mature_periods = df[periods == 'mature'].index
    decay_periods = df[periods == 'decay'].index

    dt = df.index[1] - df.index[0]

    # if there's more than one period
    if len([item for item in df['periods'].unique() if (pd.notnull(item) and item != 'residual')]) > 2:
        # Find blocks of continuous indexes for 'decay' periods
        blocks = np.split(decay_periods, np.where(np.diff(decay_periods) != dt)[0] + 1)

        # Iterate over the blocks
        for block in blocks:
            if len(block) > 0:
                first_index = block[0]

                if first_index == df.index[0]:
                    df.loc[block, 'periods'] = 'incipient'

                else:
                    prev_index = first_index - dt
                    # Check if the previous index is incipient AND before mature stage
                    if (df.loc[prev_index, 'periods'] == 'incipient' or pd.isna(df.loc[prev_index, 'periods'])) and \
                    (len(mature_periods) > 0 and prev_index < mature_periods[-1]):
                        # Set the first period of the block to incipient
                        df.loc[block, 'periods'] = 'incipient'

    
    df['periods'].fillna('incipient', inplace=True)

    # If there's more than 2 unique phases other than residual and life cycle begins with
    # incipient fill first 6 hours with incipient.
    # If the life cycle begins with intensification, incipient phase will be from the
    # beginning of it, until 2/5 to the next dz_valley
    if len([item for item in df['periods'].unique() if (pd.notnull(item) and item != 'residual')]) > 2:
        phase_order = [item for item in df['periods'].unique() if pd.notnull(item)]
        if phase_order[0] in ['incipient', 'intensification'] or (phase_order[0] == 'incipient' and phase_order[1] == 'intensification'):
            start_time = df.iloc[0].name
            # Check if there's a dz valley before the next mature stage
            next_dz_valley = df[1:][df[1:]['dz_peaks_valleys'] == 'valley'].index.min()
            next_mature = df[df['periods'] == 'mature'].index.min()
            if next_dz_valley < next_mature:
                time_range = start_time + 2 * (next_dz_valley - start_time) / 5
                df.loc[start_time:time_range, 'periods'] = 'incipient'

    return df

import pandas as pd

def post_process_periods(df):
    dt = df.index[1] - df.index[0]
    
    # Find consecutive blocks of intensification and decay
    intensification_blocks = np.split(df[df['periods'] == 'intensification'].index, np.where(np.diff(df[df['periods'] == 'intensification'].index) != dt)[0] + 1)
    decay_blocks = np.split(df[df['periods'] == 'decay'].index, np.where(np.diff(df[df['periods'] == 'decay'].index) != dt)[0] + 1)
    
    # Fill NaN periods between consecutive intensification or decay blocks
    for blocks in [intensification_blocks, decay_blocks]:
        if len(blocks) > 1:
            phase = df.loc[blocks[0][0], 'periods']
            for i in range(len(blocks)):
                block = blocks[i]
                if i != 0:
                    if len(block) > 0:
                        last_index_prev_block = blocks[i -1][-1]
                        first_index_current_block = block[0]
                        preiods_between = df.loc[
                            (last_index_prev_block + dt):(first_index_current_block - dt)]['periods']
                        if all(pd.isna(preiods_between.unique())):
                            df.loc[preiods_between.index, 'periods'] = phase
    
    # Replace periods of length dt with previous or next phase
    for index in df.index:
        period = df.loc[index, 'periods']
        if pd.notna(period) and len(period) == dt:
            prev_index = index - dt
            next_index = index + dt
            if prev_index in df.index and prev_index != df.index[0]:
                df.loc[index, 'periods'] = df.loc[prev_index, 'periods']
            elif next_index in df.index:
                df.loc[index, 'periods'] = df.loc[next_index, 'periods']
    
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

    zeta = df_copy['z_unfil']
    zeta_smoothed = df_copy['z']

    colors_phases = {'incipient': '#65a1e6', 'intensification': '#f7b538',
                     'mature': '#d62828', 'decay': '#9aa981', 'residual': 'gray'}

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
        ax.fill_between(df_copy.index, zeta, where=(df_copy.index >= start) &
                        (df_copy.index <= end), alpha=0.7, color=colors_phases[phase])

    ax.plot(df_copy.index, zeta, c='gray', lw=3, alpha=0.8)
    ax2 = ax.twinx()
    ax2.axis('off')
    ax2.plot(df_copy.index, zeta_smoothed, c='k')

    if show_title:
        title = ax.set_title(f'{phase}', fontweight='bold', horizontalalignment='center')
        title.set_position([0.5, 1.05])  # Adjust the title position as needed


    if ax is None:
        plt.show()

def plot_specific_peaks_valleys(df, ax, *kwargs):
    # Define the series and colors for plotting
    series_colors = {'z': 'k', 'dz': '#d62828', 'dz2': '#f7b538'}
    marker_sizes = {'z': 190, 'dz': 120, 'dz2': 50}

    zeta = df['z']

    ax2 = ax.twinx()
    ax2.axis('off')

    for key in kwargs:
        key_name = key.split('_')[0]
        peak_or_valley = key.split('_')[1][:-1]

        peaks_valleys_series = df[f"{key_name}_peaks_valleys"]

        color = series_colors[key_name]
        marker_size = marker_sizes[key_name]
        zorder = 99 if key_name == 'z' else 100 if key_name == 'dz' else 101

        mask_notna = peaks_valleys_series.notna()
        mask_peaks = peaks_valleys_series == 'peak'

        # Plot peaks
        ax2.scatter(df.index[mask_notna & mask_peaks],
                   zeta[mask_notna & mask_peaks],
                   marker='o', color=color, s=marker_size, zorder=zorder)

        # Plot valleys
        ax2.scatter(df.index[mask_notna & ~mask_peaks],
                   zeta[mask_notna & ~mask_peaks],
                   marker='o', edgecolors=color, facecolors='none',
                   s=marker_size, linewidth=2, zorder=zorder)

def plot_vorticity(ax, vorticity):

    zeta = vorticity.zeta
    zeta_smoothed = vorticity.zeta_smoothed2

    ax.axhline(0, c='gray', linewidth=0.5)
    
    ax.plot(vorticity.time, zeta, c='gray', linewidth=0.75, label='ζ')

    ax2 = ax.twinx()
    ax2.axis('off')
    ax2.plot(vorticity.time, zeta_smoothed, c='k', linewidth=2, label=r"$ζ_{filt}$")

    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))

    ax.legend(loc='best')

    ax.set_title("filter ζ", fontweight='bold', horizontalalignment='center')

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

    ax.set_title('derivate ζ', fontweight='bold', horizontalalignment='center')

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

    ax.set_title('stages centers', fontweight='bold', horizontalalignment='center')
    ax.title.set_position([0.5, 1.05])


def plot_all_periods(phases_dict, df, ax=None, vorticity=None, periods_outfile_path=None):
    colors_phases = {'incipient': '#65a1e6',
                      'intensification': '#f7b538',
                        'mature': '#d62828',
                          'decay': '#9aa981',
                          'residual': 'gray'}

    # Create a new figure if ax is not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 5))

    if vorticity is not None:
        ax.plot(vorticity.time, vorticity.zeta, linewidth=2.5, color='gray', label=r'ζ')

        ax2 = ax.twinx()
        ax2.axis('off')
        ax2.plot(vorticity.time, vorticity.zeta_filt, linewidth=2, c='#d68c45', label=r'$ζ_{f}$')
        ax2.plot(vorticity.time, vorticity.zeta_smoothed, linewidth=2, c='#1d3557', label=r'$ζ_{fs}$')
        ax2.plot(vorticity.time, vorticity.zeta_smoothed2, linewidth=2, c='#e63946', label=r'$ζ_{fs^{2}}$')

    else:
        ax.plot(df.time, df.z, linewidth=0.75, color='gray', label=r'ζ')

    legend_labels = set()  # To store unique legend labels

    # Shade the areas between the beginning and end of each period
    for phase, (start, end) in phases_dict.items():
        # Extract the base phase name (without suffix)
        base_phase = phase.split()[0]

        # Access the color based on the base phase name
        color = colors_phases[base_phase]

        # Fill between the start and end indices with the corresponding color
        ax.fill_between(vorticity.time, vorticity.zeta.values, where=(vorticity.time >= start) & (vorticity.time <= end),
                        alpha=0.4, color=color, label=base_phase)

        # Add the base phase name to the legend labels set
        legend_labels.add(base_phase)

    # Add legend labels for Vorticity and ζ
    for label in [r'ζ', r'$ζ_{f}$', r'$ζ_{fs}$', r'$ζ_{fs^{2}}$']:
        legend_labels.add(label)

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

    z = vorticity.zeta_smoothed2
    dz = vorticity.dz_dt_smoothed2
    dz2 = vorticity.dz_dt2_smoothed2

    df = z.to_dataframe().rename(columns={'zeta_smoothed2':'z'})
    df['z_unfil'] = vorticity.zeta.to_dataframe()
    df['dz'] = dz.to_dataframe()
    df['dz2'] = dz2.to_dataframe()

    df['z_peaks_valleys'] = find_peaks_valleys(df['z'])
    df['dz_peaks_valleys'] = find_peaks_valleys(df['dz'])
    df['dz2_peaks_valleys'] = find_peaks_valleys(df['dz2'])

    df['periods'] = np.nan

    df = find_intensification_period(df)

    df = find_decay_period(df)

    df = find_mature_stage(df)

    df = find_residual_period(df)

    # 1) Fill consecutive intensification or decay periods that have NaNs between them
    # 2) Remove periods that are too short and fill with the previous period
    # (or the next one if there is no previous period)
    df = post_process_periods(df)

    df = find_incipient_period(df)

    # Pass the periods to a dictionary with each period's name as key
    #  and their corresponding start and end times as values.
    # Also, add extra 6 hours to the start and end of the periods as "confidence intervals"
    periods_dict = periods_to_dict(df)

    return periods_dict, df

def plot_didactic(df, vorticity, output_directory):
    
    # First step: filter vorticity data
    fig = plt.figure(figsize=(10, 8.5))
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
    df_int = find_intensification_period(df.copy())
    ax4 = fig.add_subplot(334)
    plot_phase(df_int, "intensification", ax4)
    plot_specific_peaks_valleys(df_int, ax4, "z_peaks", "z_valleys")

    # Decay phase
    df_decay = find_decay_period(df.copy())
    ax5 = fig.add_subplot(335)
    plot_phase(df_decay, "decay", ax5)
    plot_specific_peaks_valleys(df_decay, ax5, "z_peaks", "z_valleys")

    # Mature phase
    df_mature = find_mature_stage(df.copy())
    ax6 = fig.add_subplot(336)
    plot_phase(df_mature, "mature", ax6)
    plot_specific_peaks_valleys(df_mature, ax6, "z_peaks", "z_valleys", "dz_valleys", "dz_peaks")

    # Residual stage
    ax7 = fig.add_subplot(337)
    plot_phase(df, "residual", ax7)
    plot_specific_peaks_valleys(df_decay, ax7, "z_peaks", "z_valleys")

    # Incipient phase
    ax8 = fig.add_subplot(338)
    plot_phase(df, "incipient", ax8, show_title=False)
    ax8.set_title("incipient", fontweight='bold', horizontalalignment='center')
    ax8.title.set_position([0.5, 1.05])
    plot_specific_peaks_valleys(df_decay, ax8, "z_valleys", "dz_valleys")

    # Put everything together
    ax9 = fig.add_subplot(339)
    plot_phase(df, "incipient", ax9, show_title=False)
    plot_phase(df, "intensification", ax9, show_title=False)
    plot_phase(df, "mature", ax9, show_title=False)
    plot_phase(df, "decay", ax9, show_title=False)
    plot_phase(df, "residual", ax9)
    ax9.set_title("post processing", fontweight='bold', horizontalalignment='center')
    ax9.title.set_position([0.5, 1.05])


    # Set y-axis labels in scientific notation (power notation) and change date format to "%m%d"
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
        date_format = mdates.DateFormatter("%m-%d %HZ")
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
    vorticity = array_vorticity(zeta_df.copy())

    # Determine the periods
    periods_dict, df = get_periods(vorticity.copy())

    # Create plots
    plot_all_periods(periods_dict, df, ax=None, vorticity=vorticity, periods_outfile_path=periods_outfile_path)
    plot_didactic(df, vorticity, periods_didatic_outfile_path)
    export_periods_to_csv(periods_dict, periods_outfile_path)
    

# Testing #
if __name__ == "__main__":

    track_file = '../inputs/track-test-periods'
    output_directory = './'
    determine_periods(track_file, output_directory)

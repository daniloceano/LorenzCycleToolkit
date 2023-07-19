# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    determine_periodv2.py                              :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo <danilo.oceano@gmail.com>           +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/05/19 19:06:47 by danilocs          #+#    #+#              #
#    Updated: 2023/07/13 15:57:10 by Danilo           ###   ########.fr        #
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


def check_create_folder(DirName, verbosity=False):
    if not os.path.exists(DirName):
                os.makedirs(DirName)
                print(DirName+' created')
    else:
        if verbosity:
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

    series_length = df.index[-1] - df.index[0]
    dt = df.index[1] - df.index[0]

    # Iterate over z valleys
    for z_valley in z_valleys:
        # Find the previous and next dz valleys relative to the current z valley
        previous_dz_valley = dz_valleys[dz_valleys < z_valley]

        # Check if there is a previous dz valley
        if len(previous_dz_valley) == 0:
            continue

        previous_dz_valley = previous_dz_valley[-1]

        # Find the previous and next z peaks relative to the current z valley
        previous_z_peak = df[(df.index < z_valley) & (df['z_peaks_valleys'] == 'peak')].index.max()
        next_z_peak = df[(df.index > z_valley) & (df['z_peaks_valleys'] == 'peak')].index.min()

        # Calculate the distances between z valley and the previous/next dz valleys
        distance_to_previous_dz_valley = z_valley - previous_dz_valley
        distance_to_previous_z_peak = z_valley - previous_z_peak
        distance_to_next_z_peak = next_z_peak - z_valley

        # Check if either the previous or next z peak is within 7.5% of the length of z away from the valley
        if (
            distance_to_previous_z_peak <= 0.075 * series_length) or (
            distance_to_next_z_peak <= 0.075 * series_length):
            continue

        # Calculate the 3/4 distances from z valley to the previous/next dz valleys
        three_fourth_previous = z_valley - (1 / 4) * distance_to_previous_dz_valley
        three_fourth_next = z_valley + (1 / 4) * distance_to_next_z_peak

        # Mature stage needs to be at least 7% of total length
        mature_indexes = df.loc[three_fourth_previous:three_fourth_next].index
        if mature_indexes[-1] - mature_indexes[0] > 0.07*series_length:
            # Fill the period between three_fourth_previous and three_fourth_next with 'mature'
            df.loc[three_fourth_previous:three_fourth_next, 'periods'] = 'mature'

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

        # If the gap between blocks is smaller than 7.5%, fill with decay
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
    mature_periods = df[df['periods'] == 'mature'].index
    decay_periods = df[df['periods'] == 'decay'].index
    intensification_periods = df[df['periods'] == 'intensification'].index

    # Find residual periods where there is no decay stage after the mature stage
    for mature_period in mature_periods:
        # Assume that there are more than two unique periods (excluding NaN)
       if len([item for item in df['periods'].unique() if pd.notnull(item)]) > 2:
            next_decay_period = decay_periods[decay_periods > mature_period].min()
            if next_decay_period is pd.NaT:
                df.loc[mature_period:, 'periods'] = 'residual'

    # Fills with residual period intensification stage if there isn't a mature stage after it
    # but only if there's more than two periods
    if len([item for item in df['periods'].unique() if (pd.notnull(item) and item != 'residual')]) > 2:
        mature_periods = df[df['periods'] == 'mature'].index
        for intensification_period in intensification_periods:
            next_mature_period = mature_periods[mature_periods > intensification_period].min()
            if next_mature_period is pd.NaT:
                df.loc[intensification_period:, 'periods'] = 'residual'

    # Fill NaNs after decay with residual if there is a decay, else, fill the NaNs after mature
    if len(df[df['periods'] == 'decay'].index) > 0:
        last_decay_index = df[df['periods'] == 'decay'].index[-1]
        dt = df.index[1] - df.index[0]
        df.loc[last_decay_index+dt:, 'periods'].fillna('residual', inplace=True)
    elif len(df[df['periods'] == 'mature'].index) > 0:
        last_decay_index = df[df['periods'] == 'mature'].index[-1]
        dt = df.index[1] - df.index[0]
        df.loc[last_decay_index+dt:, 'periods'].fillna('residual', inplace=True)

    return df

def find_incipient_period(df):

    periods = df['periods']
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
                    # Check if the previous index is incipient o
                    if df.loc[prev_index, 'periods'] == 'incipient' or pd.isna(df.loc[prev_index, 'periods']):
                        # Set the first period of the block to incipient
                        df.loc[block, 'periods'] = 'incipient'
    
    df['periods'].fillna('incipient', inplace=True)

    return df

import pandas as pd

def clean_periods(df):

    df['time'] = df.index

    # Calculate the duration of each period
    df['period_duration'] = df.groupby((df['periods'] != df['periods'].shift()).cumsum())['time'].transform(lambda x: x.max() - x.min())

    # Filter out periods with duration <= 3 hours
    df = df[df['period_duration'] > pd.Timedelta(hours=3)].copy()

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

        # # Add 6 hours to the beginning and to the end of the period as confidence intervals,
        # # except for the beginning and the end of the series
        # if start != df.index[0]:
        #     start -= pd.Timedelta(hours=6)
        # if end != df.index[0]:
        #     end += pd.Timedelta(hours=6)

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
    colors_phases = {'incipient': '#65a1e6',
                      'intensification': '#f7b538',
                        'mature': '#d62828',
                          'decay': '#9aa981',
                          'residual': 'gray'}

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

    df = find_intensification_period(df)

    df = find_decay_period(df)

    df = find_mature_stage(df)

    df = find_residual_period(df)

    df = find_incipient_period(df)

    # Remove noisy periods that are less than 3 hours long
    df = clean_periods(df)

    # Pass the periods to a dictionary with each period's name as key
    #  and their corresponding start and end times as values.
    # Also, add extra 6 hours to the start and end of the periods as "confidence intervals"
    periods_dict = periods_to_dict(df)

    return periods_dict, df

def plot_didactic(df, vorticity, output_directory):
    
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
    df_int = find_intensification_period(df.copy())
    ax4 = fig.add_subplot(334)
    plot_phase(df_int, "intensification", ax4)
    plot_specific_peaks_valleys(df_int, ax4, "z_valleys", "dz_valleys", "dz2_valleys")

    # Decay phase
    df_decay = find_decay_period(df.copy())
    ax5 = fig.add_subplot(335)
    plot_phase(df_decay, "decay", ax5)
    plot_specific_peaks_valleys(df_decay, ax5, "z_valleys", "dz_peaks", "dz2_valleys")

    # Mature phase
    df_mature = find_mature_stage(df.copy())
    ax6 = fig.add_subplot(336)
    plot_phase(df_mature, "mature", ax6)
    plot_specific_peaks_valleys(df_mature, ax6, "dz_valleys", "dz_peaks")

    # Residual stage
    ax7 = fig.add_subplot(337)
    plot_phase(df, "residual", ax7)

    # Incipient phase
    ax8 = fig.add_subplot(338)
    plot_phase(df, "incipient", ax8, show_title=False)
    ax8.set_title("add incipient phase")
    ax8.title.set_position([0.55, 1.05])

    # Put everything together
    ax9 = fig.add_subplot(339)
    plot_phase(df, "incipient", ax9, show_title=False)
    plot_phase(df, "intensification", ax9, show_title=False)
    plot_phase(df, "mature", ax9, show_title=False)
    plot_phase(df, "decay", ax9, show_title=False)
    plot_phase(df, "residual", ax9)
    ax9.set_title("combine everythng")
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
    plot_didactic(df, vorticity, periods_didatic_outfile_path)
    export_periods_to_csv(periods_dict, periods_outfile_path)
    

# Testing #
if __name__ == "__main__":

    track_file = '../LEC_results-10MostIntense/10MostIntense-19920334_ERA5_track-15x15/10MostIntense-19920334_ERA5_track-15x15_track'
    output_directory = './'
    determine_periods(track_file, output_directory)

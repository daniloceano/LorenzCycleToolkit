#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:32:27 2022

@author: daniloceano
"""

import pandas as pd
import argparse
import glob
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import cmocean
from plot_timeseries import check_create_folder
from LPS import LorenzPhaseSpace

def create_LPS_plots(fig_title, zoom=False, **kwargs):
        plt.close('all')
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        LorenzPhaseSpace(ax, zoom=zoom, **kwargs)
        zoom_suffix = "_zoom" if zoom else ""
        fname = f"{ResultsSubDirectory}/Figures/LPS/LPS_{fig_title}{zoom_suffix}.png"
        with plt.rc_context({'savefig.dpi': 500}):
                plt.savefig(fname)
        print(f"{fname} created!")

def smooth_data(df, period):
        smoothed = df.groupby(pd.Grouper(key="Datetime", freq=period)).mean(numeric_only=True)
        # Set datetime to the date range
        starts = pd.Series(smoothed.index).dt.strftime('%Y-%m-%d %H:%M')
        ends = pd.Series(pd.DatetimeIndex(starts) + \
                        pd.Timedelta(hours=12)).dt.strftime('%Y-%m-%d %H:%M')
        smoothed['Datetime'] = pd.DataFrame(starts.astype(str)+' - '+\
                                        ends.astype(str)).values
        return smoothed

def period_data(df):
        periods_file = glob.glob(f"{ResultsSubDirectory}/periods.csv")[0]
        if not periods_file:
            raise FileNotFoundError("Periods file not found.")
        periods = pd.read_csv(periods_file, index_col=[0])
        periods = periods.dropna()
        for i in range(len(periods)):
                start,end = periods.iloc[i]['start'],periods.iloc[i]['end']
                selected_dates = df[(df['Datetime'] >= start) & (df['Datetime'] <= end)]
                if i == 0:
                        period = selected_dates.drop(['Datetime','Date','Hour'],axis=1).mean()
                        period = period.to_frame(name=periods.iloc[i].name).transpose()
                else:
                        tmp = selected_dates.drop(['Datetime','Date','Hour'],axis=1).mean()
                        tmp = tmp.to_frame(name=periods.iloc[i].name).transpose()
                        period = pd.concat([period,tmp]) 
        # Set datetime to the period date range
        period['Datetime'] = (periods['start'].astype(str)+' - '+\
                                                periods['end'].astype(str)).values
        return period 

    
if __name__ == "__main__":
 
    parser = argparse.ArgumentParser(description = "\
Lorenz Phase Space.")
    parser.add_argument("outfile", help = "The .csv file containing the \
results from the main.py program.")

    args = parser.parse_args()
    
    outfile = args.outfile

    ResultsSubDirectory = '/'.join(outfile.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/LPS/'
    check_create_folder(FigsDir)

    system = outfile.split('/')[-1].split('_')[0]
    datasource = outfile.split('/')[-1].split('_')[1]

    df = pd.read_csv(outfile, index_col=[0])
    df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')

    # Set datetime to the date range
    start = pd.to_datetime(df['Datetime'].iloc[0]).strftime('%Y-%m-%d %H:%M')
    end = pd.to_datetime(df['Datetime'].iloc[-1]).strftime('%Y-%m-%d %H:%M')


    # Plot example
    kwargs = {'terms':[], 'title':system, 'datasource': datasource, 'start': start, 'end': end}
    terms = {'Ca': df['Ca']*0, 'Ck': df['Ck']*0,  'Ge': df['Ge']*0, 'Ke': df['Ke']*0}
    kwargs['terms'].append(terms)      
    create_LPS_plots("example", zoom=False, **kwargs)

    for period in ['1H', '6H', '12H', '24H', '48H']:
                kwargs = {'terms':[], 'title':system, 'datasource': datasource, 'start': start, 'end': end}
                smoothed = smooth_data(df, period)
                terms = {'Ca': smoothed['Ca'], 'Ck': smoothed['Ck'],  'Ge': smoothed['Ge'], 'Ke': smoothed['Ke']}
                kwargs['terms'].append(terms)                      

                create_LPS_plots(f"{period}", zoom=False, **kwargs)
                create_LPS_plots(f"{period}", zoom=True, **kwargs)

    try:
        df_periods = period_data(df)
        kwargs = {'terms':[], 'title':system, 'datasource': datasource, 'start': start, 'end': end}
        terms = {'Ca': df_periods['Ca'], 'Ck': df_periods['Ck'],  'Ge': df_periods['Ge'], 'Ke': df_periods['Ke']}
        kwargs['terms'].append(terms) 
        create_LPS_plots("periods", zoom=False, **kwargs)
        create_LPS_plots("periods", zoom=True, **kwargs)
    except FileNotFoundError as e:
        print(e)
    except Exception as e:
        print("An error occurred:", str(e))

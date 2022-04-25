#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 17:56:57 2022

@author: danilocoutodsouza
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import os
import sys

# Specs for plotting
linecolors = ['#A53860','#C9B857','#384A0F','#473BF0']
markerfacecolors = ['#A53860','w','#384A0F','w']
conversion_labels = ['Cz','Ca','Ck','Ce']
markers = ['s','s','o','o']         
linestyles = ['-','-','-','-']
energy_labels = ['Az','Ae','Kz','Ke']
linewidth = 4

def plot_timeseries(df,DataDirectory):
    # Guarantee no plots are open
    plt.close('all')
    # Tines for x axs
    date = df['Date']
    times = pd.date_range(date[0],date.iloc[-1],periods=len(date))
    # Loop through the distinct group of terms
    for labels in [energy_labels,conversion_labels]:
        print('Printing '+str(labels)+'...')
        # Get values for setting plot range
        maxval = np.amax(np.amax(df[labels]))
        minval = np.amin(np.amin(df[labels]))
        print('Data range: '+str(minval)+' to '+str(maxval))
        # Create figure
        plt.figure(figsize=(8,8))
        # Loop trhough terms that are being plotted..
        for term,i in zip(labels,range(len(labels))):
            plt.grid(b=True,c='gray',linewidth=0.25,linestyle='dashdot')
            plt.plot(times,df[term],c=linecolors[i], marker= markers[i],
                    markerfacecolor=markerfacecolors[i], label=term,
                    linewidth=linewidth,markersize=6, linestyle=linestyles[i])
            plt.tick_params(axis='x', labelrotation=20)
            plt.legend()
            plt.xlim(times[0],times[-1])
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
        # Set x labels as dates
        ax = plt.gca()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        if term in energy_labels:
            plt.ylabel('Energy '+r' $(J\,m^{-2})$',fontsize=14)
            # Saving figure
            fname = DataDirectory+'/energy_terms.png'
            plt.savefig(fname)
            print(fname+ 'created')
        elif term in conversion_labels:
            # Horizontal line for 0
            plt.axhline(y = 0, color = 'k', linestyle = '-',
                        linewidth=1, zorder=1,alpha=0.8)
            plt.ylabel('Conversion '+r' $(W\,m^{-2})$',fontsize=14)
            # Saving figure
            fname = DataDirectory+'/conversion_terms.png'
            plt.savefig(fname)
            print(fname+ 'created')

def main():
    # Open data from energy and conversion terms as input from user from command line
    data = sys.argv[1]
    print('Reading data from: '+data)
    df = pd.read_csv(data)
    print('Ok!')
    # Diectory for saving figures
    FigsDirectory = 'LEC_Results'
    DataDirectory = FigsDirectory+'/'+data.split('/')[-1].split('.csv')[0]
    # Check if the LEC_Figures directory exists. If not, creates it
    if not os.path.exists(FigsDirectory):
                os.makedirs(FigsDirectory)
                print(FigsDirectory+' created')
    else:
        print(FigsDirectory+' directory exists')
    # Check if a directory for current data exists. If not, creates it
    if not os.path.exists(DataDirectory):
                os.makedirs(DataDirectory)
                print(DataDirectory+' created')
    # Make plot for timeseries
    plot_timeseries(df,DataDirectory)
    print('All done!')

if __name__ == "__main__":
    main()
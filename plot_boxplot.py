#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 14:22:40 2022

This script reads CSV files with energy and conversion terms from the Lorenz
Energy Cycle and make boxplots.


Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import os
import sys
import glob
from datetime import datetime

# Specs for plotting
linecolors = ['#A53860','#C9B857','#384A0F','#473BF0']
markerfacecolors = ['#A53860','w','#384A0F','w']
energy_labels = ['Az','Ae','Kz','Ke']
conversion_labels = ['Cz','Ca','Ck','Ce']
markers = ['s','s','o','o']         
linestyles = ['-','-','-','-']
linewidth = 4

Directory = sys.argv[1]

def get_data_dict(list_terms):
    data = {}
    # Loop to store results in dictionary
    for term in list_terms:
        file = glob.glob(Directory+'/'+term+'_*')[0]
        data[term] = pd.read_csv(file)
        print('\nOpening '+term)
        print(data[term][::data[term].shape[0]-1])
        print('Ok!')
    return data

def boxplot_time(term_list):
    data = get_data_dict(term_list)
    term = list(data.keys())[0]
    t_id = data[term].columns[0] # Time indexer
    times = data[term][t_id]
    plt.close('all') 
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(9, 9))
    i = 0
    for row in range(2):
        for col in range(2):
            ax = axs[row,col]
            term = term_list[i]
            ntime = data[term].shape[0]
            if col == 0 and term in energy_labels:
                ax.set_ylabel('Energy '+r' $(J\,m^{-2})$',fontsize=14)
            elif col == 0 and term in conversion_labels:
                ax.set_ylabel('Conversion '+r' $(W\,m^{-2})$',fontsize=14)
            for t in range(ntime):
                time_step = (datetime.fromisoformat(times.iloc[t]))
                bplot = ax.boxplot(data[term].iloc[t].values[1:],
                                     positions=[mdates.date2num(time_step)],
                                     patch_artist=True) 
                bplot['boxes'][-1].set_facecolor(linecolors[i])
                bplot['boxes'][-1].set_alpha(0.85)
                locator = mdates.AutoDateLocator(minticks=5, maxticks=12)
                formatter = mdates.AutoDateFormatter(locator)
                ax.xaxis.set_major_locator(locator)
                ax.xaxis.set_major_formatter(formatter)
                ax.tick_params(axis='x',labelrotation=20)
                ax.tick_params(size=12)
                ax.set_title(term, fontsize=14)
                fig.autofmt_xdate()
            i += 1
    if term in energy_labels:
            fname = 'Energy'
    elif term in conversion_labels:
            fname = 'Conversion'
    outfile = Directory+'/Figures/boxplot_vertical_timeseries_'+fname+'.png'
    plt.savefig(outfile)
    print('Created '+outfile)
    
def boxplot_vertical(term_list):
    data = get_data_dict(term_list)
    term = list(data.keys())[0]
    levs = data[term].columns[1:] # Vertical levels
    plt.close('all') 
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(9, 11))
    i = 0
    for row in range(2):
        for col in range(2):
            ax = axs[row,col]
            term = term_list[i]
            if col == 0 and term in energy_labels:
                ax.set_ylabel('Energy '+r' $(J\,m^{-2})$',fontsize=14)
            elif col == 0 and term in conversion_labels:
                ax.set_ylabel('Conversion '+r' $(W\,m^{-2})$',fontsize=14)
            for lev,j in zip(levs,range(len(levs))):
                bplot = ax.boxplot(data[term][lev].values,positions=[j/3],
                                   labels=[lev], patch_artist=True)
                bplot['boxes'][-1].set_facecolor(linecolors[i])
                bplot['boxes'][-1].set_alpha(0.85)
                ax.tick_params(axis='x',labelrotation=90)
                ax.tick_params(size=12)
                ax.set_title(term, fontsize=14)
            i +=1
    if term in energy_labels:
            fname = 'Energy'
    elif term in conversion_labels:
            fname = 'Conversion'
    outfile = Directory+'/Figures/boxplot_vertical_'+fname+'.png'
    plt.savefig(outfile)
    print('Created '+outfile)
    
def main():
    for term_list in [energy_labels,conversion_labels]:  
        print('\n-------------------------------------------------------------')
        print('Creating boxplot for the temporal evolution of each term')
        boxplot_time(term_list)
        print('\n-------------------------------------------------------------')
        print('Creating boxplot for each vertical level')
        boxplot_vertical(term_list)
        print('All done')

if __name__ == "__main__":
    main()


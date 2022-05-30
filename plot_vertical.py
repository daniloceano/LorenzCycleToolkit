#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:12:34 2022

This script reads CSV files with energy and conversion terms from the Lorenz
Energy Cycle, for each vertical level and time step, and plot them for each
time step. The user needs to specify the file where the CSV are located

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
import matplotlib.gridspec as gridspec
import os
import sys
import glob
from datetime import datetime
import cmocean as cmo


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

def create_dir():
    # Diectory for saving figures
    FigsDirectory = Directory+'/Figures_vertical_profiles/'
    # Check if the directory for saving figures exists. If not, creates it
    if not os.path.exists(FigsDirectory):
                os.makedirs(FigsDirectory)
                print(FigsDirectory+' created')
    else:
        print(FigsDirectory+' directory exists')
    return FigsDirectory

def plot_vertical(list_terms):
    # Get info
    data = get_data_dict(list_terms)
    term = list(data.keys())[0]
    ntime = data[term].shape[0] # Number of time steps
    levs = data[term].columns[1:] # Vertical levels
    # Create dir for saving figs
    FigsDirectory = create_dir()
    # Loop through time to plot each vertical profile
    for t in range(ntime):
        timestep = data[term].iloc[t].values[0]
        timestring = datetime.fromisoformat(timestep).strftime('%Y%m%d'+'Z'+'%H%M')
        plt.close('all') 
        plt.figure(figsize=(8,8))
        plt.grid(b=True,c='gray',linewidth=0.25,linestyle='dashdot')
        for term,i in zip(list_terms,range(4)):
            plt.plot(data[term].iloc[t].values[1:],levs,c=linecolors[i],
                     marker= markers[i],markerfacecolor=markerfacecolors[i],
                     label=term,linewidth=linewidth,markersize=6,
                     linestyle=linestyles[i])
        plt.grid(b=True,c='gray',linewidth=0.25,linestyle='dashdot')
        plt.tick_params(axis='x', labelrotation=20)
        plt.legend()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(timestep)
        plt.ylim(levs[0],levs[-1])
        plt.vlines(0, levs[0],levs[-1], color = 'k', linestyle = '-',
                            linewidth=1, zorder=1,alpha=0.8)
        # Save figure for each time step
        if term in energy_labels:
            fname = 'Energy'
        elif term in conversion_labels:
            fname = 'Conversion'
        outfile = FigsDirectory+fname+'_'+timestring+'.png'
        plt.savefig(outfile)
        print('Created '+outfile)
    
def plot_hovmoller(list_terms):
    data = get_data_dict(list_terms)
    term = list(data.keys())[0]
    TimeName = data[term].columns[0]
    dates = data[term][TimeName]
    times = pd.date_range(dates[0],dates.iloc[-1],periods=len(dates))
    levs = data[term].columns[1:]
    fig = plt.figure(figsize=(12, 10))
    gs = gridspec.GridSpec(nrows=2, ncols=2,  hspace=0.3)
    for i,term in zip(range(len(list_terms)),list_terms):
        ax = fig.add_subplot(gs[i])
        if term in energy_labels:
            fname = 'Energy'
            cmap='cmo.amp'
            title = 'Energy '+r' $(J\,m^{-2})$'
        elif term in conversion_labels:
            fname = 'Conversion'
            cmap='cmo.tarn'
            title = 'Conversion '+r' $(W\,m^{-2})$'
        cf = ax.contourf(times,levs,data[term][levs].transpose(),
                    cmap=cmap, extend='both')
        ax.contour(times,levs,data[term][levs].transpose(),
                    colors='k', extend='both')
        ax.tick_params(axis='x',labelrotation=20)
        ax.tick_params(size=12)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        start, end = ax.get_xlim()
        ax.set_title(term,fontdict={'fontsize':14})
    # colorbar
    cb_ax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(cf, cax=cb_ax,extend='both')
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(title, rotation=270,fontsize=10)
    outfile = Directory+'/Figures/hovmoller_'+fname+'.png'
    plt.savefig(outfile)
    print('Created '+outfile)
    
def main():
    for term_list in [energy_labels,conversion_labels]: 
        # This will procude a lot of figures!! Not needed for most of uses
        # print('\n-------------------------------------------------------------')
        # print('Creating figures with vertical profiles for each model time')
        # plot_vertical(term_list)
        print('\n-------------------------------------------------------------')
        print('Creating hovmoller diagrams')
        plot_hovmoller(term_list)

if __name__ == "__main__":
    main()



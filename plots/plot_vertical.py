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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import os
import cmocean as cmo
import glob
from datetime import datetime
import argparse

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
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
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
        plt.savefig(outfile,bbox_inches='tight')
        print('Created '+outfile)
    
def plot_hovmoller(list_terms):
    data = get_data_dict(list_terms)
    term = list(data.keys())[0]
    TimeName = data[term].columns[0]
    dates = data[term][TimeName]
    times = pd.date_range(dates[0],dates.iloc[-1],periods=len(dates))
    levs = data[term].columns[1:]
    if term in generation_labels:
        fig = plt.figure(figsize=(16, 5))
        gs = gridspec.GridSpec(nrows=1, ncols=2,
                               right=0.75, wspace=0.4)
    elif term in conversion_labels:
        fig = plt.figure(figsize=(12, 10))
        gs = gridspec.GridSpec(nrows=2, ncols=2,
                               wspace=0.5,hspace=0.3, top=0.95)
    elif term in energy_labels:
        fig = plt.figure(figsize=(12, 10))
        gs = gridspec.GridSpec(nrows=2, ncols=2, 
                               wspace=0.3,hspace=0.3, top=0.95)
    for i,term in zip(range(len(list_terms)),list_terms):
        ax = fig.add_subplot(gs[i])
        imin = data[term].min(numeric_only=True).min()
        imax = data[term].max(numeric_only=True).max()
        absmax = np.amax([np.abs(imin),imax])
        if term in energy_labels:
            fname = 'Energy'
            cmap='cmo.amp'
            title = 'Energy '+r' $(J\,m^{-2})$'
            norm = colors.Normalize(vmin=imin,vmax=imax)
        elif term in conversion_labels:
            fname = 'Conversion'
            cmap='cmo.tarn'
            norm = colors.TwoSlopeNorm(vcenter=0,vmin=-absmax,
                vmax=absmax)
            title = 'Conversion '+r' $(W\,m^{-2})$'
        elif term in generation_labels:
            fname = 'Generation'
            cmap='cmo.tarn'
            norm = colors.TwoSlopeNorm(vcenter=0,vmin=-absmax,
                vmax=absmax)
            title = 'Generation '+r' $(W\,m^{-2})$'
        cf = ax.contourf(times,levs,data[term][levs].transpose(),
                    cmap=cmap, extend='neither', norm=norm)
        ax.contour(times,levs,data[term][levs].transpose(),
                    colors='k', extend='neither')
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        start, end = ax.get_xlim()
        ax.set_title(term,fontdict={'fontsize':20})
        ax.tick_params(axis='x',labelrotation=50)
        # colorbar
        cbar = fig.colorbar(cf)
        if i % 2 != 0:
            cbar.ax.set_ylabel(title, rotation=270,fontsize=16)
            cbar.ax.get_yaxis().labelpad = 20
        for t in cbar.ax.get_yticklabels():
             t.set_fontsize(14)
    outfile = Directory+'/Figures/hovmoller_'+fname+'.png'
    plt.savefig(outfile,bbox_inches='tight')
    print('Created '+outfile)
    

def main():
    for term_list in [energy_labels,conversion_labels,generation_labels]: 
        if args.vertical:
            print('\n-------------------------------------------------------------')
            print('Creating figures with vertical profiles for each model time')
            plot_vertical(term_list)
        print('\n-------------------------------------------------------------')
        print('Creating hovmoller diagrams')
        plot_hovmoller(term_list)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "\
reads CSV files with energy and conversion terms from the Lorenz \
Energy Cycle, for each vertical level and time step, and plot them for each \
time step. The user needs to specify the file where the CSV are located")
    parser.add_argument("Directory", help = "Directory containing the\
 results from the main.py program.")
    parser.add_argument("-v", "--vertical", default = False, action='store_true',
    help = "If this flag is used, it will create figures containing vertical\
  profiles for each model time. This will procude a lot of figures and is not \
  needed for most of uses")
    args = parser.parse_args()
    Directory = args.Directory
        
    # Specs for plotting
    linecolors = ['#3B95BF','#87BF4B','#BFAB37','#BF3D3B']
    markerfacecolors = ['#59c0f0','#b0fa61','#f0d643','#f75452']   
    energy_labels = ['Az','Ae','Kz','Ke']
    conversion_labels = ['Cz','Ca','Ck','Ce']
    generation_labels = ['Gz','Ge']
    markers = ['s','s','o','o']         
    linestyles = ['-','-','-','-']
    linewidth = 4        
    
    main()



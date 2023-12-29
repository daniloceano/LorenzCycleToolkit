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
import glob
from datetime import datetime
import argparse
import matplotlib.gridspec as gridspec

def get_data_dict(term_list):
    data = {}
    # Loop to store results in dictionary
    for term in term_list:
        file = glob.glob(Directory+'/'+term+'_*')[0]
        data[term] = pd.read_csv(file)
        if args.verbosity:
            print('\nOpening '+term)
            print(data[term][::data[term].shape[0]-1])
            print('Ok!')
    return data

def boxplot_time(FigsSubDir):
    labels_list = [energy_labels,conversion_labels]
    for term_list in labels_list:
        data = get_data_dict(term_list)
        term = list(data.keys())[0]
        t_id = data[term].columns[0] # Time indexer
        times = data[term][t_id]
        plt.close('all') 
        fig = plt.figure(figsize=(12, 12))
        if term in energy_labels:
            fname = 'Energy'
            gs = gridspec.GridSpec(nrows=2, ncols=2,  hspace=0.15, wspace=0.2)
        elif term in conversion_labels:
            fname = 'Conversion'
            gs = gridspec.GridSpec(nrows=2, ncols=2,  hspace=0.15, wspace=0.3)
        i = 0
        for row in range(2):
            for col in range(2):
                ax = fig.add_subplot(gs[i])
                term = term_list[i]
                ntime = data[term].shape[0]
                if col == 0 and term in energy_labels:
                    ax.set_ylabel('Energy '+r' $(J\,m^{-2})$',fontsize=18)
                elif col == 0 and term in conversion_labels:
                    ax.set_ylabel('Conversion '+r' $(W\,m^{-2})$',fontsize=18)
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
                    ax.xaxis.set_tick_params(labelsize=16)
                    ax.yaxis.set_tick_params(labelsize=16)
                    ax.set_title(term, fontsize=20)
                    fig.autofmt_xdate()
                    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
                i += 1
            
        outfile = FigsSubDir+'/boxplot_vertical_timeseries_'+fname+'.png'
        plt.savefig(outfile,bbox_inches='tight')
        print('Created '+outfile)
    
def boxplot_vertical(FigsSubDir):
    labels_list = [energy_labels,conversion_labels]
    for term_list in labels_list:
        data = get_data_dict(term_list)
        term = list(data.keys())[0]
        levs = data[term].columns[1:] # Vertical levels
        plt.close('all') 
        fig = plt.figure(figsize=(13, 10))
        if term in energy_labels:
            fname = 'Energy'
            gs = gridspec.GridSpec(nrows=2, ncols=2,  hspace=0.15, wspace=0.2)
        elif term in conversion_labels:
            fname = 'Conversion'
            gs = gridspec.GridSpec(nrows=2, ncols=2,  hspace=0.15, wspace=0.3)
        i = 0
        for row in range(2):
            for col in range(2):
                ax = fig.add_subplot(gs[i])
                term = term_list[i]
                if col == 0 and term in energy_labels:
                    ax.set_ylabel('Energy '+r' $(J\,m^{-2})$',fontsize=18)
                elif col == 0 and term in conversion_labels:
                    ax.set_ylabel('Conversion '+r' $(W\,m^{-2})$',fontsize=18)
                for lev,j in zip(levs,range(len(levs))):
                    bplot = ax.boxplot(data[term][lev].values,positions=[j/3],
                                       labels=[lev], patch_artist=True)
                    bplot['boxes'][-1].set_facecolor(linecolors[i])
                    bplot['boxes'][-1].set_alpha(0.85)
                    ax.tick_params(axis='x',labelrotation=60)
                    if row == 0:
                        ax.axes.xaxis.set_ticklabels([])
                    ax.xaxis.set_tick_params(labelsize=16)
                    ax.yaxis.set_tick_params(labelsize=16)
                    ax.set_title(term, fontsize=20)
                i +=1
        outfile = FigsSubDir+'/boxplot_vertical_'+fname+'.png'
        plt.savefig(outfile,bbox_inches='tight')
        print('Created '+outfile)
    
def plot_boxplot(df,term_list,linecolor,fname,label,FigsSubDir):
    # Guarantee no plots are open
    plt.close('all')
    plt.figure(figsize=(8,8))
    for term,i in zip(term_list,range(len(term_list))):
        bplot = plt.boxplot(df[term],positions=[i/3],vert=True,
                            patch_artist=True,notch=True,labels=[term])
        bplot['boxes'][-1].set_facecolor(linecolor[i])
        bplot['boxes'][-1].set_alpha(0.7)
    # Horizontal line for 0
    if term not in energy_labels:
        plt.axhline(y = 0, color = 'k', linestyle = '-',
                    linewidth=1, zorder=1,alpha=0.8)
        plt.ylabel(label+r' $(J\,m^{-2})$',fontsize=14)
    else:
        plt.ylabel(label+r' $(W\,m^{-2})$',fontsize=14)   
    if term in budget_diff_labels:
        plt.xticks(rotation=25)
    # Saving figure
    plt.savefig(FigsSubDir+'boxplot_'+fname+'.png')
    print('boxplot_'+fname+'.png created')
    
def main():
    
    # Diectory for saving figures
    FigsDir = Directory+'/Figures/'
    FigsSubDir = FigsDir+'/boxplot/'
    check_create_folder(FigsSubDir)
    
    print()
    print('-------------------------------------------------------------')
    print('Creating boxplot for the temporal evolution of each term')
    boxplot_time(FigsSubDir)
    print('-------------------------------------------------------------')
    print('Creating boxplot for each vertical level')
    boxplot_vertical(FigsSubDir)
    print('All done')
        
    data = Directory+'/'.join(Directory.split('/')[2:-1])+'.csv'
    df = pd.read_csv(data)
    for term_list,cols,fname,label in zip(terms_list,cols_list,
                                      fnames,labels_list):
        plot_boxplot(df,term_list,cols,fname,label,FigsSubDir)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "\
Reads an .csv file with all terms from the Lorenz Energy Cycle and make boxplots.")
    parser.add_argument("Directory", help = "Directory containing the\
 results from the main.py program.")
    parser.add_argument("-r", "--residuals", default = False, action='store_true',
    help = "If this flag is used, it will plot RGe, RGz, RKz and RKe\
 instead of Ge, Gz, Dz and De")
    parser.add_argument("-v", "--verbosity", default = False,
                        action='store_true')
 
    args = parser.parse_args()
    Directory = args.Directory
    print(Directory)
    
    ## Specs for plotting ##
    # Labels for plots
    conversion_labels = ['Cz','Ca','Ck','Ce']
    energy_labels = ['Az','Ae','Kz','Ke']
    budget_diff_labels = ['∂Az/∂t (finite diff.)', '∂Ae/∂t (finite diff.)',
                     '∂Kz/∂t (finite diff.)', '∂Ke/∂t (finite diff.)']
    residuals_labels = ['RGz', 'RKz', 'RGe', 'RKe']
    gendiss_labels = ['Gz', 'Ge', 'Dz', 'De']
    residuals_labels = ['RGz', 'RGe', 'RKz', 'RKe']
    # This is for comparing terms estimated as residuals with terms computed  
    comparingG_labels = ['RGz', 'RGe', 'Gz', 'Ge']
    comparingD_labels = ['RKz', 'Dz', 'RKe', 'De']
    # Color for plots
    cols_energy = ['#3B95BF','#87BF4B','#BFAB37','#BF3D3B']
    linecolors =  ['#3B95BF','#87BF4B','#BFAB37','#BF3D3B','#873e23','#A13BF0']
    cols_conversion = ['#3B95BF','#87BF4B','#BFAB37','#BF3D3B']
    cols_residual =  ['#3B95BF','#87BF4B','#BFAB37','#BF3D3B']
    cols_boundary =   ['#3B95BF','#87BF4B','#BFAB37','#BF3D3B']
    # Specs for the markers and lines
    # markerfacecolors = ['#A53860','w','#384A0F','w','#873e23', 'w']
    markers = ['s','o','^','v','<','>']     
    markercolor =  ['#59c0f0','#b0fa61','#f0d643','#f75452','#f07243','#bc6ff7']   
    linestyles = ['-','-','-','-','-','-']
    linewidth = 3 
    
    if args.residuals:
        boundary_labels = ['BAz','BAe','BKz','BKe']
        terms_list = [energy_labels,conversion_labels,
                                boundary_labels,residuals_labels, 
                                budget_diff_labels,comparingG_labels]
        cols_list = [cols_energy,cols_conversion,cols_boundary,cols_residual,
                     cols_energy,linecolors]
        fnames = ['energy','conversion','boundary','residuals',
                  'energy_budget','comparing_GenRes']
        labels_list = ['Energy', 'Conversion', 'Transport across boundaries',
                       'Residuals', 
                       'Enery budgets (estimated using finite diffs.)',
                       'Generation/Residual']
    else:
        boundary_labels = ['BAz','BAe','BKz','BKe','BΦZ','BΦE']
        terms_list = [energy_labels,conversion_labels]
        term_labels = ['energy_labels','conversion_labels',
                               'boundary_labels','residuals_labels', 
                                'budget_diff_labels','gendiss_labels',
                                'comparingG_labels','comparingD_labels']     
    
    
    main()


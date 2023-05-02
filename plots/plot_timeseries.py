#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 17:56:57 2022

This script reads an CSV file with energy and conversion terms from the Lorenz
Energy Cycle (as input from user) and plot a timeseries for each.

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
import argparse
import os

def check_create_folder(DirName):
    
    if not os.path.exists(DirName):
                os.makedirs(DirName)
                print(DirName+' created')
    else:
        print(DirName+' directory exists')

def plot_timeseries(df, term_list, linecolor, fname, label, outdir):
 
    date = df['Date']
    times = pd.date_range(date[0],date.iloc[-1],periods=len(date))
    
    # Get values for setting plot range
    maxval = df[term_list].to_numpy().max()
    minval = df[term_list].to_numpy().min()
    if args.verbosity == True:
        print('Plotting '+str(term_list)+'...')
        print('Data range: '+str(minval)+' to '+str(maxval))
        
    plt.close('all')    
    plt.figure(figsize=(8,8))
    ax = plt.gca()
    
    # Loop trhough terms that are being plotted..
    for term,i in zip(term_list,range(len(term_list))):
        plt.plot(times,df[term],c=linecolor[i],label=term,
                 linewidth=linewidth,
                 marker=markers[i],markeredgecolor='#383838',
                 markerfacecolor=markercolor[i])
    plt.grid(c='gray',linewidth=0.25,linestyle='dashdot')
    plt.tick_params(axis='x', labelrotation=20)
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    plt.legend(prop={'size': 18})
    plt.xlim(times[0],times[-1])
    
    ax = plt.gca()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        
    if term not in energy_labels:
        plt.axhline(y = 0, color = 'k', linestyle = '-',
                    linewidth=1, zorder=1,alpha=0.8)
        plt.ylabel(label+r' $(J\,m^{-2})$',fontsize=18)
    else:
        plt.ylabel(label+r' $(W\,m^{-2})$',fontsize=18) 
        
    plt.savefig(outdir+'timeseires_'+fname+'.png')
    print('timeseires_'+fname+'.png created')

    
def plot_min_zeta_hgt(track,FigsDir):
    
    plt.close('all')
    fig, ax1 = plt.subplots(figsize=(15,10))
    
    lns1 = ax1.plot(track.index, track['min_zeta_850'],c='#554348', marker='o',
             label= '850 hPa minimum vorticity')
    ax2 = ax1.twinx()
    lns2 = ax2.plot(track.index, track['min_hgt_850'],c='#6610F2', marker='s',
             label= '850 hPa minimum geopotential height')
    # added these three lines
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc='best',prop={'size': 18})
    ax1.grid(c='gray',linewidth=0.25,linestyle='dashdot', axis='x')
    ax1.tick_params(axis='x', labelrotation=20)
    ax1.xaxis.set_tick_params(labelsize=16)
    ax1.yaxis.set_tick_params(labelsize=16)
    ax2.yaxis.set_tick_params(labelsize=16)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %HZ'))
    
    plt.savefig(FigsDir+'timeseries-min_zeta_hgt.png',bbox_inches='tight')
    print('Created:',FigsDir+'track_boxes.png')


def main():
    print('Reading data from: '+outfile)
    df = pd.read_csv(outfile)
    if args.verbosity:
        print(df)
        print(' ')
        
    # Diectory for saving figures
    ResultsSubDirectory = '/'.join(outfile.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/'
    FigsSubDir = FigsDir+'/timeseries/'
    check_create_folder(FigsSubDir)
    
    # Plot track
    trackfile = ''.join(outfile.split('.csv'))+'_track'
    try:
        track = pd.read_csv(trackfile,parse_dates=[0],delimiter=';',index_col='time')
        plot_min_zeta_hgt(track,FigsDir)
    except:
        print('Error on tracking. Possibly, track file was not found')
    
    # Make plot for timeseries
    for term_list, linecolor, fname, label in zip(terms_list, cols_list,
                                                          fnames,labels_list):
        plot_timeseries(df, term_list, linecolor, fname, label, FigsSubDir)
        
    print('All done!')

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "\
Reads an .csv file with all terms from the Lorenz Energy Cycle and plot a \
timeseries for each.")
    parser.add_argument("outfile", help = "The .cdv file containing the \
results from the main.py program.")
    parser.add_argument("-r", "--residuals", default = False, action='store_true',
    help = "If this flag is used, it will plot RGe, RGz, RKz and RKe\
 instead of Ge, Gz, Dz and De")
    parser.add_argument("-v", "--verbosity", default = False,
                        action='store_true')
    
    args = parser.parse_args()
    outfile = args.outfile

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
        terms_list = [energy_labels,conversion_labels,
                                boundary_labels,residuals_labels, 
                                budget_diff_labels,gendiss_labels,
                                comparingG_labels,comparingD_labels]
        term_labels = ['energy_labels','conversion_labels',
                               'boundary_labels','residuals_labels', 
                                'budget_diff_labels','gendiss_labels',
                                'comparingG_labels','comparingD_labels']
 
    main()

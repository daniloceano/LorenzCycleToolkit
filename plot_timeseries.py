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
import sys



def plot_timeseries(df,term_list,linecolor,fname,label,outdir):
    # Guarantee no plots are open
    plt.close('all')
    # Times for x axs
    date = df['Date']
    times = pd.date_range(date[0],date.iloc[-1],periods=len(date))
    # Loop through the distinct group of terms
    print('Plotting '+str(term_list)+'...')
    # Get values for setting plot range
    maxval = np.amax(np.amax(df[term_list]))
    minval = np.amin(np.amin(df[term_list]))
    print('Data range: '+str(minval)+' to '+str(maxval))
    # Create figure
    plt.figure(figsize=(8,8))
    # Loop trhough terms that are being plotted..
    for term,i in zip(term_list,range(len(term_list))):
        plt.plot(times,df[term],c=linecolor[i],label=term,
                 linewidth=linewidth)
    plt.grid(b=True,c='gray',linewidth=0.25,linestyle='dashdot')
    plt.tick_params(axis='x', labelrotation=20)
    plt.legend()
    plt.xlim(times[0],times[-1])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # Set x labels as dates
    ax = plt.gca()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Horizontal line for 0
    if term not in energy_labels:
        plt.axhline(y = 0, color = 'k', linestyle = '-',
                    linewidth=1, zorder=1,alpha=0.8)
        
    # Horizontal line for 0
    if term not in energy_labels:
        plt.axhline(y = 0, color = 'k', linestyle = '-',
                    linewidth=1, zorder=1,alpha=0.8)
        plt.ylabel(label+r' $(J\,m^{-2})$',fontsize=14)
    else:
        plt.ylabel(label+r' $(W\,m^{-2})$',fontsize=14)   
    # Saving figure
    plt.savefig(outdir+'timeseires_'+fname+'.png')
    print('timeseires_'+fname+'.png created')

            
def plot_boxplot(df,term_list,linecolor,fname,label,outdir):
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
    plt.savefig(outdir+'boxplot_'+fname+'.png')
    print('boxplot_'+fname+'.png created')


def main():
    print('Reading data from: '+data)
    df = pd.read_csv(data)
    print(df)
    print(' ')
    # Diectory for saving figures
    ResultsSubDirectory = '/'.join(data.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/'
    # Make plot for timeseries
    for term_list,cols,fname,label in zip(terms_list,cols_list,
                                      fnames,labels_list):
        plot_timeseries(df,term_list,cols,fname,label,FigsDir)
        plot_boxplot(df,term_list,cols,fname,label,FigsDir)
        
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
    
    args = parser.parse_args()
    data = args.outfile

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
    markers = ['s','s','o','o','^','^']         
    linestyles = ['-','-','-','-','-','-']
    linewidth = 5 
    
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

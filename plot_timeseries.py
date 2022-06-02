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
import os
import sys

# Specs for plotting

conversion_labels = ['Cz','Ca','Ck','Ce']
energy_labels = ['Az','Ae','Kz','Ke']
boundary_labels = ['BAz','BAe','BKz','BKe','BΦZ','BΦE']
budget_diff_labels = ['∂Az/∂t (finite diff.)', '∂Ae/∂t (finite diff.)',
                 '∂Kz/∂t (finite diff.)', '∂Ke/∂t (finite diff.)']
residuals_labels = ['RGz', 'RKz', 'RGe', 'RKe']
gendiss_labels = ['Gz', 'Ge', 'Dz', 'De']

# This is for comparing terms estimated as residuals with terms computed  
comparingG_labels = ['RGz', 'Gz', 'RGe', 'Ge']
comparingD_labels = ['RKz', 'Dz', 'RKe', 'De']

linecolors = ['#A53860','#C9B857','#384A0F','#473BF0','#873e23','#A13BF0']
markerfacecolors = ['#A53860','w','#384A0F','w','#873e23', 'w']
markers = ['s','s','o','o','^','^']         
linestyles = ['-','-','-','-','-','-']
linewidth = 4

def plot_timeseries(df,outdir):
    # Guarantee no plots are open
    plt.close('all')
    # Times for x axs
    date = df['Date']
    times = pd.date_range(date[0],date.iloc[-1],periods=len(date))
    # Loop through the distinct group of terms
    for labels,termlist in zip([energy_labels,conversion_labels,
                                boundary_labels,residuals_labels, 
                                budget_diff_labels,gendiss_labels,
                                comparingG_labels,comparingD_labels],
                               ['energy_labels','conversion_labels',
                               'boundary_labels','residuals_labels', 
                                'budget_diff_labels','gendiss_labels',
                                'comparingG_labels','comparingD_labels']):
        print('Plotting '+str(labels)+'...')
        # Get values for setting plot range
        maxval = np.amax(np.amax(df[labels]))
        minval = np.amin(np.amin(df[labels]))
        print('Data range: '+str(minval)+' to '+str(maxval))
        # Create figure
        plt.figure(figsize=(8,8))
        # Loop trhough terms that are being plotted..
        for term,i in zip(labels,range(len(labels))):
            plt.plot(times,df[term],c=linecolors[i], marker= markers[i],
                    markerfacecolor=markerfacecolors[i], label=term,
                    linewidth=linewidth,markersize=6, linestyle=linestyles[i])
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
        # name y axis and save figure
        if termlist == 'energy_labels':
            plt.ylabel('Energy '+r' $(J\,m^{-2})$',fontsize=14)
            fname = outdir+'/timeseires_energy_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'conversion_labels':
            plt.ylabel('Conversion '+r' $(W\,m^{-2})$',fontsize=14)
            fname = outdir+'/timeseires_conversion_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'boundary_labels':
            plt.ylabel('Transport across boundaries '+r' $(W\,m^{-2})$',fontsize=14)
            fname = outdir+'/timeseires_boundary_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'residuals_labels':
            plt.ylabel('Residuals '+r' $(W\,m^{-2})$',fontsize=14)
            fname = outdir+'/timeseires_residuals_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'budget_diff_labels':
            plt.ylabel('Enery budgets (estimated using finite diffs. '+r' $(W\,m^{-2})$',fontsize=14)
            fname = outdir+'/timeseires_budget_diff_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'comparingG_labels':
            plt.ylabel('Energy generation '+r' $(W\,m^{-2})$',fontsize=14)
            fname = outdir+'/timeseires_generation_terms_compare.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'comparingD_labels':
            plt.ylabel('Energy dissipation '+r' $(W\,m^{-2})$',fontsize=14)
            fname = outdir+'/timeseires_dissipation_terms_compare.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'gendiss_labels':
            plt.ylabel('Energy dissipation/generation '+r' $(W\,m^{-2})$',fontsize=14)
            fname = outdir+'/timeseires_gendiss_terms.png'
            plt.savefig(fname)
            print(fname+' created')
            
            
def plot_boxplot(df,outdir):
    # Guarantee no plots are open
    plt.close('all')
    for labels,termlist in zip([energy_labels,conversion_labels,
                                boundary_labels,residuals_labels, 
                                budget_diff_labels,gendiss_labels],
                               ['energy_labels','conversion_labels',
                               'boundary_labels','residuals_labels', 
                                'budget_diff_labels','gendiss_labels']):
        plt.figure(figsize=(8,8))
        plt.grid(visible=True,c='gray',linewidth=0.25,linestyle='dashdot')
        for term,i in zip(labels,range(len(labels))):
            bplot = plt.boxplot(df[term],positions=[i/3],vert=True,
                                patch_artist=True,notch=True,labels=[term])
            bplot['boxes'][-1].set_facecolor(linecolors[i])
            bplot['boxes'][-1].set_alpha(0.7)
        plt.legend()
        # Horizontal line for 0
        if term not in energy_labels:
            plt.axhline(y = 0, color = 'k', linestyle = '-',
                        linewidth=1, zorder=1,alpha=0.8)
        if termlist == 'energy_labels':
            plt.ylabel('Energy '+r' $(J\,m^{-2})$',fontsize=14)
            # Saving figure
            fname = outdir+'/boxplot_time_energy_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'conversion_labels':
            plt.ylabel('Conversion '+r' $(W\,m^{-2})$',fontsize=14)
            # Saving figure
            fname = outdir+'/boxplot_time_conversion_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'boundary_labels':
            plt.ylabel('Transport across boundaries '+r' $(W\,m^{-2})$',fontsize=14)
            # Saving figure
            fname = outdir+'/boxplot_time_boundary_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'residuals_labels':
            plt.ylabel('Residuals '+r' $(W\,m^{-2})$',fontsize=14)
            # Saving figure
            fname = outdir+'/boxplot_time_residuals_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'budget_diff_labels':
            plt.ylabel('Enery budgets (estimated using finite diffs. '+r' $(W\,m^{-2})$',fontsize=14)
            # Saving figure
            fname = outdir+'/boxplot_time_budget_diff_terms.png'
            plt.savefig(fname)
            print(fname+' created')
        elif termlist == 'gendiss_labels':
            plt.ylabel('Generation/dissipation '+r' $(W\,m^{-2})$',fontsize=14)
            # Saving figure
            fname = outdir+'/boxplot_timegendiss_terms.png'
            plt.savefig(fname)
            print(fname+' created')

def main():
    # Open data from energy and conversion terms as input from user from command line
    data = sys.argv[1]
    print('Reading data from: '+data)
    df = pd.read_csv(data)
    print(df)
    print(' ')
    # Diectory for saving figures
    ResultsSubDirectory = '/'.join(data.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/'
    # Make plot for timeseries
    plot_timeseries(df,FigsDir)
    plot_boxplot(df,FigsDir)
    print('All done!')

if __name__ == "__main__":
    main()

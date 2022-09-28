#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 15:33:53 2022

@author: danilocoutodsouza
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def bar_plot(data1,data2,terms):
    plt.close('all')
    fig, ax = plt.subplots(figsize=(8,5))
    width = 0.2
    offset = 0.1
    idx = np.asarray([i for i in range(len(data_eu[terms]))])       
    ax.bar(idx-offset, [val for key,val in sorted(data1.items())],
           width=width, color='#3B95BF')
    ax.bar(idx+offset, [val for key,val in sorted(data2.items())], 
           width=width, color='#BF3D3B')    
    ax.set_xticks(idx)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_xticklabels(terms)
    ax.legend(['Eulerian', 'Lagrangian'],fontsize=12)
    plt.title(itime)
    if terms == terms1:
        ax.set_ylabel(r' $J\,m^{-2})$',fontsize=12)
        plt.savefig(outdir+'bars_energy-'+itime)
    else:
        ax.set_ylabel(r' $W\,m^{-2})$',fontsize=12)
        plt.savefig(outdir+'bars_'+itime)
    
def time_series(data_eu,data_lag,term):    
    times = df1['Date_Hour']
    dates = []
    for t in times:
        dates.append(t[5:]+'Z')
    plt.close('all')
    fig, ax = plt.subplots(figsize=(12,8))
    ax.plot(dates,data_eu,color='#3B95BF',label='Eulerian',linewidth=2)
    ax.plot(dates,data_lag,'--', color='#BF3D3B',label='Lagrangian',linewidth=2)
    ax.axvline(dates[i],zorder=0,c='#383838',alpha=0.4,linewidth=1,
               label='Point which data should be equal')
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=-45, ha="right" )
    ax.tick_params(axis='x',rotation=45)
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.title(term)
    ax.legend(fontsize=12)
    if terms == terms1:
        ax.set_ylabel(r' $J\,m^{-2})$',fontsize=12)
    else:
        ax.set_ylabel(r' $W\,m^{-2})$',fontsize=12)    
    plt.savefig(outdir+'timeseries_'+term)    

outdir =  '../energetica_python_etc/Compare_Lagrangian-Eulerian/'
indir = '/Users/danilocoutodsouza/Documents/USP/LEC_Results/'

file1 = indir+'/Reg1-Representative_NCEP-R2_EulereianRef/Reg1-Representative_NCEP-R2_52W37W37S22S.csv'
file2 = indir+'/Reg1-Representative_NCEP-R2_lagranigan/Reg1-Representative_NCEP-R2_lagranigan.csv'
df1 = pd.read_csv(file1, parse_dates=[['Date', 'Hour']])
df2 = pd.read_csv(file2, parse_dates=[['Date', 'Hour']])

# =============================================================================
# # Timestep where data should be equal
i = round(len(df1)/2)
# =============================================================================

terms1 = ['Az','Ae','Kz','Ke']
terms2 = ['Cz','Ca','Ck','Ce','BAz','BAe','BKz','BKe','Gz','Ge']
# Get 1 timestep before and 1 after
for t in range(-1,2):  
    data_eu = df1.iloc[i+t]
    data_lag = df2.iloc[i+t]
    itime = df1.iloc[i+t]['Date_Hour']+'Z'        
    for terms in [terms1,terms2]:
        data1, data2 = data_eu[terms], data_lag[terms]
        # bar_plot(data1,data2,terms)
        
for term in [*terms1,*terms2]:
    data_eu = df1[term]
    data_lag = df2[term]
    time_series(data_eu,data_lag,term)
        


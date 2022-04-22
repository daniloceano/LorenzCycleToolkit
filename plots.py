#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 17:56:57 2022

@author: danilocoutodsouza
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import numpy as np
from datetime import datetime
import matplotlib.dates

data = 'LEC_Reg1.nc_60W30W42S17S.csv'
df = pd.read_csv(data)

def timeseries(df):
    
    plt.close('all')
    
    linecolors = ['#A53860','#C9B857','#384A0F','#473BF0']
    markerfacecolors = ['#A53860','w','#384A0F','w']
    conversion_labels = ['Cz','Ca','Ck','Ce']
    markers = ['s','s','o','o']         
    linestyles = ['-','-','-','-']
    energy_labels = ['Az','Ae','Kz','Ke']
    linewidth = 4
    
    # file = sys.argv[1].split('/')[-1]
    date = df['Date']
    times = pd.date_range(date[0],date.iloc[-1],periods=len(date))

    
    for labels in [energy_labels,conversion_labels]:
        
        print('Printing '+str(labels)+'...')
        
        maxval = np.amax(np.amax(df[labels]))
        minval = np.amin(np.amin(df[labels]))
        print('Data range: '+str(minval)+' to '+str(maxval))
        
        plt.figure(figsize=(8,8))
        for term,i in zip(labels,range(len(labels))):
            
                                                   
            plt.grid(b=True,c='gray',linewidth=0.25,linestyle='dashdot')
            plt.plot(times,df[term],c=linecolors[i], marker= markers[i],
                    markerfacecolor=markerfacecolors[i], label=term,
                    linewidth=linewidth,markersize=6, linestyle=linestyles[i])
            plt.tick_params(axis='x', labelrotation=20)
            plt.legend()
            plt.xlim(times[0],times[-1])
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            
        ax = plt.gca()
        # ax.set_xticks(dates)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            
        if term in energy_labels:
            ax.set_title('Energy Terms',fontsize=12)
            plt.ylabel('Energy (J m^-2)')
        elif term in conversion_labels:
            plt.ylabel('Conversion (W m^-2)')
            ax.set_title('Conversion Terms',fontsize=12)

    # plt.savefig(file.split('.')[0]+'_2.png')
    # plt.savefig('figures_zonal average tests/testv3_'+file.split('.')[0]+'.png')
    # plt.savefig(file.split('/')[-1][9:].split('.')[0]+'.png')  

timeseries(df)
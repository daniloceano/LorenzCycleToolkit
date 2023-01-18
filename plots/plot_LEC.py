#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 09:29:19 2022

This script reads an CSV file with all terms from the Lorenz Energy Cycle 
(as input from user) and make the deafult figures for the Lorenz energy cycle.
The transparecy in each box is set to be proportional to the energy tendency,
as well as the arrows are set to be proportional to the conversion rates.

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
import matplotlib.gridspec as gridspec
import argparse
from plot_timeseries import check_create_folder

def Cz(ax,value,i,width,head_width):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(-0.637, 0.5, 0.78, 0, head_width = head_width, width = width,
                  length_includes_head=True,
                  fc=cols_conversion[i],ec=cols_conversion[i],
                  clip_on=False,transform=ax.transAxes)
        
    else:
        ax.arrow(0.139, 0.5, -0.775, 0, head_width = head_width, width = width,
          length_includes_head=True,
          fc=cols_conversion[i],ec=cols_conversion[i],
          clip_on=False,transform=ax.transAxes)
    ax.text(-0.25,0.57,value,fontdict={'fontsize':fs},transform=ax.transAxes,
            verticalalignment='center',horizontalalignment='center')

def Ck(ax,value,i,width,head_width):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, 0.952, 0, 0.52, head_width = head_width,width = width,
                 fc=cols_conversion[i],ec=cols_conversion[i],
                 clip_on=False,transform = ax.transAxes)
    else:
        ax.arrow(0.5, 1.55, 0, -0.6, head_width = head_width,width = width,
                  length_includes_head=True,
              fc=cols_conversion[i],ec=cols_conversion[i],
              clip_on=False,transform = ax.transAxes)
    ax.text(0.55,1.25,value,fontdict={'fontsize':fs},transform=ax.transAxes,
            verticalalignment='center',horizontalalignment='left')

def Ca(ax,value,i,width,head_width):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, 0.04, 0, -0.59, head_width = head_width,width = width,
                      fc=cols_conversion[i],ec=cols_conversion[i],
                                length_includes_head=True,
                      clip_on=False,transform = ax.transAxes)
    else:
        
        ax.arrow(0.5, -0.55, 0, 0.59, head_width = head_width,width = width,
                 length_includes_head=True,
                 fc=cols_conversion[i],ec=cols_conversion[i],
                 clip_on=False,transform = ax.transAxes)
    ax.text(0.45,-0.26,value,fontdict={'fontsize':fs},transform = ax.transAxes,
            verticalalignment='center',horizontalalignment='right')

def Ce(ax,value,i,width,head_width):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.86, 0.5, 0.78, 0, head_width = head_width, width = width,
                  length_includes_head=True,
                  fc=cols_conversion[i],ec=cols_conversion[i],
                  clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(1.64, 0.5, -0.78, 0, head_width = head_width, width = width,
              fc=cols_conversion[i],ec=cols_conversion[i],
                       length_includes_head=True,
              clip_on=False,transform=ax.transAxes)
    ax.text(1.25,0.42,value,fontdict={'fontsize':fs},transform=ax.transAxes,
            verticalalignment='center',horizontalalignment='center')

def RGz_RKz(ax,value,i,width,head_width):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, 1.23, 0, -0.275, head_width = head_width, width = width,
                 length_includes_head=True,
          fc=cols_residual[i],ec=cols_residual[i],
          clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.5, 0.955, 0, 0.275, head_width = head_width, width = width,
          fc=cols_residual[i],ec=cols_residual[i],
                    length_includes_head=True,
          clip_on=False,transform=ax.transAxes)
    ax.text(0.5,1.27,value,fontdict={'fontsize':fs},transform=ax.transAxes,
            verticalalignment='center',horizontalalignment='center')


def RGe_RKe(ax,value,i,width,head_width):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, -0.23, 0, 0.275, head_width = head_width, width = width,
                 length_includes_head=True,
      fc=cols_residual[i],ec=cols_residual[i],
      clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.5, 0.05, 0, -0.275, head_width = head_width, width = width,
                 length_includes_head=True,
      fc=cols_residual[i],ec=cols_residual[i],
      clip_on=False,transform=ax.transAxes)
    ax.text(0.5,-.28,value,fontdict={'fontsize':fs},transform=ax.transAxes,
            verticalalignment='center',horizontalalignment='center')
    

def BAz_BAe(ax,value,i,width,head_width):
     if isinstance(value, str) or (isinstance(value, float) and value > 0):
         ax.arrow(-0.135,0.5, 0.275, 0, head_width = head_width, width = width,
                  length_includes_head=True,
                  fc=cols_boundary[i],ec=cols_boundary[i],
                  clip_on=False,transform=ax.transAxes)
     else:
         ax.arrow(0.14,0.5, -0.275, 0, head_width = head_width, width = width,
                  fc=cols_boundary[i],ec=cols_boundary[i],
           length_includes_head=True,
           clip_on=False,transform=ax.transAxes)
     ax.text(-0.18,0.5,value,fontdict={'fontsize':fs},transform=ax.transAxes,
            verticalalignment='center',horizontalalignment='right')

def BKz_BKe(ax,value,i,width,head_width):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(1.14,0.5, -0.275, 0, head_width = head_width, width = width,
                 length_includes_head=True,
             fc=cols_boundary[i],ec=cols_boundary[i],
             clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.863,0.5, 0.275, 0, head_width = head_width, width = width,
             fc=cols_boundary[i],ec=cols_boundary[i],
                      length_includes_head=True,
             clip_on=False,transform=ax.transAxes)
    ax.text(1.18,0.5,value,fontdict={'fontsize':fs},transform=ax.transAxes,
            verticalalignment='center',horizontalalignment='left')

def plot_LEC(idata,flag):
    
    # Adjust box size proportionally to the energy budget
    x = idata[energys]
    energy_scaled=(x-x.min())/(x.max()-x.min())
        
    # Adjust arrow size proportionally to the conversion rate, residual and 
    # boundary
    x = np.abs(idata[conversions+residuals+boundaries])
    conversion_scaled=(x-x.min())/(x.max()-x.min())
    
    plt.close('all')
    fig = plt.figure(figsize=(14, 11))
    plt.axis('off')
    gs = gridspec.GridSpec(nrows=2, ncols=2,hspace=0.5,wspace=0.5,
                           right=(0.89))
    
    if flag == 'example':
        plt.title('(Daily mean/period)',fontsize=fs, loc='center',y=0.5,
                  fontdict={'fontweight':'bold'})
    elif flag == 'daily_mean':
        plt.title(str(idata.name.date()),fontsize=fs, loc='center',y=0.5,
              fontdict={'fontweight':'bold'})
    elif flag == 'periods':
        plt.title((idata['Period']),fontsize=fs, loc='center',y=0.5,
               fontdict={'fontweight':'bold'})
    
    i = 0
    for row in range(2):
        for col in range(2):
            ax = fig.add_subplot(gs[row,col])
            
            if flag == 'example':
                energy = energys[i][:7]
                plt.text(0.5, 0.5,energy,fontdict={'fontsize':fs},
                          transform = ax.transAxes,
                          verticalalignment='center',horizontalalignment='center')
                conversion = conversions[i]
                residual = residuals[i]
                boundary = boundaries[i]
                alpha = 1
                width_conversion,width_boundary,width_residual = 0.01, 0.01,0.01
                head_conversion,head_boundary,head_residual = 0.05,0.05,0.05
                
            else:
                energy = round(idata[energys[i]],1)
                plt.text(0.5, 0.5,energy,fontdict={'fontsize':fs},
                          transform = ax.transAxes,
                          verticalalignment='center',horizontalalignment='center')
                conversion = round(idata[conversions[i]],1)
                residual = round(idata[residuals[i]],1)
                boundary = round(idata[boundaries[i]],1)
                # alpha = 0.2+(0.8*energy_scaled[energys[i]].iloc[time])
                alpha = .9
                width_conversion = 0.01+(0.05*conversion_scaled[conversions[i]])
                width_boundary = 0.01+(0.05*conversion_scaled[boundaries[i]])
                width_residual = 0.01+(0.05*conversion_scaled[residuals[i]])
                head_conversion = 0.05+(0.07*conversion_scaled[conversions[i]])
                head_boundary = 0.05+(0.07*conversion_scaled[boundaries[i]])
                head_residual = 0.05+(0.07*conversion_scaled[residuals[i]])

            # Draw energy budgets as boxes
            square = plt.Rectangle((0, 0),5,5, fc=cols_energy[i],ec='k',
                                   alpha=alpha)
            ax.add_patch(square)
            plt.axis("equal")
            plt.axis('off')
                
            if row == 0 and col == 0: 
                Ca(ax,conversion,i,width_conversion,head_conversion)
                RGz_RKz(ax,residual,i,width_residual,head_residual)
                BAz_BAe(ax,boundary,i,width_boundary,head_boundary)
            if row == 0 and col == 1: 
                Cz(ax,conversion,i,width_conversion,head_conversion)
                RGz_RKz(ax,residual,i,width_residual,head_residual)
                BKz_BKe(ax,boundary,i,width_boundary,head_boundary)
            if row == 1 and col == 0:
                Ce(ax,conversion,i,width_conversion,head_conversion)
                RGe_RKe(ax,residual,i,width_residual,head_residual)
                BAz_BAe(ax,boundary,i,width_boundary,head_boundary)
            if row == 1 and col == 1:
                Ck(ax,conversion,i,width_conversion,head_conversion)
                RGe_RKe(ax,residual,i,width_residual,head_residual)
                BKz_BKe(ax,boundary,i,width_boundary,head_boundary)
            i+=1
            
    if flag == 'example':
        plt.savefig(FigsDir+'LEC_example.png')
        print('Created LEC example figure')
    elif flag == 'daily_mean':
        datestr = idata.name.strftime("%Y-%m-%d")
        print('Created LEC (daily mean) for: '+datestr)
        plt.savefig(FigsDir+'LEC_'+datestr+'.png')
    elif flag == 'periods':
        print('Created LEC for period: '+idata['Period'])
        plt.savefig(FigsDir+'LEC_'+idata['Period']+'.png')
    

def main():
    
    df = pd.read_csv(outfile)
    df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')
    # Get mean daily values
    data = df.groupby(pd.Grouper(key="Datetime", freq="1D")).mean()
    # plot example figure
    plot_LEC(data,'example')
    # plot each deaily mean
    for t in range(len(data)):
        idata = data.iloc[t]
        plot_LEC(idata,'daily_mean')
    # plot means for each periods of the system
    periods = pd.read_csv('../inputs/periods',sep= ';',header=0)
    for i in range(len(periods)):
        start,end = periods.iloc[i]['start'],periods.iloc[i]['end']
        selected_dates = df[(df['Datetime'] >= start) & (df['Datetime'] <= end)]
        period = selected_dates.drop(['Datetime','Date','Hour'],axis=1).mean()
        period['Period'] = periods.iloc[i]['Period']
        plot_LEC(period,'periods')


if __name__ == "__main__":
 
    parser = argparse.ArgumentParser(description = "\
reads an CSV file with all terms from the Lorenz Energy Cycle \
 (as input from user) and make the deafult figures for the Lorenz energy cycle\
 The arrows are set to be proportional to the conversion rates.")
    parser.add_argument("outfile", help = "The .csv file containing the \
results from the main.py program.")

    args = parser.parse_args()
    outfile = args.outfile
    ResultsSubDirectory = '/'.join(outfile.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/LEC/'
    check_create_folder(FigsDir)

    # Specs for plotting
    fs = 30
    energys = ['∂Az/∂t (finite diff.)','∂Kz/∂t (finite diff.)',
               '∂Ae/∂t (finite diff.)', '∂Ke/∂t (finite diff.)']
    conversions = ['Ca', 'Cz', 'Ce', 'Ck']
    residuals = ['RGz', 'RKz', 'RGe', 'RKe']
    boundaries = ['BAz','BKz','BAe','BKe'] 
    cols_energy = ['#5EA4BD','#F7EF7C','#96DB6E','#F77B59'] 
    cols_conversion = ['#5C5850','#5C5850','#5C5850','#5C5850']
    cols_residual = ['#5C5850','#5C5850','#5C5850','#5C5850']
    cols_boundary =  ['#5C5850','#5C5850','#5C5850','#5C5850']

    
    main()
    

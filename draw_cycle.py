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
    ax.text(-0.25,0.57,value,fontdict={'fontsize':18},transform=ax.transAxes,
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
    ax.text(0.55,1.25,value,fontdict={'fontsize':18},transform=ax.transAxes,
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
    ax.text(0.45,-0.26,value,fontdict={'fontsize':18},transform = ax.transAxes,
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
    ax.text(1.25,0.42,value,fontdict={'fontsize':18},transform=ax.transAxes,
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
    ax.text(0.5,1.27,value,fontdict={'fontsize':18},transform=ax.transAxes,
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
    ax.text(0.5,-.28,value,fontdict={'fontsize':18},transform=ax.transAxes,
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
     ax.text(-0.18,0.5,value,fontdict={'fontsize':18},transform=ax.transAxes,
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
    ax.text(1.18,0.5,value,fontdict={'fontsize':18},transform=ax.transAxes,
            verticalalignment='center',horizontalalignment='left')


def main(time, example=False):
    data = daily_means.iloc[time]
    plt.close('all')
    fig = plt.figure(figsize=(14, 11))
    gs = gridspec.GridSpec(nrows=2, ncols=2,hspace=0.5,wspace=0.5,
                           right=(0.89))
    i = 0
    if example == True:
        plt.title('(Daily mean)',fontsize=20, loc='center',y=0.5,
                  fontdict={'fontweight':'bold'})
    else:
        plt.title(str(data.name.date()),fontsize=20, loc='center',y=0.5,
              fontdict={'fontweight':'bold'})
    
    plt.axis('off')

    for row in range(2):
        for col in range(2):
            ax = fig.add_subplot(gs[row,col])
            
            if example == True:
                energy = energys[i][:7]
                plt.text(0.5, 0.5,energy,fontdict={'fontsize':25},
                          transform = ax.transAxes,
                          verticalalignment='center',horizontalalignment='center')
                conversion = conversions[i]
                residual = residuals[i]
                boundary = boundaries[i]
                alpha = 1
                width_conversion,width_boundary,width_residual = 0.01, 0.01,0.01
                head_conversion,head_boundary,head_residual = 0.05,0.05,0.05
                
            else:
                energy = round(data[energys[i]],1)
                plt.text(0.5, 0.5,energy,fontdict={'fontsize':20},
                          transform = ax.transAxes,
                          verticalalignment='center',horizontalalignment='center')
                conversion = round(data[conversions[i]],1)
                residual = round(data[residuals[i]],1)
                boundary = round(data[boundaries[i]],1)
                # alpha = 0.2+(0.8*energy_scaled[energys[i]].iloc[time])
                alpha = .9
                width_conversion = 0.01+(0.05*conversion_scaled[conversions[i]].iloc[time])
                width_boundary = 0.01+(0.05*conversion_scaled[boundaries[i]].iloc[time])
                width_residual = 0.01+(0.05*conversion_scaled[residuals[i]].iloc[time])
                head_conversion = 0.05+(0.07*conversion_scaled[conversions[i]].iloc[time])
                head_boundary = 0.05+(0.07*conversion_scaled[boundaries[i]].iloc[time])
                head_residual = 0.05+(0.07*conversion_scaled[residuals[i]].iloc[time])

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
            
    if example == True:
        plt.savefig(FigsDir+'LEC_example.png')
    else:
        print(df.iloc[time]['Date'])
        plt.savefig(FigsDir+'LEC_'+str(daily_means.index[time])+'.png')


if __name__ == "__main__":
 
    parser = argparse.ArgumentParser(description = "\
reads an CSV file with all terms from the Lorenz Energy Cycle \
 (as input from user) and make the deafult figures for the Lorenz energy cycle\
The transparecy in each box is set to be proportional to the energy tendency,\
 as well as the arrows are set to be proportional to the conversion rates.")
    parser.add_argument("outfile", help = "The .csv file containing the \
results from the main.py program.")

    args = parser.parse_args()
    outfile = args.outfile
    ResultsSubDirectory = '/'.join(outfile.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/'

    df = pd.read_csv(outfile)
    df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')
    # Get mean daily values
    daily_means = df.groupby(pd.Grouper(key="Datetime", freq="1D")).mean()
    
    
    # Specs for plotting
    energys = ['∂Az/∂t (finite diff.)','∂Kz/∂t (finite diff.)',
               '∂Ae/∂t (finite diff.)', '∂Ke/∂t (finite diff.)']
    conversions = ['Ca', 'Cz', 'Ce', 'Ck']
    residuals = ['RGz', 'RKz', 'RGe', 'RKe']
    boundaries = ['BAz','BKz','BAe','BKe'] 
    cols_energy = ['#5EA4BD','#F7EF7C','#96DB6E','#F77B59'] 
    cols_conversion = ['#5C5850','#5C5850','#5C5850','#5C5850']
    cols_residual = ['#5C5850','#5C5850','#5C5850','#5C5850']
    cols_boundary =  ['#5C5850','#5C5850','#5C5850','#5C5850']
    # cols_conversion = ['#6785A6','#254A73','#A8893C','#DE8268']
    # cols_residual = ['#406F80','#807B40','#578040','#80402E']
    # cols_boundary = ['#1C4080','#8A8A00','#169A00','#AE241C']

    
    # Adjust circle size proportionally to the energy budget
    x = daily_means[energys]
    energy_scaled=(x-x.min().min())/(x.max().max()-x.min().min())
        
    # Adjust arrow size proportionally to the conversion rate, residual and 
    # boundary
    x = np.abs(daily_means[conversions+residuals+boundaries])
    conversion_scaled=(x-x.min().min())/(x.max().max()-x.min().min())
    
    main(0,example=True)
    # plot each deaily mean
    for t in range(len(daily_means)):
        main(t)
    

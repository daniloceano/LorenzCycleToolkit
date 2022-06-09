#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 09:29:19 2022

@author: danilocoutodsouza
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



def Cz(ax,value):
    # Az -> Ae
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.95, 0.5, 0.8, 0, head_width=0.05, width=0.002,
                  fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    # Ae -> Az
    else:
        ax.arrow(1.95, 0.5, -0.92, 0, head_width=0.05, width=0.002,
                  fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(1.38,0.55,value,fontdict={'fontsize':18},transform=ax.transAxes)

def Ca(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, 1.55, 0, -0.525, head_width=0.05,
             fc='k',ec='k',clip_on=False,transform = ax.transAxes)
    else:
        ax.arrow(0.5, 0.95, 0, 0.525, head_width=0.05,
                 fc='k',ec='k',clip_on=False,transform = ax.transAxes)
    ax.text(0.35,1.15,value,fontdict={'fontsize':18},transform=ax.transAxes)

def Ck(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, -0.95, 0, 0.92, head_width=0.05,
                      fc='k',ec='k',clip_on=False,transform = ax.transAxes)
    else:
        ax.arrow(0.5, 0.05, 0, -0.52, head_width=0.05,
                      fc='k',ec='k',clip_on=False,transform = ax.transAxes)
    ax.text(0.55,-0.35,value,fontdict={'fontsize':18},transform = ax.transAxes)

def Ce(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(-0.83, 0.5, 0.8, 0, head_width=0.05, width=0.002,
                  fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.05, 0.5, -0.8, 0, head_width=0.05, width=0.002,
              fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(-0.4,0.4,value,fontdict={'fontsize':18},transform=ax.transAxes)

def RGz(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, 1.23, 0, -0.2, head_width=0.05, width=0.002,
          fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.5, 0.95, 0, 0.2, head_width=0.05, width=0.002,
          fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(0.42,1.25,value,fontdict={'fontsize':18},transform=ax.transAxes)
    
def RKz(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, 1.23, 0, -0.2, head_width=0.05, width=0.002,
          fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.5, 0.95, 0, 0.2, head_width=0.05, width=0.002,
          fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(0.42,1.25,value,fontdict={'fontsize':18},transform=ax.transAxes)

def RGe(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, -0.23, 0, 0.2, head_width=0.05, width=0.002,
      fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.5, 0.05, 0, -0.2, head_width=0.05, width=0.002,
      fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(0.42,-.31,value,fontdict={'fontsize':18},transform=ax.transAxes)
    
def RKe(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(0.5, -0.23, 0, 0.2, head_width=0.05, width=0.002,
      fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.5, 0.05, 0, -0.2, head_width=0.05, width=0.002,
      fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(0.42,-.31,value,fontdict={'fontsize':18},transform=ax.transAxes)

def BAz(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(-0.23,0.5, 0.2, 0, head_width=0.05, width=0.002,
   fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.045,0.5, -0.2, 0, head_width=0.05, width=0.002,
  fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(-0.4,0.48,value,fontdict={'fontsize':18},transform=ax.transAxes)

def BKz(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(1.23,0.5, -0.2, 0, head_width=0.05, width=0.002,
             fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.955,0.5, 0.2, 0, head_width=0.05, width=0.002,
             fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(1.25,0.48,value,fontdict={'fontsize':18},transform=ax.transAxes)
    
def BAe(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
      ax.arrow(-0.23,0.5, 0.2, 0, head_width=0.05, width=0.002,
               fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    else:
      ax.arrow(0.045,0.5, -0.2, 0, head_width=0.05, width=0.002,
               fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(-0.4,0.48,value,fontdict={'fontsize':18},transform=ax.transAxes)
    
def BKe(ax,value):
    if isinstance(value, str) or (isinstance(value, float) and value > 0):
        ax.arrow(1.23,0.5, -0.2, 0, head_width=0.05, width=0.002,
             fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    else:
        ax.arrow(0.955,0.5, 0.2, 0, head_width=0.05, width=0.002,
             fc='k',ec='k',clip_on=False,transform=ax.transAxes)
    ax.text(1.25,0.48,value,fontdict={'fontsize':18},transform=ax.transAxes)
    


def example():
    plt.close('all')
    fig = plt.figure(figsize=(13, 11))
    gs = gridspec.GridSpec(nrows=2, ncols=2,hspace=0.5,wspace=0.5)
    i = 0
    plt.title('(Daily mean)',fontsize=20, loc='center',y=0.5)
    plt.axis('off')
    for col in range(2):
        for row in range(2):
            ax = fig.add_subplot(gs[row,col])
            # Circles for energy terms
            circle = plt.Circle((0, 0), radius=0.5, fc='gray')
            plt.gca().add_patch(circle)
            plt.axis('scaled')
            plt.axis('off')
            # Text in the circles
            energy = energys[i]
            plt.text(0.45, 0.5,energy,fontdict={'fontsize':18},
                     transform = ax.transAxes)
            conversion = conversions[i]
            residual = residuals[i]
            boundary = boundaries[i]
            if row == 0 and col == 0: 
                Cz(ax,conversion)
                RGz(ax,residual)
                BAz(ax,boundary)
            if row == 0 and col == 1: 
                Ck(ax,conversion)
                RKz(ax,residual)
                BKz(ax,boundary)
            if row == 1 and col == 0:
                Ca(ax,conversion)
                RGe(ax,residual)
                BAe(ax,boundary)
            if row == 1 and col == 1:
                Ce(ax,conversion)
                RKe(ax,residual)
                BKe(ax,boundary)
            i+=1
    plt.savefig('../LEC_Results/test.png')
        
def plot_LEC(time):
    data = daily_means.iloc[time]
    plt.close('all')
    fig = plt.figure(figsize=(13, 11))
    gs = gridspec.GridSpec(nrows=2, ncols=2,hspace=0.5,wspace=0.5)
    i = 0
    plt.title(str(data.name.date()),fontsize=20, loc='center',y=0.5)
    plt.axis('off')
    for col in range(2):
        for row in range(2):
            ax = fig.add_subplot(gs[row,col])
            # Circles for energy terms
            circle = plt.Circle((0, 0), radius=0.5, fc='gray')
            plt.gca().add_patch(circle)
            plt.axis('scaled')
            plt.axis('off')
            # Text in the circles
            energy = round(data[energys[i]]/1e5,2)
            plt.text(0.45, 0.5,energy,fontdict={'fontsize':18},
                     transform = ax.transAxes)
            conversion = round(data[conversions[i]],2)
            residual = round(data[residuals[i]],2)
            boundary = round(data[boundaries[i]],2)
            if row == 0 and col == 0: 
                Cz(ax,conversion)
                RGz(ax,residual)
                BAz(ax,boundary)
            if row == 0 and col == 1: 
                Ck(ax,conversion)
                RKz(ax,residual)
                BKz(ax,boundary)
            if row == 1 and col == 0:
                Ca(ax,conversion)
                RGe(ax,residual)
                BAe(ax,boundary)
            if row == 1 and col == 1:
                Ce(ax,conversion)
                RKe(ax,residual)
                BKe(ax,boundary)
            i+=1
    plt.savefig('../LEC_Results/test.png')
            
file = '../LEC_Results/Reg1_wo_stress_60W30W42S17S/Reg1_wo_stress_60W30W42S17S.csv'
df = pd.read_csv(file)

df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')

daily_means = df.groupby(pd.Grouper(key="Datetime", freq="1D")).mean()

energys = ['Az','Ae','Kz','Ke']
conversions = ['Cz', 'Ca', 'Ck', 'Ce']
residuals = ['RGz', 'RGe', 'RKz', 'RKe']
boundaries = ['BAz','BAe','BKz','BKe'] 

plot_LEC(0)
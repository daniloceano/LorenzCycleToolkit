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


file = '/Users/danilocoutodsouza/Documents/USP/LEC_Results/Reg1_ERA5_60W30W42S17S/Reg1_ERA5_60W30W42S17S.csv'
df = pd.read_csv(file)

df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')

daily_means = df.groupby(pd.Grouper(key="Datetime", freq="1D")).mean()

energys = ['Az','Ke','Ae','Kz']

plt.close('all')
fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(nrows=2, ncols=2)
i = 0
plt.title(daily_means.index[i].date())
plt.axis('off')
for row in range(2):
    for col in range(2):
        energy = energys[i]
        energy_value = str(round(daily_means[energy][0]/1e5,2))
        ax = fig.add_subplot(gs[row,col])
        circle = plt.Circle((0, 0), radius=0.5, fc='gray')
        plt.gca().add_patch(circle)
        plt.axis('scaled')
        plt.axis('off')
        plt.text(0.37, 0.5,energy+' = '+energy_value,
                 transform = ax.transAxes)
        i+=1
        
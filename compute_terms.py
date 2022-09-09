#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:56:32 2022

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
    
"""

from metpy.units import units
import numpy as np


# Compute the budget equation for the energy terms (Az, Ae, Kz and Ke) using
# finite differences method (used for estimating generation, disspation and
# boundary work terms as residuals)
def calc_budget_diff(df,time):
    # get time delta in seconds
    dt = float((time[1]-time[0]) / np.timedelta64(1, 's'))
    # Estimate budget values for all energy terms
    for term in ['Az','Ae','Kz','Ke']:
        name = '∂'+term+'/∂t (finite diff.)'
        print('\nEstimating '+name)
        # forward finite difference for the first value
        forward = (df[term].iloc[1]-df[term].iloc[0])/dt
        # central finite differentes for the second and the penultimate values
        central_second = (df[term].iloc[2]-df[term].iloc[0])/dt
        central_penultimate = (df[term].iloc[-1]-df[term].iloc[-3])/dt
        # fourth order for the third to antepenultimate value
        fourthorder1 = (4/3)*(
            df[term].iloc[3:-1].values-df[term].iloc[1:-3].values)/(2*dt)
        fourthorder2 = (1/3)*(
            df[term].iloc[4:].values-df[term].iloc[1:-3].values)/(4*dt)
        fourthorder = fourthorder1-fourthorder2
        # backward finite difference for the last value
        backward = (df[term].iloc[-1]-df[term].iloc[-2])/dt
        # put all values together
        df[name] = [forward,central_second,*fourthorder,
                    central_penultimate,backward]
        print(df[name].values*units('W/ m **2'))
    return df

# Compute the residuals RGz, RKz, RGe and RKe using the budget terms estimated
# via finite differences
def calc_residuals(df):
    print('\nResiduals ('+str((1*units('W/ m **2')).units)+'):')
    df['RGz'] = df['∂Az/∂t (finite diff.)'] + df['Cz'] + df['Ca'] - df['BAz']
    df['RGe'] = df['∂Ae/∂t (finite diff.)'] - df['Ca'] + df['Ce'] - df['BAe']
    df['RKz'] = df['∂Kz/∂t (finite diff.)'] - df['Cz'] - df['Ck'] - df['BKz']
    df['RKe'] = df['∂Ke/∂t (finite diff.)'] - df['Ce'] + df['Ck'] - df['BKe']
    print(df[['RGz', 'RKz', 'RGe', 'RKe']])
    return df
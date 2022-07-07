#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 18:36:50 2022

@author: daniloceano
"""

import pandas as pd
import os
import glob

## GET STATS FOR ALL KE VALUES
files = glob.glob('../LEC_Results/*NCEP-R2*/*NCEP-R2*')
kes, cas, cks, ges = [],[],[],[]
for outfile in files:
    df = pd.read_csv(outfile)
    for ke,ca,ck,ge in zip(df['Ke'].values, df['Ca'].values,
                           df['Ck'].values,df['Ge'].values):
        kes.append(ke)
        cas.append(ca)
        cks.append(ck)
        ges.append(ge)

for data,term in zip([kes,cas,cks,ges],['ke','ca','ck','ge']):
    print('\n--------')
    print(term)
    print(pd.DataFrame(data).describe())

for outfile in files:
    os.system("python LorenzPhaseSpace.py "+outfile)
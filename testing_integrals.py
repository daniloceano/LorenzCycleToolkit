#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 31 13:00:20 2021

@author: danilocoutodsouza
"""

import aux
# from conversion_terms import (calc_cz, calc_ca, calc_ce, calc_ck)
# from energy_contents import (calc_AZ, calc_AE, calc_KZ, calc_KE)
# from boundary_terms import (calc_BAz, calc_BAe, calc_BKz, calc_BKe,
#                             calc_BOz, calc_BOe)
from calc import Differentiate

from calc import (HorizontalTrazpezoidalIntegration,
                       VerticalTrazpezoidalIntegration,
                       CalcZonalAverage, CalcAreaAverage,
                       StaticStability)
from EnergyContents import EnergyContents
from ConversionTerms import ConversionTerms

from metpy.units import units

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import numpy as np

from BoxData_v3 import BoxData

import glob
    
def Plotter(EnergyTerms,ConversionTerms,TimeName,file):
    
    time = EnergyTerms[0][TimeName]
    
    cs = ['r','y','g','k']
    markers = ['s','s','o','o']
    markerfacecolors = ['r','w','g','w']
    elabels = ['Az','Ae','Kz','Ke']
    clabels = ['Cz','Ca','Ck','Ce']

    plt.close('all')
    fig = plt.figure(figsize=(8,4))
    gs = gridspec.GridSpec(1, 2)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    axs = [ax0,ax1]
    
    for j in range(2):
        for i in range(4):
            if j==0:
                label=elabels
                axs[j].set_ylabel('Energy (10^5 J m^-2)')
                Terms = EnergyTerms
                label = elabels
                
            elif j==1:
                label=clabels
                axs[j].set_ylabel('Conversion (W m^-2)')
                Terms = ConversionTerms
                label = clabels
                   
            axs[j].grid(b=True,c='gray',linewidth=0.25,linestyle='dashdot')
            axs[j].plot(time,Terms[i],c= cs[i], marker= markers[i],
                    markerfacecolor=markerfacecolors[i], label = label[i])
            axs[j].tick_params(axis='x', labelrotation=20)
            axs[j].legend()
            axs[j].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            axs[j].set_xlim(time[0].values,time[-1].values)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
     
    
    axs[0].set_yticks(range(0,round(np.amax(EnergyTerms))+2,2))
    axs[1].set_yticks(range(round(np.amin(ConversionTerms)-2),
                            round(np.amax(ConversionTerms))+2,2))        
    axs[0].set_title('Energy Terms',fontsize=12)
    axs[1].set_title('Conversion Terms',fontsize=12)
    plt.suptitle(file.split('.')[0],fontsize=14)

    # plt.savefig(file.split('.')[0]+'_2.png')
    plt.savefig('figures_zonal average tests/testv3_'+file.split('.')[0]+'.png')
    # plt.savefig(file.split('/')[-1][9:].split('.')[0]+'.png')  

    
# =============================================================================
# # TESTING ENERGY CONTENTS
# =============================================================================

plt.close('all')

files = ['Reg1.nc','Reg2.nc','Reg3.nc','Catarina.nc']
minlons = [-60,-65,-82.5,-60]
maxlons = [-30,-20,-15,-30]
minlats = [-42,-65,-60,-35]
maxlats = [-17.5,-20,-30,-20]
fvar = 'FVars/NCEP_R2_RDA.csv'

# files = glob.glob('/Users/danilocoutodsouza/Documents/USP/Master/MPAS/Catarina/*/isobaric_*')
# minlons = np.ones(7)*[-60]
# maxlons = np.ones(7)*[-30]
# minlats = np.ones(7)*[-35]
# maxlats = np.ones(7)*[-20]
# fvar = 'FVars/MPAS.csv'

# files = ['Reg1.nc']
# minlons = [-60,]
# maxlons = [-30,]
# minlats = [-42,]
# maxlats = [-17.5,]
# fvar = 'FVars/NCEP_R2_RDA.csv'


for file,min_lon, max_lon, min_lat, max_lat in zip(files,minlons,
                                                    maxlons,minlats,maxlats):
    data = aux.get_data(file,fvar,min_lon, max_lon, min_lat, max_lat)
    
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = data[0],data[1],data[2], data[3]
        
    
    tair = data[4]*units(data[4].units).to('K')
    hgt = data[5]*units(data[5].units).to('gpm')
    rhum = data[6]*units(data[6].units)
    omega = data[7]*units(data[7].units).to('Pa/s')
    u = data[8]*units(data[8].units).to('m/s')
    v = data[9]*units(data[9].units).to('m/s')
    slp = data[10]*units(data[10].units).to('Pa')
    pres = tair[VerticalCoordIndexer]*units(tair[VerticalCoordIndexer].units).to('Pa')
    
    box_obj = BoxData(TemperatureData=tair, PressureData=pres,
                     UWindComponentData=u, VWindComponentData=v,
                     OmegaData=omega,
                     LonIndexer=LonIndexer, LatIndexer=LatIndexer, TimeName=TimeName,
                     VerticalCoordIndexer=VerticalCoordIndexer,
                     western_limit=min_lon, eastern_limit=max_lon,
                     southern_limit=min_lat, northern_limit=max_lat)
    
    ec_obj = EnergyContents(box_obj)
    EnergyList = [ec_obj.calc_az()/1e5, ec_obj.calc_ae()/1e5,
                      ec_obj.calc_kz()/1e5,ec_obj.calc_ke()/1e5]
    
    ct_obj = ConversionTerms(box_obj)
    ConversionList = [ct_obj.calc_cz(),ct_obj.calc_ca(),ct_obj.calc_ck(),ct_obj.calc_ce()]
    
    # EnergyTerms,ConversionTerms = [AZ,AE,KZ,KE], [CZ,CA,CK,CE]
    
    Plotter(EnergyList,ConversionList,TimeName,file)


    
# # =============================================================================
# # Test reduced area
# # =============================================================================
# data = aux.get_data(
#         'Reg1.nc','FVars/NCEP_R2_RDA.csv',-30, -25, -20, -15
#         ) 
# LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = data[0],data[1],data[2], data[3]
# u = data[8]*units(data[8].units).to('m/s')
# v = data[9]*units(data[9].units).to('m/s')
# pres = u[VerticalCoordIndexer]*units(tair[VerticalCoordIndexer].units).to('Pa')

# tair,u,v = tair[0][:3],u[0][:3],v[0][:3]

# KZ = Calc_Kz(u,v,pres,LonIndexer,LatIndexer,VerticalCoordIndexer)

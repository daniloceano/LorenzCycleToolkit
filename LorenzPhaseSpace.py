#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:32:27 2022

@author: daniloceano
"""

import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import cmocean
import glob


def plot3D():

    plt.close('all')
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    # Data for a three-dimensional line
    zline = df['Ca']
    xline = df['Ck']
    yline = df['RGe']
    # ax.plot(xline, yline, zline)
    
    # generate a list of (x,y,z) points
    points = np.array([xline,yline,zline]).transpose().reshape(-1,1,3)
    print(points.shape)  # Out: (len(x),1,3)
    
    # set up a list of segments
    segs = np.concatenate([points[:-1],points[1:]],axis=1)
    print(segs.shape)  # Out: ( len(x)-1, 2, 3 )
                      # see what we've done here -- we've mapped our (x,y,z)
                      # points to an array of segment start/end coordinates.
                      # segs[i,0,:] == segs[i-1,1,:]
    lc = Line3DCollection(segs, cmap=cmocean.cm.curl,linewidth=5)
    lc.set_array(yline) # color the segments by our parameter
    ax.add_collection3d(lc)
    
    ax.set_xlabel('Ck')
    ax.set_ylabel('RGe')
    ax.set_zlabel('Ca')
    
    ax.set_xlim(-20,10)
    ax.set_ylim(-10,10)
    ax.set_zlim(-3,10)

    # ax.set_box_aspect([np.ptp(i) for i in df])  # equal aspect ratio
    c = np.arange(len(yline)) / len(yline)  # create some colours
    norm = plt.Normalize(-5,5)
    # lc = LineCollection(yline,cmap='RdBu', norm=norm)
    # c = plt.cm.RdBu(np.linspace(0,max(yline),len(yline)))
    dots = ax.scatter(xline, yline, zline, cmap = cmocean.cm.curl,
                     c=yline, s=100, vmax=5,vmin=-5)
    plt.colorbar(dots)
    
    # set angle
    ax.view_init(25,-75)
    plt.savefig(FigsDir+"test3D.png")
    
if __name__ == "__main__":
 
#     parser = argparse.ArgumentParser(description = "\
# reads an CSV file with all terms from the Lorenz Energy Cycle \
#  (as input from user) and make the deafult figures for the Lorenz energy cycle\
# The transparecy in each box is set to be proportional to the energy tendency,\
#  as well as the arrows are set to be proportional to the conversion rates.")
#     parser.add_argument("outfile", help = "The .csv file containing the \
# results from the main.py program.")

#     args = parser.parse_args()
#     outfile = args.outfile
    outfile = '../LEC_Results/Reg1_NCEP-R2_60W30W42S17S/Reg1_NCEP-R2_60W30W42S17S.csv'
    ResultsSubDirectory = '/'.join(outfile.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/'
    
    # df = pd.read_csv(outfile)
    # df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')
    
    AllResults = glob.glob('/'.join(outfile.split('/')[:-2])+'/*')
    data,names = [],[]
    for result in AllResults:
        file = result +'/'+''.join(result.split('/')[-1])+'.csv'
        data.append(pd.read_csv(file))
        names.append(file.split('/')[-1].split('.')[0].split('_')[0])

    plt.close('all')
    fig = plt.figure()
    ax = plt.axes()
    for df,name in zip(data,names):

        x = df['Ca']
        y = df['Ck']
        z = df['RGe']
        
        ax.set_xlabel('Ck')
        ax.set_ylabel('Ca')
        ax.set_xlim(-5,10)
        ax.set_ylim(-3,10)
    
        plt.plot(x,y,'-o',label=name)
        plt.legend()
        
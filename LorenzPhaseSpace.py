#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:32:27 2022

@author: daniloceano
"""

import pandas as pd
import argparse
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import cmocean
import glob

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def plot3D():

    orig_cmap = cmocean.cm.curl

    plt.close('all')
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    # Data for a three-dimensional line
    xline = df['Ck']
    yline = df['Ca']
    zline = df['RGe']
    
    # generate a list of (x,y,z) points
    points = np.array([xline,yline,zline]).transpose().reshape(-1,1,3)
    print(points.shape)  # Out: (len(x),1,3)
    
    # set up a list of segments
    segs = np.concatenate([points[:-1],points[1:]],axis=1)
    print(segs.shape)  # Out: ( len(x)-1, 2, 3 )
                      # see what we've done here -- we've mapped our (x,y,z)
                      # points to an array of segment start/end coordinates.
                      # segs[i,0,:] == segs[i-1,1,:]
    lc = Line3DCollection(segs, cmap=orig_cmap,linewidth=5)
    lc.set_array(yline) # color the segments by our parameter
    # ax.add_collection3d(lc)
    ax.plot3D(xline, yline, zline,c='gray')
    
    ax.set_xlabel('Ck')
    ax.set_ylabel('Ca')
    ax.set_zlabel('RGe')
    
    ax.set_xlim(-10,5)
    ax.set_ylim(-3,3)
    ax.set_zlim(-5,10)

    c = np.arange(len(yline)) / len(yline)  # create some colours
    # norm = plt.Normalize(-3,10)
    dots = ax.scatter(xline, yline, zline, cmap = orig_cmap,
                      norm=MidpointNormalize(midpoint=0),
                     c=zline, s=100)
    ax.text(xline[0], yline[0], zline[0],'A')
    ax.text(xline.iloc[-1], yline.iloc[-1], zline.iloc[-1], 'Z')
    
    cbar = plt.colorbar(dots)
    cbar.ax.set_ylabel('RGe', rotation=270)

    
    # set angle
    ax.view_init(25,-75)
    plt.savefig(FigsDir+"test3D.png")
    
def plotSurface():
    
    plt.close('all')
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    # Data for a three-dimensional line
    x = df['Ck']
    y = df['Ca']
    z = df['RGe']
    
    ax.set_xlim(-10,5)
    ax.set_ylim(-3,5)
    ax.set_zlim(-3,10)

    norm = plt.Normalize(-3,10)
    ax.set_xlabel('Ck')
    ax.set_ylabel('Ca')
    ax.set_zlabel('RGe')
    
    ax.plot_trisurf(x, y, z, cmap = cmocean.cm.curl, norm=norm)

def plot2D():
    
    # AllResults = glob.glob('/'.join(outfile.split('/')[:-2])+'/*')
    # data,names = [],[]
    # for result in AllResults:
    #     file = result +'/'+''.join(result.split('/')[-1])+'.csv'
    #     data.append(pd.read_csv(file))
    #     names.append(file.split('/')[-1].split('.')[0].split('_')[0])


    x = df['Ca']
    y = df['Ck']
    z = df['RGe']
    
    plt.close('all')
    fig = plt.plot(fig_size=(10,10))
    ax = plt.gca()

    ax.plot(x,y,'-',c='gray',zorder=2,linewidth=3)
    dots = ax.scatter(x,y,c=z,cmap=cmocean.cm.curl,s=100,zorder=100)
    ax.set_xlabel('Ck',fontsize=12)
    ax.set_ylabel('Ca',fontsize=12)
    plt.tick_params(labelsize=8)
    ax.set_xlim(-4,4)
    ax.set_ylim(-10,5)
    
    for i in range(10):
        ax.axhline(y=0+(i/10),zorder=0+(i/5),linewidth=5.3,
                   alpha=0.3-(i/30),c='#65b6fc')
        ax.axhline(y=0-(i/10),zorder=0+(i/5),linewidth=5.3,
                   alpha=0.3-(i/30),c='#65b6fc')
        ax.axvline(x=0+(i/20),zorder=0+(i/5),linewidth=5.3,
               alpha=0.3-(i/30),c='#65b6fc')
        ax.axvline(x=0-(i/20),zorder=0+(i/5),linewidth=5.3,
               alpha=0.3-(i/30),c='#65b6fc')
    
    ax.text(x[0], y[0],'A',zorder=101,fontsize=12)
    ax.text(x.iloc[-1], y.iloc[-1], 'Z',zorder=101,fontsize=12)
    
    cbar = plt.colorbar(dots)
    cbar.ax.set_ylabel('RGe', rotation=270,fontsize=12)
    for t in cbar.ax.get_yticklabels():
         t.set_fontsize(8)
        
    
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
    # outfile = '../LEC_Results/Reg1_NCEP-R2_60W30W42S17S/Reg1_NCEP-R2_60W30W42S17S.csv'
    outfile = '../LEC_Results/Catarina_NCEP-R2_55W36W35S20S/Catarina_NCEP-R2_55W36W35S20S.csv' 
    ResultsSubDirectory = '/'.join(outfile.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/'
    
    df = pd.read_csv(outfile)
    df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')
    

    plot2D()
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

def plot3D(bw=False):

    orig_cmap = cmocean.cm.curl

    plt.close('all')
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    # Data for a three-dimensional line
    xline = smoothed['Ck']
    yline = smoothed['Ca']
    zline = smoothed['RGe']
    
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
    ax.text(xline[0], yline[0], zline[0],'A')
    ax.text(xline.iloc[-1], yline.iloc[-1], zline.iloc[-1], 'Z')
    
    ax.set_xlabel('Ck')
    ax.set_ylabel('Ca')
    ax.set_zlabel('RGe')
    
    ax.set_xlim(-10,5)
    ax.set_ylim(-3,3)
    ax.set_zlim(-5,10)
    
    # set angle
    ax.view_init(25,-75)

    if bw == False:
        ax.plot3D(xline, yline, zline,c='gray')
        divnorm = colors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=8)
        dots = ax.scatter(xline, yline, zline, cmap = orig_cmap,
                          norm=divnorm,
                          c=zline, s=100)
        cbar = plt.colorbar(dots)
        cbar.ax.set_ylabel('RGe', rotation=270)
        plt.savefig(FigsDir+"test3D.png")
        
    else:
        ax.plot3D(xline, yline, zline,c='gray',linewidth=3)
        # ax.add_collection3d(lc)
        plt.savefig(FigsDir+"test3D_bw.png")
    
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

def MarkerSizeKe(df):
    
    msizes = [60,100,140,200,240]
    intervals = [3e5,3.5e5,4.5e5,5e5]


    sizes = []
    for val in df['Ke']:
        if val <= intervals[0]:
            sizes.append(60)
        elif val > intervals[0] and val <= intervals[1]:
            sizes.append(100)
        elif val > intervals[1] and val <= intervals[2]:
            sizes.append(140)
        elif val > intervals[2] and val <= intervals[3]:
            sizes.append(200)
        else:
            sizes.append(240)
    df['sizes'] = sizes
    
    # Plot legend
    labels = ['< '+str(intervals[0]),
              '< '+str(intervals[1]),
              '< '+str(intervals[2]),
              '< '+str(intervals[3]),
              '> '+str(intervals[3])]
    l1 = plt.scatter([],[],c='k', s=msizes[0],label=labels[0])
    l2 = plt.scatter([],[], c='k', s=msizes[1],label=labels[1])
    l3 = plt.scatter([],[],c='k', s=msizes[2],label=labels[2])
    l4 = plt.scatter([],[],c='k', s=msizes[3],label=labels[3])
    l5 = plt.scatter([],[],c='k', s=msizes[4],label=labels[4])
    leg = plt.legend([l1, l2, l3, l4, l5], labels, ncol=1, frameon=True,
                     fontsize=6, handlelength = 0.5, handleheight = 3,
                     borderpad = 1.5, scatteryoffsets = [0.1],
                handletextpad=1.5, title='Eddy Kinect Engery (Ke)', scatterpoints = 1, loc = 1,
                bbox_to_anchor=(-0.22, -0.56, 0.5, 1),labelcolor = '#383838')
    plt.setp(leg.get_title(), color='#383838')
    plt.setp(leg.get_title(),fontsize=5)
    for i in range(len(leg.legendHandles)):
        leg.legendHandles[i].set_color('#383838')
        leg.legendHandles[i].set_edgecolor('gray')
    
    return df

def plot2D():
    
    x = smoothed['Ca']
    y = smoothed['Ck']
    z = smoothed['Ge']
    w = smoothed['Ke']
    
    plt.close('all')
    plt.plot(fig_size=(30,30))
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(right=0.98)
    plt.gcf().subplots_adjust(left=0.135)
    ax = plt.gca()
    
    # Line plot
    ax.plot(x,y,'-',c='gray',zorder=2,linewidth=3)
    
    # Scatter plot
    s = MarkerSizeKe(smoothed)['sizes']
    divnorm = colors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=8)
    dots = ax.scatter(x,y,c=z,cmap=cmocean.cm.curl,s=s,zorder=100,
                      edgecolors='grey', norm=divnorm)
    
    
    # Plot limits
    ax.set_xlim(-4,8)
    ax.set_ylim(-23,8)
    
    # Gradient lines in the center of the plot
    for i in range(7):
        alpha, offsetalpha = 0.3, 20
        offsetx,offsety = 8.8,17.6
        # c,lw = '#65b6fc',2
        c,lw = 'grey',2
        ax.axhline(y=0+(i/offsetx),zorder=0+(i/5),linewidth=lw,
                   alpha=alpha-(i/offsetalpha),c=c)
        ax.axhline(y=0-(i/offsetx),zorder=0+(i/5),linewidth=lw,
                   alpha=alpha-(i/offsetalpha),c=c)
        ax.axvline(x=0+(i/offsety),zorder=0+(i/5),linewidth=lw,
               alpha=alpha-(i/offsetalpha),c=c)
        ax.axvline(x=0-(i/offsety),zorder=0+(i/5),linewidth=lw,
               alpha=alpha-(i/offsetalpha),c=c)
        
    # Labels
    ax.set_xlabel('Conversion from zonal to eddy Kinetic Energy (Ck)',
                  fontsize=10,labelpad=18,c='#383838')
    ax.set_ylabel('Conversion from zonal to eddy Potential Energy (Ca)',
                  fontsize=10,labelpad=20,c='#383838')
    plt.tick_params(labelsize=8)
    
    # Colorbar
    cbar = plt.colorbar(dots, extend='both')
    cbar.ax.set_ylabel('Generation of eddy Potential Energy (Ge)',
                       rotation=270,fontsize=10,verticalalignment='bottom',
                       c='#383838',labelpad=20)
    for t in cbar.ax.get_yticklabels():
         t.set_fontsize(8)

    # Annotate plot
    system = outfile.split('/')[-1].split('_')[0]
    datasource = outfile.split('/')[-1].split('_')[1]
    start, end = str(df['Datetime'][0]),str(df['Datetime'].iloc[-1]) 
    ax.text(0,1.1,'System: '+system+' - Data from: '+datasource,
            fontsize=12,c='#242424',horizontalalignment='left',
            transform=ax.transAxes)
    ax.text(0,1.06,'Start (A):',fontsize=8,c='#242424',
            horizontalalignment='left',transform=ax.transAxes)
    ax.text(0,1.025,'End (Z):',fontsize=8,c='#242424',
            horizontalalignment='left',transform=ax.transAxes)
    ax.text(0.14,1.06,start,fontsize=8,c='#242424',
            horizontalalignment='left',transform=ax.transAxes)
    ax.text(0.14,1.025,end,fontsize=8,c='#242424',
            horizontalalignment='left',transform=ax.transAxes)
    ax.text(-0.11,0.1,'Eddy is providing potential energy \n to the mean flow',
            rotation=90,fontsize=6,horizontalalignment='center',c='#19616C',
            transform=ax.transAxes)
    ax.text(-0.11,0.6,'Eddy is gaining potential energy \n from the mean flow',
            rotation=90,fontsize=6,horizontalalignment='center',c='#CF6D66',
            transform=ax.transAxes)
    ax.text(0.17,-0.11,'Eddy is gaining kinetic energy \n from the mean flow',
            fontsize=6,horizontalalignment='center',c='#CF6D66',
            transform=ax.transAxes)
    ax.text(0.7,-0.11,'Eddy is providing kinetic energy \n to the mean flow',
            fontsize=6,horizontalalignment='center',c='#19616C',
            transform=ax.transAxes)
    ax.text(1.2,0.2,'Subisidence decreases \n eddy potential energy',
            rotation=270,fontsize=6,horizontalalignment='center',c='#19616C'
            ,transform=ax.transAxes)
    ax.text(1.2,0.65,'Latent heat release feeds \n eddy potential energy',
            rotation=270,fontsize=6,horizontalalignment='center',c='#CF6D66',
            transform=ax.transAxes)
    
    # Marking start and end of the system
    ax.text(x[0], y[0],'A',
            zorder=101,fontsize=12,horizontalalignment='center',
            verticalalignment='center')
    ax.text(x.iloc[-1], y.iloc[-1], 'Z',
            zorder=101,fontsize=12,horizontalalignment='center',
            verticalalignment='center')
        
    fname = FigsDir+"LPS.png"
    plt.savefig(fname,dpi=500)
    print(fname+' created!')
    
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
    # outfile = '../LEC_Results/Reg1_NCEP-R2_60W30W42S17S/Reg1_NCEP-R2_60W30W42S17S.csv'
    # outfile = '../LEC_Results/Catarina_NCEP-R2_55W36W35S20S/Catarina_NCEP-R2_55W36W35S20S.csv' 
    ResultsSubDirectory = '/'.join(outfile.split('/')[:-1])
    FigsDir = ResultsSubDirectory+'/Figures/'
    
    df = pd.read_csv(outfile)
    df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')
    smoothed = df.groupby(pd.Grouper(key="Datetime", freq="12H")).mean()

    plot2D()
    # plot3D()
    # plot3D(bw=True)
    
    
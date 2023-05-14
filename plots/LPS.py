'''
 # @ Author: Dnailo Couto de Souza
 # @ Create Time: 2022-06-14 16:32:27
 # @ Modified by: Danilo Couto de Souza
 # @ Modified time: 2023-05-14 19:08:49
 # @ Description: Lorenz Phase Space main function
 '''

import pandas as pd
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import cmocean
import glob
import numpy as np


def MarkerSizeKe(Ke, ke_label, labelsize):

    msizes = [200,400,600,800,1000]
    
    intervals = [3e5,4e5,5e5,6e5]

    sizes = []
    for val in Ke:
        if val <= intervals[0]:
            sizes.append(msizes[0])
        elif val > intervals[0] and val <= intervals[1]:
            sizes.append(msizes[1])
        elif val > intervals[1] and val <= intervals[2]:
            sizes.append(msizes[2])
        elif val > intervals[2] and val <= intervals[3]:
            sizes.append(msizes[3])
        else:
            sizes.append(msizes[4])
    Ke = pd.DataFrame(Ke)
    Ke['sizes'] = sizes
    
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
    leg = plt.legend([l1, l2, l3, l4, l5], labels, ncol=1, frameon=False,
                     fontsize = 10, handlelength = 0.3, handleheight = 4,
                     borderpad = 1.5, scatteryoffsets = [0.1], framealpha = 1,
                handletextpad = 1.5, title = ke_label,
                scatterpoints = 1, loc = 1,
                bbox_to_anchor=(0.73, -0.57, 0.5, 1),labelcolor = '#383838')
    leg._legend_box.align = "center"
    plt.setp(leg.get_title(), color = '#383838')
    plt.setp(leg.get_title(),fontsize = labelsize)
    for i in range(len(leg.legendHandles)):
        leg.legendHandles[i].set_color('#383838')
        leg.legendHandles[i].set_edgecolor('gray')
    
    return Ke

def annotate_plot(ax, **kwargs):

    fontsize = kwargs.get('fontsize', 10)
    title = kwargs.get('title', '')
    datasource = kwargs.get('datasource', '')

    ax.text(0,1.12,'System: '+title+' - Data from: '+datasource,
            fontsize=16,c='#242424',horizontalalignment='left',
            transform=ax.transAxes)
    ax.text(0,1.07,'Start (A):',fontsize=14,c='#242424',
            horizontalalignment='left',transform=ax.transAxes)
    ax.text(0,1.025,'End (Z):',fontsize=14,c='#242424',
            horizontalalignment='left',transform=ax.transAxes)
    ax.text(0.14,1.07,str(kwargs['start']),fontsize=14,c='#242424',
            horizontalalignment='left',transform=ax.transAxes)
    ax.text(0.14,1.025,str(kwargs['end']),fontsize=14,c='#242424',
            horizontalalignment='left',transform=ax.transAxes)
    y_upper = 'Eddy is gaining potential energy \n from the mean flow'
    y_lower = 'Eddy is providing potential energy \n to the mean flow'
    x_left = 'Eddy is gaining kinetic energy \n from the mean flow'
    x_right = 'Eddy is providing kinetic energy \n to the mean flow'
    col_lower = 'Subsidence decreases \n eddy potential energy'
    col_upper = 'Latent heat release feeds \n eddy potential energy'
    lower_left = 'Barotropic instability'
    upper_left = 'Barotropic and baroclinic instabilities'
    lower_right = 'Eddy is feeding the local atmospheric circulation'
    upper_right = 'Baroclinic instability'
        
    ax.text(-0.07,-0.02,y_lower,
            rotation=90,fontsize=fontsize,horizontalalignment='center',c='#19616C',
            transform=ax.transAxes)
    ax.text(-0.07,0.5,y_upper,
            rotation=90,fontsize=fontsize,horizontalalignment='center',c='#CF6D66',
            transform=ax.transAxes)
    ax.text(0.22,-0.07,x_left,
            fontsize=fontsize,horizontalalignment='center',c='#CF6D66',
            transform=ax.transAxes)
    ax.text(0.75,-0.07,x_right,
            fontsize=fontsize,horizontalalignment='center',c='#19616C',
            transform=ax.transAxes)
    ax.text(1.13,0.49,col_lower,
            rotation=270,fontsize=fontsize,horizontalalignment='center',c='#19616C'
            ,transform=ax.transAxes)
    ax.text(1.13,0.75,col_upper,
            rotation=270,fontsize=fontsize,horizontalalignment='center',c='#CF6D66',
            transform=ax.transAxes)
    ax.text(0.22,0.03,lower_left,
            fontsize=fontsize,horizontalalignment='center',c='#660066',
            verticalalignment='center', transform=ax.transAxes)
    ax.text(0.22,0.97,upper_left,
            fontsize=fontsize,horizontalalignment='center',c='#800000',
            verticalalignment='center', transform=ax.transAxes)
    ax.text(0.75,0.03,lower_right,
            fontsize=fontsize,horizontalalignment='center',c='#000066',
            verticalalignment='center', transform=ax.transAxes)
    ax.text(0.75,0.97,upper_right,
            fontsize=fontsize,horizontalalignment='center',c='#660066',
            verticalalignment='center', transform=ax.transAxes)
    
def gradient_lines(ax):
    alpha, offsetalpha = 0.3, 20
    lw, c = 2.5, '#383838'
    offsetx, offsety = 18, 6
    for i in range(7):
        ax.axhline(y=0+(i/offsetx),zorder=0+(i/5),linewidth=lw,
                   alpha=alpha-(i/offsetalpha),c=c)
        ax.axhline(y=0-(i/offsetx),zorder=0+(i/5),linewidth=lw,
                   alpha=alpha-(i/offsetalpha),c=c)
        ax.axvline(x=0+(i/offsety),zorder=0+(i/5),linewidth=lw,
               alpha=alpha-(i/offsetalpha),c=c)
        ax.axvline(x=0-(i/offsety),zorder=0+(i/5),linewidth=lw,
               alpha=alpha-(i/offsetalpha),c=c)
        # Vertical line showing when Ca is more important than Ck
        n = 15
        plt.plot(np.arange(0-(i/offsety),-n-(i/offsety),-1),
                 np.arange(0,n), c=c,zorder=1, 
                 alpha=0.2-(i/offsetalpha*.5))
        plt.plot(np.arange(0+(i/offsety),-n+(i/offsety),-1),
                 np.arange(0,n), c=c,zorder=1, linewidth=lw,
                 alpha=0.2-(i/offsetalpha*.5))
        
def get_max_vals(term, **kwargs):
    max_values = []
    for term_list in kwargs['terms']:
        max_term = term_list[term].max()
        max_values.append(max_term)
    return max(max_values)

def get_min_vals(term, **kwargs):
    min_values = []
    for term_list in kwargs['terms']:
        min_term = term_list[term].min()
        min_values.append(min_term)
    return min(min_values)

def limits_zoomed(ax, **kwargs):
    # Get limits
    minCk, maxCk =  get_min_vals('Ck', **kwargs), get_max_vals('Ck', **kwargs)
    minCa, maxCa = get_min_vals('Ca', **kwargs), get_max_vals('Ca', **kwargs)
    # Plot limits for Ck
    if minCk < -1:
        minLimitCk = minCk*1.3
    else:
        minLimitCk = -1
    if maxCk > 3:
        maxLimitCk = maxCk*1.3
    else:
        maxLimitCk = 3
    ax.set_xlim(minLimitCk,maxLimitCk)
    # Plot limits for Ca
    if minCa < -0.5:
        minLimitCa =minCa*1.3
    else:
        minLimitCa = -0.5
    if maxCa > 2:
        maxLimitCa = maxCa*1.3
    else:
        maxLimitCa = 2.5
    ax.set_ylim(minLimitCa,maxLimitCa)

def LorenzPhaseSpace(ax, zoom=False, example=False, **kwargs):
    
    if zoom == False:
        # Labels
        ck_label = 'Conversion from zonal to eddy Kinetic Energy (Ck - '+r' $W\,m^{-2})$'
        ca_label = 'Conversion from zonal to eddy Potential Energy (Ca - '+r' $W\,m^{-2})$'
        ge_label = 'Generation of eddy Potential Energy (Ge - '+r' $W\,m^{-2})$'
        ke_label = 'Eddy Kinect\n    Energy\n(Ke - '+r' $J\,m^{-2})$'
        ax.set_xlim(-30,30)
        ax.set_ylim(-6,12)
        # Write physical meaning of each quadrant
        plt.tick_params(labelsize=kwargs.get('fontsize', 10))
        annotate_plot(ax, **kwargs)
        # Gradient lines in the center of the plot
        gradient_lines(ax)
        # limits for Ge
        norm = colors.TwoSlopeNorm(vmin=-7, vcenter=0, vmax=15)
        # pad for labels
        labelpad = 38
        
    else:
        # Labels
        ck_label = 'Ck - '+r' $W\,m^{-2}$'
        ca_label = 'Ca - '+r' $W\,m^{-2}$'
        ge_label = 'Ge - '+r' $W\,m^{-2}$'
        ke_label = 'Ke - '+r' $J\,m^{-2}$'
        limits_zoomed(ax, **kwargs)
        minGe, maxGe =  get_min_vals('Ge', **kwargs), get_max_vals('Ge', **kwargs)
        if minGe > 0:
            minGe = -1
        elif maxGe < 0:
            maxGe = 1
        norm =  colors.TwoSlopeNorm(vmin=minGe, vcenter=0, vmax=maxGe)
        # pad for labels
        labelpad = 5
        # Lines in the center of the plot
        c,lw,alpha = '#383838',20,0.2
        ax.axhline(y=0,linewidth=lw,c=c,alpha=alpha,zorder=1)
        ax.axvline(x=0,linewidth=lw,c=c,alpha=alpha,zorder=1)
        ax.plot(range(0,-40,-1),range(0,40,1),
                linewidth=lw/3,c=c,alpha=alpha,zorder=1)

    labelsize = kwargs.get('labelsize', 14)

    # Loop through all list of terms in kwargs
    for term_list in kwargs['terms']:    

        Ca = term_list['Ca']
        Ck = term_list['Ck']
        Ge = term_list['Ge']
        Ke = term_list['Ke']

        # Line plot
        ax.plot(Ck,Ca,'-',c='gray',linewidth=3)
        
        # Label for eddy kinectinc energy (Ke)
        s = MarkerSizeKe(Ke, ke_label, labelsize)['sizes']
    
        # arrows connecting dots
        ax.quiver(Ck[:-1], Ca[:-1],
                (Ck[1:].values-Ck[:-1].values)*.97,
                (Ca[1:].values-Ca[:-1].values)*.97,
                angles='xy', scale_units='xy',
                scale=1, color='k')

        # plot the moment of maximum intensity
        ax.scatter(Ck.loc[s.idxmax()],Ca.loc[s.idxmax()],
                c='None',s=s.loc[s.idxmax()]*1.1,
                zorder=100,edgecolors='k', linewidth=3)
        
        # Circles representing Ck on x-axis and Ca on y-axis, while the
        # colors represent Ge and the circle sizes, Ke.
        dots = ax.scatter(Ck,Ca,c=Ge,cmap=cmocean.cm.curl,s=s,zorder=100,
                        edgecolors='grey', norm=norm)
        
        # Marking start and end of the system
        ax.text(Ck[0], Ca[0],'A',
                zorder=101,fontsize=22,horizontalalignment='center',
                verticalalignment='center')
        ax.text(Ck.iloc[-1], Ca.iloc[-1], 'Z',
                zorder=101,fontsize=22,horizontalalignment='center',
                verticalalignment='center')

    # Colorbar
    cax = ax.inset_axes([ax.get_position().x1+0.12,
                    ax.get_position().y0+0.35,0.02, ax.get_position().height/1.5])
    cbar = plt.colorbar(dots, extend='both',cax=cax)
    
    # Write labels
    ax.set_xlabel(ck_label, fontsize=labelsize,labelpad=labelpad,c='#383838')
    ax.set_ylabel(ca_label, fontsize=labelsize,labelpad=labelpad,c='#383838')
    cbar.ax.set_ylabel(ge_label, rotation=270,fontsize=labelsize,
                       verticalalignment='bottom', c='#383838',
                       labelpad=labelpad, y=0.59)
    for t in cbar.ax.get_yticklabels():
         t.set_fontsize(10) 
    
    plt.subplots_adjust(right=0.84, bottom=0.1)
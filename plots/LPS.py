# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    LPS.py                                             :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo <danilo.oceano@gmail.com>           +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/06/14 16:32:27 by Danilo            #+#    #+#              #
#    Updated: 2023/07/20 17:41:50 by Danilo           ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import pandas as pd
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import cmocean
import numpy as np


def calculate_marker_size(term):
    msizes = [200, 400, 600, 800, 1000]
    intervals = [3e5, 4e5, 5e5, 6e5]
    sizes = []
    for val in term:
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
    return pd.Series(sizes)

def plot_legend(ax, cmap):
    msizes = [200, 400, 600, 800, 1000]
    intervals = [3e5, 4e5, 5e5, 6e5]
    labels = ['< ' + str(intervals[0]),
              '< ' + str(intervals[1]),
              '< ' + str(intervals[2]),
              '< ' + str(intervals[3]),
              '> ' + str(intervals[3])]
    for i, label in enumerate(labels):
        ax.scatter([], [], c=cmap(msizes[i]), s=msizes[i], label=label)

    ax.legend(title='Eddy Potential Energy (Ge - $W\,m^{-2}$)',
              fontsize=10, loc='lower left', bbox_to_anchor=(0.73, -0.57, 0.5, 1),
              labelcolor='#383838', frameon=False, handlelength=0.3, handleheight=4,
              borderpad=1.5, scatteryoffsets=[0.1], framealpha=1,
              handletextpad=1.5, scatterpoints=1)

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

def get_labels(label_type, zoom=False):
    labels_dict = {}

    if label_type == 'mixed':
        labels_dict['y_upper'] = 'Eddy is gaining potential energy \n from the mean flow'
        labels_dict['y_lower'] = 'Eddy is providing potential energy \n to the mean flow'
        labels_dict['x_left'] = 'Eddy is gaining kinetic energy \n from the mean flow'
        labels_dict['x_right'] = 'Eddy is providing kinetic energy \n to the mean flow'
        labels_dict['col_lower'] = 'Subsidence decreases \n eddy potential energy'
        labels_dict['col_upper'] = 'Latent heat release feeds \n eddy potential energy'
        labels_dict['lower_left'] = 'Barotropic instability'
        labels_dict['upper_left'] = 'Barotropic and baroclinic instabilities'
        labels_dict['lower_right'] = 'Eddy is feeding the local atmospheric circulation'
        labels_dict['upper_right'] = 'Baroclinic instability'

        if zoom == False:
            labels_dict['ck_label'] = 'Conversion from zonal to eddy Kinetic Energy (Ck - $W\,m^{-2})$'
            labels_dict['ca_label'] = 'Conversion from zonal to eddy Potential Energy (Ca - $W\,m^{-2})$'
            labels_dict['ge_label'] = 'Generation of eddy Potential Energy (Ge - $W\,m^{-2})$'
            labels_dict['ke_label'] = 'Eddy Kinect\n    Energy\n(Ke - $J\,m^{-2})$'
        elif zoom == True:
            labels_dict['ck_label'] = 'Ck - $W\,m^{-2}$'
            labels_dict['ca_label'] = 'Ca - $W\,m^{-2}$'
            labels_dict['ge_label'] = 'Ge - $W\,m^{-2}$'
            labels_dict['ke_label'] = 'Ke - $J\,m^{-2}$'

    if label_type == 'baroclinic':
        labels_dict['y_upper'] = ''
        labels_dict['y_lower'] = ''
        labels_dict['x_left'] = ''
        labels_dict['x_right'] = ''
        labels_dict['col_lower'] = ''
        labels_dict['col_upper'] = ''
        labels_dict['lower_left'] = ''
        labels_dict['upper_left'] = ''
        labels_dict['lower_right'] = ''
        labels_dict['upper_right'] = ''
        
        if zoom == False:
            labels_dict['ck_label'] = 'Conversion from zonal to eddy Kinetic Energy (Ck - $W\,m^{-2})$'
            labels_dict['ca_label'] = 'Conversion from zonal to eddy Potential Energy (Ca - $W\,m^{-2})$'
            labels_dict['ge_label'] = 'Generation of eddy Potential Energy (Ge - $W\,m^{-2})$'
            labels_dict['ke_label'] = 'Eddy Kinect\n    Energy\n(Ke - $J\,m^{-2})$'
        elif zoom == True:
            labels_dict['ck_label'] = 'Ck - $W\,m^{-2}$'
            labels_dict['ca_label'] = 'Ca - $W\,m^{-2}$'
            labels_dict['ge_label'] = 'Ge - $W\,m^{-2}$'
            labels_dict['ke_label'] = 'Ke - $J\,m^{-2}$'

    return labels_dict


def LorenzPhaseSpace(ax, type, zoom=False, example=False, **kwargs):
    
    if zoom == False:
        #Limits
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

        # whether to extend cbar
        extend = 'both'
        
    else:
        # Limits
        limits_zoomed(ax, **kwargs)
        minGe, maxGe =  get_min_vals('Ge', **kwargs), get_max_vals('Ge', **kwargs)
        if abs(maxGe) > abs(minGe):
            minGe = -maxGe
        else:
            maxGe = -minGe
        norm =  colors.TwoSlopeNorm(vmin=minGe, vcenter=0, vmax=maxGe)

        # pad for labels
        labelpad = 5

        # Lines in the center of the plot
        c,lw,alpha = '#383838',20,0.2
        ax.axhline(y=0,linewidth=lw,c=c,alpha=alpha,zorder=1)
        ax.axvline(x=0,linewidth=lw,c=c,alpha=alpha,zorder=1)
        ax.plot(range(0,-40,-1),range(0,40,1),
                linewidth=lw/3,c=c,alpha=alpha,zorder=1)
        extend = 'neither'
    
    labels = get_labels(label_type=type, zoom=zoom)
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
        marker_sizes = calculate_marker_size(Ke)
        # s = MarkerSizeKe(Ke, ke_label, labelsize)['sizes']
    
        # arrows connecting dots
        ax.quiver(Ck[:-1], Ca[:-1],
                (Ck[1:].values-Ck[:-1].values)*.97,
                (Ca[1:].values-Ca[:-1].values)*.97,
                angles='xy', scale_units='xy',
                scale=1, color='k')

        # plot the moment of maximum intensity
        ax.scatter(Ck.loc[marker_sizes.idxmax()],Ca.loc[marker_sizes.idxmax()],
                c='None',s=marker_sizes.loc[marker_sizes.idxmax()]*1.1,
                zorder=100,edgecolors='k', linewidth=3)
        
        # Circles representing Ck on x-axis and Ca on y-axis, while the
        # colors represent Ge and the circle sizes, Ke.
        dots = ax.scatter(Ck,Ca,c=Ge,cmap=cmocean.cm.curl,s=marker_sizes,zorder=100,
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
    cbar = plt.colorbar(dots, extend=extend,cax=cax)
    
    # Write labels
    ax.set_xlabel(labels['ck_label'], fontsize=labelsize,labelpad=labelpad,c='#383838')
    ax.set_ylabel(labels['ca_label'], fontsize=labelsize,labelpad=labelpad,c='#383838')
    cbar.ax.set_ylabel(labels['ge_label'], rotation=270,fontsize=labelsize,
                       verticalalignment='bottom', c='#383838',
                       labelpad=labelpad, y=0.59)
    
    for t in cbar.ax.get_yticklabels():
         t.set_fontsize(10) 
    
    plt.subplots_adjust(right=0.84, bottom=0.1)

if __name__ == '__main__':

    outfile = '../inputs/sample_results.csv'

    df = pd.read_csv(outfile, index_col=[0])
    df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')

    # Set datetime to the date range
    start = pd.to_datetime(df['Datetime'].iloc[0]).strftime('%Y-%m-%d %H:%M')
    end = pd.to_datetime(df['Datetime'].iloc[-1]).strftime('%Y-%m-%d %H:%M')

    # Plot example
    kwargs = {'terms':[], 'title': 'sample', 'datasource': 'sample', 'start': start, 'end': end}
    terms = {'Ca': df['Ca'], 'Ck': df['Ck'],  'Ge': df['Ge'], 'Ke': df['Ke']}
    kwargs['terms'].append(terms) 

    for type in ['mixed', 'baroclinic', 'barotropic']:
        for zoom in [False, True]:
            plt.close('all')
            plt.figure(figsize=(10,10))
            ax = plt.gca()
            LorenzPhaseSpace(ax, type, zoom=zoom, **kwargs)
            zoom_suffix = "_zoom" if zoom else ""
            fname = f"./LPS_test_{type}{zoom_suffix}.png"
            with plt.rc_context({'savefig.dpi': 500}):
                    plt.savefig(fname)
            print(f"{fname} created!")

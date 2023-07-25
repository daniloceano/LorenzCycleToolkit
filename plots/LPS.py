# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    LPS.py                                             :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo <danilo.oceano@gmail.com>           +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/06/14 16:32:27 by Danilo            #+#    #+#              #
#    Updated: 2023/07/24 20:38:53 by Danilo           ###   ########.fr        #
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

def plot_legend(ax):
    msizes = [200, 400, 600, 800, 1000]
    intervals = [3e5, 4e5, 5e5, 6e5]
    labels = ['< ' + str(intervals[0]),
              '< ' + str(intervals[1]),
              '< ' + str(intervals[2]),
              '< ' + str(intervals[3]),
              '> ' + str(intervals[3])]

    # Create separate scatter plots for each size category
    for i in range(len(msizes)):
        ax.scatter([], [], c='#383838', s=msizes[i], label=labels[i])

    ax.legend(title='Eddy Kinetic Energy \n(Ke - $W\,m^{-2}$)',
              fontsize=10, loc='lower left', bbox_to_anchor=(0.97, 0, 0.5, 1),
              labelcolor='#383838', frameon=False, handlelength=0.3, handleheight=4,
              borderpad=1.5, scatteryoffsets=[0.1], framealpha=1,
              handletextpad=1.5, scatterpoints=1)

def annotate_plot(ax, LPS_type, zoom, **kwargs):
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
    
    labels = get_labels(LPS_type, zoom)
        
    # Centering text annotations on y-axis
    yticks, xticks = ax.get_yticks(), ax.get_xticks()
    y_tick_0 = len(yticks) // 2
    y_offset = 0.5 * (yticks[y_tick_0] - yticks[-1])  # Half the distance between two consecutive y-ticks
    x_tick_pos = xticks[0] - ((xticks[1] - xticks[0])/12)

    ax.text(x_tick_pos, yticks[0] - y_offset, labels['y_lower'], rotation=90, fontsize=fontsize,
            horizontalalignment='center', c='#19616C', verticalalignment='center')
    ax.text(x_tick_pos, yticks[-1] + y_offset, labels['y_upper'], rotation=90, fontsize=fontsize,
            horizontalalignment='center', c='#CF6D66', verticalalignment='center')
    
    ax.text(0.22,-0.07, labels['x_left'], fontsize=fontsize,
             horizontalalignment='center', c='#CF6D66', transform=ax.transAxes)
    ax.text(0.75,-0.07,labels['x_right'], fontsize=fontsize,
            horizontalalignment='center', c='#19616C', transform=ax.transAxes)
    
    ax.text(1.13,0.49, labels['col_lower'], rotation=270, fontsize=fontsize, 
            horizontalalignment='center', c='#19616C', transform=ax.transAxes)
    ax.text(1.13,0.75, labels['col_upper'], rotation=270,fontsize=fontsize,
            horizontalalignment='center', c='#CF6D66', transform=ax.transAxes)
    
    ax.text(0.22,0.03, labels['lower_left'], fontsize=fontsize, horizontalalignment='center',
            c='#660066', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.22,0.97, labels['upper_left'], fontsize=fontsize,horizontalalignment='center',
            c='#800000', verticalalignment='center', transform=ax.transAxes)
    
    ax.text(0.75,0.03, labels['lower_right'], fontsize=fontsize,horizontalalignment='center',
            c='#000066', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.75,0.97,labels['upper_right'], fontsize=fontsize,horizontalalignment='center',
            c='#660066', verticalalignment='center', transform=ax.transAxes)
    
def gradient_lines(ax, LPS_type):
    lw, c = 0.5, '#383838'
    num_lines = 20
    x_ticks = ax.get_xticks()
    y_ticks = ax.get_yticks()

    x_previous0 = x_ticks[int((len(x_ticks))/2)-1] * 0.17
    y_previous0 = y_ticks[int((len(y_ticks))/2)-1] * 0.17

    x_offsets = np.linspace(x_previous0, 0, num_lines)
    y_offsets = np.linspace(y_previous0, 0, num_lines)

    alpha_values = np.linspace(0, 0.6, num_lines)

    for i, alpha in enumerate(alpha_values):
        ax.axhline(y=0 + y_offsets[i], linewidth=lw, alpha=alpha, c=c)
        ax.axhline(y=0 - y_offsets[i], linewidth=lw, alpha=alpha, c=c)
        ax.axvline(x=0 + x_offsets[i], linewidth=lw, alpha=alpha, c=c)
        ax.axvline(x=0 - x_offsets[i], linewidth=lw, alpha=alpha, c=c)

    # Diagonal line
    if LPS_type == 'mixed':
        y_ticks = -x_ticks
        for i, alpha in enumerate(alpha_values):
            x, y = x_offsets[i], y_offsets[i]
            ax.plot([x, -x_ticks[-1] + x], [y, -y_ticks[-1] + y], linewidth=lw, alpha=alpha, c=c)
            ax.plot([-x, -x_ticks[-1] - x], [-y, -y_ticks[-1] - y], linewidth=lw, alpha=alpha, c=c)
    
        
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
    min_x, max_x =  get_min_vals('x_axis', **kwargs), get_max_vals('x_axis', **kwargs)
    min_y, max_y = get_min_vals('y_axis', **kwargs), get_max_vals('y_axis', **kwargs)

    # Plot limits for x_label
    if min_x < -1:
        xlabel_min = min_x*1.3
    else:
        xlabel_min = -1
    if max_x > 3:
        xlabel_max = max_x*1.3
    else:
        xlabel_max = 3
    ax.set_xlim(xlabel_min,xlabel_max)

    # Plot limits for y-axis
    if min_y < -0.5:
        ylabel_min = min_y*1.3
    else:
        ylabel_min = -0.5
    if max_y > 2:
        ylabel_max = max_y*1.3
    else:
        ylabel_max = 2.5
    ax.set_ylim(ylabel_min,ylabel_max)

def get_labels(LPS_type, zoom=False):
    labels_dict = {}

    if LPS_type == 'mixed':
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
            labels_dict['x_label'] = 'Conversion from zonal to eddy Kinetic Energy (Ck - $W\,m^{-2})$'
            labels_dict['y_label'] = 'Conversion from zonal to eddy Potential Energy (Ca - $W\,m^{-2})$'
            labels_dict['color_label'] = 'Generation of eddy Potential Energy (Ge - $W\,m^{-2})$'
            labels_dict['size_label'] = 'Eddy Kinect\n    Energy\n(Ke - $J\,m^{-2})$'
        elif zoom == True:
            labels_dict['x_label'] = 'Ck - $W\,m^{-2}$'
            labels_dict['y_label'] = 'Ca - $W\,m^{-2}$'
            labels_dict['color_label'] = 'Ge - $W\,m^{-2}$'
            labels_dict['size_label'] = 'Ke - $J\,m^{-2}$'

    elif LPS_type == 'baroclinic':
        labels_dict['y_upper'] = 'Zonal temperature gradient feeds \n eddy potential energy'
        labels_dict['y_lower'] = 'Eddy potential energy feeds \n zonal temperature gradient'
        labels_dict['x_left'] = 'Meridional temperature gradient feeds \n eddy kinetic energy'
        labels_dict['x_right'] = 'Eddy kinetic energy consumes \n meridional temperature gradient'
        labels_dict['col_lower'] = 'Subsidence decreases \n eddy potential energy'
        labels_dict['col_upper'] = 'Latent heat release feeds \n eddy potential energy'
        labels_dict['lower_left'] = 'Baroclinic stability'
        labels_dict['upper_left'] = ''
        labels_dict['lower_right'] = ''
        labels_dict['upper_right'] = 'Baroclinic instability'
        
        if zoom == False:
            labels_dict['x_label'] = 'Conversion from zonal to eddy Kinetic Energy (Ce - $W\,m^{-2})$'
            labels_dict['y_label'] = 'Conversion from zonal to eddy Potential Energy (Ca - $W\,m^{-2})$'
            labels_dict['color_label'] = 'Generation of eddy Potential Energy (Ge - $W\,m^{-2})$'
            labels_dict['size_label'] = 'Eddy Kinect\n    Energy\n(Ke - $J\,m^{-2})$'
        elif zoom == True:
            labels_dict['x_label'] = 'Ce - $W\,m^{-2}$'
            labels_dict['y_label'] = 'Ca - $W\,m^{-2}$'
            labels_dict['color_label'] = 'Ge - $W\,m^{-2}$'
            labels_dict['size_label'] = 'Ke - $J\,m^{-2}$'

    elif LPS_type == 'barotropic':
        labels_dict['y_upper'] = 'Importation of Kinectic Energy'
        labels_dict['y_lower'] = 'Exportation of Kinectic Energy'
        labels_dict['x_left'] = 'Eddy is gaining kinetic energy \n from the mean flow'
        labels_dict['x_right'] = 'Eddy is providing kinetic energy \n to the mean flow'
        labels_dict['col_lower'] = 'Subsidence decreases \n eddy potential energy'
        labels_dict['col_upper'] = 'Latent heat release feeds \n eddy potential energy'
        labels_dict['lower_left'] = 'Barotropic instability wihtout \n downstream development'
        labels_dict['upper_left'] = 'Barotropic instability and \n downstream development'
        labels_dict['lower_right'] = 'Barotropic stability without \n downstream development'
        labels_dict['upper_right'] = 'Barotropic stability and \n downstream development'
        
        if zoom == False:
            labels_dict['x_label'] = 'Conversion from zonal to eddy Kinetic Energy (Ck - $W\,m^{-2})$'
            labels_dict['y_label'] = ' Kinetic Energy transport across boundaries (BKz - $W\,m^{-2})$'
            labels_dict['color_label'] = 'Generation of eddy Potential Energy (Ge - $W\,m^{-2})$'
            labels_dict['size_label'] = 'Eddy Kinect\n    Energy\n(Ke - $J\,m^{-2})$'
        elif zoom == True:
            labels_dict['x_label'] = 'Ck - $W\,m^{-2}$'
            labels_dict['y_label'] = 'Bkz - $W\,m^{-2}$'
            labels_dict['color_label'] = 'Ge - $W\,m^{-2}$'
            labels_dict['size_label'] = 'Ke - $J\,m^{-2}$'

    return labels_dict


def LorenzPhaseSpace(ax, LPS_type, zoom=False, example=False, **kwargs):
    if zoom == False:
        #Limits
        if LPS_type == 'mixed':
            ax.set_xlim(-70,70)
            ax.set_ylim(-20,20)
        elif LPS_type == 'baroclinic':
            ax.set_xlim(-70,70)
            ax.set_ylim(-20,20)
        elif LPS_type == 'barotropic':
            ax.set_xlim(-70,70)
            ax.set_ylim(-200,200)
        
        # Write physical meaning of each quadrant
        plt.tick_params(labelsize=kwargs.get('fontsize', 10))
        annotate_plot(ax, LPS_type, zoom, **kwargs)

        # Gradient lines in the center of the plot
        gradient_lines(ax, LPS_type)

        # limits for coors
        norm = colors.TwoSlopeNorm(vmin=-30, vcenter=0, vmax=30)

        # pad for labels
        labelpad = 38

        # whether to extend cbar
        extend = 'both'
        
    else:
        # Limits
        limits_zoomed(ax, **kwargs)
        min_colors, max_colors =  get_min_vals('circles_colors', **kwargs), get_max_vals('circles_colors', **kwargs)
        if abs(max_colors) > abs(min_colors):
            min_colors = -max_colors
        else:
            max_colors = -min_colors
        norm =  colors.TwoSlopeNorm(vmin=min_colors, vcenter=0, vmax=max_colors)

        # pad for labels
        labelpad = 5

        # Lines in the center of the plot
        c,lw,alpha = '#383838',20,0.2
        ax.axhline(y=0,linewidth=lw,c=c,alpha=alpha,zorder=1)
        ax.axvline(x=0,linewidth=lw,c=c,alpha=alpha,zorder=1)
        extend = 'neither'

        # Vertical lines for mixed LPS
        if LPS_type == 'mixed':
            min_x = int(round(ax.get_xlim()[0],-1)+5)
            max_y = int(round(ax.get_ylim()[1],-1)-5)
            if abs(min_x) > abs(max_y):
                max_y = -min_x
            else:
                min_x = -max_y
            ax.plot(range(0, min_x,-1),range(0, max_y, 1), linewidth=lw/3,
                    c=c, alpha=alpha,zorder=1)
    
    labels = get_labels(LPS_type, zoom)
    labelsize = kwargs.get('labelsize', 14)

    # Loop through all list of terms in kwargs: this allows to a plotting multiple systems at once
    for term_list in kwargs['terms']:  
        y_axis = term_list['y_axis']
        x_axis = term_list['x_axis']
        circles_colors = term_list['circles_colors']
        circles_sizes = term_list['circles_size']

        # Line plot
        ax.plot(x_axis, y_axis,'-',c='gray',linewidth=3)
        
        # Label for eddy kinectinc energy (Ke)
        marker_sizes = calculate_marker_size(circles_sizes)

        # arrows connecting dots
        ax.quiver(x_axis[:-1], y_axis[:-1],
                (x_axis[1:].values - x_axis[:-1].values)*.97,
                (y_axis[1:].values-y_axis[:-1].values)*.97,
                angles='xy', scale_units='xy',
                scale=1, color='k')

        # plot the moment of maximum intensity
        ax.scatter(x_axis.loc[marker_sizes.idxmax()],y_axis.loc[marker_sizes.idxmax()],
                c='None',s=marker_sizes.loc[marker_sizes.idxmax()]*1.1,
                zorder=100,edgecolors='k', linewidth=3)
        
        # Circles colors represent Ge and the circle sizes represent Ke.
        dots = ax.scatter(x_axis, y_axis, c=circles_colors, cmap=cmocean.cm.curl,s=marker_sizes,zorder=100,
                        edgecolors='grey', norm=norm)
        
    # Marking start and end of the system
    ax.text(x_axis[0], y_axis[0],'A',
            zorder=101,fontsize=22,horizontalalignment='center',
            verticalalignment='center')
    ax.text(x_axis.iloc[-1], y_axis.iloc[-1], 'Z',
            zorder=101,fontsize=22,horizontalalignment='center',
            verticalalignment='center')

    # Colorbar
    cax = ax.inset_axes([ax.get_position().x1+0.12,
                    ax.get_position().y0+0.35,0.02, ax.get_position().height/1.5])
    cbar = plt.colorbar(dots, extend=extend,cax=cax)
    
    # Write labels
    ax.set_xlabel(labels['x_label'], fontsize=labelsize,labelpad=labelpad,c='#383838')
    ax.set_ylabel(labels['y_label'], fontsize=labelsize,labelpad=labelpad,c='#383838')
    cbar.ax.set_ylabel(labels['color_label'], rotation=270,fontsize=labelsize,
                    verticalalignment='bottom', c='#383838',
                    labelpad=labelpad, y=0.59)
    
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(10) 

    plot_legend(ax)
    
    plt.subplots_adjust(right=0.84, bottom=0.1)

if __name__ == '__main__':

    outfile = '../inputs/sample_results.csv'

    df = pd.read_csv(outfile, index_col=[0])
    df['Datetime'] = pd.to_datetime(df.Date) + pd.to_timedelta(df.Hour, unit='h')

    # Set datetime to the date range
    start = pd.to_datetime(df['Datetime'].iloc[0]).strftime('%Y-%m-%d %H:%M')
    end = pd.to_datetime(df['Datetime'].iloc[-1]).strftime('%Y-%m-%d %H:%M')

    # Kwargs for plotting LPS example
    kwargs = {'terms':[], 'title': 'sample', 'datasource': 'sample', 'start': start, 'end': end}

    # Kwargs for each kind of plot
    terms_mixed = {'y_axis': df['Ca'], 'x_axis': df['Ck'],
                     'circles_colors': df['Ge'], 'circles_size': df['Ke']}
    term_baroclinic = {'y_axis': df['Ca'], 'x_axis': df['Ce'],
                     'circles_colors': df['Ge'], 'circles_size': df['Ke']}
    terms_barotropic = {'y_axis': df['BKz'], 'x_axis': df['Ck'],
                     'circles_colors': df['Ge'], 'circles_size': df['Ke']}

    for LPS_type, term_list in zip(['mixed', 'baroclinic', 'barotropic'],
                                [terms_mixed, term_baroclinic, terms_barotropic]):
        kwargs['terms'] = []
        kwargs['terms'].append(term_list)

        print("\n--------------------------")
        print(f"Plotting for {LPS_type} terms\n")
        for term in kwargs['terms'][0]:
            print(f"axis: {term} term: {kwargs['terms'][0][term].name}")
        
        for zoom in [False, True]:
            plt.close('all')
            plt.figure(figsize=(10,10))
            ax = plt.gca()
            LorenzPhaseSpace(ax, LPS_type, zoom=zoom, **kwargs)
            zoom_suffix = "_zoom" if zoom else ""
            fname = f"./LPS_test_{LPS_type}{zoom_suffix}.png"
            with plt.rc_context({'savefig.dpi': 500}):
                    plt.savefig(fname)
            print(f"{fname} created!")


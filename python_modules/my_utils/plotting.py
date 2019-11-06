#!/usr/bin/env python3


#-------------------------------------------
# All matplotlib plotting related stuff
#-------------------------------------------

import matplotlib as mpl


#====================================================
def setplotparams_single_plot(  left=0.07, 
                                right=0.97, 
                                bottom=0.10, 
                                top=0.92, 
                                wspace=0.15, 
                                hspace=0.12,
                                tinyfont=10,
                                smallfont=12,
                                bigfont=16,
                                for_presentation=False
                            ):
#====================================================
    """
    Set rcParams for nice looking plots.
    Optimized for a figure with a single subplot.
    For presentation: make all text bigger
    """
    # Plot parameters
    params = {
        'axes.labelsize': smallfont,
        'axes.titlesize': smallfont*1.1,
        'font.size': smallfont,
        'font.family': 'serif', 
        'font.serif': 'Computer Modern',
        'legend.fontsize': smallfont,
        'xtick.labelsize': tinyfont,
        'ytick.labelsize': tinyfont,
        'text.usetex': True,
        'figure.subplot.left'    : left,
        'figure.subplot.right'   : right,
        'figure.subplot.bottom'  : bottom,
        'figure.subplot.top'     : top,
        'figure.subplot.wspace'  : wspace,
        'figure.subplot.hspace'  : hspace,
        'figure.figsize' : (5, 5.5), 
        'figure.dpi' : 300,
        'lines.markersize' : 6,
        'lines.linewidth' : 3.
    }

    if for_presentation:
        params['axes.labelsize'] *= 1.5
        params['axes.titlesize'] *= 1.5
        params['font.size'] *= 1.5
        params['legend.fontsize'] *= 1.4
        params['xtick.labelsize'] *= 1.2
        params['ytick.labelsize'] *= 1.2

    mpl.rcParams.update(params)

    return






#======================================================
def setplotparams_multiple_plots(
                                    left=0.05, 
                                    right=0.97, 
                                    bottom=0.10, 
                                    top=0.92, 
                                    wspace=0.15, 
                                    hspace=0.15,
                                    tinyfont=10,
                                    smallfont=12,
                                    bigfont=16,
                                    for_presentation=False
                                ):
#======================================================
    """
    Set rcParams for nice looking plots. Intended for about figsize 5 x 5.5 per 
    subplot Set those when calling plt.figure()
    For presentation: make all text bigger
    """
    # Plot parameters
    params = {
        'axes.labelsize': smallfont,
        'axes.titlesize': smallfont*1.1,
        'font.size': smallfont,
        'font.family': 'serif', 
        'font.serif': 'Computer Modern',
        'legend.fontsize': smallfont,
        'xtick.labelsize': tinyfont,
        'ytick.labelsize': tinyfont,
        'text.usetex': True,
        'figure.subplot.left'    : left,
        'figure.subplot.right'   : right,
        'figure.subplot.bottom'  : bottom,
        'figure.subplot.top'     : top,
        'figure.subplot.wspace'  : wspace,
        'figure.subplot.hspace'  : hspace,
        'figure.dpi' : 300,
        'lines.markersize' : 6,
        'lines.linewidth' : 3.
    }

    if for_presentation:
        params['axes.labelsize'] *= 1.5
        params['axes.titlesize'] *= 1.5
        params['font.size'] *= 1.5
        params['legend.fontsize'] *= 1.4
        params['xtick.labelsize'] *= 1.2
        params['ytick.labelsize'] *= 1.2


    mpl.rcParams.update(params)

    return



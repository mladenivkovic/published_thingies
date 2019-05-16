#!/usr/bin/env python3

#---------------------------------------------
# Compare the solutions of the exact Riemann
# solver to the results of the Godunov scheme
# using the exact Riemann solver.
#
# This script assumes that there are outputs
# in the following directories:
#  ./bin/EXACT for the Godunov solutions
#  ../riemann/exact-full/bin for the Riemann
#
# It will only plot the last output, so make
# sure by yourself that the last output is
# correct!
#
# It will create output pngs in 
#  ./png_compare_exact
#---------------------------------------------


from matplotlib import pyplot as plt
import numpy as np
import os, sys


usage="""

Compare the solutions of the exact Riemann
solver to the results of the Godunov scheme
using the exact Riemann solver.

This script assumes that there are outputs
in the following directories:
    ./bin/EXACT for the Godunov solutions
    ../riemann/exact-full/bin for the Riemann

It will create output pngs in
    ./png_compare_exact

usage:
    presentation_compare_exact.py
"""


# globals

srcdir = ''
filebase = ''
noutputs = 0

t = 0

fig = None




#=====================
def get_ics():
#=====================
    """
    Read in what initial conditions have been solved
    for in with the godunov solver
    """

    global srcdir
    srcdir = os.path.join(os.getcwd(), 'bin/EXACT')

    allfiles = os.listdir(srcdir)
    
    ics = []
    for f in allfiles:
        if os.path.isdir(os.path.join(srcdir, f)):
            ics.append(f)

    return ics





#================================
def read_and_plot_godunov():
#================================
    """
    Read in data of godunov solution
    """ 

    global fig

    godunovdir = os.path.join(os.getcwd(), 'bin/EXACT/')
    godunovdir = os.path.join(godunovdir, filebase)
    print("Plotting Godunov solution")

    allfiles = os.listdir(godunovdir)
    outfiles = []
    for f in allfiles:
        if f[-4:] == '.out':
            outfiles.append(f)

    outfiles.sort()

    fname = outfiles[-1] # only take last file
    fname = os.path.join(godunovdir, fname)

    x, rho, u, p = np.loadtxt(fname, dtype=np.float, skiprows=2, unpack = True)


    ax1, ax2, ax3, ax4, ax5, ax6 = fig.axes


    label = 'Godunov method'
    ax4.plot(x, rho, lw=1, label=label)
    ax5.plot(x, u,   lw=1, label=label)
    ax6.plot(x, p,   lw=1, label=label)

    return




#=============================
def read_and_plot_riemann():
#=============================
    """
    Read in data of exact Riemann solver and plot them
    """
    global t, fig
    
    riemanndir = os.path.join(os.getcwd(), '../riemann/exact-full/bin/')
    riemanndir = os.path.join(riemanndir, filebase)

    allfiles = os.listdir(riemanndir)
    outfiles = []
    for f in allfiles:
        if f[-4:] == '.out':
            outfiles.append(f)

    outfiles.sort()

    firstfile = outfiles[0]
    lastfile = outfiles[-1]

    print("Plotting exact Riemann solution")

    # read in first file (which are ICs)
    fname = os.path.join(riemanndir, firstfile)
    xf, rhof, uf, pf = np.loadtxt(fname, dtype=np.float, skiprows=2, unpack = True)

    # read in last file 
    fname = os.path.join(riemanndir, lastfile)
    xl, rhol, ul, pl = np.loadtxt(fname, dtype=np.float, skiprows=2, unpack = True)

    f = open(fname, 'r')
    tline = f.readline()
    f.close()

    tstr, equal, t = tline.partition("=")
    t = float(t)

    ax1, ax2, ax3, ax4, ax5, ax6 = fig.axes


    label = 'Exact solution'
    ax1.plot(xf, rhof, lw=1, zorder=0, label=label)
    ax2.plot(xf, uf,   lw=1, zorder=0, label=label)
    ax3.plot(xf, pf,   lw=1, zorder=0, label=label)
    ax4.plot(xl, rhol, lw=1, zorder=0, label=label)
    ax5.plot(xl, ul,   lw=1, zorder=0, label=label)
    ax6.plot(xl, pl,   lw=1, zorder=0, label=label)


    return




#=======================
def setup_figure():
#=======================
    """
    Sets up the figure for all the plots
    """

    global fig

    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    return





#================================
def tweak_and_save_plot():
#================================
    """
    Make plot prettier, save figure
    """

    ax1, ax2, ax3, ax4, ax5, ax6 = fig.axes

    ax1.set_xlabel('x')
    ax2.set_xlabel('x')
    ax3.set_xlabel('x')
    ax4.set_xlabel('x')
    ax5.set_xlabel('x')
    ax6.set_xlabel('x')
    ax1.set_ylabel('density')
    ax2.set_ylabel('velocity')
    ax3.set_ylabel('pressure')
    ax4.set_ylabel('density')
    ax5.set_ylabel('velocity')
    ax6.set_ylabel('pressure')


    ax1.set_xlim(-1.01, 1.01)
    ax2.set_xlim(-1.01, 1.01)
    ax3.set_xlim(-1.01, 1.01)
    ax4.set_xlim(-1.01, 1.01)
    ax5.set_xlim(-1.01, 1.01)
    ax6.set_xlim(-1.01, 1.01)

    ax2.set_title(filebase+" at t = 0.000")
    ax5.set_title(filebase+" at t = %5.3f" % t)

    ax6.legend()


    #  plt.show()
    outdir = os.path.join(os.getcwd(), 'png_compare_exact')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    fname = os.path.join(outdir, filebase + '.png')

    plt.tight_layout(rect=(0,0,1,0.95))
    print("saving figure", fname)
    plt.savefig(fname, dpi=300)
    plt.close()


    return
    



#==========================
def main():
#==========================
    global filebase
    ics = get_ics()

    for ic in ics:
        print()
        print("Working for", ic)
        filebase = ic
        setup_figure()
        read_and_plot_riemann()
        read_and_plot_godunov()
        tweak_and_save_plot()

    return



#=======================================
if __name__ == '__main__':
#=======================================

    main()


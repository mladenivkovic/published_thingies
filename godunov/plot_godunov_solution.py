#!/usr/bin/env python3


from matplotlib import pyplot as plt
import numpy as np
import os, sys


usage="""
    plot the results made with the riemann code.
    usage:
        plot_riemann_solution.py <output_directory>
"""


# globals

srcdir = ''
filebase = ''
noutputs = 0

t = 0
x = 0
rho = 0
u = 0
p = 0


#=============================
def read_cmdlineargs():
#=============================
    """
    Read in source directory from cmdline
    """

    global srcdir, filebase
    try:
        srcdir = sys.argv[1].strip();
    except IndexError:
        print(usage)
        quit()
    if srcdir.strip() == '-h' or srcdir.strip() == '--help':
        print(usage)
        quit()
    else:
        if not os.path.exists(srcdir):
            print("Directory '", srcdir, "doesn't exist. Can't work like this.")
            print(usage)
            quit()

        else:
            # remove whitespaces and trailing slash
            if srcdir[-1]=='/':
                srcdir = os.path.dirname(srcdir)
            junk, filebase = os.path.split(srcdir)


    return



#========================
def get_noutputs():
#========================
    """
    Get how many output files you have to work with
    """
    global noutputs
    noutputs = 0
    allfiles = os.listdir(srcdir)
    for f in allfiles:
        # look for '.out' in filename
        if f.endswith('.out'):
            if (f.startswith(filebase)):
                noutputs += 1

    return



#========================
def read_data(out):
#========================
    """
    Read in data of output out
    """ 

    global t, x, rho, u, p

    fname = os.path.join(srcdir, filebase + "-" + str(out).zfill(2) + '.out')
    print("Plotting", fname)

    x, rho, u, p = np.loadtxt(fname, dtype=np.float, skiprows=2, unpack = True)

    f = open(fname, 'r')
    tline = f.readline()
    f.close()

    tstr, equal, t = tline.partition("=")
    t = float(t)

    return t, x, rho, u, p



#===========================
def plot_results(out):
#===========================
    """
    Plot and save the results of output out
    """

    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    ax1.plot(x, rho, c='C0', lw=1, zorder=0)
    ax2.plot(x, u, c='C0', lw=1, zorder=0)
    ax3.plot(x, p, c='C0', lw=1, zorder=0)
    ax1.scatter(x, rho, c='C1', s=5, alpha=0.5, zorder=1)
    ax2.scatter(x, u, c='C1', s=5, alpha=0.5, zorder=1)
    ax3.scatter(x, p, c='C1', s=5, alpha=0.5, zorder=1)

    ax1.set_xlabel('x')
    ax2.set_xlabel('x')
    ax3.set_xlabel('x')
    ax1.set_ylabel('density')
    ax2.set_ylabel('velocity')
    ax3.set_ylabel('pressure')

    ax1.set_xlim(-1.01, 1.01)
    ax2.set_xlim(-1.01, 1.01)
    ax3.set_xlim(-1.01, 1.01)

    fig.suptitle(srcdir+" at t = %5.3f" % t)


    fname = os.path.join(srcdir, filebase + "-" + str(out).zfill(2) + '.png')
    plt.tight_layout(rect=(0,0,1,0.95))
    plt.savefig(fname, dpi=72)
    plt.close()





#=======================================
if __name__ == '__main__':
#=======================================

    read_cmdlineargs()
    get_noutputs()

    for i in range(noutputs):
        read_data(i)
        plot_results(i)

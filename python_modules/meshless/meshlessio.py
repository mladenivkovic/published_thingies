#!/usr/bin/env python3

#======================
# Contains IO routines
#======================



#==========================================
def read_file(srcfile, ptype):
#==========================================
    """
    Read swift output hdf5 file.
    """

    import h5py

    f = h5py.File(srcfile)

    x = f[ptype]['Coordinates'][:,0]
    y = f[ptype]['Coordinates'][:,1]
    h = f[ptype]['SmoothingLength'][:]
    rho = f[ptype]['Density'][:]
    m = f[ptype]['Masses'][:]
    ids = f[ptype]['ParticleIDs'][:]

    npart = x.shape[0]

    f.close()

    return x, y, h, rho, m, ids, npart





#===========================
def get_sample_size():
#===========================
    """
    Count how many files we're dealing with
    Assumes snapshots start with "snapshot-" string and contain
    two numbers: snpshot-XXX-YYY_ZZZZ.hdf5, where both XXX and YYY
    are integers, have the same minimal, maximal value and same
    difference between two consecutive numbers.

    this is intended for numbered output.
    Returns:
        nx : number of files (in one direction)
        filenummax: highest XXX
        fileskip: integer difference between two XXX or YYY
    """

    import os
    import numpy as np

    filelist = os.listdir()
    snaplist = [ ]
    for f in filelist:
        if f.startswith('snapshot-'):
            snaplist.append(f)

    snaplist.sort()
    first = snaplist[0]
    s, dash, rest = first.partition("-")
    num, dash, junk = rest.partition("-")
    lowest = int(num)
    
    finalsnap = snaplist[-1]
    s, dash, rest = finalsnap.partition("-")
    num, dash, junk = rest.partition("-")

    highest = int(num)

    steps = int(np.sqrt(len(snaplist)))

    nx = steps
    filenummax = highest
    fileskip = int((highest - lowest)/(steps - 1))

    return nx, filenummax, fileskip 




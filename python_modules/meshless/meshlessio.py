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
    m = f[ptype]['Masses'][:]
    ids = f[ptype]['ParticleIDs'][:]

    try:
        # old SWIFT header versions
        h = f[ptype]['SmoothingLength'][:]
        rho = f[ptype]['Densities'][:]
    except KeyError:
        # new SWIFT header versions
        h = f[ptype]['SmoothingLengths'][:]
        rho = f[ptype]['Densities'][:]
    npart = x.shape[0]

    f.close()

    return x, y, h, rho, m, ids, npart





#====================================
def get_sample_size(prefix=None):
#====================================
    """
    Count how many files we're dealing with
    Assumes snapshots start with "snapshot-" string and contain
    two numbers: snashot-XXX-YYY_ZZZZ.hdf5, where both XXX and YYY
    are integers, have the same minimal, maximal value and same
    difference between two consecutive numbers.

    if prefix is given, it will prepend it to snapshots.

    this is intended for numbered output.
    Returns:
        nx : number of files (in one direction)
        filenummax: highest XXX
        fileskip: integer difference between two XXX or YYY
    """

    import os
    import numpy as np

    if prefix is not None:
        filelist = os.listdir(prefix)
    else:
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




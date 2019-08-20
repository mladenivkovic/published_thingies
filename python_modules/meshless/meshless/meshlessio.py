#!/usr/bin/env python3

#======================
# Contains IO routines
#======================



#=========================================================
def read_file(srcfile, ptype='PartType0', sort=False):
#=========================================================
    """
    Read swift output hdf5 file.
    srcfile:    string of file to be read in
    ptype:      which particle type to work with
    sort:       whether to sort read in arrays by particle ID

    returns:
    x, y, h, rho, m, ids: numpy arrays of x, y position, 
        smoothing length, density, mass, particle ID
    npart: Number of particles
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
        rho = f[ptype]['Density'][:]
    except KeyError:
        # new SWIFT header versions
        h = f[ptype]['SmoothingLengths'][:]
        rho = f[ptype]['Densities'][:]
    npart = x.shape[0]

    f.close()

    if sort:
        from numpy import argsort
        inds = argsort(ids)
        x = x[inds]
        y = y[inds]
        h = h[inds]
        rho = rho[inds]
        m = m[inds]
        ids = ids[inds]

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




#!/usr/bin/env python3

#===============================
# Particle related methods
#===============================

import numpy as np




#===============================
def find_index(x, y, h, pcoord):
#===============================
    """
    Find the index in the read-in arrays where
    the particle with coordinates of your choice is
    """

    tolerance = 1e-6

    for i in range(x.shape[0]):
        if abs(x[i]-pcoord[0]) < tolerance and abs(y[i] - pcoord[1]) < tolerance:
            pind = i
            break

    print("got pind", pind, 'with x=', x[pind], 'y=', y[pind], 'h=', h[pind])
 

    return pind





#============================================
def find_neighbours(ind, x, y, h, fact=2):
#============================================
    """
    Find indices of all neighbours of a particle with index ind
    within fact*h (where kernel != 0)
    x, y, h: arrays of positions/h of all particles
    returns list of neighbour indices
    """

    x0 = x[ind]
    y0 = y[ind]
    fhsq = h[ind]*h[ind]*fact*fact
    neigh = [None for i in x]

    j = 0
    for i in range(x.shape[0]):
        if i==ind:
            continue

        dist = (x[i]-x0)**2 + (y[i]-y0)**2
        if dist < fhsq:
            neigh[j] = i
            j+=1

    return neigh[:j]




#===========================================================
def find_neighbours_arbitrary_x(x0, y0, x, y, h, fact=2):
#===========================================================
    """
    Find indices of all neighbours around position x0, y0
    within fact*h (where kernel != 0)
    x, y, h: arrays of positions/h of all particles
    returns list of neighbour indices
    """


    neigh = [None for i in x]
    j = 0


    if isinstance(h, np.ndarray):
        fsq = fact*fact

        for i in range(x.shape[0]):
            dist = (x[i]-x0)**2 + (y[i]-y0)**2
            fhsq = h[i]*h[i]*fsq
            if dist < fhsq:
                neigh[j]=i
                j+=1

    else:
        fhsq = fact*fact*h*h
        for i in range(x.shape[0]):
            dist = (x[i]-x0)**2 + (y[i]-y0)**2
            if dist < fhsq:
                neigh[j] = i
                j+=1


    return neigh[:j]






#===================
def V(ind, m, rho):
#===================
    """
    Volume estimate for particle with index ind
    """

    return m[ind]/rho[ind]



#======================================
def find_central_particle(L, ids):
#======================================
    """
    Find the index of the central particle at (0.5, 0.5)
    """

    i = L//2-1
    cid = i*L + i + 1
    cind = np.asscalar(np.where(ids==cid)[0])

    return cind


#======================================
def find_added_particle(ids):
#======================================
    """
    Find the index of the added particle (has highest ID)
    """

    pid = ids.shape[0]
    pind = np.asscalar(np.where(ids==pid)[0])

    return pind

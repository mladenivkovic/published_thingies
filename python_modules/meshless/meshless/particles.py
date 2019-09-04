#!/usr/bin/env python3

#===============================
# Particle related methods
#===============================

import numpy as np




#===============================================
def find_index(x, y, pcoord, tolerance=1e-3):
#===============================================
    """
    Find the index in the read-in arrays where
    the particle with coordinates of your choice is

    x, y:       arrays of x, y positions of all particles
    pcoors:     array/list of x,y position of particle to look for
    tolerance:  how much of a tolerance to use to identify particle
                useful when you have random perturbations of a
                uniform grid with unknown exact positions
    """

    for i in range(x.shape[0]):
        if abs(x[i]-pcoord[0]) < tolerance and abs(y[i] - pcoord[1]) < tolerance:
            pind = i
            break

    return pind




#===============================================
def find_index_by_id( ids, id_to_look_for ):
#===============================================
    """
    Find the index in the read-in arrays where
    the particle with id_to_look_for is

        ids:    numpy array of particle IDs
        id_to_look_for : which ID to find

    returns:
        pind:  index of particle with id_to_look_for

    """

    pind = np.asscalar(np.where(ids==id_to_look_for)[0])

    return pind





#================================================================
def find_neighbours(ind, x, y, h, fact=1, L=1, periodic=True):
#================================================================
    """
    Find indices of all neighbours of a particle with index ind
    within fact*h (where kernel != 0)
    x, y, h:    arrays of positions/h of all particles
    fact:       kernel support radius factor: W = 0 for r > fact*h
    L:          boxsize
    periodic:   Whether you assume periodic boundary conditions

    returns list of neighbour indices in x,y,h array
    """
    

    # None for Gaussian
    if fact is not None:

        x0 = x[ind]
        y0 = y[ind]
        fhsq = h[ind]*h[ind]*fact*fact
        neigh = [None for i in x]

        j = 0
        for i in range(x.shape[0]):
            if i==ind:
                continue

            dx, dy = get_dx(x0, x[i], y0, y[i], L=L, periodic=periodic)

            dist = dx**2 + dy**2 

            if dist < fhsq:
                neigh[j] = i
                j+=1

        return neigh[:j]

    else:
        neigh = [i for i in range(x.shape[0])] 
        neigh.remove(ind)
        return neigh






#=================================================================================
def find_neighbours_arbitrary_x(x0, y0, x, y, h, fact=1, L=1, periodic=True):
#=================================================================================
    """
    Find indices of all neighbours around position x0, y0
    within fact*h (where kernel != 0)
    x, y, h:    arrays of positions/h of all particles
    fact:       kernel support radius factor: W = 0 for r > fact*h
    L:          boxsize
    periodic:   Whether you assume periodic boundary conditions

    returns list of neighbour indices
    """


    # None for Gaussian
    if fact is not None:
        neigh = [None for i in x]
        j = 0


        if isinstance(h, np.ndarray):
            fsq = fact*fact

            for i in range(x.shape[0]):

                dx, dy = get_dx(x0, x[i], y0, y[i], L=L, periodic=periodic)

                dist = dx**2 + dy**2

                fhsq = h[i]*h[i]*fsq
                if dist < fhsq:
                    neigh[j]=i
                    j+=1

        else:
            fhsq = fact*fact*h*h
            for i in range(x.shape[0]):

                dx, dy = get_dx(x0, x[i], y0, y[i], L=L, periodic=periodic)

                dist = dx**2 + dy**2

                if dist < fhsq:
                    neigh[j] = i
                    j+=1


        return neigh[:j]

    else:
        neigh = [i for i in range(x.shape[0])] 
        return neigh








#===================
def V(ind, m, rho):
#===================
    """
    Volume estimate for particle with index ind
    """
    V = m[ind]/rho[ind]
    if V > 1:
        print("Got particle volume V=", v, ". Did you put the arguments in the correct places?")
    return V







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




#=====================================================
def get_dx(x1, x2, y1, y2, L=1, periodic=True):
#=====================================================
    """
    Compute difference of vectors [x1 - x2, y1 - y2] while
    checking for periodicity if necessary
    L:          boxsize
    periodic:   whether to assume periodic boundaries
    """

    dx = x1 - x2
    dy = y1 - y2

    if periodic:

        Lhalf = 0.5*L

        if dx > Lhalf:
            dx -= L
        elif dx < -Lhalf:
            dx += L

        if dy > Lhalf:
            dy -= L
        elif dy < -Lhalf:
            dy += L


    return dx, dy

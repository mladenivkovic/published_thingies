#!/usr/bin/env python3

###########################################################################################
#  package:   astro-meshless-surfaces
#  file:      particles.py
#  brief:     particle related functions
#  copyright: GPLv3
#             Copyright (C) 2019 EPFL (Ecole Polytechnique Federale de Lausanne)
#             LASTRO - Laboratory of Astrophysics of EPFL
#  author:    Mladen Ivkovic <mladen.ivkovic@epfl.ch>
#
# This file is part of astro-meshless-surfaces.
###########################################################################################

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
def find_neighbours(ind, x, y, h, fact=1.0, L=1.0, periodic=True):
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
            if i == ind:
                continue

            dx, dy = get_dx(x0, x[i], y0, y[i], L=L, periodic=periodic)

            dist = dx**2 + dy**2 

            if dist < fhsq:
                neigh[j] = i
                j += 1

        return neigh[:j]

    else:
        neigh = [i for i in range(x.shape[0])] 
        neigh.remove(ind)
        return neigh









#=================================================================================
def find_neighbours_arbitrary_x(x0, y0, x, y, h, fact=1.0, L=1.0, periodic=True):
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
        print("Got particle volume V=", V, ". Did you put the arguments in the correct places?")
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
def get_dx(x1, x2, y1, y2, L=1.0, periodic=True):
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








#========================================================================
def get_neighbour_data_for_all(x, y, h, fact=1.0, L=1.0, periodic=True):
#========================================================================
    """
    Gets all the neighbour data for all particles ready.
    x, y, h:    arrays of positions/h of all particles
    fact:       kernel support radius factor: W = 0 for r > fact*h
    L:          boxsize
    periodic:   Whether you assume periodic boundary conditions

    returns neighbour_data object:
        self.neighbours :   List of lists of every neighbour of every particle
        self.maxneigh :     Highest number of neighbours
        self.nneigh:        integer array of number of neighbours for every particle
        self.iinds:         iinds[i, j] = which index does particle i have in the neighbour
                            list of particle j, where j is the j-th neighbour of i
                            Due to different smoothing lengths, particle j can be the
                            neighbour of i, but i not the neighbour of j. 
                            In that case, the particles will be assigned indices j > nneigh[i]

    """


    npart = x.shape[0]

    # find and store all neighbours;
    neighbours = [[] for i in x]
    for i in range(npart):
        neighbours[i] = find_neighbours(i, x, y, h, fact=fact, L=L, periodic=periodic)


    # get neighbour counts array
    nneigh = np.zeros((npart), dtype=np.int)
    for i in range(npart):
        nneigh[i] = len(neighbours[i])

    
    # max number of neighbours; needed for array allocation
    maxneigh = nneigh.max()


    # store the index of particle i when required as the neighbour of particle j in arrays[npart, maxneigh]
    # i.e. find index 0 <= i < maxneigh for ever j
    iinds = np.zeros((npart, 2*maxneigh), dtype=np.int) 
    current_count = nneigh[:]

    for i in range(npart):
        for jc,j in enumerate(neighbours[i]):

            try:
                iinds[i, jc] = (neighbours[j]).index(i)
            except ValueError:
                # it is possible that j is a neighbour for i, but i is not a neighbour
                # for j depending on their respective smoothing lengths
                dx, dy = get_dx(x[i], x[j], y[i], y[j], L=L, periodic=periodic)
                r = np.sqrt(dx**2 + dy**2)
                if r/h[j] < 1:
                    print("something went wrong when computing gradients.")
                    print("i=", i, "j=", j, "r=", r, "H=", h[j], "r/H=", r/h[j])
                    print("neighbours i:", neighbours[i])
                    print("neighbours j:", neighbours[j])
                    print("couldn't find i as neighbour of j")
                    print("exiting")
                    quit()
                else:
                    # append after nneigh[j]
                    iinds[i, jc] = current_count[j]
                    current_count[j] += 1


    class neighbour_data:
        def __init__(self, neighbours=None, maxneigh=None, nneigh=None, iinds=None):
            self.neighbours = neighbours
            self.maxneigh = maxneigh
            self.nneigh = nneigh
            self.iinds = iinds

    nd = neighbour_data(neighbours=neighbours,
                        maxneigh=maxneigh,
                        nneigh=nneigh,
                        iinds=iinds) 


    return nd

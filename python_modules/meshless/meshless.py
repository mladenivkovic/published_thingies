#!/usr/bin/env python3

#===========================================
#
# A module containing common routines for 
# the meshless effective area visualisiation
# with 2d datasets
#
#===========================================


from meshlessio import read_file, get_sample_size
from kernels import W, kernels
from particles import find_index, find_neighbours, V, find_central_particle, find_added_particle


import numpy as np


#=================================================================
def Aij_Hopkins(pind, x, y, h, m, rho, kernel='cubic_spline', fact=2):
#=================================================================
    """
    Compute A_ij as defined by Hopkins 2015
    i, j: particle indices for x, y, h, m arrays
    x, y, h, m: full data arrays as read in from hdf5 file
    kernel: which kernel to use
    fact:   factor for h for limit of neighbour search; neighbours are closer than fact*h
    """

    debug = False

    nbors = find_neighbours(pind, x, y, h, fact)

    xj = x[nbors]
    yj = y[nbors]
    hj = h[nbors]

    #-------------------------------------------------------
    # Part 1: For particle at x_i (Our chosen particle)
    #-------------------------------------------------------

    # save and store psi(dx=0, dy=0)
    psi_null = psi(0,0,0,0,1,kernel)

    # compute psi_j(x_i)
    psi_j = compute_psi(x[pind], y[pind], xj, yj, h[pind], kernel)

    if debug:
        print("unnormalized psi_j:", psi_j)

    # normalize psi_j
    omega_xi =  (np.sum(psi_j) + psi_null)
    psi_j /= omega_xi
    psi_j = np.float64(psi_j)

    # compute B_i
    B_i = get_matrix(x[pind], y[pind], xj, yj, psi_j)

    # compute psi_tilde_j(x_i)
    psi_tilde_j = np.empty((len(nbors), 2), dtype=np.float)
    for i, n in enumerate(nbors):
        dx = np.array([xj[i]-x[pind], yj[i]-y[pind]])
        psi_tilde_j[i] = np.dot(B_i, dx) * psi_j[i]

    if debug:
        print('psi_tilde_j', psi_tilde_j)

    #---------------------------------------------------------------------------
    # Part 2: values of psi/psi_tilde of particle i at neighbour positions x_j
    #---------------------------------------------------------------------------

    psi_i = np.zeros(len(nbors), dtype=np.float128)         # psi_i(xj)
    psi_tilde_i = np.empty((len(nbors), 2), dtype=np.float)  # psi_tilde_i(x_j)

    for i, n in enumerate(nbors):
        # first compute all psi(xj) from neighbour's neighbours to get weight omega
        nneigh = find_neighbours(n, x, y, h, fact)
        xk = x[nneigh]
        yk = y[nneigh]
        for j, nn in enumerate(nneigh):
            psi_k = compute_psi(x[n], y[n], xk, yk, h[n], kernel)
            if nn == pind: # store psi_i, which is the psi for the particle whe chose at position xj; psi_i(xj)
                psi_i[i] = psi_k[j]
    
        omega_xj = (np.sum(psi_k) + psi_null)

        psi_i[i]/= omega_xj
        psi_k /= omega_xj
        psi_k = np.float64(psi_k)


        # now compute B_j^{\alpha \beta}
        B_j = get_matrix(x[n], y[n], xk, yk, psi_k)

        # get psi_i_tilde(x = x_j)
        dx = np.array([x[pind]-x[n], y[pind]-y[n]])
        psi_tilde_i[i] = np.dot(B_j, dx) * np.float64(psi_i[i])


    if debug:
        print("unnormalized psi_i:", psi_i)
        print("psi_tilde_i", psi_tilde_i)


    #-------------------------------
    # Part 3: Compute A_ij    #-------------------------------

    A_ij = np.empty((len(nbors),2), dtype = np.float)

    for i,n in enumerate(nbors):
        A_ij[i] = V(pind, m, rho)*psi_tilde_j[i] - V(n, m, rho)*psi_tilde_i[i]

    if debug:
        print("neighbour volumes")
        for n in nbors:
            print(V(n, m, rho), end=' ')
        print()
        print("particle volume estimate:", 1.0/(np.sum(psi_j)+psi_null))


    return A_ij








#====================================================
def x_ij(pind, x, y, h, nbors=None, which=None):
#====================================================
    """
    compute x_ij for all neighbours of particle with index pind
    if which=integer is given, instead compute only for specific particle
    where which= that particle's ID
    """


    if which is not None:
        hfact = h[pind]/(h[pind]+h[which])
        x_ij = np.array([x[pind]-hfact*(x[pind]-x[which]), y[pind]-hfact*(y[pind]-y[which])])
        return x_ij

    elif nbors is not None:
        x_ij = np.empty((len(nbors),2), dtype = np.float)
        for i,n in enumerate(nbors):
            hfact = h[pind]/(h[pind]+h[n])
            x_ij[i] = np.array([x[pind]-hfact*(x[pind]-x[n]), y[pind]-hfact*(y[pind]-y[n])])

        return x_ij

    else:
        print("Gotta give me a list of neighbours or a single particle info for x_ij")
        quit()

    return







#==================================================
def A_ij_Ivanova(x, y, h, rho, m, nbors):
#==================================================
    """
    Compute the effective area A_ij(x)

    CURRENTLY COMPLETELY WRONG LOL
    """

    xj = x[nbors]
    yj = y[nbors]
    hj = h[nbors]

    #-------------------------------------------------------
    # Part 1: For particle at x_i (Our chosen particle)
    #-------------------------------------------------------

    # compute psi_j(x_i)
    psi_j = compute_psi(x[pind], y[pind], xj, yj, h[pind], kernel)

    # normalize psi_j. Don't forget to add the self-contributing value!
    omega_xi =  (np.sum(psi_j) + psi(x[pind], y[pind], x[pind], y[pind], h[pind]))
    psi_j /= omega_xi
    psi_j = np.float64(psi_j)

    # compute B_i
    B_i = get_matrix(x[pind], y[pind], xj, yj, psi_j)

    # get index of central particle in nbors list
    cind_n = nbors.index(cind)

    # compute grad_psi_j(x_i)
    grad_psi_j = np.empty((1, 2), dtype=np.float)
    dx = np.array([x[cind]-x[pind], y[cind]-y[pind]])
    grad_psi_j = np.dot(B_i, dx) * psi_j[cind_n]



    #---------------------------------------------------------------------------
    # Part 2: values of psi/grad_psi of particle i at neighbour positions x_j
    #---------------------------------------------------------------------------

    psi_i = 0.0                                    # psi_i(xj)
    grad_psi_i = np.empty((1, 2), dtype=np.float)  # grad_psi_i(x_j)

    # first compute all psi(xj) from central's neighbours to get weight omega
    nneigh = find_neighbours(cind, x, y, h)
    xk = x[nneigh]
    yk = y[nneigh]
    for j, nn in enumerate(nneigh):
        psi_k = compute_psi(x[cind], y[cind], xk, yk, h[cind])
        if nn == pind: # store psi_i, which is the psi for the particle whe chose at position xj; psi_i(xj)
            psi_i = psi_k[j]

    omega_xj = (np.sum(psi_k) + psi(x[cind], y[cind], x[cind], y[cind], h[cind]))

    psi_i/= omega_xj
    psi_i = np.float64(psi_i)

    # now compute B_j^{\alpha \beta}
    B_j = get_matrix(x[cind], y[cind], xk, yk, h[cind])

    # get gradient
    dx = np.array([x[pind]-x[cind], y[pind]-y[cind]])
    grad_psi_i = np.dot(B_j, dx) * psi_i



    #-------------------------------
    # Part 3: Compute A_ij, x_ij
    #-------------------------------

    V = m/rho

    A_ij = None

    A_ij = V[pind]*grad_psi_j - V[cind]*grad_psi_i

    if A_ij is None:
        print("PROBLEM: A_IJ IS NONE")
        raise ValueError
    else:
        return A_ij












#=============================================================
def compute_psi(xi, yi, xj, yj, hi, kernel='cubic_spline'):
#=============================================================
    """
    Compute all psi_j(x_i)
    xi, yi: floats; position for which to compute psi's
    xj, yj: arrays of neighbour's positions
    hi:     float; smoothing length at position xi, yi
    """

    psi_j = np.zeros(xj.shape[0], dtype=np.float128)

    for i in range(xj.shape[0]):
        psi_j[i] = psi(xi, yi, xj[i], yj[i], hi, kernel)

    return psi_j







#====================================================
def psi(x, y, xi, yi, hi, kernel='cubic_spline'):
#====================================================
    """
    UNNORMALIZED Volume fraction at position x of some particle
    with coordinates xi, yi
    ind: neighbour index in x/y/h array
    """
    q = np.float128(np.sqrt((x - xi)**2 + (y - yi)**2)/hi)

    return W(q, hi, kernel)








#=============================================
def get_matrix(xi, yi, xj, yj, psi_j):
#=============================================
    """
    Get B_i ^{alpha beta}
    xi, yi: floats; Evaluate B at this position
    xj, yj: arrays; Neighbouring points
    psi_j:  array;  volume fraction of neighbours at position x_i; psi_j(x_i)
    """

    E00 = np.sum((xj-xi)**2 * psi_j)
    E01 = np.sum((xj-xi)*(yj-yi) * psi_j)
    E11 = np.sum((yj-yi)**2 * psi_j)
          
    E = np.matrix([[E00, E01], [E01, E11]])

    B = E.getI()
    return B



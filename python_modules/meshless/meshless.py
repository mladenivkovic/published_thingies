#!/usr/bin/env python3

#===========================================
#
# A module containing common routines for 
# the meshless effective area visualisiation
# with 2d datasets
#
#===========================================


from meshlessio import *
from kernels import *
from particles import *


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







#===========================================================================================
def Integrand_Aij_Ivanova(iind, jind, xx, yy, hh, x, y, h, m, rho, kernel='cubic_spline', fact=2):
#===========================================================================================
    """
    Compute the effective area integrand for the particles iind jind at
    the positions xx, yy

    (Note that this should be integrated to get the proper effective surface)

    A_ij    = psi_j(x) \nabla psi_i(x) - psi_i (x) \nabla psi_j(x)
            = sum_k [ psi_j(x_k) psi_i(x) - psi_i(x_k) psi_j(x) ] * psi_tilde_k(x)
            = psi_i(x) * sum_k psi_j(x_k) * psi_tilde_k(x) - psi_j(x) * sum_k psi_i(x_k) * psi_tilde_k(x)
    
    The last line is what is actually computed here, with the expression for the gradient
    inserted.
    """

    nbors = find_neighbours_arbitrary_x(xx, yy, x, y, h, fact)

    xk = x[nbors]
    yk = y[nbors]
    hk = h[nbors]



    #----------------------
    # compute psi_i/j(x)
    #----------------------

    # compute all psi(x)
    psi_x = compute_psi(xx, yy, xk, yk, hh, kernel)

    # normalize psis
    omega = np.sum(psi_x)
    psi_x /= omega
    psi_x = np.float64(psi_x)



    # find where psi_i and psi_j are in that array
    inb = nbors.index(iind)
    jnb = nbors.index(jind)

    psi_i_of_x = psi_x[inb] # psi_i(xx, yy)
    psi_j_of_x = psi_x[jnb] # psi_j(xx, yy)




    #------------------------------------------------
    # Compute psi_i/j(x_k) at neighbouring positions
    #------------------------------------------------

    psi_i_xk = [None for n in nbors]
    psi_j_xk = [None for n in nbors]

    psi_null = psi(0, 0, 0, 0, 1, kernel)

    omegas = [0, 0]
    
    for i, n in enumerate([iind, jind]):
        # first compute all psi(xl) from neighbour's neighbours to get weights omega
        nneigh = find_neighbours(n, x, y, h, fact)

        xl = x[nneigh]
        yl = y[nneigh]

        for j, nn in enumerate(nneigh):
            psi_l = compute_psi(x[n], y[n], xl, yl, h[n], kernel)

        omegas[i] = np.sum(psi_l) + psi_null


    # now compute psi_i/j(x_k)
    for i, n in enumerate(nbors):
        psi_i_xk[i] = psi(xk[i], yk[i], x[iind], y[iind], h[iind], kernel) / omegas[0]
        psi_j_xk[i] = psi(xk[i], yk[i], x[jind], y[jind], h[jind], kernel) / omegas[1]





    #---------------------------------------------
    # Compute psi_tilde_k(x)
    #---------------------------------------------

    # compute matrix B
    B = get_matrix(xx, yy, xk, yk, psi_x)

    # compute psi_tilde_k(xx)
    psi_tilde_k = np.empty((2, xk.shape[0]))
    for i in range(xk.shape[0]):
        dx = np.array([xk[i]-xx, yk[i]-yy])
        psi_tilde_k[:, i] = np.multiply(np.dot(B, dx), psi_x[i])




    #----------------------------------
    # Compute A_ij
    #----------------------------------

    sum_i = np.sum(np.multiply(psi_tilde_k, psi_i_xk ), axis=1)
    sum_j = np.sum(np.multiply(psi_tilde_k, psi_j_xk ), axis=1)


    A_ij = psi_i_of_x * sum_j - psi_j_of_x * sum_i 
    
    return A_ij















#=============================================================
def compute_psi(xi, yi, xj, yj, h, kernel='cubic_spline'):
#=============================================================
    """
    Compute all psi_j(x_i)
    xi, yi: floats; position for which to compute psi's
    xj, yj: arrays of neighbour's positions
    h:      float; smoothing length at position xi, yi
            or array of h for xj, yj [used to compute h(x)]

    return numpy array of psi_j(x) 
    """

    psi_j = np.zeros(xj.shape[0], dtype=np.float128)

    if isinstance(h, np.ndarray):
        for i in range(xj.shape[0]):
            psi_j[i] = psi(xi, yi, xj[i], yj[i], h[i], kernel)

    else:
        for i in range(xj.shape[0]):
            psi_j[i] = psi(xi, yi, xj[i], yj[i], h, kernel)

    return psi_j







#====================================================
def psi(x, y, xi, yi, hi, kernel='cubic_spline'):
#====================================================
    """
    UNNORMALIZED Volume fraction at position x of some particle
    with coordinates xi, yi, smoothing length hi
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





#========================================================
def h_of_x(xx, yy, x, y, h, m, rho, kernel='cubic_spline', fact=2):
#========================================================
    """
    Compute h(x) at position (xx, yy), where there is 
    not necessariliy a particle
    x, y, h : full particle arrays
    """

    nbors = find_neighbours_arbitrary_x(xx, yy, x, y, h, fact)

    xj = x[nbors]
    yj = y[nbors]
    hj = h[nbors]

    psi_j = compute_psi(xx, yy, xj, yj, hj, kernel)
    psi_j /= np.sum(psi_j)

    hh = np.sum(hj*psi_j)

    return hh




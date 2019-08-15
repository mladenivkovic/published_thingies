#!/usr/bin/env python3

#===========================================
#
# A module containing common routines for 
# the meshless effective area visualisiation
# with 2d datasets
#
#===========================================

try:
    from .meshlessio import *
    from .kernels import *
    from .particles import *
except ImportError:
    # in case you're not using it as a package, but directly in the pythonpath
    from meshlessio import *
    from kernels import *
    from particles import *


import numpy as np


#===========================================================================================
def Aij_Hopkins(pind, x, y, h, m, rho, kernel='cubic_spline', fact=1, L=1, periodic=True):
#===========================================================================================
    """
    Compute A_ij as defined by Hopkins 2015
    pind:           particle index for which to work with. (The i in A_ij) 
    x, y, m, rho:   full data arrays as read in from hdf5 file
    h:              kernel support radius array
    kernel:         which kernel to use
    fact:           factor for h for limit of neighbour search; neighbours are closer than fact*h
    L:              boxsize
    periodic:       whether to assume periodic boundaries

    returns:
        A_ij: array of A_ij, containing x and y component for every neighbour j of particle i
    """

    debug = False

    nbors = find_neighbours(pind, x, y, h, fact=fact, L=L, periodic=periodic)

    xj = x[nbors]
    yj = y[nbors]
    hj = h[nbors]

    #-------------------------------------------------------
    # Part 1: For particle at x_i (Our chosen particle)
    #-------------------------------------------------------

    # compute psi_j(x_i)
    psi_j = compute_psi(x[pind], y[pind], xj, yj, h[pind], kernel, fact=fact, L=L, periodic=periodic)


    # normalize psi_j
    omega_xi =  (np.sum(psi_j) + psi(0,0,0,0,h[pind],kernel))
    psi_j /= omega_xi
    psi_j = np.atleast_1d(np.float64(psi_j))


    # compute B_i
    B_i = get_matrix(x[pind], y[pind], xj, yj, psi_j, L=L, periodic=periodic)

    # compute psi_tilde_j(x_i)
    psi_tilde_j = np.empty((len(nbors), 2), dtype=np.float)
    for i, n in enumerate(nbors):
        dx = np.array([xj[i]-x[pind], yj[i]-y[pind]])
        psi_tilde_j[i] = np.dot(B_i, dx) * psi_j[i]

    #---------------------------------------------------------------------------
    # Part 2: values of psi/psi_tilde of particle i at neighbour positions x_j
    #---------------------------------------------------------------------------

    psi_i = np.zeros(len(nbors), dtype=np.float128)          # psi_i(xj)
    psi_tilde_i = np.empty((len(nbors), 2), dtype=np.float)  # psi_tilde_i(x_j)

    for i, n in enumerate(nbors):
        # first compute all psi(xj) from neighbour's neighbours to get weight omega
        nneigh = find_neighbours(n, x, y, h, fact=fact, L=L, periodic=periodic)
        xk = x[nneigh]
        yk = y[nneigh]
        for j, nn in enumerate(nneigh):
            psi_k = compute_psi(x[n], y[n], xk, yk, h[n], kernel, fact=fact, L=L, periodic=periodic)
            if nn == pind: # store psi_i, which is the psi for the particle whe chose at position xj; psi_i(xj)
                psi_i[i] = psi_k[j]
    
        omega_xj = (np.sum(psi_k) + psi(0,0,0,0,h[nn],kernel))

        psi_i[i]/= omega_xj
        psi_k /= omega_xj
        psi_k = np.float64(psi_k)

        # now compute B_j^{\alpha \beta}
        B_j = get_matrix(x[n], y[n], xk, yk, psi_k, L=L, periodic=periodic)

        # get psi_i_tilde(x = x_j)
        dx = np.array([x[pind]-x[n], y[pind]-y[n]])
        psi_tilde_i[i] = np.dot(B_j, dx) * np.float64(psi_i[i])


    #-------------------------------
    # Part 3: Compute A_ij    
    #-------------------------------

    A_ij = np.empty((len(nbors),2), dtype = np.float)

    

    for i,n in enumerate(nbors):
        A_ij[i] = V(pind, m, rho)*psi_tilde_j[i] - V(n, m, rho)*psi_tilde_i[i]

    return A_ij







#===================================================================================================
def Aij_Hopkins_v2(pind, x, y, h, m, rho, kernel='cubic_spline', fact=1, L=1, periodic=True):
#===================================================================================================
    """
    Compute A_ij as defined by Hopkins 2015, second version
    pind:           particle index for which to work for (The i in A_ij)
    x, y, m, rho:   full data arrays as read in from hdf5 file
    h:              kernel support radius array
    kernel:         which kernel to use
    fact:           factor for h for limit of neighbour search; neighbours are closer than fact*h
    L:          boxsize
    periodic:   whether to assume periodic boundaries

    returns:
        A_ij: array of A_ij, containing x and y component for every neighbour j of particle i
    """


    npart = x.shape[0]


    # compute all psi_k(x_l) for all l, k
    # first index: index k of psi: psi_k(x)
    # second index: index of x_l: psi(x_l)
    psi_k_at_l = np.zeros((npart, npart), dtype=np.float128)

    for k in range(npart):
        for l in range(npart):
            # kernels are symmetric in x_i, x_j, but h can vary!!!!
            psi_k_at_l[k,l] = psi(x[l], y[l], x[k], y[k], h[l], kernel=kernel, fact=fact, L=L, periodic=periodic)

    neighbours = [[] for i in x]
    omega = np.zeros(npart, dtype=np.float128)

    for l in range(npart):

        # find and store all neighbours;
        neighbours[l] = find_neighbours(l, x, y, h, fact=fact, L=L, periodic=periodic)

        # compute normalisation omega for all particles
        # needs psi_k_at_l to be computed already
        omega[l] = np.sum(psi_k_at_l[:, l])
        # omega_k = sum_l W(x_k - x_l, h_k) = sum_l psi_l(x_k) as it is currently stored in memory



    # normalize psi's and convert to float64 for linalg module
    for k in range(npart):
        psi_k_at_l[:, k] /= omega[k]
    psi_k_at_l = np.float64(psi_k_at_l)


    # compute all matrices B_k
    B_k = np.zeros((npart), dtype=np.matrix)
    for k in range(npart):
        nbors = neighbours[k]
        # nbors now contains all neighbours l
        B_k[k] = get_matrix(x[k], y[k], x[nbors], y[nbors], psi_k_at_l[nbors, k], L=L, periodic=periodic)



    # compute all psi_tilde_k at every l
    psi_tilde_k_at_l = np.zeros((npart, npart, 2))
    for k in range(npart):
        for l in range(npart):

            dx = np.array([x[k]-x[l], y[k]-y[l]])
            psi_tilde_k_at_l[k,l] = np.dot(B_k[l], dx) * psi_k_at_l[k,l]

 

    # now compute A_ij for all neighbours j of i
    nbors = neighbours[pind]

    A_ij = np.zeros((len(nbors), 2), dtype=np.float64)

    for i,j in enumerate(nbors): 

        A_ij[i] = V(pind, m, rho) * psi_tilde_k_at_l[j, pind] - V(j, m, rho) * psi_tilde_k_at_l[pind, j]

 
    return A_ij









#=====================================================================================================================
def Aij_Ivanova_approximate_gradients(pind, x, y, h, m, rho, kernel='cubic_spline', fact=1, L=1, periodic=True):
#=====================================================================================================================
    """
    Compute A_ij as defined by Ivanova 2013
    pind:           particle index for which to work for (The i in A_ij)
    x, y, m, rho:   full data arrays as read in from hdf5 file
    h:              kernel support radius array
    kernel:         which kernel to use
    fact:           factor for h for limit of neighbour search; neighbours are closer than fact*h
    L:          boxsize
    periodic:   whether to assume periodic boundaries

    returns:
        A_ij: array of A_ij, containing x and y component for every neighbour j of particle i
    """


    npart = x.shape[0]


    # compute all psi_k(x_l) for all l, k
    # first index: index k of psi: psi_k(x)
    # second index: index of x_l: psi(x_l)
    psi_k_at_l = np.zeros((npart, npart), dtype=np.float128)

    for k in range(npart):
        for l in range(npart):
            # kernels are symmetric in x_i, x_j, but h can vary!!!!
            psi_k_at_l[k,l] = psi(x[l], y[l], x[k], y[k], h[l], kernel=kernel, fact=fact, L=L, periodic=periodic)


    neighbours = [[] for i in x]
    omega = np.zeros(npart, dtype=np.float128)

    for l in range(npart):

        # find and store all neighbours;
        neighbours[l] = find_neighbours(l, x, y, h, fact=fact, L=L, periodic=periodic)

        # compute normalisation omega for all particles
        # needs psi_k_at_l to be computed already
        omega[l] =  np.sum(psi_k_at_l[:, l])
        # omega_k = sum_l W(x_k - x_l) = sum_l psi_l(x_k) as it is currently stored in memory




    # normalize psi's and convert to float64 for linalg module
    for k in range(npart):
        psi_k_at_l[:, k] /= omega[k]
    psi_k_at_l = np.float64(psi_k_at_l)


    # compute all matrices B_k
    B_k = np.zeros((npart), dtype=np.matrix)
    for k in range(npart):
        nbors = neighbours[k]
        # nbors now contains all neighbours l
        B_k[k] = get_matrix(x[k], y[k], x[nbors], y[nbors], psi_k_at_l[nbors, k], L=L, periodic=periodic)



    # compute all psi_tilde_k at every l
    psi_tilde_k_at_l = np.zeros((npart, npart, 2))
    for k in range(npart):
        for l in range(npart):

            dx = np.array([x[k]-x[l], y[k]-y[l]])
            psi_tilde_k_at_l[k,l] = np.dot(B_k[l], dx) * psi_k_at_l[k,l]


    # now compute A_ij for all neighbours j of i
    nbors = neighbours[pind]

    A_ij = np.zeros((len(nbors), 2), dtype=np.float64)

    for i,j in enumerate(nbors): 
        
        A = np.array([0.0,0.0])
        for k in range(npart): 
            psi_i_xk = psi_k_at_l[pind, k]
            psi_j_xk = psi_k_at_l[j, k]
            Vk = V(k, m, rho)
            temp = np.array([0.0,0.0])
            for l in range(npart):
                psi_i_xl = psi_k_at_l[pind, l]
                psi_j_xl = psi_k_at_l[j, l]
                psi_tilde_l = psi_tilde_k_at_l[l, k]

                temp += (psi_j_xk * psi_i_xl - psi_i_xk * psi_j_xl) * psi_tilde_l
            
            temp *= Vk
            A += temp
    
        A_ij[i] = A

 
    # return -A_ij: You will actually use A_ji . F in the formula
    # for the hydrodynamics, not A_ij . F
    return -A_ij








#==================================================================================================
def Aij_Ivanova(pind, x, y, h, m, rho, kernel='cubic_spline', fact=1, L=1, periodic=True):
#==================================================================================================
    """
    Compute A_ij as defined by Ivanova 2013, using the discretization by Taylor 
    expansion as Hopkins does it. Use analytical expressions for the 
    gradient of the kernels instead of the matrix representation.

    pind:           particle index for which to work for (The i in A_ij)
    x, y, m, rho:   full data arrays as read in from hdf5 file
    h:              kernel support radius array
    kernel:         which kernel to use
    fact:           factor for h for limit of neighbour search; neighbours are closer than fact*h
    L:              boxsize
    periodic:       whether to assume periodic boundaries

    returns:
        A_ij: array of A_ij, containing x and y component for every neighbour j of particle i
    """


    npart = x.shape[0]


    # compute all psi_k(x_l) for all l, k
    # first index: index k of psi: psi_k(x)
    # second index: index of x_l: psi(x_l)
    psi_k_at_l = np.zeros((npart, npart), dtype=np.float128)

    for k in range(npart):
        for l in range(npart):
            # kernels are symmetric in x_i, x_j, but h can vary!!!!
            psi_k_at_l[k,l] = psi(x[l], y[l], x[k], y[k], h[l], kernel=kernel, fact=fact, L=L, periodic=periodic)


    neighbours = [[] for i in x]
    omega = np.zeros(npart, dtype=np.float128)

    for l in range(npart):

        # find and store all neighbours;
        neighbours[l] = find_neighbours(l, x, y, h, fact=fact, L=L, periodic=periodic)

        # compute normalisation omega for all particles
        # needs psi_k_at_l to be computed already
        omega[l] =  np.sum(psi_k_at_l[:, l])
        # omega_k = sum_l W(x_k - x_l) = sum_l psi_l(x_k) as it is currently stored in memory

    grad_psi_k_at_l = get_grad_psi_k_at_l_analytical(x, y, h, omega, psi_k_at_l, 
            kernel=kernel, fact=fact)


    # normalize psi's and convert to float64 for linalg module
    for k in range(npart):
        psi_k_at_l[:, k] /= omega[k]
    psi_k_at_l = np.float64(psi_k_at_l)

    print("psi_k_at_l", psi_k_at_l)

    # now compute A_ij for all neighbours j of i
    nbors = neighbours[pind]

    A_ij = np.zeros((len(nbors), 2), dtype=np.float64)

    V_i = 1/omega[pind]

    for i,j in enumerate(nbors): 
        
        grad_psi_i_xj = grad_psi_k_at_l[pind, j]
        grad_psi_j_xi = grad_psi_k_at_l[j, pind]
        V_j = 1/omega[j]
    
        A_ij[i] = V_j * grad_psi_i_xj - V_i * grad_psi_j_xi
 
    # return -A_ij: You will actually use A_ji . F in the formula
    # for the hydrodynamics, not A_ij . F
    return -A_ij










#==================================================================================================================
def Aij_Ivanova_analytical_gradients(pind, x, y, h, m, rho, kernel='cubic_spline', fact=1, L=1, periodic=True):
#==================================================================================================================
    """
    Compute A_ij as defined by Ivanova 2013. Use analytical expressions for the 
    gradient of the kernels instead of the matrix representation.

    pind:           particle index for which to work for (The i in A_ij)
    x, y, m, rho:   full data arrays as read in from hdf5 file
    h:              kernel support radius array
    kernel:         which kernel to use
    fact:           factor for h for limit of neighbour search; neighbours are closer than fact*h
    L:              boxsize
    periodic:       whether to assume periodic boundaries

    returns:
        A_ij: array of A_ij, containing x and y component for every neighbour j of particle i
    """


    npart = x.shape[0]


    # compute all psi_k(x_l) for all l, k
    # first index: index k of psi: psi_k(x)
    # second index: index of x_l: psi(x_l)
    psi_k_at_l = np.zeros((npart, npart), dtype=np.float128)

    for k in range(npart):
        for l in range(npart):
            # kernels are symmetric in x_i, x_j, but h can vary!!!!
            psi_k_at_l[k,l] = psi(x[l], y[l], x[k], y[k], h[l], kernel=kernel, fact=fact, L=L, periodic=periodic)

 
    neighbours = [[] for i in x]
    omega = np.zeros(npart, dtype=np.float128)


    for l in range(npart):

        # find and store all neighbours;
        neighbours[l] = find_neighbours(l, x, y, h, fact=fact, L=L, periodic=periodic)

        # compute normalisation omega for all particles
        # needs psi_k_at_l to be computed already
        omega[l] =  np.sum(psi_k_at_l[:, l])
        # omega_k = sum_l W(x_k - x_l) = sum_l psi_l(x_k) as it is currently stored in memory


    grad_psi_k_at_l = get_grad_psi_k_at_l_analytical(x, y, h, omega, psi_k_at_l, 
            kernel=kernel, fact=fact)



    # normalize psi's and convert to float64 for linalg module
    for k in range(npart):
        psi_k_at_l[:, k] /= omega[k]
    psi_k_at_l = np.float64(psi_k_at_l)


    # now compute A_ij for all neighbours j of i
    nbors = neighbours[pind]

    A_ij = np.zeros((len(nbors), 2), dtype=np.float64)

    for i,j in enumerate(nbors): 
        
        A = np.array([0.0,0.0], dtype=np.float64)
        for k in range(npart):
        #  for
            psi_i_xk = psi_k_at_l[pind, k]
            psi_j_xk = psi_k_at_l[j, k]
            grad_psi_i_xk = grad_psi_k_at_l[pind, k]
            grad_psi_j_xk = grad_psi_k_at_l[j, k]
            V_k = 1/omega[k]

            A += (psi_j_xk * grad_psi_i_xk - psi_i_xk*grad_psi_j_xk)*V_k
    
        A_ij[i] = A

 
    # return -A_ij: You will actually use A_ji . F in the formula
    # for the hydrodynamics, not A_ij . F
    return -A_ij











#=========================================================================================================================
def Integrand_Aij_Ivanova(iind, jind, xx, yy, hh, x, y, h, m, rho, kernel='cubic_spline', fact=1, L=1, periodic=True):
#=========================================================================================================================
    """
    Compute the effective area integrand for the particles iind jind at
    the positions xx, yy

    (Note that this should be integrated to get the proper effective surface)

    integrand A_ij  = psi_j(x) \nabla psi_i(x) - psi_i (x) \nabla psi_j(x)
                    = sum_k [ psi_j(x_k) psi_i(x) - psi_i(x_k) psi_j(x) ] * psi_tilde_k(x)
                    = psi_i(x) * sum_k psi_j(x_k) * psi_tilde_k(x) - psi_j(x) * sum_k psi_i(x_k) * psi_tilde_k(x)
    
    The last line is what is actually computed here, with the expression for the gradient
    inserted.


    iind, jind:     particle index for which to work for (The i and j in A_ij)
    xx, yy:         position at which to evaluate
    hh:             kernel support radius at xx, yy
    x, y, m, rho:   full data arrays as read in from hdf5 file
    h:              kernel support radius array
    kernel:         which kernel to use
    fact:           factor for h for limit of neighbour search; neighbours are closer than fact*h
    L:              boxsize
    periodic:       whether to assume periodic boundaries

    returns:
        A_ij: array of integrands A_ij, containing x and y component for every neighbour j of particle i
  


    """

    nbors = find_neighbours_arbitrary_x(xx, yy, x, y, h, fact=fact, L=L, periodic=periodic)

    xk = x[nbors]
    yk = y[nbors]
    hk = h[nbors]



    #----------------------
    # compute psi_i/j(x)
    #----------------------

    # compute all psi(x)
    psi_x = compute_psi(xx, yy, xk, yk, hh, kernel=kernel, fact=fact, L=L, periodic=periodic)

    # normalize psis
    omega = np.sum(psi_x)
    psi_x /= omega
    psi_x = np.float64(psi_x)

    # find where psi_i and psi_j are in that array
    try:
        inb = nbors.index(iind)
        psi_i_of_x = psi_x[inb] # psi_i(xx, yy)
    except ValueError:
        psi_i_of_x = 0 # can happen for too small smoothing lengths
        print("Exception in psi_i_of_x: iind not found in neighbour list")

    try:
        jnb = nbors.index(jind)
        psi_j_of_x = psi_x[jnb] # psi_j(xx, yy)
    except ValueError:
        psi_j_of_x = 0 # can happen for too small smoothing lengths
        print("Exception in psi_j_of_x: jind not found in neighbour list")




    #------------------------------------------------
    # Compute psi_i/j(x_k) at neighbouring positions
    #------------------------------------------------

    psi_i_xk = [None for n in nbors]
    psi_j_xk = [None for n in nbors]


    omegas = [0, 0]
    
    for i, n in enumerate([iind, jind]):
        # first compute all psi(xl) from neighbour's neighbours to get weights omega
        nneigh = find_neighbours(n, x, y, h, fact=fact, L=L, periodic=periodic)

        xl = x[nneigh]
        yl = y[nneigh]

        for j, nn in enumerate(nneigh):
            psi_l = compute_psi(x[n], y[n], xl, yl, h[n], kernel=kernel, fact=fact, L=L, periodic=periodic)

        omegas[i] = np.sum(psi_l) + psi(0, 0, 0, 0, h[iind], kernel=kernel, fact=fact, L=L, periodic=periodic)


    # now compute psi_i/j(x_k)
    for i, n in enumerate(nbors):
        psi_i_xk[i] = psi(xk[i], yk[i], x[iind], y[iind], h[iind], kernel=kernel, fact=fact, L=L, periodic=periodic) / omegas[0]
        psi_j_xk[i] = psi(xk[i], yk[i], x[jind], y[jind], h[jind], kernel=kernel, fact=fact, L=L, periodic=periodic) / omegas[1]





    #---------------------------------------------
    # Compute psi_tilde_k(x)
    #---------------------------------------------

    # compute matrix B
    B = get_matrix(xx, yy, xk, yk, psi_x, L=L, periodic=periodic)

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








#========================================================================================================================
def get_grad_psi_k_at_l_analytical(x, y, h, omega, psi_k_at_l, kernel='cubic_spline', fact=1, L=1, periodic=True):
#========================================================================================================================
    """
    Compute \nabla \psi_k (x_l) for all particles k and l 
    x, y, h:    arrays of positions and compact support radius of all particles
    omega:      weights; sum_j W(x - xj) for all particles x=x_k
    psi_k_at_l: UNNORMED psi_k(x_l) npart x npart array for all k, l
    kernel:     which kernel to use
    fact:       factor for h for limit of neighbour search; neighbours are closer than fact*h
    L:          boxsize
    periodic:   whether to assume periodic boundaries

    returns:

        grad_psi_k_at_l: npart x npart x 2 array; grad psi_k (x_l) for all k,l for both x and y direction
    """



    npart = x.shape[0]

    grad_psi_k_at_l = np.zeros((npart, npart, 2), dtype=np.float128)
    grad_W_k_at_l = np.zeros((npart, npart, 2), dtype=np.float128)


    for k in range(npart):
        for l in range(npart):
            # get kernel gradients

            dx, dy = get_dx(x[l], x[k], y[l], y[k], L=L, periodic=periodic)

            r = np.sqrt(dx**2 + dy**2)
            if r == 0:
                grad_W_k_at_l[k, l, 0] = 0
                grad_W_k_at_l[k, l, 1] = 0
            else:
                grad_W_k_at_l[k, l, 0] = dWdr(r/h[l], h[l], kernel) * dx / r
                grad_W_k_at_l[k, l, 1] = dWdr(r/h[l], h[l], kernel) * dy / r



    sum_grad_W = np.zeros((npart, 2), dtype=np.float128)

    for l in range(npart):
        sum_grad_W[l] = np.sum(grad_W_k_at_l[:, l], axis=0)


    # first finish computing the gradients: Need W(r, h), which is currently stored as psi
    for k in range(npart):
        for l in range(npart):
            grad_psi_k_at_l[k, l, 0] = grad_W_k_at_l[k, l, 0]/omega[l] - psi_k_at_l[k, l] * sum_grad_W[l, 0]/omega[l]**2
            grad_psi_k_at_l[k, l, 1] = grad_W_k_at_l[k, l, 1]/omega[l] - psi_k_at_l[k, l] * sum_grad_W[l, 1]/omega[l]**2


    return grad_psi_k_at_l















#==========================================================================================
def compute_psi(xi, yi, xj, yj, h, kernel='cubic_spline', fact=1, L=1, periodic=True):
#==========================================================================================
    """
    Compute all psi_j(x_i)
    xi, yi:     floats; position for which to compute psi's
    xj, yj:     arrays of neighbour's positions
    h:          float; smoothing length at position xi, yi
            or array of h for xj, yj [used to compute h(x)]
    kernel:     which kernel to use
    fact:       factor for h for limit of neighbour search; neighbours are closer than fact*h
    L:          boxsize
    periodic:   whether to assume periodic boundaries

    return numpy array of psi_j(x) 
    """

    psi_j = np.zeros(xj.shape[0], dtype=np.float128)

    if isinstance(h, np.ndarray):
        for i in range(xj.shape[0]):
            psi_j[i] = psi(xi, yi, xj[i], yj[i], h[i], kernel=kernel, fact=fact, L=L, periodic=periodic)

    else:
        for i in range(xj.shape[0]):
            psi_j[i] = psi(xi, yi, xj[i], yj[i], h, kernel=kernel, fact=fact, L=L, periodic=periodic)

    return psi_j










#============================================================================
def psi(x, y, xi, yi, h, kernel='cubic_spline', fact=1, L=1, periodic=True):
#============================================================================
    """
    UNNORMALIZED Volume fraction at position x of some particle
    with coordinates xi, yi, smoothing length h(x)

    i.e. psi_i(x) = W([x - xi, y - yi], h(x))

    kernel:     which kernel to use
    fact:       factor to increase h with
    L:          boxsize
    periodic:   Whether you assume periodic boundary conditions

    !!!! returns type np.float128! 
    Needed to prevent precision errors for normalisation
    """


    dx, dy = get_dx(x, xi, y, yi, L=L, periodic=periodic)

    q = np.float128(np.sqrt(dx**2 + dy**2)/(fact*h))

    return W(q, h, kernel)








#============================================================
def get_matrix(xi, yi, xj, yj, psi_j, L=1, periodic=True):
#============================================================
    """
    Get B_i ^{alpha beta}

    xi, yi:         floats; Evaluate B at this position
    xj, yj:         arrays; Neighbouring points
    psi_j:          array;  volume fraction of neighbours at position x_i; psi_j(x_i)
    L:              boxsize
    periodic:       whether to assume periodic boundaries
    """

    dx = np.zeros(xj.shape[0])
    dy = np.zeros(xj.shape[0])

    for i in range(xj.shape[0]):
        dx[i], dy[i] = get_dx(xj[i], xi, yj[i], yi, L=L, periodic=periodic)

    E00 = np.sum(dx * dx * psi_j)
    E01 = np.sum(dx * dy * psi_j)
    E11 = np.sum(dy * dy * psi_j)
          
    E = np.matrix([[E00, E01], [E01, E11]])

    try:
        return E.getI()
    except np.linalg.LinAlgError:
        print("Exception: Singular Matrix")
        print("E:", E)
        print("dx:", xj - xi)
        print("dy:", yj - yi)
        print("psi:", psi_j)
        quit(2)

    return






#=============================================================================================
def h_of_x(xx, yy, x, y, h, m, rho, kernel='cubic_spline', fact=1, L=1, periodic=True):
#=============================================================================================
    """
    Compute h(x) at position (xx, yy), where there is 
    not necessariliy a particle
    by approximating it as h(x) = sum_j h_j * psi_j(x)

    x, y, h :   full particle arrays
    fact:       factor to increase h with
    L:          boxsize
    periodic:   whether to assume periodic boundaries
    """

    nbors = find_neighbours_arbitrary_x(xx, yy, x, y, h, fact=fact, L=L, periodic=periodic)

    xj = x[nbors]
    yj = y[nbors]
    hj = h[nbors]

    psi_j = compute_psi(xx, yy, xj, yj, hj, kernel=kernel, fact=fact, L=L, periodic=periodic)
    psi_j /= np.sum(psi_j)

    hh = np.sum(hj*psi_j)

    return hh




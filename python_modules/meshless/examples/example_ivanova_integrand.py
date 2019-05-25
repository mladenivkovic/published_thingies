#!/usr/bin/env python3

#===============================================================
# Compute A(x) between two specified particles at various
# positions x
#===============================================================


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size 

try:
    import meshless as ms
except ImportError:
    print("Didn't find 'meshless' module... be sure to add it to your pythonpath!")
    quit(2)



srcfile = './snapshot_perturbed.hdf5'
ptype = 'PartType0'             # for which particle type to look for


# border limits for plots
lowlim = 0.37
uplim = 0.53
# number of points to plot in each direction
nx = 50
# tolerance for particle finding
tol = 5e-2






#========================
def main():
#========================
    

    #-----------------------------
    # Part1 : compute all A
    #-----------------------------

    print("Computing effective surfaces")

    # read in file
    x, y, h, rho, m, ids, npart = ms.read_file(srcfile, ptype)

    # transform smoothing length to compact support radius
    H = ms.get_H(h)

    # find where particles i (0.4, 0.4) and j (0.5, 0.5) are
    iind = None
    jind = None

    for i in range(x.shape[0]):
        if abs(x[i] - 0.4) < tol and abs(y[i] - 0.4) < tol:
            iind = i
        if abs(x[i] - 0.5) < tol and abs(y[i] - 0.5) < tol:
            jind = i


    if iind is None or jind is None:
        raise ValueError("iind=", iind, "jind=", jind)

    

    A = np.zeros((nx, nx, 2), dtype=np.float) # storing computed effective surfaces
    dx = (uplim - lowlim)/nx


    for i in range(nx):
        xx = lowlim + dx * i

        print("i = ", i, "/", nx)

        for j in range(nx):
            yy = lowlim + dx * j

            hh = ms.h_of_x(xx, yy, x, y, H, m, rho)

            A[j, i] = ms.Integrand_Aij_Ivanova(iind, jind, xx, yy, hh, x, y, H, m, rho) # not a typo: need A[j,i] for imshow






    #-----------------------------
    # Part2: Plot results
    #-----------------------------

    print("Plotting")

    fig = plt.figure(figsize=(14,5))
    ax1 = fig.add_subplot(131, aspect='equal')
    ax2 = fig.add_subplot(132, aspect='equal')
    ax3 = fig.add_subplot(133, aspect='equal')


    Ax = A[:,:,0] 
    Ay = A[:,:,1]
    Anorm = np.sqrt(Ax**2 + Ay**2)
    xmin = Ax.min()
    xmax = Ax.max()
    ymin = Ay.min()
    ymax = Ay.max()
    normmin = Anorm.min()
    normmax = Anorm.max()


    # reset lowlim and maxlim so cells are centered around the point they represent
    dx = (uplim - lowlim) / A.shape[0]

    cmap = 'YlGnBu_r'



    im1 = ax1.imshow(Ax, origin='lower', 
            vmin=xmin, vmax=xmax, cmap=cmap,
            extent=(lowlim, uplim, lowlim, uplim))
    im2 = ax2.imshow(Ay, origin='lower', 
            vmin=ymin, vmax=ymax, cmap=cmap,
            extent=(lowlim, uplim, lowlim, uplim))
    im3 = ax3.imshow(Anorm, origin='lower', 
            vmin=normmin, vmax=normmax, cmap=cmap,
            extent=(lowlim, uplim, lowlim, uplim))

    for ax, im in [(ax1, im1), (ax2, im2), (ax3, im3)]:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.05)
        fig.colorbar(im, cax=cax)




    # superpose particles

    inds = np.argsort(ids)

    mask = np.logical_and(x>=lowlim-tol, x<=uplim+tol)
    mask = np.logical_and(mask, y>=lowlim-tol)
    mask = np.logical_and(mask, y<=uplim+tol)



    for ax in [ax1, ax2, ax3]:

        ps = 50
        fc = 'grey'
        ec = 'black'
        lw = 2

        # plot neighbours (and the ones you drew anyway)
        ax.scatter(x[mask], y[mask], s=ps, lw=lw,
                facecolor=fc, edgecolor=ec)

        # plot the chosen one
        ps = 100
        fc = 'white'
        ax.scatter(x[iind], y[iind], s=ps, lw=lw,
                facecolor=fc, edgecolor=ec)


        # plot central (and the ones you drew anyway)
        fc = 'black'
        ax.scatter(x[jind], y[jind], s=ps, lw=lw,
                facecolor=fc, edgecolor=ec)


        ax.set_xlim((lowlim,uplim))
        ax.set_ylim((lowlim,uplim))

        ax.set_xlabel('x')
        ax.set_ylabel('y')

        ax1.set_title(r'$x$ component of $\mathbf{A}_{ij}$')
        ax2.set_title(r'$y$ component of $\mathbf{A}_{ij}$')
        ax3.set_title(r'$|\mathbf{A}_{ij}|$')



    fig.suptitle(r'Effective Area $\mathbf{A}_{ij}(\mathbf{x}) = \psi_i(\mathbf{x}) \nabla \psi_j(\mathbf{x}) - \psi_j(\mathbf{x}) \nabla \psi_i(\mathbf{x})$ of a particle (white) w.r.t. the central particle (black) in a uniform distribution')
    plt.tight_layout()
    plt.show()
    #  plt.savefig('effective_area_A_of_x.png', dpi=300)

    print('finished.')

    return





if __name__ == '__main__':
    main()


#!/usr/bin/env python3

#=============================
# Contains kernels
#=============================

import numpy as np


# Names of all available kernels
kernels = ['cubic_spline', 'quintic_spline', 
        'gaussian', 'gaussian_compact', 'supergaussian',
        'wendland_C2', 'wendland_C4', 'wendland_C6']

# factors of all kernels for which fact*H = 0
kernelfacts = [2, 2,
        1000, 2, 2,
        2, 2, 2]

kernels_shortlist = ['cubic_spline', 'quintic_spline',
        'wendland_C2', 'wendland_C4', 'wendland_C6']

kernel_shortlist_facts = [2, 2,
        2, 2, 2]

kernel_shortlist_H_over_h = [1.778002, 2.158131,
        1.897367, 2.171239, 2.415230]







#=====================================
def W(q, h, kernel='cubic_spline'):
#=====================================
    """
    Various kernels

    Currently implemented:
        cubic_spline, 
        quintic_spline, 
        gaussian, 
        gaussian_compact, 
        supergaussian,
        wendland_C2, 
        wendland_C4, 
        wendland_C6
    """ 
    #  https://pysph.readthedocs.io/en/latest/reference/kernels.html#liu2010

    if kernel == 'cubic_spline': 

        sigma = 10./(7*np.pi*h**2)
        if q < 1:
            return 1. - q*q * (1.5 - 0.75*q) 
        elif q < 2:
            return 0.25*(2-q)**3
        else:
            return 0




    elif kernel == 'quintic_spline':
        # re-scale to H = 3
        q = 1.5*q
        sigma = 9.0/4.0 * 7.0/(478*np.pi*h*h)

        if q < 1:
            return sigma * ((3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5)
        elif q<2:
            return sigma * ((3-q)**5 - 6*(2-q)**5)
        elif q<3:
            return sigma * ((3-q)**5)
        else:
            return 0





    elif kernel == 'gaussian':
        # gaussian without compact support
        return 1./(np.sqrt(0.5*np.pi)*h)**3*np.exp(-2*q**2)

   



    elif kernel == 'gaussian_compact':
        # gaussian with compact support
        # re-scale to H = 3
        q = 1.5*q

        sigma = 1./(np.pi*h*h)

        if q <= 3:
            return sigma * np.exp(-q**2)
        else:
            return 0




    elif kernel == 'supergaussian':
        # re-scale to H = 3
        q = 1.5*q

        if q <= 3:
            sigma = 1./(np.sqrt(np.pi)*h)**3
            return sigma * np.exp(-q**2)*(2.5 - q**2)
        else:
            return 0





    elif kernel == 'wendland_C2':

        if q <= 2:
            sigma = 7/(4*np.pi * h**2)
            return sigma * (1 - 0.5*q)**4*(2*q+1)
        else:
            return 0





    elif kernel == 'wendland_C4':

        if q <= 2:
            sigma = 9/(4*np.pi*h**2)
            return sigma*(1-0.5*q)**6 * (35/12*q**2 + 3*q + 1)
        else:
            return 0






    elif kernel == 'wendland_C6':
        
        if q <= 2:
            sigma = 78/(28*np.pi*h**2)
            return sigma * (1 - 0.5*q)**8*(4*q**3 + 6.25*q**2 + 4*q + 1)
        else:
            return 0





    else:
        raise ValueError("Didn't find kernel", kernel)


    return

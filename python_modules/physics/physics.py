#!/usr/bin/env python3


def get_help():
    help_message =   """
        ==========================
        === THE PHYSICS MODULE ===
        ==========================

        
        1) CONSTANTS
        -------------------

        The physics module contains physical constants, available 
        in various units.

        Possible choices for units:

            1) SI
            -------------------------------
            Mass:           kg
            Length:         m
            Time:           s
            Temperature:    K
            Charge:         C
            Energy:         J       [1J = 1kg * m^2/s^2]


            2) cgs (centimetre-gram-second)
            --------------------------------
            Mass:           g
            Length:         cm
            Time:           s
            Temperature:    K
            Charge:         esu     [electrostatic unit]
            Energy:         erg     [1erg = 10^-7 J]


            3) eV
            --------------------------------
            Mass:           eV/c^2
            Length:         m
            Time:           s
            Temperature:    eV/k_B
            Charge:         C
            Energy:         eV      [1eV = J]

         
        To see the constants you are currently using, call 'physics.get_const()'.
        To change the units of the constants, call 'physics.set_const(str units),
        where str units can be ['SI', 'cgs', 'eV'].
        Default units are SI.


        2) FUNCTIONS
        -----------------
        Following physical functions are ready to use:
            gamma(v)            gives Lorentz factor
        

        """

    print(help_message)
    return



#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================

#============
# CONSTANTS
#============




value_dict = dict()

names =       ['pi',    'c',                'e',                    'm_e',              'm_p',          'm_n',          'm_mu',         'u',                'k_B',                  'sigma',                        'N_A',                        'epsilon_0',                        'mu_0',                      'h',               'hq',                       'a_0',          'alpha',                    'G',                      'r_e',           'M_e',         'r_Sol',      'M_Sol',    'au',                'pc',      'sigma_T',                      'H_0',                          "Joule_per_eV",   "kCal_per_Joule",    'day',                    'year'                  ]
description = ['pi',    'speed of light',   'elementary charge',    'elektron mass',    'proton mass',  'neutron mass', 'muon mass',    'atomic mass unit', 'Boltzmann constant',   'Stefan-Boltzmann constant',    "Avogadro's number [mol^-1]", 'electric constant [SI:As/(Vm)]',   'magnetic constant [Vs/Am]', 'Planck constant', 'reduced Planck constant',  'Bohr radius',  'fine structure constant',  'gravitational constant', 'Earth radius',  'Earth mass',  'Sun radius', 'Sun mass', 'astronomical unit', 'parsec',  'Thomson cross section [L^2]',  'Hubble constant [1/years]',    'J/eV',           'kCal/J',            'seconds in a day',       "seconds in a year"     ]
value_SI =    [3.14159, 2.998e8,            1.602e-19,              9.109e-31,          1.672e-27,      1.674e-27,      1.884e-28,      1.660538e-27,       1.381e-23,              5.671e-8,                        6.022e+23,                   8.854e-12,                          4*3.14159*1e-7 ,             6.626e-34,         1.055e-34,                  5.3e-11,        7.297e-3,                    6.673e-11,               6371e3,          5.972e24 ,     695700e3,     1.988e30,   1.496e11,            3.0857e16,  6.65e-29,                      1/13.7e9,                       1.602e-19,        4184,                24*3600,                  24*3600*365             ]
value_eV =    [3.14159, 2.998e8,            1.602e-19,              0.511e6,            939.565e6,      939.565e6,      931.988e6,      105.658e6,          8.617e-5,               5.671e-8,                        6.022e+23,                   8.854e-12,                          4*3.14159*1e-7 ,             4.136e-15,         6.583e-16,                  5.3e-11,        7.297e-3,                    6.673e-11,               6371e3,          3.359e60,      695700e3,     1.115e66,   1.496e11,            3.0857e16,  6.65e-29,                      1/13.7e9,                       1.602e-19,        4184,                24*3600,                  24*3600*365             ]
value_cgs =   [3.14159, 2.998e10,           4.803e-10,              9.109e-28,          1.672e-24,      1.674e-24,      1.884e-25,      1.660538e-24,       1.381e-16,              5.671e-5,                        6.022e+23,                   1.0/(4.0*3.14159) ,                 4*3.14159/(2.998e10)**2 ,    6.626e-27,         1.055e-27,                  5.3e-9,         7.297e-3,                    6.673e-8,                6371e5,          5.972e27,      695700e5,     1.988e33,   1.496e13,            3.0857e18,  6.65e-25,                      1/13.7e9,                       1.602e-19,        4184,                24*3600,                  24*3600*365             ]





def set_const(units='SI'):
 
    if units=='SI':
        constant_list=value_SI
    elif units=='cgs':
        constant_list=value_cgs
    elif units=='eV':
        constant_list=value_eV
    else:
        print("Something went wrong.")
        return
        
    # go through list of constants; assign name that is stored in list the corresponding value
    for x,y in zip(names,constant_list):
        globals()[x] = y

    # write the corresponding dict for chosen units so it will be printed correctly when using get_const!
    global value_dict
    value_dict = dict(zip(names, constant_list))


    return










def get_const():

    #print header
    print("="*64)
    print("="*18, " "*3 , "DEFINED CONSTANTS",  " "*4 , "="*18)
    print("="*64)
    print("")
    print('{0:20}{1:30}{2:14}'.format("constant name", "description","value" ))
    print("-"*64)


    manual_breaks=['pi', 'u', 'N_A', 'mu_0', 'alpha', 'H_0','kCal_per_Joule']

    #print constants
    for i,const in enumerate(names):
        print('{0:20}{1:30}{2:14.4E}'.format(const, description[i], value_dict[const] ))
        if const in manual_breaks:
            print()
    
    return





#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================

#============
# FUNCTIONS 
#============

def gamma(v):
    import math
    return math.sqrt(1.0-v**2/c**2)




#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================

#=================
# INITIAL SET UP
#=================

#Call initial set-up
if __name__ != "__main__": #if not executed directly, but for example imported

    #Greeting message
    print("Welcome to the physics module!\nIf you need help, call get_help()")
    
    #Set up constants
    set_const()



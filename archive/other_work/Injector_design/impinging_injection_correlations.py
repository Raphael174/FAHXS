""" 
@ author : Raphael Aubry



"""

import numpy as np 

def jet_breakup_length_single_jet_GrantMiddleman(Weber_jet, do):
    """
    Break up length of single turbulent jet
    Grant, R. P. and Middleman, S., “Newtonian Jet Stability,” AIChE Journal ,
    Vol. 12, No. 4, 1966, pp. 669–678.

    ratio of length break up (single jet) / impingingeement length -> used by Sweenie
    """
    return do*8.51*Weber_jet**(0.32)

def length_impingement(alpha, distance_jet_orifice):
    """ 
    Compute length from orifice exit to impingement point
    distance_jet_orifice is either the hole-hole distance of classic like-doublet OR
                                    the gas-diameter (hole-hole) distance
    alpha is the doublet impingement angle (not the half angle)
    """
    return (distance_jet_orifice/2)/np.sin(alpha*np.pi/180/2)

def Webber_jet(density_L, velocity_jet, l, surface_tension_L):
    return density_L*velocity_jet**2*l/surface_tension_L

def Webber_transition_breakup_mode_like_doublet_Sweenie(theta, l_b, l_i):
    """ 
    Transition Webber that defines the transition point between the ruffled-sheet 
    and fully-developed breakup modes

    Valid for for lb/li>=1
    """
    if l_b/l_i<1:
        raise Warning("Transition Webber correlation not built for lb/li<1")
    
    return 455*(np.sin(theta*np.pi/180))**(-2.5)


def sheet_breakup_length_like_doublet_Sweenie(Weber_jet, Weber_transition, do, theta, l_b, l_i):
    """ 
    Sheet break up length for either Ruffled-sheet (RS) mode and Fully-developped (FD) 
    """
    if l_b/l_i<1:
        raise Warning("Sheet break up length correlation not built for lb/li<1")
    
    if Weber_jet<Weber_transition:
        print("Ruffled mode like_doublet break up")
        return do * 3.47*Weber_jet**(0.22) * np.sin(theta*np.pi/180)**(-0.42)
    
    else:
        print("Fully-developped mode like_doublet break up")
        return do * 51.28*Weber_jet**(-0.22) * np.sin(theta*np.pi/180)**(-1.52)


    
def d10_like_doublet_Sweenie(Weber_jet, theta, l_b, l_i):
    """ 

    lb= break-up length of liquid jet
    li=length of impingement 

    Good breakup and atomization for lb/li>=1
    """

    if l_b/l_i<1:
        raise Warning("Like doublet d10 correlation not built for lb/li<1")
    return 4751e-6 * Weber_jet**(-0.38)*(np.sin(theta*np.pi/180))**(-0.46)
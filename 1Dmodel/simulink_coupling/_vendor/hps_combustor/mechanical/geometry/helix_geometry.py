""" 
@ author : Raphaël Aubry
"""

import numpy as np
from scipy.integrate import simpson
# from scipy.interpolate import interp2d #interp1d
from scipy.optimize import curve_fit
from scipy.optimize import fsolve

def compute_Dh_shell(D_coil, d_coil_outer, shell_diameter, coil_pitch):
    """ 
    Shell with coil tube hydraulic diameter according to Salimpour 
    10.1016/j.expthermflusci.2008.07.015
    """
    R_coil = D_coil/2
    gamma = coil_pitch/(2*np.pi*R_coil)
    return (shell_diameter**2 - 2*np.pi*R_coil*d_coil_outer**2/gamma) / (shell_diameter + 2*np.pi*R_coil*d_coil_outer/gamma)

def HelixGeometryRadiusCST (coil_pitch, D_coil, L_coil):

    """ 
    Function to produce coordinates and path length data for a given helical
    HX around the combustor

    This enables direct computation of local nozzle area for each location 
    along HX channel

    pitch_length : distance between 2 consecutive coils

    Equations of helix coiling around x axis : 
        x = h*t 
        y = R_coil*cos(t)
        z = R_coil*cos(t)

        where pitch length = 2*pi*h

    Equation for arc length s of any curve with certain polar coordinates: 

        s = integral|theta0-->theta1 (sqrt(R^2 + (dr/dtheta)^2)) d theta

        where R is the helix local radius as a function of theta
        and dr/dtheta is its first derivative with theta
        theta1 can be much above 2*pi since more than one turn is possible

    returns (x,y,z,s,x_to_R_HX,x_to_R_engine)

            x,y,z : cartesian coordinates of HX coil centerline 
                    where coil wraps around x axis
            s : path length (arc length)

    
    """
    h = coil_pitch/(2*np.pi)
    theta_max = L_coil/h

    # run through combustor for N_turns
    angular_pos = np.linspace(0, theta_max, 10000, endpoint=True) # 10 turns 
    # first extract the axial (x) coordinate of the helix (needed for plotting Rcc engine vs x)
    # x is HX pipe centerline position, relative to its distance from the nozzle exit
    x = h*angular_pos 

    # convert angular positions to remaining axial coordinates of helix
    y = (D_coil/2)*np.cos(angular_pos)
    z = (D_coil/2)*np.sin(angular_pos)

    # arc length derivatives method 2 
    dx__dt = np.gradient(x, angular_pos)
    dy__dt = np.gradient(y, angular_pos)
    dz__dt = np.gradient(z, angular_pos)

    # building pipe path length array
    s = np.zeros(len(angular_pos))

    #! pre-compute arc length cummulation to link it to helix coordinates
    for i in range(len(angular_pos)) : 
        # initially 0m arc length
        if i == 0 :
            s[i] = 0
        # numerical integration according to arc length equation 
        else:
            # term inside arc length integral up to current index
            #term = np.sqrt((x_to_R_HX[:i])**2 + (dr__dtheta[:i])**2 + h**2)
            integrand = np.sqrt(dx__dt[:i]**2 + dy__dt[:i]**2 + dz__dt[:i]**2)
            s[i] = simpson(y=integrand, x=angular_pos[:i])

    def poly1(x, a, b):
        return a * x 

    def poly2(x, a, b, c):
        return a * x**2 + b * x + c

    params_x, _ = curve_fit(poly1, s, x)
    params_theta, _ = curve_fit(poly1, s, angular_pos)

    # ready the functions for quick use in solver
    def func_s_to_x (s):
        return poly1(s, *params_x)
    def func_s_to_theta (s):
        return poly1(s, *params_theta)
    
    L_max_pipe = s[-1]

    return (func_s_to_x, func_s_to_theta, L_max_pipe)


# def Lcoil_to_pipeLength (func_s_to_x, Dh_ch, s_w, L_coil):

#     def f_x (s):

#         return  - func_s_to_x(s)

# self.T_wg, self.T_wc, self.T_c_check_f = fsolve(func=f_x, x0=[self.T_wg_0, self.T_wc_0, self.T_c_check_0], xtol=1e-8)
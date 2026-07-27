""" 
@author : Raphaël Aubry

Injector sizing for combustor-HX design & manufacturing
"""

#%%
import numpy as np 
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

#%%
def Cd_swirl_simplex_RizkLefebvre1985 (Ap, Ds, do):
    """ 
    Rizk, N.K. and Lefebvre, A.H., "Internal Flow Characteristics
    of Simplex Swirl Atomizers," Journal of Propulsion and Power, Vol.
    1, May-June 1985, pp. 193-199.

    Ap : total inlet ports area, m²
    Ds : swirl chamber diameter, m
    do : discharge orifice diameter, m
    """
    return 0.35*(Ap/(Ds*do))**0.5 * (Ds/do)**0.25 




def FlowNumber_swirl_simplex_SuyariLefebvre1986 (mass_flow_L, density_L, dP_L):
    """
    Flow number of simplex swirl injector
    M. Suyari, A.H. Lefebvre, Film thickness measurements in a simplex swirl atomizer, 
    J. Propuls. Power 2 (1986) 528–533.
    """
    return mass_flow_L/(dP_L*density_L)**0.5

def u_axial_swirl_simplex (mass_flow_L, density_L, Ao, Aa):
    """
    Atomization and Sprays, 2nd edition, p 113 eq. 5.32

    Ao : total orifice area, m²
    Aa : air core area, m² 

    """
    return mass_flow_L/(density_L*(Ao - Aa))

def film_thickness_swirl_simplex_SuyariLefebvre1986 (do, FN, viscosity_L, dP_L, density_L):
    """ 
    Film thickness of simplex swirl injector
    M. Suyari, A.H. Lefebvre, Film thickness measurements in a simplex swirl atomizer, 
    J. Propuls. Power 2 (1986) 528–533.

    The value of the constant is unclear. 3.66 in the original work and 2.7 in other cited works, probably a unit thing
    """
    beta = 3.66 #
    return beta*(do*FN*viscosity_L/(dP_L*density_L)**0.5)**0.25 

def film_thickness_swirl_simplex_RizkLefebvre (FN, density_L, viscosity_L, dP_L, do):
    """
    N.K. Rizk and A.H. Lefebvre, 1985
    """

    t0 = do/10

    def f_x(t_):

        X = (do - 2*t_)**2/do**2

        fluid_term = viscosity_L/(density_L**0.5*dP_L**0.5)
        X_term = (1 + X)/(1 - X)**2   

        return 1560*FN*fluid_term/do * X_term - t_**2
    
        
    t = fsolve(func=f_x, x0=[t0], xtol=1e-8)


    return t


def VelocityCoefficient_swirl_simplex_RizkLefebvreDATE(K, dP_L, density_L, viscosity_L):
    """ 
    Ratio of the actual discharge velocity to the theoretical velocity corresponding to the total
    pressure differential across the nozzle
    """
    return 0.00367* K**0.29 * (dP_L*density_L/viscosity_L)**0.2

def DischargeVelocity_swirl_simplex(Kv, dP_L, density_L):
    """
    Kv : Velocity coefficient 
    """
    return Kv*(2*dP_L/density_L)**0.5

def SMD_swirl_simplex_WangLefebvreDATE (surface_tension_L, viscosity_L, density_gas, dP_L, t, theta):
    """
    Sauter Mean Diameter based on fluide properties, film thickness and spray cone angle
    Wang and Lefebvre, Atomization and Sprays, Eq. 6.31 

    theta : spray cone anle, rad
    t : film thickness, m

    """

    term1 = surface_tension_L/density_gas
    term2 = t*np.cos(theta)

    return  4.52*(term1*viscosity_L**2/dP_L**2)**0.25 * term2**0.25 \
            + 0.39*(term1*density_gas/dP_L)**0.75 * term2**0.75


def spray_angle_swirl_simplex_RizkLefebvre(mass_flow_L, viscosity_L, density_L, do, t, dP_L, X, cst=400):
    """
    N.K. Rizk and A.H. Lefebvre, 1985
    """
    costheta = np.sqrt(12*mass_flow_L*viscosity_L*cst/(np.pi*density_L*do*dP_inj*t**2*(1-X)))
    return np.arccos(costheta)

def K_swirl_simplex(Ap, Ds, do):
    return Ap/(Ds*do)


#%%

m_LOX = 19
N_inj =28
dP_inj = 10e5 
p0 = 53e5
T_inj = 115 

m_LOX_unit = m_LOX/N_inj

N_p = 3
L_char = 12e-3 # characteristic diameter of injector element
e_inj = 2e-3 # thickness of injector wall
K = 2.5 # set swirl number according to good practice
lo__do = 3 # seems legit

# starting data
Cd = 0.7

density_LOx = PropsSI('D','T', T_inj,'P',p0+dP_inj,"OXYGEN")
viscosity_LOx = PropsSI('V','T', T_inj,'P',p0+dP_inj,"OXYGEN")
tension_LOx = PropsSI('I','Q', 0,'P',50e5,"OXYGEN")

#0.35*(Ap/(Ds*do))**0.5 * (Ds/do)**0.25 

# swirl geometry computations
Ao = m_LOX_unit/(Cd*(2*density_LOx*(dP_inj))**0.5) # swirl injector mass flow
# Ao = m_LOX_unit/(Cd*(2*density_LOx*(p0+dP_inj))**0.5) # swirl injector mass flow USED FOR RESULT TO VISHESH

do = np.sqrt(4*Ao/np.pi) # area to diameter
Ds = do*(Cd/(0.35*K**0.5))**4 # Rizk, N.K. and Lefebvre, A.H. 1985
Ap = K*Ds*do
Dp = np.sqrt(4*Ap/N_p/np.pi)

FN = FlowNumber_swirl_simplex_SuyariLefebvre1986(mass_flow_L=m_LOX_unit, dP_L=dP_inj, density_L=density_LOx)
#film_thickness = film_thickness_swirl_simplex_SuyariLefebvre1986(do=do, FN=FN, viscosity_L=viscosity_LOx, dP_L=dP_inj, density_L=density_LOx)
film_thickness = film_thickness_swirl_simplex_RizkLefebvre(FN=FN, density_L=density_LOx, viscosity_L=viscosity_LOx, dP_L=dP_inj, do=do)[0]
X = (do - 2*film_thickness)**2/do**2

spray_half_angle = spray_angle_swirl_simplex_RizkLefebvre(mass_flow_L=m_LOX_unit, viscosity_L=viscosity_LOx, density_L=density_LOx, do=do, t=film_thickness, dP_L=dP_inj, X=X)
SMD = SMD_swirl_simplex_WangLefebvreDATE(surface_tension_L=tension_LOx, viscosity_L=viscosity_LOx, density_gas=0.3, dP_L=dP_inj, t=film_thickness, theta=spray_half_angle)

theta_2 =np.arccos(np.sqrt((1-X)/(1+X)))


print("First iteration of swirl design : ")
print("GEOMETRY")
print(f"Ds={Ds*1e3} mm")
print(f"do={do*1e3} mm")
print(f"Dp={Dp*1e3} mm")
print(f"X={X}")
print(f"N_ports={N_p}")
print(f"N_inj={N_inj}")

print("FLOW")
print(f"Cd={Cd}")
print(f"mass_flow_unit={np.round(m_LOX_unit,4)} kg/s")
print(f"dP={np.round(dP_inj/1e5,4)} bar")
print(f"exit_film_thickness={np.round(film_thickness*1e3,4)} mm")
print(f"theta_spray1={spray_half_angle*180/np.pi}° | theta_spray2={theta_2*180/np.pi}°")
print(f"SMD={SMD*1e6} µm")
print(" ")
print("Length guidelines : ")
print("port length = (3-6)*Rp")
print("length swirl chamber > 2(Rs - Rp)")
print("Length nozzle ie. orifice = (0.5-2)*ro")


#%%
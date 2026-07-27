####### TEMPERATURE
"""
@ author : Raphaël Aubry
"""

"""
COMPRESSIBLE IDEAL GAS FLOW
valid for :
    - single-component gas (fluid name/molar mass supplied by the caller via
      coolantProp; no gas is hardcoded here)
    - ideal gas (no real-gas compressibility correction; see main_solve.py's
      Z-correction TODO for the He supercritical deviation, Z ~ 1.04-1.06)
    - quasi-1D
"""


def dT__dx_IdealGas (q_w, P_w, m_dot, U, dU__dx, cp):
    """
    Coolant temperature gradient along channel in counterflow configuration

    UP : conductance according to the total resistance across wall in 1D, 1/R, where R is
            computed with heat transfer coefficients and perimeter

    C_c : coolant capacitance, m_dot * cp
    """

    return (q_w*P_w/m_dot - U * dU__dx)/cp


####### VELOCITY

def dU__dx_IdealGas (U, A, p, q_w, P_w, T, cp, m_dot, dA__dx, f, Dh) :
    """
    Velocity gradient through channel

    """
    num = (- A*p*q_w*P_w/(T*cp*m_dot) + p*dA__dx - (0.5*f/Dh)*m_dot*U)
    den = m_dot - A*p/U - A*p*U/(T*cp)

    return num/den

####### PRESSURE

def dp__dx_IdealGas (p, T, dT__dx, U, dU__dx, A, dA__dx):
    """
    Reurns the pressure gradient along the channel
    """
    return p*(dT__dx/T - dU__dx/U - dA__dx/A)

def dp__dx_IdealGas_logical (p, T, dT__dx, rho, drho__dx):
    """
    Reurns the pressure gradient along the channel
    """
    return p*(dT__dx/T + drho__dx/rho)

####### DENSITY

def drho__dx_IdealGas (rho, p, dp__dx, T, dT__dx) :
    """
    Returns density gradient of coolant along channel
    """
    return rho*(dp__dx/p - dT__dx/T)

def drho__dx_IdealGas_logical (rho, U, dU__dx, A, dA__dx) :
    """
    Returns density gradient of coolant along channel
    """
    return -rho*(dU__dx/U - dA__dx/A)



def dT_g__dx (UP, T_g, T_c, mass_flow_g, cp_g):
    """
    Hot gas temperature evolution assuming no pressure loss, according to total heat convection
    """
    return -(T_g - T_c)*UP/(mass_flow_g*cp_g)

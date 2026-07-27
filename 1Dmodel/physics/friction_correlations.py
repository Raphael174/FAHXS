import numpy as np 
from scipy.optimize import fsolve



#%% EMPIRICAL LAWS FOR FRICTION


def getFrictionColebrook1939 (Re, e__Dh):

    # Colebrook (1939) presents an implicit expression for the fully developed, turbulent 
    # friction factor in a duct with surface roughness, e
    #p 654 Nellis 2009

    def f_x (f):
        return (-2*np.log10(e__Dh/3.7 + 2.51/(Re*np.sqrt(f))))**(-2) - f

    friction_coef = fsolve(func=f_x, x0=0.01)[0]

    return friction_coef 

def getFrictionCurvedPipeAli2024(Re, Dh, Rc, *,
                                  c_lo: float = 0.316, c_hi: float = 0.325,
                                  I_split: float = 0.868):
    """
    Universal skin friction laws for turbulent flow in curved tube.
    https://doi.org/10.1063/5.0222083

    Rc      : pipe radius of curvature
    Dh      : pipe hydraulic diameter
    Re      : Reynolds number
    c_lo    : prefactor for I <= I_split  (weak curvature branch)
    c_hi    : prefactor for I >  I_split  (strong curvature branch, active at design)
    I_split : Dean-group threshold; at design (Re~50k, alpha~0.05) I~3 >> I_split
    """
    alpha = Dh / (2 * Rc)
    I = (Re * alpha**2)**(1/4)

    if I <= I_split:
        return c_lo * I * alpha**(1/2)
    else:
        return c_hi * I**(-4/5) * alpha**(1/2)
    

def friction_factor_curved_tube_Mishra1979 (Re, d, Rc):
    """
    Compute friction factor f from Eq. (14.52):

        f = 0.079 / Re^0.25 + 0.0075 * sqrt(d / (2 Rc))

    Parameters
    ----------
    Re : float or array-like
        Reynolds number.
    d : float or array-like
        Tube (hydraulic) diameter.
    Rc : float or array-like
        Radius of curvature (same length unit as d).

    Returns
    -------
    f : float or np.ndarray
        Darcy friction factor (same shape as numpy-broadcasted inputs).
    """
    Re = np.asarray(Re, dtype=float)
    d = np.asarray(d, dtype=float)
    Rc = np.asarray(Rc, dtype=float)

    if np.any(Re <= 0) or np.any(d <= 0) or np.any(Rc <= 0):
        raise ValueError("Re, d, and Rc must all be > 0.")

    return 0.079 / Re**0.25 + 0.0075 * np.sqrt(d / (2.0 * Rc))




def getFrictionZygSyl1982 (Re, e_Dh):

    #Zigrang and Sylvester (1982) present an explicit and therefore more convenient 
    #correlation for the fully developed, turbulent friction factor in a rough duct:
    #p 654 Nellis 2009
        
    e__Dh_term = e_Dh/3.7

    return (-2*np.log10(e__Dh_term - 5.02/Re*np.log10(e__Dh_term - 5.02/Re*np.log10(e__Dh_term + 13/Re))))**(-2)

def getFrictionSwameeJain1976 (Re, e_Dh):

    """ 
    Used in :
    HELIUM GAS EVACUATION IN SUPERCONDUCTING RFQ STRUCTURE
    A. Lombardi, G. Bisoffi, F. Chiurlotto, E. Tovo, A.M. Porcellato, L. Badan
    1999
    https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://cds.cern.ch/record/553463/files/tua107.pdf&ved=2ahUKEwimktHl0qyLAxWmK_sDHRMTKAsQFnoECBYQAQ&usg=AOvVaw0v0nRCJ82VH2J4nrO1mPi5
    
    Originates from :
    O.G. Swamee Jain, 1976 
    https://doi.org/10.1061/JYCEAJ.0004542
    """

    num = (np.log10(5.74/Re**0.9 + e_Dh/3.7))**2

    return 0.25/num

def getFrictionPetukhov1970 (Re):

    #The friction factor for fully developed turbulent ﬂow in an aerodynamically 
    # smooth duct is provided by Petukhov (1970):
    # 3000 < Re < 5e6
    #p 654 Nellis 2009

    return 1/(0.79*np.log(Re) - 1.64)**2

def Filonenko1954 (Re_D):
    """ 
    Friction factor smooth pipe
    G.K. Filoneko. Hydraulic resistance in pipes. Teploenergetika, 1(4) 40-44, 1954
    """
    return 1/(1.82*np.log10(Re_D) - 1.64)**2

def getFrictionDeveloping (f_fd, Dh, x_duct):

    # entrance effect on turbulent friction factor
    # correction to give developing friction
    #p 655 Nellis 2009

    if x_duct <=0:
        x_duct = 1e-6

    return f_fd*(1 + (Dh/x_duct)**0.7)



#%% HEAT TRANSFER COOLANT 
def getNusseltGnielinski1976 (f_fd, Re, Pr):
    """ 
    The Nusselt number for fully developed turbulent ﬂow is provided by Gnielinski (1976)
    """
    num = f_fd/8 * (Re - 1000)*Pr 
    den = 1 + 12.7*(Pr**(2/3) - 1)* np.sqrt(f_fd/8)
    
    return num/den

def getNusseltDeveloping (Nu_fd, Dh, x_duct):

    """
    The average Nusselt number, [Handbook of Heat Transfer, 1998, p329]:
    """
    if x_duct <=0:
        x_duct = 1e-6
    return Nu_fd*(1 + 0.9756*(Dh/x_duct)**(0.76))


# ---------------------------------------------------------------------------
# Dispatcher — select friction correlation from combustorProp.friction_coil
# ---------------------------------------------------------------------------

def dispatch_friction_coil(selector: str, Re: float, Dh: float, Rc: float,
                            roughness: float, x: float, error_factor: float,
                            corrCoeffs) -> float:
    """Select and evaluate the coil-side friction correlation.

    Returns the DARCY-Weisbach friction factor f_c (fully-developed + entry
    correction).  Pass x=10e10 to obtain the fully-developed limit.

    Convention note: Ali & Dey 2024 (10.1063/5.0222083) call f the "skin friction
    coefficient" but use the Blasius prefactor 0.316 (= Darcy; Fanning Blasius is
    0.079), so the returned f is Darcy.  This is consistent with the momentum
    closures in governing_equations.py (0.5*f/Dh term) and the hot-gas dp form
    (-f*rho*U^2/(2*Dh)).  Do NOT apply a Fanning->Darcy x4.

    error_factor = 1 + numericalProp.artificial_error_friction_cold
    """
    import warnings
    if selector == "CurvedPipeAli2024":
        f_fd = getFrictionCurvedPipeAli2024(
            Re=Re, Dh=Dh, Rc=Rc,
            c_lo=corrCoeffs.ali_c_lo,
            c_hi=corrCoeffs.ali_c_hi,
            I_split=corrCoeffs.ali_I_split,
        )
    else:
        if selector not in ("Colebrook1939",):
            warnings.warn(f"friction_coil selector '{selector}' not recognised — falling back to Colebrook1939")
        f_fd = getFrictionColebrook1939(Re, roughness / Dh)
    return getFrictionDeveloping(f_fd * error_factor, Dh, x)


def dispatch_friction_tube_straight(Re, roughness, Dh, x=10e10,
                                    Re_lo=2300.0, Re_hi=4000.0, error_factor=1.0):
    """
    Tube-side (hot combustion gas) DARCY friction factor for a straight circular
    tube, blended laminar/transitional/turbulent to match the tube-side Nusselt
    (dispatch_nu_tube_straight). See DESIGN_PLAN_shellntube_transient.md section 2.1.

      Re < Re_lo : laminar Hagen-Poiseuille  f = 64/Re
      Re > Re_hi : Colebrook (rough turbulent)
      between    : linear blend (same gamma as the Nu blend — no discontinuity)
    """
    def _lam(Re_):
        return 64.0 / max(Re_, 1e-6)

    if Re <= Re_lo:
        f_fd = _lam(Re)
    elif Re >= Re_hi:
        f_fd = getFrictionColebrook1939(Re, roughness / Dh)
    else:
        gamma = (Re - Re_lo) / (Re_hi - Re_lo)
        f_fd = (1 - gamma) * _lam(Re_lo) + gamma * getFrictionColebrook1939(Re_hi, roughness / Dh)
    return getFrictionDeveloping(f_fd * error_factor, Dh, x)


def friction_corrugated_tube_vicente(Re, phi, Re_lo=2000.0, Re_hi=4000.0,
                                     error_factor=1.0):
    """Darcy friction factor for helically corrugated tubes.

    Uses the Vicente et al. helical-corrugated-tube form summarized by Cruz et
    al. (2021), with severity index phi = e^2/(p*D_i). The published equations
    are Darcy-factor forms here:

      laminar: f = 119.6 * phi^0.11 * Re^-0.97
      low-Re turbulent: f = 6.12 * phi^0.46 * Re^-0.16

    The turbulent branch was reported for soft corrugations and roughly
    2000 < Re < 8000. This function blends through the transition band rather
    than introducing a discontinuity in the transient march.
    """
    Re = max(float(Re), 1e-6)
    phi = max(float(phi), 1e-12)

    def _lam(Re_):
        return 119.6 * phi ** 0.11 * max(Re_, 1e-6) ** (-0.97)

    def _turb(Re_):
        return 6.12 * phi ** 0.46 * max(Re_, 1e-6) ** (-0.16)

    if Re <= Re_lo:
        f_fd = _lam(Re)
    elif Re >= Re_hi:
        f_fd = _turb(Re)
    else:
        gamma = (Re - Re_lo) / (Re_hi - Re_lo)
        f_fd = (1.0 - gamma) * _lam(Re_lo) + gamma * _turb(Re_hi)
    return float(f_fd) * error_factor

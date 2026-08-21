""" 
@author : Raphaël Aubry
"""
#%%
import numpy as np 
import math
from typing import Optional
from scipy.optimize import fsolve

####### COOLANT

def getNusseltDeveloping (Nu_fd, Dh, x_duct):

    """ 
    The average Nusselt number, [Handbook of Heat Transfer, 1998, p329]:
    """
    if x_duct <=0:
        x_duct = 1e-6
    return Nu_fd*(1 + 0.9756*(Dh/x_duct)**(0.76))

def compute_Nusselt_Gnielinski (f, Re_c, Pr_c):
    """ 
    Returns Nusselt number of fluid according to turbulent pipe flow, equation of Gnielinski
    2300 <=Re_D<=5e6 
    0.6<=Pr<=1e5
    accuracy = +-6% to +-10%
    """
    num = f/8 * (Re_c - 1000)*Pr_c
    den = 1 + 12.7*np.sqrt(f/8) * (Pr_c**(2/3) - 1)
    
    return num/den 



def shell_side_nusselt_geometryterm(
    D_c,       # coil diameter [m]
    d_v,       # shell inlet/outlet diameter
    D_cc,      # shell diameter [m]
    L_c,       # coil height [m]
    L_cc,      # shell height [m]
    f,         # distance from inlet to outlet of shell
    coil_pitch,         # pitch [m]
    d_coil_outer,      # outer tube diameter [m]
    # Re_sh=1,     # Reynolds number in shell
    # Pr_sh=1,     # Prandtl number in shell
):
    """
    Calculate shell-side Nusselt number for a shell-and-tube heat exchanger.

    Alimoradi2016

    Parameters:
        Re_sh : float - Reynolds number in the shell
        Pr_sh : float - Prandtl number in the shell
        D_c   : float - coil diameter
        d_v   : float - some flow diameter (possibly equivalent or bypass gap)
        D_cc  : float - shell diameter
        L_c   : float - coil height
        L_cc  : float - shell height
        f     : float - distance, possibly fin length or baffle spacing
        p     : float - pitch
        d_coil_outer : float - outer diameter of the tube

    Returns:
        Nu_t_sh : float - Nusselt number on the shell side

        https://doi.org/10.1016/j.ijthermalsci.2016.04.010 
    """

    
    Nu_t_sh_geometryterm = 0.247  \
        * (D_c / d_coil_outer)**0.378 \
        * (d_v / d_coil_outer)**0.556 \
        * (D_cc / d_coil_outer)**(-0.82) \
        * (L_c / d_coil_outer)**0.043 \
        * (L_cc / d_coil_outer)**(-1.03) \
        * (f / d_coil_outer)**0.561 \
        * (coil_pitch / d_coil_outer)**0.138 
    
    return Nu_t_sh_geometryterm

def shell_side_nusselt(
    geometry_term,
    Re_sh,     # Reynolds number in shell
    Pr_sh,     # Prandtl number in shell
    ):

    """
    Alimoradi2016
    https://doi.org/10.1016/j.ijthermalsci.2016.04.010 
    """

    if Re_sh > 50e3:
        print("Reynolds of shell combustor above 49e3 limit for correlation, by ", (Re_sh-49e3)/49e3*100, "%")
    if Pr_sh > 7.1:
        print("Prandlt of shell combustor above 7.1 limit for correlation, by ", (Pr_sh-7.1)/7.1*100, "%")


    Nu_t_sh = geometry_term * Re_sh**0.723 * Pr_sh**0.717

    return Nu_t_sh


def nu_churchill_bernstein(Re: float, Pr: float, Pr_s: Optional[float] = None) -> float:
    """
    Churchill–Bernstein Nusselt number for a circular cylinder in crossflow.

    Nu_D = 0.3
           + [0.62 * Re^(1/2) * Pr^(1/3)] / [1 + (0.4/Pr)^(2/3)]^(1/4)
             * [1 + (Re/282000)^(5/8)]^(4/5)
           * (Pr/Pr_s)^(1/4)   # optional wall-property correction

    Parameters
    ----------
    Re : Reynolds number based on cylinder diameter (use film properties).
    Pr : Prandtl number at film temperature.
    Pr_s : (optional) Prandtl number at the wall temperature; if given,
           applies the (Pr/Pr_s)**0.25 correction.

    Returns
    -------
    Nu : Nusselt number based on cylinder diameter.

    Notes
    -----
    Valid for ~0.2 ≲ Pr ≲ O(10^3) and wide Re range in external crossflow.
    For gases, keep Mach ≲ 0.35 unless you add compressibility corrections.
    """
    # guard against zero/negative inputs
    Re = max(Re, 1e-16)
    Pr = max(Pr, 1e-16)

    term1 = 0.3
    term2 = (0.62 * math.sqrt(Re) * (Pr ** (1.0/3.0))) / ((1.0 + (0.4/Pr) ** (2.0/3.0)) ** 0.25)
    term3 = (1.0 + (Re / 282000.0) ** (5.0/8.0)) ** (4.0/5.0)

    Nu = term1 + term2 * term3

    if Pr_s is not None:
        Pr_s = max(Pr_s, 1e-16)
        Nu *= (Pr / Pr_s) ** 0.25

    return Nu


def phi_sl_multiplier(SL_over_D: float,
                      ref_SL_over_D: float = 3.0,
                      cap: Optional[float] = 1.4,
                      enforce_range: bool = True) -> float:
    """
    KCY-derived streamwise-spacing multiplier φ_SL(S_L/D) to apply on top of a
    single-cylinder baseline h for consecutive coil turns (wide S_T/D).

    φ_SL(S_L) = ([0.25 + exp(-0.55 S_L)] * S_L^0.212) /
                ([0.25 + exp(-0.55 S_ref)] * S_ref^0.212)

    • Validated range (from KCY in-line bank fit): 1.05 ≤ S_L/D ≤ 3.
    • For S_L/D ≥ S_ref (default 3), returns 1.0 (no added enhancement).
    • Optional 'cap' (default 1.4) limits enhancement pending CFD/experiments.
    • Set enforce_range=False to evaluate outside KCY’s range (not recommended).

    Parameters
    ----------
    SL_over_D : float
        Streamwise pitch ratio S_L/D (center-to-center spacing between turns / tube OD).
    ref_SL_over_D : float, optional
        Reference spacing where φ_SL is defined to be 1.0 (default 3.0).
    cap : float or None, optional
        Maximum allowed multiplier; set None to disable capping (default 1.4).
    enforce_range : bool, optional
        If True, S_L/D is clamped to [1.05, ref_SL_over_D]. (default True)

    Returns
    -------
    φ_SL : float
        Multiplier to apply to the baseline h (e.g., Churchill–Bernstein at gap Re).

    Notes
    -----
    Source trend from: Khan–Culham–Yovanovich (2006), IJHMT, in-line bank fit.
    Use with single-cylinder baseline when S_T/D ≫ 3; compute Re with maximum
    (gap) velocity to capture blockage.
    """
    # If user asks beyond reference, no enhancement
    if SL_over_D >= ref_SL_over_D:
        return 1.0

    SL = SL_over_D
    if enforce_range:
        SL = max(1.05, min(SL, ref_SL_over_D))

    num = (0.25 + math.exp(-0.55 * SL)) * (SL ** 0.212)
    den = (0.25 + math.exp(-0.55 * ref_SL_over_D)) * (ref_SL_over_D ** 0.212)
    phi = num / den

    if cap is not None:
        phi = min(phi, cap)

    # Guard against numerical underflow
    return max(phi, 1.0)

def nusselt_shell_salimpour2008(
    Re_s: float,
    Pr_s: float,
    coil_pitch: float,
    D_o: float,
    T_bulk: float = None,
    T_wall: float = None,
    *,
    a: float = 0.317,
    b: float = 0.643,
    c: float = -0.215,
    n: float = 0.25,
) -> float:
    """
    Shell-side Nusselt number for a helical-coil-in-cylindrical-shell heat exchanger.

    Source
    ------
    Salimpour, M.R. (2008). "Heat transfer coefficients of shell and coiled tube heat
    exchangers." Experimental Thermal and Fluid Science 33(2): 203–207.
    https://doi.org/10.1016/j.expthermflusci.2008.07.015

    Correlation
    -----------
    Nu_s = 0.317 · Re_s^0.643 · Pr_s^(1/3) · (p / D_o)^(−0.215)

    where:
      Re_s  — Reynolds number based on the shell hydraulic diameter (Dh_shell)
      Pr_s  — bulk Prandtl number of the shell-side fluid
      p     — coil pitch [m] (= coil_gap + D_o, axial centre-to-centre spacing)
      D_o   — outer tube diameter [m] (= Dh_coil + 2·t_wall)

    Characteristic length
    ---------------------
    h_s = Nu_s · k / Dh_shell  (use the same Dh_shell that enters Re_s)

    Applicability
    -------------
    - Multi-turn helical coil inside a cylindrical shell with axial shell-side flow.
    - Turbulent regime: Re_s ≳ 10 000.
    - Original experimental data: water and ethylene-glycol/water mixtures, Pr ≈ 4–15.
    - The correlation has been validated for the specific shell-and-helical-coil geometry
      (unlike Ahmed-1997 which was derived for a toroidal cavity, or Churchill-Bernstein
      which assumes an isolated cylinder in infinite unbounded crossflow).
    - The Dh_shell formula used elsewhere in this code (compute_Dh_shell) originates
      from the same Salimpour paper, ensuring dimensional consistency.

    Property correction for gases
    ------------------------------
    Salimpour's data used liquids (Pr ≫ 1). When applied to combustion gases (Pr ≈ 0.65),
    the Pr^(1/3) term handles the thermal boundary layer thickness difference. The dominant
    remaining effect is the large bulk-to-wall temperature ratio typical in combustors.

    A Sieder-Tate viscosity-ratio correction (μ_b/μ_w)^0.14 is NOT appropriate for gases:
    for gases μ increases with T (Sutherland), so near a cooled wall μ_wall < μ_bulk, but
    the physically dominant effect is the density change (ρ ∝ 1/T), not the viscosity change.

    Instead, a Kays & Crawford temperature-ratio correction is applied when T_bulk and T_wall
    are provided (recommended for combustor applications):
        Nu_corrected = Nu_Salimpour · (T_bulk / T_wall)^0.25
    This correction is non-trivial: at T_bulk ≈ 2000 K and T_wall ≈ 800 K the factor is ~1.26.

    Parameters
    ----------
    Re_s    : Reynolds number (use shell hydraulic diameter Dh_shell)
    Pr_s    : Prandtl number at bulk gas temperature
    coil_pitch : axial centre-to-centre pitch [m]
    D_o     : outer coil tube diameter [m]
    T_bulk  : bulk gas temperature [K]  (optional — enables gas correction)
    T_wall  : gas-side wall temperature [K]  (optional — enables gas correction)

    Returns
    -------
    Nu_s : float — shell-side Nusselt number
    """
    Nu = a * Re_s**b * Pr_s**(1.0/3.0) * (coil_pitch / D_o)**c
    if T_bulk is not None and T_wall is not None and T_wall > 0:
        Nu *= (T_bulk / T_wall)**n
    return Nu


def nusselt_toroid_Ahmed1997 (U_g, rho_g, mu_g, Pr_g, Asqrt_toroid):

    
    def f_x(gamma):
        Re_gamA = gamma*Asqrt_toroid*U_g*rho_g/mu_g
        return 1/(1 + 0.49*Re_gamA**1.25)**0.2 - gamma

    gamma_A = fsolve(func=f_x, x0=0.2, xtol=1e-8)[0]
    func_Pr_gamA = Pr_g**(1/3) / ((2*gamma_A + 1)**3 + 1/Pr_g)**(1/6)

    Re_A_sqrt = Asqrt_toroid*U_g*rho_g/mu_g
    return 3.41 + 1.58*Re_A_sqrt**0.5*func_Pr_gamA
    

def nusselt_inner_curved_tube_mori1967(Re, Pr, d, R, *,
                                        pr_switch: float = 1.0,
                                        a_lo: float = 26.2, b_lo: float = 0.074, c_lo: float = 0.098,
                                        a_hi: float = 41.0, c_hi: float = 0.061):
    """
    Compute Nusselt number Nu from the correlation shown in the screenshot.

    Parameters
    ----------
    Re : float or array-like
        Reynolds number.
    Pr : float or array-like
        Prandtl number.
    d : float or array-like
        Tube (hydraulic) diameter.
    R : float or array-like
        Radius of curvature (same length unit as d).
    pr_switch : float, optional
        Switch between the two expressions:
        - use Eq. (14.53a) when Pr <= pr_switch (Pr ≈ 1 region)
        - use Eq. (14.53b) when Pr >  pr_switch

    Returns
    -------
    Nu : float or np.ndarray
        Nusselt number (same shape as numpy-broadcasted inputs).
    """
    Re = np.asarray(Re, dtype=float)
    Pr = np.asarray(Pr, dtype=float)
    d = np.asarray(d, dtype=float)
    R = np.asarray(R, dtype=float)

    if np.any(Re <= 0) or np.any(Pr <= 0) or np.any(d <= 0) or np.any(R <= 0):
        raise ValueError("Re, Pr, d, and R must all be > 0.")

    dr = d / (2.0 * R)  # (d / 2R)

    # Eq. (14.53a): for Pr ≈ 1  [Mori & Nakayama 1967]
    Nu_a = (
        (Pr / (a_lo * (Pr ** (2.0 / 3.0) - b_lo)))
        * Re ** (4.0 / 5.0)
        * dr ** (1.0 / 10.0)
        * (1.0 + c_lo * (Re * dr ** 2.0)) ** (1.0 / 5.0)
    )

    # Eq. (14.53b): for Pr > 1  [Mori & Nakayama 1967]
    Nu_b = (
        (Pr ** 0.4 / a_hi)
        * Re ** (5.0 / 6.0)
        * dr ** (1.0 / 12.0)
        * (1.0 + c_hi * (Re * dr ** 2.5)) ** (1.0 / 6.0)
    )

    use_a = Pr <= pr_switch
    return np.where(use_a, Nu_a, Nu_b)


# Example:
# Nu = nusselt_curved_tube(Re=2e4, Pr=0.9, d=0.01, R=0.05)
# Nu = nusselt_curved_tube(Re=2e4, Pr=3.0, d=0.01, R=0.05)




"""
Natural convection around a horizontal cylinder — water side (film-property version)

This module provides:
  1) nu_horizontal_cylinder_filmprops(): Churchill–Chu Nusselt using *film* properties and ΔT only.
  2) get_water_film_properties_coolprop(): convenience helper to obtain film properties from CoolProp.

Notes
-----
• Evaluate properties at film temperature: T_film = 0.5*(T_bulk + T_surface).  You still need ΔT = |T_bulk − T_surface|.
• Valid for single-phase, quiescent water (no boiling at the tube). Radiation in water is negligible.
• Correlation: Churchill–Chu (1975) for a horizontal cylinder; wide Ra coverage (≈ 1e-10 to 1e12).
"""

from typing import Dict, Optional

try:
    from CoolProp.CoolProp import PropsSI  # type: ignore
    _HAS_COOLPROP = True
except Exception:  # CoolProp not available
    PropsSI = None
    _HAS_COOLPROP = False


def Nusselt_horizontal_cylinder_filmprops(
    deltaT: float,     # [K] absolute temperature difference |T_bulk − T_surface|
    D: float,          # [m] cylinder/tube outside diameter
    k: float,          # [W/m/K] thermal conductivity at T_film
    nu: float,         # [m^2/s] kinematic viscosity ν at T_film
    alpha: float,      # [m^2/s] thermal diffusivity α at T_film
    Pr: float,         # [-] Prandtl number at T_film
    beta: Optional[float] = None,  # [1/K] isobaric thermal expansion coefficient at T_film; if None uses 1/T_film (approx)
    T_film: Optional[float] = None,  # [K] only needed if beta is None (to use beta ≈ 1/T_film)
    g: float = 9.80665  # [m/s^2]
) -> Dict[str, float]:
    """
    Churchill–Chu Nusselt for a horizontal cylinder using *film* properties.

    Returns dict with Nu, h, Ra.
    """
    if D <= 0:
        raise ValueError("D must be > 0")
    if deltaT <= 0:
        raise ValueError("deltaT must be > 0 (K)")
    if min(k, nu, alpha, Pr) <= 0:
        raise ValueError("k, nu, alpha, Pr must be > 0")

    if beta is None:
        if T_film is None or T_film <= 0:
            raise ValueError("Provide beta or T_film (>0 K) to use beta ≈ 1/T_film")
        beta_used = 1.0 / T_film  # quick approximation for liquids; prefer tabulated β
    else:
        beta_used = beta

    # Rayleigh number (based on D)
    Ra = g * beta_used * deltaT * (D**3) / (nu * alpha)

    # Churchill–Chu (horizontal cylinder)
    term = 0.387 * (Ra ** (1.0/6.0)) / ((1.0 + (0.559/Pr) ** (9.0/16.0)) ** (8.0/27.0))
    Nu = (0.60 + term) ** 2

    h = Nu * k / D

    return {"Nu": Nu, "h": h, "Ra": Ra, "beta_used": beta_used}


# ---------------------------------------------------------------------------
# Dispatchers — select Nu correlation from combustorProp.Nusselt_coil / Nusselt_shell
# ---------------------------------------------------------------------------

def dispatch_nu_coil(selector: str, Re: float, Pr: float, d: float, R: float,
                     f_fd: float, x: float, error_factor: float, corrCoeffs) -> float:
    """Select and evaluate the coil-side (helium) Nusselt correlation.

    Returns Nu_c (dimensionless).
    f_fd : fully-developed Fanning friction factor (needed for Gnielinski).
    x    : axial duct position for entry-length correction; pass 10e10 for fully developed.
    error_factor = 1 + numericalProp.artificial_error_Nu_cold
    """
    import warnings
    if selector == "mori1967":
        Nu_fd = nusselt_inner_curved_tube_mori1967(
            Re=Re, Pr=Pr, d=d, R=R,
            a_lo=corrCoeffs.mori_a_lo,
            b_lo=corrCoeffs.mori_b_lo,
            c_lo=corrCoeffs.mori_c_lo,
        )
        return float(getNusseltDeveloping(Nu_fd, d, x)) * error_factor
    else:
        if selector not in ("Gnielinski",):
            warnings.warn(f"Nusselt_coil selector '{selector}' not recognised — falling back to Gnielinski")
        Nu_fd = compute_Nusselt_Gnielinski(f=f_fd, Re_c=Re, Pr_c=Pr)
        return float(getNusseltDeveloping(Nu_fd, d, x)) * error_factor


def nu_tube_laminar_entrance(Re, Pr, d, x):
    """
    Laminar Nusselt in a circular tube with combined entrance region, constant-Tw
    (VDI Heat Atlas G1 / Gnielinski laminar composite):

        Nu = [ 3.66^3 + 0.7^3 + (1.615*(Re*Pr*d/L)^(1/3) - 0.7)^3 ]^(1/3)

    At L/D~67, Re*Pr*d/L ~ 20-25 -> Nu ~ 5-6, i.e. +40-60% over the fully-developed
    3.66 (see DESIGN_PLAN_shellntube_transient.md section 2.1). `x` is the local
    axial position [m]; a small floor avoids the singularity at x->0.
    """
    L = max(x, 1e-4)
    gz = Re * Pr * d / L                      # inverse Graetz-type group
    return (3.66 ** 3 + 0.7 ** 3 + (1.615 * gz ** (1.0 / 3.0) - 0.7) ** 3) ** (1.0 / 3.0)


def dispatch_nu_tube_straight(selector, Re, Pr, d, x, f_fd,
                              T_bulk=None, T_wall=None,
                              error_factor=1.0, corrCoeffs=None):
    """
    Tube-side (hot combustion gas) Nusselt for a straight circular tube, blended
    smoothly across laminar / transitional / turbulent regimes.

    Regimes (section 2.1):
      Re < Re_lo (2300)      : laminar entrance composite (nu_tube_laminar_entrance)
      Re > Re_hi (4000)      : Gnielinski turbulent (compute_Nusselt_Gnielinski)
      Re_lo <= Re <= Re_hi   : linear blend (Gnielinski 2013) — no discontinuity as
                               the gas cools and its Re marches down through the band.

    Variable-property correction for a gas being COOLED: Nu *= (T_wall/T_bulk)^n,
    n = corrCoeffs.n_tube_gas (default 0.0 — the hot-end T_b/T_w~2.5 is far outside
    correlation databases, so n is a prime calibration knob, like Nusselt_correction).
    """
    import warnings
    Re_lo = getattr(corrCoeffs, "Re_transition_lo", 2300.0) if corrCoeffs else 2300.0
    Re_hi = getattr(corrCoeffs, "Re_transition_hi", 4000.0) if corrCoeffs else 4000.0
    n = getattr(corrCoeffs, "n_tube_gas", 0.0) if corrCoeffs else 0.0

    if selector not in ("gnielinski_blended",):
        warnings.warn(f"Nusselt_tube selector '{selector}' not recognised — using gnielinski_blended")

    def _turb(Re_):
        f = f_fd if f_fd is not None else (0.790 * np.log(Re_) - 1.64) ** (-2)  # Konakov
        return compute_Nusselt_Gnielinski(f=f, Re_c=Re_, Pr_c=Pr)

    if Re <= Re_lo:
        Nu = nu_tube_laminar_entrance(Re, Pr, d, x)
    elif Re >= Re_hi:
        Nu = _turb(Re)
    else:  # transitional linear blend
        gamma = (Re - Re_lo) / (Re_hi - Re_lo)
        Nu = (1 - gamma) * nu_tube_laminar_entrance(Re_lo, Pr, d, x) + gamma * _turb(Re_hi)

    if T_bulk is not None and T_wall is not None and T_bulk > 0 and n != 0.0:
        Nu *= (T_wall / T_bulk) ** n
    return float(Nu) * error_factor


def nu_corrugated_tube_vicente(Re, Pr, phi, D_i, x, D_h=None,
                               Re_lo=2000.0, Re_hi=4000.0,
                               error_factor=1.0):
    """Nusselt number for helically corrugated tubes.

    Uses the Vicente et al. turbulent correlation summarized by Cruz et al.
    (2021):

      Nu_tilde = 0.3741 * phi^0.25 * (Re_tilde - 1500)^0.74 * Pr^0.44

    where Re_tilde and Nu_tilde are based on the smooth inner diameter D_i.
    The solver uses D_i as its current hydraulic diameter for the representative
    tube; D_h is accepted so a future exact corrugated hydraulic diameter can be
    converted with Nu = Nu_tilde * D_h/D_i.

    For Re below the turbulent correlation's lower limit, this falls back to the
    existing laminar entrance model and blends through the transition band.
    """
    Re = max(float(Re), 1e-6)
    Pr = max(float(Pr), 1e-12)
    phi = max(float(phi), 1e-12)
    D_i = max(float(D_i), 1e-12)
    D_h = D_i if D_h is None else max(float(D_h), 1e-12)

    def _lam(Re_):
        return nu_tube_laminar_entrance(Re_, Pr, D_i, x)

    def _turb(Re_):
        Re_tilde = (D_i / D_h) * max(Re_, 1500.0 + 1e-6)
        Nu_tilde = 0.3741 * phi ** 0.25 * (Re_tilde - 1500.0) ** 0.74 * Pr ** 0.44
        return Nu_tilde * (D_h / D_i)

    if Re <= Re_lo:
        Nu = _lam(Re)
    elif Re >= Re_hi:
        Nu = _turb(Re)
    else:
        gamma = (Re - Re_lo) / (Re_hi - Re_lo)
        Nu = (1.0 - gamma) * _lam(Re_lo) + gamma * _turb(Re_hi)
    return float(Nu) * error_factor


def dispatch_nu_shell(
    selector: str,
    Re_sh: float, Re_g: float, Pr_g: float, k_g: float,
    U_g: float, rho_g: float, mu_g: float,
    coil_pitch: float, Dh_cc: float,
    Dh_ch: float, D_coil: float, thickness_wall: float,
    T_bulk: float, T_wall: float,
    Nusselt_correction: float, error_factor: float,
    corrCoeffs,
) -> tuple:
    """Select and evaluate the shell-side (hot gas) Nusselt correlation.

    Returns (Nu_g, h_g [W/m²K]).
    Each branch uses a different characteristic length for h_g — the dispatcher
    handles this internally so the caller always receives the correct h_g.

    error_factor = 1 + numericalProp.artificial_error_Nu_hot
    """
    import warnings
    D_tube_outer = Dh_ch + 2.0 * thickness_wall

    if selector == "salimpour2008":
        Nu_g = nusselt_shell_salimpour2008(
            Re_s=Re_g, Pr_s=Pr_g,
            coil_pitch=coil_pitch, D_o=D_tube_outer,
            T_bulk=T_bulk, T_wall=T_wall,
            a=corrCoeffs.salimpour_a, b=corrCoeffs.salimpour_b,
            c=corrCoeffs.salimpour_c, n=corrCoeffs.kays_crawford_n,
        ) * Nusselt_correction * error_factor
        return float(Nu_g), float(Nu_g) * k_g / Dh_cc

    elif selector == "churchill_bernstein_tightcoil":
        phi = phi_sl_multiplier(SL_over_D=coil_pitch / D_tube_outer)
        Nu_g = nu_churchill_bernstein(Re=Re_sh, Pr=Pr_g) * phi * Nusselt_correction * error_factor
        return float(Nu_g), float(Nu_g) * k_g / D_tube_outer

    elif selector == "churchill_bernstein":
        Nu_g = nu_churchill_bernstein(Re=Re_sh, Pr=Pr_g) * error_factor
        return float(Nu_g), float(Nu_g) * k_g / D_tube_outer

    else:
        if selector not in ("ahmed_toroid",):
            warnings.warn(f"Nusselt_shell selector '{selector}' not recognised — falling back to ahmed_toroid")
        Asqrt = np.sqrt(D_tube_outer * D_coil * np.pi**2)
        Nu_g = nusselt_toroid_Ahmed1997(
            U_g=U_g, rho_g=rho_g, mu_g=mu_g, Pr_g=Pr_g, Asqrt_toroid=Asqrt,
        ) * error_factor
        return float(Nu_g), float(Nu_g) * k_g / Asqrt


def get_water_film_properties_coolprop(T_film: float, P: float = 101325.0) -> Dict[str, float]:
    """
    Convenience helper to obtain *film* properties of water from CoolProp at (T_film, P).
    Returns a dict with: rho, mu, k, cp, Pr, nu, alpha, beta.

    Requires CoolProp. If unavailable, raises RuntimeError.
    """
    if not _HAS_COOLPROP:
        raise RuntimeError("CoolProp not available. Install with `pip install CoolProp`. "
                           "Docs: https://coolprop.org/")
    if T_film <= 0 or P <= 0:
        raise ValueError("T_film and P must be > 0")

    rho = PropsSI('D', 'T', T_film, 'P', P, 'Water')                 # [kg/m^3]
    mu  = PropsSI('VISCOSITY', 'T', T_film, 'P', P, 'Water')         # [Pa·s]
    k   = PropsSI('CONDUCTIVITY','T', T_film, 'P', P, 'Water')       # [W/m/K]
    cp  = PropsSI('CPMASS', 'T', T_film, 'P', P, 'Water')            # [J/kg/K]
    # CoolProp provides Pr directly; otherwise compute from mu, cp, k
    try:
        Pr = PropsSI('Prandtl', 'T', T_film, 'P', P, 'Water')        # [-]
    except Exception:
        Pr = mu * cp / k
    # Isobaric expansion coefficient β
    try:
        beta = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','T',T_film,'P',P,'Water')  # [1/K]
    except Exception:
        beta = 1.0 / T_film  # last-resort approximation

    nu = mu / rho                                              # [m^2/s]
    alpha = k / (rho * cp)                                     # [m^2/s]

    return {"rho": rho, "mu": mu, "k": k, "cp": cp, "Pr": Pr, "nu": nu, "alpha": alpha, "beta": beta}


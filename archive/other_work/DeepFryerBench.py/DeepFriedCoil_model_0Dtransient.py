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


def nu_horizontal_cylinder_filmprops(
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


if __name__ == "__main__":
    # Example usage
    Tf = 358.0  # K (~85 °C)
    props = {}
    if _HAS_COOLPROP:
        props = get_water_film_properties_coolprop(Tf, 101325.0)
        vals = nu_horizontal_cylinder_filmprops(
            deltaT=20.0,           # K, example |T_bulk − T_surface|
            D=0.006,               # m
            k=props['k'],
            nu=props['nu'],
            alpha=props['alpha'],
            Pr=props['Pr'],
            beta=props['beta'],
        )
        print({**vals, **{k: props[k] for k in ('nu','alpha','Pr')}})
    else:
        print("CoolProp not installed; only the correlation function is available.")


# =============================
# Transient bath–coil simulator
# =============================

from dataclasses import dataclass
import math
from typing import Callable, List, Tuple, Optional

# -- Optional CoolProp-backed property helpers (water + helium + generic incompressible bath)

def _water_film_props(T_film: float, P: float = 101325.0):
    d = get_water_film_properties_coolprop(T_film, P)
    return d['k'], d['nu'], d['alpha'], d['Pr'], d['beta']

def _incomp_film_props(fluid: str, T_film: float, P: float = 101325.0):
    """Generic incompressible bath property helper using CoolProp INCOMP backend.
    Returns (k, nu, alpha, Pr, beta) at (T_film, P).
    Fluid examples: 'INCOMP::Water', 'INCOMP::Therminol66'.
    """
    if not _HAS_COOLPROP:
        raise RuntimeError("CoolProp not available. Install with `pip install CoolProp`.")
    from CoolProp.CoolProp import PropsSI
    rho = PropsSI('D','T',T_film,'P',P,fluid)
    mu  = PropsSI('VISCOSITY','T',T_film,'P',P,fluid)
    k   = PropsSI('CONDUCTIVITY','T',T_film,'P',P,fluid)
    cp  = PropsSI('CPMASS','T',T_film,'P',P,fluid)
    try:
        Pr  = PropsSI('Prandtl','T',T_film,'P',P,fluid)
    except Exception:
        Pr = mu*cp/k
    try:
        beta = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','T',T_film,'P',P,fluid)
    except Exception:
        beta = 1.0/T_film
    nu = mu/rho
    alpha = k/(rho*cp)
    return k, nu, alpha, Pr, beta

def make_tabulated_bath_props(T_points: List[float], rho: List[float], mu: List[float], k: List[float], cp: List[float], beta: Optional[List[float]] = None):
    """Create a bath film-property provider from tabulated data (e.g., frying oil).
    Linear interpolation in temperature; returns (k, nu, alpha, Pr, beta).
    Units: T [K]; rho [kg/m^3]; mu [Pa·s]; k [W/m/K]; cp [J/kg/K]; beta [1/K] optional.
    """
    import bisect
    if not (len(T_points)==len(rho)==len(mu)==len(k)==len(cp)):
        raise ValueError("All arrays must have same length")
    if beta is not None and len(beta)!=len(T_points):
        raise ValueError("beta array must match length of T_points")
    T = T_points
    def _interp(arr, Tq):
        i = max(0, min(len(T)-2, bisect.bisect_left(T, Tq)-1))
        x0,x1 = T[i], T[i+1]
        y0,y1 = arr[i], arr[i+1]
        w = 0.0 if x1==x0 else (Tq - x0)/(x1 - x0)
        return y0 + w*(y1 - y0)
    def provider(T_film: float, P_unused: float = 101325.0):
        _rho = _interp(rho, T_film)
        _mu  = _interp(mu,  T_film)
        _k   = _interp(k,   T_film)
        _cp  = _interp(cp,  T_film)
        _Pr  = _mu*_cp/_k
        _nu  = _mu/_rho
        _alpha = _k/(_rho*_cp)
        _beta = _interp(beta, T_film) if beta is not None else 1.0/T_film
        return _k, _nu, _alpha, _Pr, _beta
    return provider

def _bath_film_props(bath, T_film: float):
    if getattr(bath, 'props_provider', None) is not None:
        return bath.props_provider(T_film, bath.P)
    if getattr(bath, 'fluid_name', None):
        name = bath.fluid_name
        if name.lower() in ('water','incomp::water'):
            return _water_film_props(T_film, bath.P)
        return _incomp_film_props(name, T_film, bath.P)
    # Default to water if nothing specified
    return _water_film_props(T_film, bath.P)


def _helium_bulk_props(T: float, P: float = 101325.0):
    """Return (rho, mu, k, cp, Pr, nu) for Helium at (T, P). If CoolProp
    is not installed, use ideal-gas-ish constants as a fallback.
    """
    if _HAS_COOLPROP:
        from CoolProp.CoolProp import PropsSI
        rho = PropsSI('D','T',T,'P',P,'Helium')
        mu  = PropsSI('VISCOSITY','T',T,'P',P,'Helium')
        k   = PropsSI('CONDUCTIVITY','T',T,'P',P,'Helium')
        cp  = PropsSI('CPMASS','T',T,'P',P,'Helium')
        try:
            Pr = PropsSI('Prandtl','T',T,'P',P,'Helium')
        except Exception:
            Pr = mu*cp/k
    else:
        # Rough constants (STP-order) — for preliminary runs only
        rho = P/(2077.0*T)             # ideal gas with R=2077 J/kg/K
        mu  = 1.96e-5                  # Pa·s
        k   = 0.15                     # W/m/K
        cp  = 5193.0                   # J/kg/K
        Pr  = mu*cp/k
    nu = mu/rho
    return rho, mu, k, cp, Pr, nu


# -- Internal convection (helium) — default: Gnielinski with optional helical correction

def hi_gnielinski_helical(
    m_dot: float,      # kg/s
    D_i: float,        # m (tube inner diameter)
    R_coil: float,     # m (coil centerline radius); set to math.inf for straight tube
    T_g: float,        # K (bulk helium T)
    P_g: float = 101325.0,
    roughness: float = 0.0  # m (Darcy friction not using roughness explicitly here; kept for future)
) -> Tuple[float, float]:
    """Return (h_i, Nu_i) for turbulent helium using Gnielinski on a smooth tube
    plus a simple curvature multiplier for helical coils. Use with care for
    Dean-number effects; user may replace this with their in-house correlation.
    """
    rho, mu, k, cp, Pr, nu = _helium_bulk_props(T_g, P_g)
    # Flow quantities
    A = 0.25*math.pi*D_i**2
    v = m_dot/(rho*A)
    Re = v*D_i/nu
    Re = max(Re, 1.0)
    # Darcy friction factor (smooth-tube Petukhov fit)
    if Re < 3000:
        # Transitional; revert to Dittus-Boelter style to avoid singularities
        f = 64.0/max(Re,1.0)
    else:
        f = (0.79*math.log(Re)-1.64)**-2
    Nu_gn = (f/8.0)*(Re-1000.0)*Pr/(1.0+12.7*math.sqrt(f/8.0)*(Pr**(2.0/3.0)-1.0))
    # Helical curvature correction (Mori–Nakayama style multiplier)
    if math.isfinite(R_coil) and R_coil > 0.0:
        delta = D_i/(2.0*R_coil)       # curvature ratio
        phi = 1.0 + 3.5*delta + 0.7*delta**2
    else:
        phi = 1.0
    Nu = Nu_gn*phi
    h_i = Nu*k/D_i
    return h_i, Nu


@dataclass
class CoilGeom:
    D_o: float      # m outer diameter
    D_i: float      # m inner diameter
    k_wall: float   # W/m/K tube wall conductivity
    R_coil: float   # m coil centerline radius
    N_turn: int = 8 # number of turns

    @property
    def length(self) -> float:
        """Total wetted length in bath, assuming N_turn full circles at R_coil."""
        return 2.0 * math.pi * self.R_coil * max(self.N_turn, 0)




@dataclass
class Bath:
    m_fluid: float = 30.0      # kg total bath mass (water or oil)
    T0: float = 293.0          # K initial bath temperature
    P: float = 101325.0        # Pa bath pressure (assume ~1 atm)
    UA_loss: float = 0.0       # W/K parasitic loss to ambient (set >0 if uninsulated)
    T_amb: float = 293.0       # K ambient temp for loss calc
    fluid_name: Optional[str] = 'Water'  # default to water properties
    props_provider: Optional[Callable[[float,float],Tuple[float,float,float,float,float]]] = None
    cp_const: Optional[float] = None  # If set, use this constant cp for bath [J/kg/K]


@dataclass
class Heater:
    mode: str = "off"            # 'off' | 'fixed' | 'thermostat'
    P_max: float = 0.0           # [W] nameplate power
    eta: float = 0.98            # [-] electrical→thermal efficiency (immersion ≈ 1)
    P_fixed: float = 0.0         # [W] used if mode == 'fixed'
    T_set: float = 368.15        # [K] thermostat setpoint (e.g., 95 °C)
    deadband: float = 2.0        # [K] ON at T_set - DB/2, OFF at T_set + DB/2


@dataclass
class HeliumInlet:
    m_dot: float = 0.015    # kg/s (15 g/s)
    T_in: float = 293.0     # K
    P: float = 30e5         # Pa (30 bar)


def _solve_segment(
    T_w: float,
    T_g: float,
    p_g: float,
    geom: CoilGeom,
    bath: Bath,
    m_dot: float,
    hi_func: Callable[[float,float,float,float,float], Tuple[float,float]],
    max_iter: int = 12,
    tol: float = 1e-3,
) -> dict:
    """Solve one axial segment (per-unit-length) thermal network consistently to get
    q' (W/m), local wall temps, h_ext, h_i, and helium properties at (T_g, p_g).

    Iterates because h_ext depends on T_wall,o via film properties and ΔT.
    """
    r_o = 0.5*geom.D_o
    r_i = 0.5*geom.D_i

    # Helium properties and internal h_i (based on gas bulk at segment inlet)
    rho_h, mu_h, k_h, cp_h, Pr_h, nu_h = _helium_bulk_props(T_g, p_g)
    h_i, Nu_i = hi_func(m_dot, geom.D_i, geom.R_coil, T_g, p_g)

    # Initial guess for outer wall temperature
    Tso = T_w - 0.2*abs(T_w - T_g)
    Tso = min(T_w-1e-6, max(T_g+1e-6, Tso))  # keep between fluid temps

    h_ext = None
    Nu_ext = None
    Ra_ext = None
    Tf_w = None

    for _ in range(max_iter):
        Tf_w = 0.5*(T_w + Tso)
        k_w, nu_w, alpha_w, Pr_w, beta_w = _bath_film_props(bath, Tf_w)
        vals = nu_horizontal_cylinder_filmprops(
            deltaT=abs(T_w - Tso),
            D=geom.D_o,
            k=k_w,
            nu=nu_w,
            alpha=alpha_w,
            Pr=Pr_w,
            beta=beta_w
        )
        h_ext = vals['h']
        Nu_ext = vals['Nu']
        Ra_ext = vals['Ra']

        R_ext = 1.0/(h_ext*2.0*math.pi*r_o)
        R_cond = math.log(r_o/r_i)/(2.0*math.pi*geom.k_wall)
        R_int = 1.0/(h_i*2.0*math.pi*r_i)
        R_tot = R_ext + R_cond + R_int
        qprime = (T_w - T_g)/R_tot  # W per meter

        Tso_new = T_w - qprime*R_ext
        if abs(Tso_new - Tso) < tol:
            Tso = Tso_new
            break
        Tso = 0.5*(Tso + Tso_new)

    Tsi = Tso - qprime*math.log(r_o/r_i)/(2.0*math.pi*geom.k_wall)

    # Reynolds and friction / pressure drop for this segment (use Petukhov + curvature multiplier)
    A = 0.25*math.pi*geom.D_i**2
    v = m_dot/(rho_h*A)
    Re = max(rho_h*v*geom.D_i/mu_h, 1.0)
    if Re < 3000:
        f = 64.0/Re
    else:
        f = (0.79*math.log(Re)-1.64)**-2
    # Curvature multiplier similar to Mori–Nakayama; treat as model parameter
    delta = geom.D_i/(2.0*geom.R_coil) if math.isfinite(geom.R_coil) and geom.R_coil>0 else 0.0
    phi = 1.0 + 3.5*delta + 0.7*delta**2
    f_helix = f #*phi

    return {
        'qprime': qprime,
        'h_ext': h_ext,
        'Nu_ext': Nu_ext,
        'Ra_ext': Ra_ext,
        'Tso': Tso,
        'Tsi': Tsi,
        'Tf_w': Tf_w,
        'h_i': h_i,
        'Nu_i': Nu_i,
        'rho_h': rho_h,
        'mu_h': mu_h,
        'k_h': k_h,
        'cp_h': cp_h,
        'Pr_h': Pr_h,
        'nu_h': nu_h,
        'Re': Re,
        'f_helix': f_helix,
        'v': v,
    }


def simulate_bath_coil(
    bath: Bath,
    geom: CoilGeom,
    He: HeliumInlet,
    t_end: float,           # s total simulation time
    dt: float,              # s time step
    Nseg: int = 50,
    hi_func: Callable[[float,float,float,float,float], Tuple[float,float]] = hi_gnielinski_helical,
    callback: Optional[Callable[[float, dict], None]] = None,
    heater: Optional[Heater] = None,
    debug: bool = False,
) -> dict:
    """Explicit time-marching model for N_turn-turn coil dipped in hot, nearly
    quiescent water. Water is a single, well-mixed node. The coil is discretized
    axially; each segment solves a consistent wall energy balance (outer/inner
    convection + wall conduction). External h uses Churchill–Chu (film properties).

    Returns a dict with time histories: water T, helium T_out, Q_dot, and per-time-step
    coil-averaged diagnostics (wall T, film T, helium bulk T, p, Re, Nu_i, rho, mu, k).
    """
    """Explicit time-marching model for N_turn-turn coil dipped in hot, nearly
    quiescent water. Water is a single, well-mixed node. The coil is discretized
    axially; each segment solves a consistent wall energy balance (outer/inner
    convection + wall conduction). External h uses Churchill–Chu (film properties).

    Returns a dict with time histories: water T, helium T_out, Q_dot, and per-time-step
    coil-averaged diagnostics (wall T, film T, helium bulk T, p, Re, Nu_i, rho, mu, k).
    """
    # Histories
    times: List[float] = []
    T_w_hist: List[float] = []
    T_out_hist: List[float] = []
    Q_hist: List[float] = []
    p_out_hist: List[float] = []
    P_in_hist: List[float] = []
    heater_on_hist: List[bool] = []
    avg_hist = {
        'T_wall_avg': [],
        'T_film_avg': [],
        'T_He_avg': [],
        'p_avg': [],
        'Re_avg': [],
        'Nu_He_avg': [],
        'rho_He_avg': [],
        'mu_He_avg': [],
        'lambda_He_avg': []
    }

    T_w = bath.T0
    L = geom.length
    dx = max(geom.length/float(Nseg), 1e-9)

    # --- Preflight sanity and thermostat guardrails ---
    if heater is not None and heater.mode in ('fixed','thermostat') and heater.P_max <= 0.0:
        raise ValueError("Heater.P_max must be > 0 W when mode is 'fixed' or 'thermostat'.")
    if heater is not None and heater.mode == 'fixed' and heater.P_fixed < 0.0:
        raise ValueError("Heater.P_fixed must be >= 0 W for fixed mode.")

    # Estimate initial UA (W/K) and helium m*c_p (W/K) for expectations
    try:
        seg0 = _solve_segment(T_w, He.T_in, He.P, geom, bath, He.m_dot, hi_func)
        r_o = 0.5*geom.D_o; r_i = 0.5*geom.D_i
        R_ext0 = 1.0/(seg0['h_ext']*2.0*math.pi*r_o)
        R_cond0 = math.log(r_o/r_i)/(2.0*math.pi*geom.k_wall)
        R_int0 = 1.0/(seg0['h_i']*2.0*math.pi*r_i)
        UA_prime0 = 1.0/(R_ext0 + R_cond0 + R_int0)                 # W/K per meter
        UA_total0 = UA_prime0 * L                                   # W/K total
        mcp_he0 = He.m_dot * seg0['cp_h']                           # W/K
    except Exception:
        UA_prime0, UA_total0, mcp_he0 = float('nan'), float('nan'), He.m_dot*5200.0
    if debug:
        print(f"[preflight] L={L:.3f} m, UA'≈{UA_prime0:.1f} W/K/m, UA≈{UA_total0:.1f} W/K, m·c_p(He)≈{mcp_he0:.1f} W/K")

    t = 0.0
    # initialize thermostat state based on initial T
    heater_on = False
    if heater is not None and heater.mode == 'thermostat':
        T_on0 = heater.T_set - 0.5*heater.deadband
        heater_on = (T_w <= T_on0)

    while t <= t_end + 1e-12:
        T_g = He.T_in
        p_g = He.P
        Q_dot_total = 0.0

        # Averages along the coil this time step
        acc = {k:0.0 for k in ['Tso','Tf_w','Tg','p','Re','Nu_i','rho','mu','k']}

        for i in range(Nseg):
            seg = _solve_segment(T_w, T_g, p_g, geom, bath, He.m_dot, hi_func)
            # Heat into helium over this segment
            # Use q' per m and local cp for consistency
            dQ = seg['qprime'] * dx
            dTg = dQ/(He.m_dot*seg['cp_h'])
            T_g_next = T_g + dTg

            # Pressure drop over dx: Δp = f*(dx/D)*(ρ v^2/2)
            dp = seg['f_helix'] * (dx/geom.D_i) * 0.5*seg['rho_h']*(seg['v']**2)
            p_next = max(p_g - dp, 1.0)  # keep positive

            # Accumulate
            Q_dot_total += dQ
            acc['Tso'] += seg['Tso']
            acc['Tf_w'] += seg['Tf_w']
            acc['Tg'] += T_g
            acc['p'] += p_g
            acc['Re'] += seg['Re']
            acc['Nu_i'] += seg['Nu_i']
            acc['rho'] += seg['rho_h']
            acc['mu'] += seg['mu_h']
            acc['k'] += seg['k_h']

            print(f"Re={seg['Re']/1e3}e3")

            # March
            T_g = T_g_next
            p_g = p_next

        # Averages along coil
        for key in acc:
            acc[key] /= Nseg

        T_out = T_g
        p_out = p_g

        # Bath cp and energy balance with parasitic loss UA_loss
        if bath.cp_const is not None:
            cp_w = bath.cp_const
        elif _HAS_COOLPROP and bath.fluid_name is not None:
            from CoolProp.CoolProp import PropsSI
            cp_w = PropsSI('CPMASS','T',T_w,'P',bath.P,bath.fluid_name)
        else:
            cp_w = 2000.0  # fallback typical for oils; override via cp_const

        # Heater power (W)
        P_in = 0.0
        if heater is not None and heater.mode != 'off':
            if heater.mode == 'fixed':
                P_in = max(0.0, min(heater.P_fixed, heater.P_max)) * heater.eta
            elif heater.mode == 'thermostat':
                T_on = heater.T_set - 0.5*heater.deadband
                T_off = heater.T_set + 0.5*heater.deadband
                if heater_on and T_w >= T_off:
                    heater_on = False
                elif (not heater_on) and T_w <= T_on:
                    heater_on = True
                P_in = (heater.P_max * heater.eta) if heater_on else 0.0

        Q_loss = bath.UA_loss*(T_w - bath.T_amb)  # W
        dT_w_dt = (P_in - Q_dot_total - Q_loss)/(bath.m_fluid*cp_w)
        T_w_next = T_w + dt*dT_w_dt

        # Save histories
        times.append(t)
        T_w_hist.append(T_w)
        T_out_hist.append(T_out)
        Q_hist.append(Q_dot_total)  # W
        avg_hist['T_wall_avg'].append(acc['Tso'])
        avg_hist['T_film_avg'].append(acc['Tf_w'])
        avg_hist['T_He_avg'].append(acc['Tg'])
        avg_hist['p_avg'].append(acc['p'])
        avg_hist['Re_avg'].append(acc['Re'])
        avg_hist['Nu_He_avg'].append(acc['Nu_i'])
        avg_hist['rho_He_avg'].append(acc['rho'])
        avg_hist['mu_He_avg'].append(acc['mu'])
        avg_hist['lambda_He_avg'].append(acc['k'])
        p_out_hist.append(p_out)
        P_in_hist.append(P_in)
        heater_on_hist.append(bool(heater_on))

        if callback is not None:
            callback(t, {
                'T_w': T_w,
                'T_out': T_out,
                'Q_dot': Q_dot_total,
                'p_out': p_out,
                'P_in': P_in,
                'heater_on': heater_on,
                'cp_w': cp_w,
                'averages': {k: avg_hist[k][-1] for k in avg_hist}
            })

        T_w = T_w_next
        t += dt

    return {
        't': times,
        'T_w': T_w_hist,
        'T_out': T_out_hist,
        'Q_dot': Q_hist,
        'p_out': p_out_hist,
        'averages': avg_hist,
        'P_in': P_in_hist,
        'heater_on': heater_on_hist,
        'meta': {
            'Nseg': Nseg,
            'dt': dt,
            'geom': geom,
            'bath': bath,
            'He': He,
            'UA_prime_W_per_mK_initial': UA_prime0,
            'UA_total_W_per_K_initial': UA_total0,
            'mcp_He_W_per_K_initial': mcp_he0,
        }
    }

def make_live_dashboard(He_inlet: HeliumInlet, show_every: int = 1, heater: Optional[Heater] = None,
                        figsize=(9.5, 6.0), dpi: int = 110, font_size: int = 9):
    """Live dashboard (single figure mosaic):
      • Top-left: Nu_He_avg (purple, left) and Re_He_avg (red, right)
      • Bottom-left: Bath/He outlet temperatures in °C (left) and He outlet pressure in bar (right)
      • Right column (spans rows): Heater power in W (brown)
    Usage: cb = make_live_dashboard(He, heater=heater); simulate_bath_coil(..., callback=cb, heater=heater)
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.gridspec import GridSpec
    mpl.rcParams.update({'figure.dpi': dpi, 'savefig.dpi': dpi, 'font.size': font_size, 'text.antialiased': True})

    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[3.0, 1.6], height_ratios=[1.0, 1.0])

    # --- Top-left: Nu/Re ---
    axNu = fig.add_subplot(gs[0, 0])
    axRe = axNu.twinx()
    axNu.set_ylabel('Nu_He_avg [-]')
    axRe.set_ylabel('Re_He_avg [-]')
    axNu.set_xlabel('Time [s]')
    ln_Nu, = axNu.plot([], [], label='Nu_He_avg', color='tab:purple', linewidth=1.6)
    ln_Re, = axRe.plot([], [], label='Re_He_avg', color='tab:red', linewidth=1.6)
    h_Nu0 = None
    h_Re0 = None
    xs_NR, yNu, yRe = [], [], []

    # --- Bottom-left: Temps (°C) + Pressure (bar) ---
    axT = fig.add_subplot(gs[1, 0])
    axP = axT.twinx()
    axT.set_ylabel('Temperature [°C]')
    axP.set_ylabel('Outlet pressure [bar]')
    axT.set_xlabel('Time [s]')
    ln_Tw, = axT.plot([], [], label='Bath bulk T', color='tab:blue', linewidth=1.6)
    ln_Tout, = axT.plot([], [], label='He outlet T', color='tab:orange', linewidth=1.6)
    ln_Pout, = axP.plot([], [], label='He outlet p', color='tab:green', linewidth=1.6)
    # Baselines (dash-dot) in matching colors (temps converted to °C)
    T_in_C = He_inlet.T_in - 273.15
    h_Tin = axT.axhline(T_in_C, linestyle='-.', linewidth=1.2, color='tab:orange', label='He T_in (start)')
    h_Pin = axP.axhline(He_inlet.P/1e5, linestyle='-.', linewidth=1.2, color='tab:green', label='He p_in (start)')
    xs_TP, yTwC, yToutC, yPbar = [], [], [], []

    # --- Right column: Heater power ---
    axPow = fig.add_subplot(gs[:, 1])
    axPow.set_ylabel('Heater power [W]')
    axPow.set_xlabel('Time [s]')
    ln_Pin, = axPow.plot([], [], label='Heater power', color='tab:brown', linewidth=1.6)
    h_Pmax = None
    if heater is not None and getattr(heater, 'P_max', 0.0) > 0:
        h_Pmax = axPow.axhline(heater.P_max*heater.eta, linestyle='-.', linewidth=1.2,
                               color='tab:brown', label='Heater max')
    xs_P, yPin = [], []

    # Simple, non-overlapping legends (inside axes to avoid layout fights)
    axNu.legend(loc='upper left', framealpha=0.25)
    axRe.legend(loc='upper right', framealpha=0.25)
    axT.legend([ln_Tw, ln_Tout, h_Tin], [l.get_label() for l in [ln_Tw, ln_Tout, h_Tin]],
               loc='upper left', framealpha=0.25)
    axP.legend([ln_Pout, h_Pin], [l.get_label() for l in [ln_Pout, h_Pin]],
               loc='upper right', framealpha=0.25)
    handlesPow = [ln_Pin] + ([h_Pmax] if h_Pmax is not None else [])
    axPow.legend(handlesPow, [h.get_label() for h in handlesPow], loc='upper left', framealpha=0.25)

    # Heater-on shading tracks (on temps and power)
    def shade_if_on(t_now: float, heater_on: bool):
        if heater_on:
            axT.axvspan(t_now - max(1e-9, show_every*1e-6), t_now, color='tab:brown', alpha=0.08)
            axPow.axvspan(t_now - max(1e-9, show_every*1e-6), t_now, color='tab:brown', alpha=0.08)

    def cb(t, data):
        nonlocal h_Nu0, h_Re0
        # Shade if heater is ON
        shade_if_on(t, bool(data.get('heater_on', False)))

        # --- Nu/Re ---
        xs_NR.append(t)
        yNu.append(data['averages']['Nu_He_avg'])
        yRe.append(data['averages']['Re_avg'])
        if len(xs_NR) == 1:
            h_Nu0 = axNu.axhline(yNu[0], linestyle='-.', linewidth=1.2, color='tab:purple', label='Nu start')
            h_Re0 = axRe.axhline(yRe[0], linestyle='-.', linewidth=1.2, color='tab:red', label='Re start')
        if (len(xs_NR) % max(1, show_every)) == 0:
            ln_Nu.set_data(xs_NR, yNu)
            ln_Re.set_data(xs_NR, yRe)
            axNu.relim(); axNu.autoscale_view()
            axRe.relim(); axRe.autoscale_view()

        # --- Temps (°C) + Pressure (bar) ---
        xs_TP.append(t)
        yTwC.append(data['T_w'] - 273.15)
        yToutC.append(data['T_out'] - 273.15)
        yPbar.append(data['p_out']/1e5)
        if (len(xs_TP) % max(1, show_every)) == 0:
            ln_Tw.set_data(xs_TP, yTwC)
            ln_Tout.set_data(xs_TP, yToutC)
            ln_Pout.set_data(xs_TP, yPbar)
            axT.relim(); axT.autoscale_view()
            axP.relim(); axP.autoscale_view()

        # --- Heater power ---
        xs_P.append(t)
        yPin.append(data.get('P_in', 0.0))
        if (len(xs_P) % max(1, show_every)) == 0:
            ln_Pin.set_data(xs_P, yPin)
            axPow.relim(); axPow.autoscale_view()

        plt.pause(0.001)

    return cb

# --- Minimal usage example (requires CoolProp) ---

def make_live_plots(He_inlet: HeliumInlet, show_every: int = 1, heater: Optional[Heater] = None,
                    figsize_tp=(6, 3.6), figsize_pow=(6, 2.6), figsize_nr=(6, 3.6),
                    font_size: int = 9, dpi: int = 110):
    """Return a callback that live-plots with separate figures:
      Fig T&P (left y): Bath bulk T [K] (blue), He outlet T [K] (orange)
               (right y): He outlet p [bar] (green)
               Dash-dot baselines: He T_in, He p_in.
      Fig Power: Heater power [W] (brown) with dash-dot baseline at P_max*eta if provided.
      Fig Nu/Re: Nu_He_avg (purple) + Re_He_avg (red) on twin axes with dash-dot initial baselines.
    Use: cb = make_live_plots(He); simulate_bath_coil(..., callback=cb)
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams.update({'figure.dpi': dpi, 'savefig.dpi': dpi, 'font.size': font_size, 'text.antialiased': True})

    # --- Figure 1: Temps + outlet pressure (twin axes) ---
    figTP, axT = plt.subplots(figsize=figsize_tp, dpi=dpi)
    axP = axT.twinx()
    figTP.subplots_adjust(bottom=0.30)
    axT.set_xlabel('Time [s]')
    axT.set_ylabel('Temperature [K]')
    axP.set_ylabel('Outlet pressure [bar]')

    # Primary lines with distinct colors
    ln_Tw, = axT.plot([], [], label='Bath bulk T', color='tab:blue', linewidth=1.6)
    ln_Tout, = axT.plot([], [], label='He outlet T', color='tab:orange', linewidth=1.6)
    ln_Pout, = axP.plot([], [], label='He outlet p', color='tab:green', linewidth=1.6)

    # Baselines (dashdot)
    h_Tin = axT.axhline(He_inlet.T_in-273, linestyle='-.', linewidth=1.2, color='tab:orange', label='He T_in (start)')
    h_Pin = axP.axhline(He_inlet.P/1e5, linestyle='-.', linewidth=1.2, color='tab:green', label='He p_in (start)')

    xs1, yTw, yTout, yPout = [], [], [], []
    legT_obj = None
    legP_obj = None

    def _update_legends_figTP():
        nonlocal legT_obj, legP_obj
        if legT_obj is None:
            legT_obj = axT.legend([ln_Tw, ln_Tout, h_Tin],
                                  [l.get_label() for l in [ln_Tw, ln_Tout, h_Tin]],
                                  loc='upper center', bbox_to_anchor=(0.5, -0.20),
                                  ncol=3, frameon=False)
        if legP_obj is None:
            legP_obj = axP.legend([ln_Pout, h_Pin],
                                  [l.get_label() for l in [ln_Pout, h_Pin]],
                                  loc='upper center', bbox_to_anchor=(0.5, -0.30),
                                  ncol=2, frameon=False)

    # --- Figure 2: Heater power (separate figure) ---
    figPow, axPow = plt.subplots(figsize=figsize_pow, dpi=dpi)
    figPow.subplots_adjust(bottom=0.18)
    axPow.set_xlabel('Time [s]')
    axPow.set_ylabel('Heater power [W]')
    ln_Pin, = axPow.plot([], [], label='Heater power', color='tab:brown', linewidth=1.6)
    h_Pmax = None
    if heater is not None and heater.P_max > 0:
        h_Pmax = axPow.axhline(heater.P_max*heater.eta, linestyle='-.', linewidth=1.2,
                               color='tab:brown', label='Heater max')
    axPow.grid(True, alpha=0.3)

    xsPow, yPin = [], []
    legPow_obj = None

    def _update_legends_figPow():
        nonlocal legPow_obj
        if legPow_obj is None:
            handlesPow = [ln_Pin] + ([h_Pmax] if h_Pmax is not None else [])
            labelsPow = [h.get_label() for h in handlesPow]
            legPow_obj = axPow.legend(handlesPow, labelsPow, loc='upper center', bbox_to_anchor=(0.5, -0.15),
                                      ncol=len(handlesPow), frameon=False)

    # --- Figure 3: Nu_avg + Re_avg (with heater-on shading applied to Fig T&P and Power) ---
    figNR, axNu = plt.subplots(figsize=figsize_nr, dpi=dpi)
    axRe = axNu.twinx()
    figNR.subplots_adjust(bottom=0.30)
    axNu.set_xlabel('Time [s]')
    axNu.set_ylabel('Nu_He_avg [-]')
    axRe.set_ylabel('Re_He_avg [-]')

    ln_Nu, = axNu.plot([], [], label='Nu_He_avg', color='tab:purple', linewidth=1.6)
    ln_Re, = axRe.plot([], [], label='Re_He_avg', color='tab:red', linewidth=1.6)

    # To be created on first data point
    h_Nu0 = None
    h_Re0 = None
    legNu_obj = None
    legRe_obj = None

    xs2, yNu, yRe = [], [], []

    def _update_legends_figNR():
        nonlocal legNu_obj, legRe_obj
        if legNu_obj is None:
            handlesL = [ln_Nu] + ([h_Nu0] if h_Nu0 is not None else [])
            labelsL = [h.get_label() for h in handlesL]
            legNu_obj = axNu.legend(handlesL, labelsL, loc='upper center', bbox_to_anchor=(0.5, -0.20),
                                    ncol=len(handlesL), frameon=False)
        if legRe_obj is None:
            handlesR = [ln_Re] + ([h_Re0] if h_Re0 is not None else [])
            labelsR = [h.get_label() for h in handlesR]
            legRe_obj = axRe.legend(handlesR, labelsR, loc='upper center', bbox_to_anchor=(0.5, -0.30),
                                    ncol=len(handlesR), frameon=False)

    def cb(t, data):
        nonlocal h_Nu0, h_Re0
        # optional heater-on shading on temperature/pressure and power figures
        if data.get('heater_on', False):
            axT.axvspan(t - max(1e-9, show_every*1e-6), t, color='tab:brown', alpha=0.08)
            axPow.axvspan(t - max(1e-9, show_every*1e-6), t, color='tab:brown', alpha=0.08)

        # --- Fig T&P updates ---
        xs1.append(t)
        yTw.append(data['T_w']-273)
        yTout.append(data['T_out']-273)
        yPout.append(data['p_out']/1e5)
        if (len(xs1) % max(1, show_every)) == 0:
            ln_Tw.set_data(xs1, yTw)
            ln_Tout.set_data(xs1, yTout)
            ln_Pout.set_data(xs1, yPout)
            axT.relim(); axT.autoscale_view()
            axP.relim(); axP.autoscale_view()
            _update_legends_figTP()
            plt.pause(0.001)

        # --- Fig Power updates ---
        xsPow.append(t)
        yPin.append(data.get('P_in', 0.0))
        if (len(xsPow) % max(1, show_every)) == 0:
            ln_Pin.set_data(xsPow, yPin)
            axPow.relim(); axPow.autoscale_view()
            _update_legends_figPow()
            plt.pause(0.001)

        # --- Fig Nu/Re updates ---
        xs2.append(t)
        yNu.append(data['averages']['Nu_He_avg'])
        yRe.append(data['averages']['Re_avg'])
        if len(xs2) == 1:
            h_Nu0 = axNu.axhline(yNu[0], linestyle='-.', linewidth=1.2, color='tab:purple', label='Nu start')
            h_Re0 = axRe.axhline(yRe[0], linestyle='-.', linewidth=1.2, color='tab:red', label='Re start')
            _update_legends_figNR()
        if (len(xs2) % max(1, show_every)) == 0:
            ln_Nu.set_data(xs2, yNu)
            ln_Re.set_data(xs2, yRe)
            axNu.relim(); axNu.autoscale_view()
            axRe.relim(); axRe.autoscale_view()
            _update_legends_figNR()
            plt.pause(0.001)

    return cb

if __name__ == "__main__":
    if _HAS_COOLPROP:
        import matplotlib.pyplot as plt
        # Example A: Water bath via CoolProp with live plotting
        bath = Bath(m_fluid=30, T0=95+273, UA_loss=50, T_amb=20+273, fluid_name='Water')
        geom = CoilGeom(D_o=9e-3, D_i=7e-3, k_wall=20, R_coil=35e-3, N_turn=5)
        He   = HeliumInlet(m_dot=0.015, T_in=20+273, P=25e5)
        cb = make_live_plots(He, show_every=4, heater=Heater(mode='off'))
        out = simulate_bath_coil(bath, geom, He, t_end=10, dt=0.25, Nseg=100, callback=cb, heater=Heater(mode='off'))
        print(f"[Water] Final bath T: {out['T_w'][-1]:.2f} K, He Tout: {out['T_out'][-1]:.2f} K, Re_He=")
        plt.show()
        print("curvature effect on dP OFF for now, to be verified")



    #     # Example B: Thermal oil via CoolProp INCOMP (demo; not edible oil)
    #     oil_bath = Bath(m_fluid=12.0, T0=453.0, UA_loss=12.0, T_amb=295.0, fluid_name='INCOMP::Therminol66')
    #     cb2 = make_live_plots(He, show_every=2, heater=Heater(mode='thermostat', P_max=6000.0, T_set=368.15, deadband=5.0))
    #     out_oil = simulate_bath_coil(oil_bath, geom, He, t_end=300.0, dt=0.5, Nseg=100, callback=cb2,
    #                                  heater=Heater(mode='thermostat', P_max=6000.0, T_set=368.15, deadband=5.0))
    #     print(f"[Therminol66] Final bath T: {out_oil['T_w'][-1]:.2f} K, He Tout: {out_oil['T_out'][-1]:.2f} K")
    #     plt.show()

    #     # Example C: Frying oil using your tabulated data (replace with measured data)
    #     Ttab = [373.0, 423.0, 473.0]  # 100, 150, 200 °C
    #     rho  = [900.0, 880.0, 860.0]  # kg/m^3
    #     mu   = [0.010, 0.0045, 0.0030]  # Pa·s (example)
    #     k    = [0.17, 0.16, 0.15]      # W/m/K
    #     cp   = [2200.0, 2300.0, 2400.0]  # J/kg/K
    #     oil_provider = make_tabulated_bath_props(Ttab, rho, mu, k, cp)
    #     fry_bath = Bath(m_fluid=12.0, T0=443.0, UA_loss=12.0, T_amb=295.0, props_provider=oil_provider, cp_const=2300.0)
    #     cb3 = make_live_plots(He, show_every=2, heater=Heater(mode='fixed', P_max=6000.0, P_fixed=3000.0))
    #     out_fr = simulate_bath_coil(fry_bath, geom, He, t_end=300.0, dt=0.5, Nseg=100, callback=cb3,
    #                                 heater=Heater(mode='fixed', P_max=6000.0, P_fixed=3000.0))
    #     print(f"[Fryer oil (tabulated)] Final bath T: {out_fr['T_w'][-1]:.2f} K, He Tout: {out_fr['T_out'][-1]:.2f} K")
    #     plt.show()
    # else:
    #     print("Install CoolProp to run the example. Core correlation functions still usable.")




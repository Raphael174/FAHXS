"""
pressure_drop_prescribedT.py
─────────────────────────────────────────────────────────────────────────────
Standalone helium-side pressure-drop march under a *prescribed* coolant
temperature profile.

Motivation
----------
The hot-gas Nusselt number is not yet validated against experiment, so any
result that depends on it is untrustworthy.  But we are confident the design
can deliver *enough* energy to hit the target temperature rise.  If we simply
**prescribe** the helium temperature along the coil (a linear rise from inlet
to outlet), the entire heat-transfer side drops out — no Cantera, no
conduction, no hot-gas Nusselt.  What remains is a pure compressible
pressure-drop march: as helium heats, its density collapses, it accelerates,
and friction + acceleration eat static pressure.

Physics
-------
Prescribed temperature along coil arc length s in [0, L_pipe]:

    T(s) = T_in + (T_out - T_in) * s / L_pipe        (linear)

Mass flux is constant (continuity, constant area, single channel):

    G = mdot / (A_ch * N_ch) = rho * U = const,   A_ch = pi * Dh**2 / 4

Per cell i -> i+1 (length ds) we integrate the steady 1-D compressible
momentum equation with real-gas density from CoolProp, splitting the
static-pressure change into an acceleration term and a friction term
(Darcy convention, matching the solver's hot-gas form in main_solve.py):

    p_{i+1} - p_i = - G * (U_{i+1} - U_i)             <- acceleration (momentum) loss
                    - (f/Dh) * (rho*U**2/2)_avg * ds   <- wall friction loss

The equation is implicit in p_{i+1} (rho depends on it); solved with a short
fixed-point iteration (dp-per-cell is tiny so it converges in a few sweeps).

This module is fully standalone — it does NOT touch the coupled solver and
leaves the original code base usable exactly as-is.  It only *reads* the
shared friction correlation `dispatch_friction_coil`.

@ author : (study addition)
"""

import numpy as np
from CoolProp.CoolProp import PropsSI

from .friction_correlations import dispatch_friction_coil, getFrictionColebrook1939


def roughness_multiplier(Re, Dh, eps):
    """Moody rough/smooth Darcy ratio at fixed Re, used to add a roughness
    increment on top of the (hydraulically smooth) Ali 2024 curved-tube law.

    Justified by Ali & Dey (2024): roughness and curvature enter the curved-tube
    skin-friction law as *separable* factors (roughness via the roughness-curvature
    number, curvature via the alpha^(1/2) prefactor), and the law recovers the
    straight-tube Strickler [rough] / Blasius [smooth] limits.  So the *relative*
    roughness effect carries over from the straight-tube Colebrook ratio.  Returns
    1.0 for eps <= 0 (hydraulically smooth).
    """
    if eps <= 0:
        return 1.0
    return getFrictionColebrook1939(Re, eps / Dh) / getFrictionColebrook1939(Re, 0.0)


def march_pressure_prescribedT(
    Dh, mdot, N_ch, T_in, T_out, p_in,
    *,
    D_coil, Rc, L_pipe,
    roughness,
    friction_selector,
    corrCoeffs,
    fluid="Helium",
    n_steps=None,
    p_floor=2e5,
    friction_error_factor=1.0,
    roughness_height=0.0,
):
    """March static pressure along a helical coil under a prescribed linear T(s).

    Parameters
    ----------
    Dh : float        coil inner hydraulic diameter [m]
    mdot : float      total coolant mass flow [kg/s]
    N_ch : int        number of parallel coil channels (mass flow splits evenly)
    T_in, T_out : float   prescribed inlet / outlet temperatures [K]
    p_in : float      inlet static pressure [Pa]
    D_coil : float    coil centre-to-centre diameter [m] (for friction curvature)
    Rc : float        coil radius of curvature [m]
    L_pipe : float    true coil arc length [m]
    roughness : float channel inner roughness [m]
    friction_selector : str   "CurvedPipeAli2024" | "Colebrook1939" (dispatch_friction_coil)
    corrCoeffs : CorrelationCoefficients   supplies ali_c_* etc.
    fluid : str       CoolProp fluid name
    n_steps : int|None  number of march cells; default ~50 per coil turn (solver resolution)

    Returns
    -------
    dict with scalar summaries (dp_total, dp_friction, dp_accel, p_out, U_in,
    U_out, Mach_out, Re_in, Re_out, L_pipe, n_turns, f_in, f_out, rho_in,
    rho_out) and profile arrays (s, T, p, rho, U, Re, f).
    """
    A_ch = np.pi * Dh ** 2 / 4.0
    G = mdot / (A_ch * N_ch)          # mass flux per channel [kg/m^2/s] = rho*U

    if n_steps is None:
        # match the solver's arc resolution: ~50 sub-steps per coil turn
        n_steps = int(np.ceil(L_pipe / (np.pi * D_coil / 50.0)))
        n_steps = max(n_steps, 200)

    s = np.linspace(0.0, L_pipe, n_steps + 1)
    ds = s[1] - s[0]
    T = T_in + (T_out - T_in) * s / L_pipe   # prescribed linear profile

    p = np.empty_like(s)
    rho = np.empty_like(s)
    U = np.empty_like(s)
    Re = np.empty_like(s)
    f = np.empty_like(s)

    # The Ali 2024 curved-tube law is hydraulically smooth (no eps term); add the
    # roughness increment via the Moody ratio.  Colebrook1939 already carries eps,
    # so skip the multiplier there to avoid double-counting.
    apply_rough = (friction_selector == "CurvedPipeAli2024") and (roughness_height > 0)

    def _friction(Re_local):
        f_base = dispatch_friction_coil(
            friction_selector,
            Re=Re_local, Dh=Dh, Rc=D_coil / 2.0,
            roughness=roughness, x=10e10,   # fully developed
            error_factor=friction_error_factor,
            corrCoeffs=corrCoeffs,
        )
        if apply_rough:
            f_base *= roughness_multiplier(Re_local, Dh, roughness_height)
        return f_base

    # ---- node 0 (inlet) ----
    p[0] = p_in
    rho[0] = PropsSI('D', 'T', T[0], 'P', p[0], fluid)
    U[0] = G / rho[0]
    mu0 = PropsSI('V', 'T', T[0], 'P', p[0], fluid)
    Re[0] = G * Dh / mu0
    f[0] = _friction(Re[0])

    dp_friction = 0.0
    dp_accel = 0.0
    valid = True
    last = n_steps            # index of last successfully solved node

    # ---- march ----
    for i in range(n_steps):
        # fixed-point on p[i+1] (rho_{i+1} depends on it)
        p_next = p[i]
        blew = False
        for _ in range(8):
            if p_next < p_floor:
                blew = True
                break
            rho_next = PropsSI('D', 'T', T[i + 1], 'P', p_next, fluid)
            U_next = G / rho_next
            mu_next = PropsSI('V', 'T', T[i + 1], 'P', p_next, fluid)
            Re_next = G * Dh / mu_next
            f_next = _friction(Re_next)

            f_avg = 0.5 * (f[i] + f_next)
            rhoU2_avg = 0.5 * (rho[i] * U[i] ** 2 + rho_next * U_next ** 2)

            d_accel = -G * (U_next - U[i])
            d_fric = -(f_avg / Dh) * (0.5 * rhoU2_avg) * ds
            p_new = p[i] + d_accel + d_fric

            if abs(p_new - p_next) < 1.0:   # 1 Pa tolerance
                p_next = p_new
                break
            p_next = p_new

        if blew or p_next < p_floor:
            # static pressure collapsed before reaching the outlet: this coil
            # cannot pass the flow within the available supply pressure
            # (effectively choked / dp exceeds p_in). Flag and stop.
            valid = False
            last = i
            # truncate profile arrays to what was actually solved
            s, T, p = s[:i + 1], T[:i + 1], p[:i + 1]
            rho, U, Re, f = rho[:i + 1], U[:i + 1], Re[:i + 1], f[:i + 1]
            break

        p[i + 1] = p_next
        rho[i + 1] = rho_next
        U[i + 1] = U_next
        Re[i + 1] = Re_next
        f[i + 1] = f_next
        dp_accel += d_accel
        dp_friction += d_fric

    n_turns = L_pipe / (np.pi * D_coil)
    if valid:
        a_out = PropsSI('A', 'T', T[-1], 'P', p[-1], fluid)   # sound speed
        mach_out = U[-1] / a_out
        dp_total = p[0] - p[-1]
    else:
        # did not reach the outlet — report as exceeding the full supply head
        mach_out = float("nan")
        dp_total = p_in - p_floor   # lower bound on the (infeasible) loss

    return {
        "valid": valid,
        "frac_reached": last / n_steps,
        "dp_total": dp_total,
        "dp_friction": -dp_friction,     # report as positive loss
        "dp_accel": -dp_accel,
        "p_out": p[-1],
        "U_in": U[0], "U_out": U[-1],
        "Mach_out": mach_out,
        "Re_in": Re[0], "Re_out": Re[-1],
        "f_in": f[0], "f_out": f[-1],
        "rho_in": rho[0], "rho_out": rho[-1],
        "L_pipe": L_pipe, "n_turns": n_turns,
        # profiles
        "s": s, "T": T, "p": p, "rho": rho, "U": U, "Re": Re, "f": f,
    }

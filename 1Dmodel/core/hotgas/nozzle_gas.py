"""Quasi-1D nozzle gas expansion + Bartz hot-side HTC — Stage F groundwork.

See docs/solver_design/FV_CORE_REWORK_PLAN.md section 5.2. Standalone module
(no core/mesh or core/residual dependency yet) so it can be built and tested
now, ahead of Stages B-E.

**Chemistry mode implemented here: "frozen" with a SINGLE constant ratio of
specific heats** (evaluated once at the chamber equilibrium state), matching
the design doc's "frozen -- gamma fixed at chamber composition. Cheap, the
validation default" and the user's "keep it simple for now" scope
(2026-07-31). This is the standard first-pass rocket-nozzle sizing
approximation. Transport properties for the Bartz correlation (mu, cp, k, Pr)
are NOT held constant -- they are evaluated from the frozen-composition gas
object at each station's actual local static temperature, which is a real
(if commonly made) internal inconsistency worth naming: the T(M) relation
assumes constant cp/gamma, but mu(T)/cp(T)/k(T) then use the real
temperature-dependent NASA-polynomial values at that T. This does NOT
involve per-node Cantera ``equilibrate()`` calls (forbidden on a repeated
march hot path per CLAUDE.md) -- composition is frozen, so each station is
just a property lookup at fixed Y, cheap and appropriate for this one-off
standalone sizing script. A shifting-equilibrium mode (tabulated, matching
the finite-rate/equilibrium manifold discipline used elsewhere in this repo)
is the natural next step, not implemented here.

**Bartz form implemented: Cornelisse (1979) Eq. 8.3-3, film-property
variant** -- per the user's explicit choice, 2026-07-31 ("we'll use film
properties"), over the RPE Eq. 8-22 stagnation-referenced variant. See the
design doc for why these are NOT interchangeable (different correction-term
reference state and exponent).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from hps_combustor.core.geometry.nozzle_contour import NozzleContour


def area_mach_ratio(M: float, gamma: float) -> float:
    """A/A* for a given Mach number and (constant) ratio of specific heats."""
    if M <= 0.0:
        raise ValueError("M must be positive")
    g = gamma
    term = (2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M**2)
    return (1.0 / M) * term ** ((g + 1.0) / (2.0 * (g - 1.0)))


def mach_from_area_ratio(area_ratio: float, gamma: float, supersonic: bool) -> float:
    """Invert area_mach_ratio for M, on the subsonic (M<1) or supersonic
    (M>1) branch. area_ratio must be >= 1 (== 1 only exactly at the throat)."""
    if area_ratio < 1.0 - 1e-12:
        raise ValueError(f"area_ratio must be >= 1, got {area_ratio}")
    if area_ratio <= 1.0 + 1e-9:
        return 1.0
    f = lambda M: area_mach_ratio(M, gamma) - area_ratio
    if supersonic:
        M_hi = 2.0
        while f(M_hi) < 0.0:
            M_hi *= 2.0
            if M_hi > 1.0e4:
                raise RuntimeError("failed to bracket supersonic Mach root")
        return brentq(f, 1.0 + 1e-9, M_hi)
    return brentq(f, 1.0e-6, 1.0 - 1e-9)


@dataclass(frozen=True)
class ChamberState:
    """Stagnation (chamber) conditions feeding the nozzle expansion."""

    T0_K: float
    p0_Pa: float
    gamma: float
    R_specific_J_kgK: float  # = R_universal / mean molecular weight

    @property
    def rho0_kg_m3(self) -> float:
        return self.p0_Pa / (self.R_specific_J_kgK * self.T0_K)


def isentropic_static_state(M: float, chamber: ChamberState) -> tuple[float, float, float, float]:
    """Constant-gamma isentropic (T, p, rho, V) at Mach M from chamber
    stagnation conditions."""
    g = chamber.gamma
    T0_over_T = 1.0 + 0.5 * (g - 1.0) * M**2
    T = chamber.T0_K / T0_over_T
    p = chamber.p0_Pa * T0_over_T ** (-g / (g - 1.0))
    rho = chamber.rho0_kg_m3 * T0_over_T ** (-1.0 / (g - 1.0))
    V = M * np.sqrt(g * chamber.R_specific_J_kgK * T)
    return T, p, rho, V


def choked_mass_flux(chamber: ChamberState) -> float:
    """Throat mass flux [kg/m2/s] at M=1 for these chamber conditions.

    For a given chamber pressure/temperature/gamma, throat area and total
    mass flow are NOT independent -- choking fixes ``mdot = G* * A_t``. Use
    this (or ``throat_diameter_for_mass_flow``) to check that a stated
    throat diameter and a stated total mass flow are mutually consistent
    before trusting a Bartz result computed from one of them; a
    ~48% mismatch was found and corrected in the user's first C2H4/O2 design
    point, 2026-07-31 (120 mm throat implied ~30.5 kg/s at 50 bar, not the
    stated ~45 kg/s -- see docs/solver_design/FV_CORE_REWORK_PLAN.md).
    """
    T_t, p_t, rho_t, V_t = isentropic_static_state(1.0, chamber)
    return rho_t * V_t


def throat_diameter_for_mass_flow(mdot_kg_s: float, chamber: ChamberState) -> float:
    """Throat diameter [m] consistent with a target total mass flow at these
    chamber conditions (inverse of the choked-flow relation)."""
    G_t = choked_mass_flux(chamber)
    A_t = mdot_kg_s / G_t
    return float(np.sqrt(4.0 * A_t / np.pi))


def adiabatic_wall_temperature(
    T_static_K: float, chamber: ChamberState, M: float, Pr: float, turbulent: bool = True
) -> float:
    """Recovery-corrected driving temperature for convective heat transfer.

    T_aw = T * (1 + r*(gamma-1)/2*M^2), r = Pr^(1/3) turbulent, Pr^(1/2)
    laminar -- NOT T0 and NOT T_static; using either of those instead is the
    most common modelling error in regen-cooling heat-flux estimates (see
    design doc section 5.2).
    """
    r = Pr ** (1.0 / 3.0) if turbulent else Pr**0.5
    g = chamber.gamma
    return T_static_K * (1.0 + r * 0.5 * (g - 1.0) * M**2)


def bartz_cornelisse_htc(
    *,
    mu_Pa_s: float,
    cp_J_kgK: float,
    Pr: float,
    rho_kg_m3: float,
    V_m_s: float,
    D_m: float,
    rho_film_kg_m3: float,
    mu_film_Pa_s: float,
) -> float:
    """Bartz hot-gas HTC, Cornelisse (1979) Eq. 8.3-3, film-property form.

    h_c = 0.026 * (mu^0.2 * cp / Pr^0.6) * (rho*V)^0.8 / D^0.2
          * (rho_f/rho) * (mu_f/mu)

    All of mu, cp, Pr, rho, V, D are LOCAL free-stream (static) values at the
    station of interest; rho_film/mu_film are evaluated at the film
    temperature T_f = 0.5*(T_wall + T_static) (arithmetic mean, per
    Cornelisse's definition). Coefficients verified from a rendered page
    image supplied by the user, 2026-07-31 -- cross-checked against RPE
    (Sutton & Biblarz) Eq. 8-22, which shares the same core group
    (0.026, mu^0.2, Pr^-0.6, (rho*V)^0.8, D^-0.2) but a different correction
    term (see design doc section 5.2 for the discrepancy).
    """
    core = 0.026 * (mu_Pa_s**0.2 * cp_J_kgK / Pr**0.6) * (rho_kg_m3 * V_m_s) ** 0.8 / D_m**0.2
    return core * (rho_film_kg_m3 / rho_kg_m3) * (mu_film_Pa_s / mu_Pa_s)


@dataclass(frozen=True)
class NozzleGasStation:
    z_m: float
    D_m: float
    M: float
    T_K: float
    p_Pa: float
    rho_kg_m3: float
    V_m_s: float
    T_aw_K: float
    h_g_W_m2K: float
    q_w_W_m2: float


def solve_frozen_expansion(
    contour: NozzleContour,
    chamber: ChamberState,
    frozen_gas,  # cantera.Solution with composition already set (frozen)
    *,
    T_wall_guess_K: float,
    turbulent_recovery: bool = True,
) -> list[NozzleGasStation]:
    """March the frozen-composition constant-gamma expansion along the
    contour and evaluate Bartz (Cornelisse, film-property) at each station.

    ``T_wall_guess_K`` is a placeholder, uniform assumed hot-wall temperature
    used only to evaluate Bartz's film properties and a first-pass q_w. This
    is NOT a coupled wall/coolant solve (that is core/wall.py +
    core/residual.py, not built yet) -- it is a deliberate first-order
    estimate, exactly the "keep it simple for now" scope requested
    2026-07-31. Revisit once the coupled solve exists.
    """
    stations: list[NozzleGasStation] = []
    for z, D, A in zip(contour.z_m, contour.D_m2, contour.A_m2):
        area_ratio = A / contour.A_t_m2
        supersonic = z > 0.0
        M = mach_from_area_ratio(area_ratio, chamber.gamma, supersonic=supersonic)
        T, p, rho, V = isentropic_static_state(M, chamber)

        frozen_gas.TP = T, p
        mu = frozen_gas.viscosity
        cp = frozen_gas.cp_mass
        k = frozen_gas.thermal_conductivity
        Pr = mu * cp / k

        T_aw = adiabatic_wall_temperature(T, chamber, M, Pr, turbulent=turbulent_recovery)

        T_film = 0.5 * (T_wall_guess_K + T)
        frozen_gas.TP = T_film, p
        rho_film = frozen_gas.density
        mu_film = frozen_gas.viscosity

        h_g = bartz_cornelisse_htc(
            mu_Pa_s=mu,
            cp_J_kgK=cp,
            Pr=Pr,
            rho_kg_m3=rho,
            V_m_s=V,
            D_m=D,
            rho_film_kg_m3=rho_film,
            mu_film_Pa_s=mu_film,
        )
        q_w = h_g * (T_aw - T_wall_guess_K)

        stations.append(
            NozzleGasStation(
                z_m=float(z), D_m=float(D), M=float(M), T_K=float(T), p_Pa=float(p),
                rho_kg_m3=float(rho), V_m_s=float(V), T_aw_K=float(T_aw),
                h_g_W_m2K=float(h_g), q_w_W_m2=float(q_w),
            )
        )
    return stations

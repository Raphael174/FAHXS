"""Thermodynamic and heat-transfer regime detection for the coolant closure.

This module is the fluid-agnostic front door to state closure and regime
classification. It answers three questions the downstream correlation registry
needs, for ANY CoolProp fluid at any (p, h):

1. Which thermodynamic regime is this state in? -- subcooled liquid, two-phase,
   superheated vapor (all subcritical), OR supercritical liquid-like /
   pseudo-critical / gas-like (p >= p_crit, no phase dome at all).
2. What is the (p, h) state closure, WITHOUT assuming a saturation dome exists?
   ``equilibrium_state_ph`` in correlations.py raises for p >= p_crit by design
   (it is a subcritical, dome-based closure); ``real_fluid_state_ph`` here
   dispatches to it below the critical pressure and to a direct single-phase
   real-EOS flash above it.
3. What are the buoyancy / property-variation indicators (buoyancy parameter,
   Eckert number) that supercritical correlations and heat-transfer-
   deterioration (HTD) detection depend on?

Design note: regime labels are used to SELECT a correlation and to raise honest
diagnostics, not as hard physical switches with discontinuous behavior at the
boundary. The pseudo-critical band in particular is a labeling aid, not a sharp
transition -- property variation through the pseudo-critical point is smooth.
"""

from __future__ import annotations

from functools import lru_cache

import CoolProp.CoolProp as CP
from scipy.optimize import minimize_scalar

from hps_combustor.physics.liquid_flow.coolprop_state_cache import get_cached_state
from hps_combustor.physics.liquid_flow.correlations import (
    EquilibriumState,
    equilibrium_state_ph,
)

# Regime labels. The three subcritical ones mirror EquilibriumState.phase
# values produced by equilibrium_state_ph; the three supercritical ones are new.
SUBCOOLED_LIQUID = "subcooled_liquid"
TWO_PHASE = "two_phase"
SUPERHEATED_VAPOR = "superheated_vapor"
SUPERCRITICAL_LIQUID_LIKE = "supercritical_liquid_like"
PSEUDO_CRITICAL = "pseudo_critical"
SUPERCRITICAL_GAS_LIKE = "supercritical_gas_like"

SUPERCRITICAL_REGIMES = frozenset(
    {SUPERCRITICAL_LIQUID_LIKE, PSEUDO_CRITICAL, SUPERCRITICAL_GAS_LIKE}
)

# Half-width of the pseudo-critical band, as a fraction of the pseudo-critical
# temperature. |T - T_pc| < band * T_pc is labeled ``pseudo_critical``. Purely a
# labeling aid for correlation selection/diagnostics (property variation is
# smooth through T_pc); not a physical discontinuity. For N2 at 80 bar
# (T_pc ~ 145.7 K) this is a +/- ~7 K window.
PSEUDO_CRITICAL_BAND_FRACTION = 0.05


def is_supercritical(fluid: str, p_Pa: float) -> bool:
    """True when pressure is at or above the fluid's critical pressure."""
    return p_Pa >= get_cached_state(fluid).p_crit_Pa


def pseudo_critical_temperature(fluid: str, p_Pa: float) -> float:
    """Pseudo-critical temperature T_pc(p): the temperature of maximum isobaric
    cp at a supercritical pressure (the smeared-out remnant of the latent-heat
    spike). Cached per (fluid, pressure rounded to 100 Pa) since it drives only
    regime labeling and the Eckert number, never a per-node hot path.

    Raises ValueError for subcritical pressure (no pseudo-critical point below
    p_crit -- there is a real saturation temperature there instead).
    """
    if not is_supercritical(fluid, p_Pa):
        raise ValueError("pseudo_critical_temperature requires supercritical pressure")
    return _pseudo_critical_temperature_cached(fluid, round(float(p_Pa), -2))


@lru_cache(maxsize=256)
def _pseudo_critical_temperature_cached(fluid: str, p_Pa: float) -> float:
    cached = get_cached_state(fluid)
    T_crit = cached.T_crit_K
    # The cp peak sits just above T_crit at p just above p_crit and rises with
    # pressure; [T_crit, 2.5*T_crit] brackets it for the pressures of interest
    # (e.g. N2 80 bar: T_crit=126.2 K, T_pc~145.7 K, well inside [126, 315]).
    result = minimize_scalar(
        lambda T: -cached.cp_at_tp(T, p_Pa),
        bounds=(T_crit * 1.0001, T_crit * 2.5),
        method="bounded",
    )
    return float(result.x)


def classify(fluid: str, p_Pa: float, h_J_kg: float) -> str:
    """Return the regime label for a (fluid, p, h) state.

    Subcritical: delegates to the dome-based EquilibriumState phase. Supercritical:
    liquid-like / pseudo-critical / gas-like relative to T_pc(p).
    """
    if not is_supercritical(fluid, p_Pa):
        eq = equilibrium_state_ph(p_Pa, h_J_kg, fluid)
        if eq.phase == "two_phase":
            return TWO_PHASE
        return SUBCOOLED_LIQUID if eq.quality < 0.0 else SUPERHEATED_VAPOR
    T = get_cached_state(fluid).flash_ph(p_Pa, h_J_kg).T()
    return classify_supercritical(fluid, p_Pa, T)


def classify_supercritical(fluid: str, p_Pa: float, T_K: float) -> str:
    """Supercritical regime label from temperature relative to T_pc(p)."""
    T_pc = pseudo_critical_temperature(fluid, p_Pa)
    band = PSEUDO_CRITICAL_BAND_FRACTION * T_pc
    if T_K < T_pc - band:
        return SUPERCRITICAL_LIQUID_LIKE
    if T_K > T_pc + band:
        return SUPERCRITICAL_GAS_LIKE
    return PSEUDO_CRITICAL


def supercritical_state_ph(fluid: str, p_Pa: float, h_J_kg: float) -> EquilibriumState:
    """Single-phase real-EOS state closure at a supercritical (p, h).

    There is no dome, so quality and void fraction are undefined: quality is set
    to NaN and void to 0.0 (never raises, unlike the subcritical closure). The
    ``phase`` field carries the supercritical regime label so downstream
    consumers can branch on it.
    """
    st = get_cached_state(fluid).flash_ph(p_Pa, h_J_kg)
    T = float(st.T())
    return EquilibriumState(
        fluid=fluid,
        p_Pa=float(p_Pa),
        h_J_kg=float(h_J_kg),
        T_K=T,
        quality=float("nan"),
        void_fraction=0.0,
        rho_kg_m3=float(st.rhomass()),
        phase=classify_supercritical(fluid, p_Pa, T),
    )


def real_fluid_state_ph(fluid: str, p_Pa: float, h_J_kg: float) -> EquilibriumState:
    """State closure valid at any pressure: subcritical dome-based below p_crit,
    single-phase real-EOS above it. Replaces bare ``equilibrium_state_ph`` calls
    that would otherwise crash on supercritical inputs.
    """
    if is_supercritical(fluid, p_Pa):
        return supercritical_state_ph(fluid, p_Pa, h_J_kg)
    return equilibrium_state_ph(p_Pa, h_J_kg, fluid)


def eckert_number(*, T_pc_K: float, T_bulk_K: float, T_wall_K: float) -> float:
    """Eckert number E = (T_pc - T_b) / (T_w - T_b), the phase-condition
    criterion used by the Yamagata/Wang2023 regime split: E > 1 liquid-like,
    E < 0 gas-like, 0 <= E <= 1 straddling the pseudo-critical point. Returns
    NaN when wall and bulk temperatures coincide (criterion undefined)."""
    denom = T_wall_K - T_bulk_K
    if denom == 0.0:
        return float("nan")
    return (T_pc_K - T_bulk_K) / denom


def buoyancy_parameter(
    *,
    beta_bulk_1_K: float,
    heat_flux_W_m2: float,
    k_bulk: float,
    rho_bulk: float,
    mu_bulk: float,
    diameter_m: float,
    Re_bulk: float,
    Pr_bulk: float,
) -> float:
    """Jackson (McEligot-Jackson) buoyancy parameter Bo* for supercritical
    channel flow, exactly as used by Urbano & Nasuti (2013), verified against
    docs/reference/JTHT2013.pdf (nomenclature + Sec. IV):

        Gr* = g * beta_b * q_w * D^4 / (k_b * nu_b^2),   nu_b = mu_b / rho_b
        Bo* = Gr* / (Re_b^3.425 * Pr_b^0.8)

    where beta_b = -(1/rho)(d rho/dT)_p is the isobaric thermal expansion
    coefficient at bulk (CoolProp isobaric_expansion_coefficient). This is a
    HEAT-FLUX Grashof (buoyancy driven by the wall-heat-flux-induced density
    defect), NOT the wall-bulk density-difference form used previously. Buoyancy
    (mixed convection) is significant, i.e. the pure-forced-convection
    property-ratio correlations are suspect, when Bo* > 6e-7 (McEligot-Jackson).

    This is a DETECTION flag only: there is no validated deteriorated-HTC
    magnitude model for nitrogen at these conditions (Urbano-Nasuti's q_w/G
    threshold correlation is methane-specific). Flagged nodes are treated as
    outside the validated forced-convection envelope, not silently corrected.
    """
    nu_b = mu_bulk / rho_bulk
    gr_star = 9.80665 * beta_bulk_1_K * heat_flux_W_m2 * diameter_m**4 / (k_bulk * nu_b**2)
    return gr_star / (Re_bulk**3.425 * Pr_bulk**0.8)


def acceleration_parameter(
    *,
    beta_bulk_1_K: float,
    heat_flux_W_m2: float,
    mass_flux_kg_m2_s: float,
    cp_bulk: float,
    Re_bulk: float,
    Pr_bulk: float,
) -> float:
    """McEligot flow-acceleration parameter (Urbano-Nasuti 2013, Sec. IV):

        q+ = (q_w / G) / (cp_b / beta_b) = q_w * beta_b / (G * cp_b)
        Kv = 4 q+ / (Re_b^0.625 * Pr_b^0.4)

    Thermal-acceleration effects (bulk speed-up as the fluid heats and expands)
    are significant, again violating the pure-forced-convection assumption, when
    Kv > 2.9e-5.
    """
    G = max(float(mass_flux_kg_m2_s), 1.0e-30)
    q_plus = heat_flux_W_m2 * beta_bulk_1_K / (G * cp_bulk)
    return 4.0 * q_plus / (Re_bulk**0.625 * Pr_bulk**0.4)


BUOYANCY_THRESHOLD = 6.0e-7      # Bo* above this => buoyancy/mixed convection significant
ACCELERATION_THRESHOLD = 2.9e-5  # Kv above this => flow-acceleration effects significant


def htd_risk(buoyancy_param: float, acceleration_param: float = 0.0) -> bool:
    """True when buoyancy OR flow-acceleration is significant -- i.e. the
    forced-convection property-ratio correlations (McCarthy-Wolf/Taylor) assume
    negligible buoyancy and acceleration, and here that assumption is violated,
    so the closure result should be treated with caution. Named ``htd_risk``
    because these are the classic precursors to heat-transfer deterioration; it
    is an assumption-validity flag, not a deteriorated-magnitude prediction."""
    return abs(buoyancy_param) > BUOYANCY_THRESHOLD or abs(acceleration_param) > ACCELERATION_THRESHOLD

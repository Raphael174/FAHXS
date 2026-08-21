"""Liquid and equilibrium two-phase coolant correlations.

The functions in this module are intentionally geometry-light. They provide the
state closure and straight-pipe correlations needed for the first liquid
coolant proof of concept before wiring the same ideas into HX-specific solvers.

CHF/dryout lookup lives in ``physics/liquid_flow/chf.py`` — it is independent
of the property/HTC/friction closures here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import CoolProp.CoolProp as CP
import numpy as np

from hps_combustor.physics.liquid_flow.coolprop_state_cache import get_cached_state


@dataclass(frozen=True)
class SaturationState:
    fluid: str
    p_Pa: float
    T_sat_K: float
    h_l_J_kg: float
    h_v_J_kg: float
    rho_l_kg_m3: float
    rho_v_kg_m3: float
    mu_l_Pa_s: float
    mu_v_Pa_s: float
    k_l_W_m_K: float
    k_v_W_m_K: float
    cp_l_J_kg_K: float
    cp_v_J_kg_K: float
    sigma_N_m: float
    molar_mass_kg_mol: float
    p_crit_Pa: float
    # Saturated-phase sound speed [m/s], stored alongside the rest of the
    # saturation state so two_phase_sound_speed() doesn't need two more
    # redundant flashes at Q=0/Q=1 on top of the ones already done here.
    c_l_m_s: float = float("nan")
    c_v_m_s: float = float("nan")

    @property
    def h_fg_J_kg(self) -> float:
        return self.h_v_J_kg - self.h_l_J_kg

    @property
    def pr_l(self) -> float:
        return self.cp_l_J_kg_K * self.mu_l_Pa_s / self.k_l_W_m_K

    @property
    def pr_v(self) -> float:
        return self.cp_v_J_kg_K * self.mu_v_Pa_s / self.k_v_W_m_K


@dataclass(frozen=True)
class EquilibriumState:
    fluid: str
    p_Pa: float
    h_J_kg: float
    T_K: float
    quality: float
    void_fraction: float
    rho_kg_m3: float
    phase: str


def saturation_state(fluid: str, p_Pa: float) -> SaturationState:
    """Return saturated liquid/vapor properties at pressure ``p_Pa``.

    Routed through a cached low-level CoolProp AbstractState (see
    coolprop_state_cache.py) instead of ~10 independent high-level PropsSI
    calls - a performance refactor only for the default "HEOS" backend (same
    equation of state, same answer to floating-point precision). ``fluid``
    may be backend-tagged (e.g. "BICUBIC::Water", via
    coolprop_fluid_string()) to opt into a faster, interpolated - and
    therefore approximate - property backend; see
    validation/liquid_ttse_backend_validation.py for measured error bounds.
    """
    if p_Pa <= 0.0:
        raise ValueError("pressure must be positive")
    cached = get_cached_state(fluid)
    if p_Pa >= cached.p_crit_Pa:
        raise ValueError("saturation_state requires subcritical pressure")
    liq = cached.saturated_liquid(p_Pa)
    vap = cached.saturated_vapor(p_Pa)
    return SaturationState(
        fluid=fluid,
        p_Pa=p_Pa,
        T_sat_K=liq.T(),
        h_l_J_kg=liq.hmass(),
        h_v_J_kg=vap.hmass(),
        rho_l_kg_m3=liq.rhomass(),
        rho_v_kg_m3=vap.rhomass(),
        mu_l_Pa_s=liq.viscosity(),
        mu_v_Pa_s=vap.viscosity(),
        k_l_W_m_K=liq.conductivity(),
        k_v_W_m_K=vap.conductivity(),
        cp_l_J_kg_K=liq.cpmass(),
        cp_v_J_kg_K=vap.cpmass(),
        sigma_N_m=cached.surface_tension_liquid(p_Pa),
        molar_mass_kg_mol=cached.molar_mass_kg_mol,
        p_crit_Pa=cached.p_crit_Pa,
        c_l_m_s=liq.speed_sound(),
        c_v_m_s=vap.speed_sound(),
    )


def two_phase_sound_speed(
    *, p_Pa: float, void_fraction: float, rho_mix_kg_m3: float, fluid: str
) -> float:
    """Wood's equation homogeneous-mixture sound speed for a saturated
    liquid-vapor mixture (Wood 1930; see e.g. Brennen, "Fundamentals of
    Multiphase Flow", Eq. 4.20-4.22).

    Uses volumetric VOID fraction (not mass quality) - the two are related
    but distinct for a density-mismatched pair like water/steam. Only
    meaningful for 0 <= void_fraction <= 1 (inside the two-phase dome):
    the characteristic collapse this formula predicts (mixture sound speed
    can drop to tens of m/s, far below either pure-phase sound speed) is a
    real two-phase compressibility effect and does not apply to a
    single-phase state - use CoolProp's real-EOS SPEED_OF_SOUND directly
    outside the dome instead.
    """
    sat = saturation_state(fluid, p_Pa)
    c_l = sat.c_l_m_s
    c_v = sat.c_v_m_s
    alpha = min(max(float(void_fraction), 0.0), 1.0)
    inv_rho_c2 = (
        alpha / (sat.rho_v_kg_m3 * c_v**2)
        + (1.0 - alpha) / (sat.rho_l_kg_m3 * c_l**2)
    )
    return 1.0 / math.sqrt(max(float(rho_mix_kg_m3), 1.0e-9) * inv_rho_c2)


def thermodynamic_quality(p_Pa: float, h_J_kg: float, fluid: str = "Water") -> float:
    """Return equilibrium thermodynamic quality, unclipped outside [0, 1]."""
    sat = saturation_state(fluid, p_Pa)
    return (h_J_kg - sat.h_l_J_kg) / sat.h_fg_J_kg


def homogeneous_void_fraction(x: float, rho_l: float, rho_v: float) -> float:
    """Homogeneous equilibrium void fraction from mass quality and densities."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    vapor_volume = x / rho_v
    liquid_volume = (1.0 - x) / rho_l
    return vapor_volume / (vapor_volume + liquid_volume)


def equilibrium_state_ph(
    p_Pa: float,
    h_J_kg: float,
    fluid: str = "Water",
) -> EquilibriumState:
    """Close a pure-fluid state from pressure and enthalpy.

    Two-phase states use homogeneous equilibrium density. Single-phase states
    are reconstructed directly with CoolProp.
    """
    sat = saturation_state(fluid, p_Pa)
    x = (h_J_kg - sat.h_l_J_kg) / sat.h_fg_J_kg
    if 0.0 <= x <= 1.0:
        alpha = homogeneous_void_fraction(x, sat.rho_l_kg_m3, sat.rho_v_kg_m3)
        rho = 1.0 / (x / sat.rho_v_kg_m3 + (1.0 - x) / sat.rho_l_kg_m3)
        return EquilibriumState(
            fluid=fluid,
            p_Pa=p_Pa,
            h_J_kg=h_J_kg,
            T_K=sat.T_sat_K,
            quality=x,
            void_fraction=alpha,
            rho_kg_m3=rho,
            phase="two_phase",
        )
    flashed = get_cached_state(fluid).flash_ph(p_Pa, h_J_kg)
    rho = flashed.rhomass()
    T = flashed.T()
    return EquilibriumState(
        fluid=fluid,
        p_Pa=p_Pa,
        h_J_kg=h_J_kg,
        T_K=T,
        quality=x,
        void_fraction=0.0 if x < 0.0 else 1.0,
        rho_kg_m3=rho,
        phase="liquid" if x < 0.0 else "vapor",
    )


def darcy_friction_smooth_pipe(Re: float) -> float:
    """Darcy friction factor for a smooth circular tube."""
    Re = max(float(Re), 1.0e-12)
    if Re < 2300.0:
        return 64.0 / Re
    if Re > 4000.0:
        return 1.0 / (0.79 * math.log(Re) - 1.64) ** 2
    gamma = (Re - 2300.0) / 1700.0
    f_lam = 64.0 / 2300.0
    f_turb = 1.0 / (0.79 * math.log(4000.0) - 1.64) ** 2
    return (1.0 - gamma) * f_lam + gamma * f_turb


def liquid_single_phase_nusselt(Re: float, Pr: float, f_darcy: float | None = None) -> float:
    """Circular-tube single-phase liquid Nusselt number.

    Uses a conservative fully developed laminar value below Re=2300 and a
    Gnielinski turbulent value above Re=4000, with linear transition blending.
    """
    Re = max(float(Re), 1.0e-12)
    Pr = max(float(Pr), 1.0e-12)
    if Re < 2300.0:
        return 3.66
    if f_darcy is None:
        f_darcy = darcy_friction_smooth_pipe(Re)
    def _gnielinski(Re_value: float) -> float:
        f_value = darcy_friction_smooth_pipe(Re_value)
        num = (f_value / 8.0) * (Re_value - 1000.0) * Pr
        den = 1.0 + 12.7 * math.sqrt(f_value / 8.0) * (Pr ** (2.0 / 3.0) - 1.0)
        return max(num / den, 3.66)
    if Re > 4000.0:
        num = (f_darcy / 8.0) * (Re - 1000.0) * Pr
        den = 1.0 + 12.7 * math.sqrt(f_darcy / 8.0) * (Pr ** (2.0 / 3.0) - 1.0)
        return max(num / den, 3.66)
    gamma = (Re - 2300.0) / 1700.0
    return (1.0 - gamma) * 3.66 + gamma * _gnielinski(4000.0)


def post_chf_dispersed_flow_htc(
    *, p_Pa: float, mass_flux_kg_m2_s: float, diameter_m: float, fluid: str
) -> float:
    """Simplified post-CHF (post-dryout) coolant-side heat transfer coefficient.

    Conservative simplification: treats the coolant as single-phase saturated
    vapor carrying the FULL mass flux (ignoring the heat-transfer-enhancing
    effect of entrained liquid droplets that a real dispersed-flow/mist-flow
    correlation, e.g. Groeneveld 5.7 or Bromley, would capture). This
    deliberately under-predicts post-CHF HTC relative to reality, so it is a
    defensible bound but not a validated post-dryout closure — if tighter
    fidelity is needed, replace with a literature dispersed-flow correlation.
    """
    sat = saturation_state(fluid, p_Pa)
    Re_v = mass_flux_kg_m2_s * diameter_m / sat.mu_v_Pa_s
    f_v = darcy_friction_smooth_pipe(Re_v)
    Nu_v = liquid_single_phase_nusselt(Re_v, sat.pr_v, f_v)
    return Nu_v * sat.k_v_W_m_K / diameter_m


def bergles_rohsenow_onb_wall_superheat(*, p_Pa: float, heat_flux_W_m2: float) -> float:
    """Bergles-Rohsenow (1964) onset-of-nucleate-boiling (ONB) wall superheat.

    Returns the wall superheat DeltaT_ONB = T_wall - T_sat [K] at which
    nucleate boiling is predicted to begin for a given wall heat flux,
    independent of bulk quality - subcooled boiling can start well before the
    bulk fluid reaches saturation if the local wall is hot enough, which the
    quality-only blend window in ``dispatch.py`` cannot represent. SI form
    (Collier & Thome, "Convective Boiling and Condensation", 3rd ed.):

        q_ONB [W/m^2] = 1082 * p^1.156 * (1.8 * DeltaT_ONB)^(2.16 / p^0.0234)

    with p in bar. Inverted here for DeltaT_ONB given q''. Returns +inf for
    zero or negative heat flux (no boiling can be triggered by cooling).
    """
    if heat_flux_W_m2 <= 0.0:
        return float("inf")
    p_bar = max(p_Pa, 1.0) / 1.0e5
    exponent = p_bar ** 0.0234 / 2.16
    return (heat_flux_W_m2 / (1082.0 * p_bar ** 1.156)) ** exponent / 1.8


def martinelli_parameter_tt(x: float, rho_l: float, rho_v: float, mu_l: float, mu_v: float) -> float:
    """Turbulent-turbulent Martinelli parameter X_tt."""
    x = min(max(float(x), 1.0e-9), 1.0 - 1.0e-9)
    return ((1.0 - x) / x) ** 0.9 * (rho_v / rho_l) ** 0.5 * (mu_l / mu_v) ** 0.1


def martinelli_parameter_laminar_liquid_turbulent_vapor(
    *,
    quality: float,
    rho_l: float,
    rho_v: float,
    Re_l: float,
    Re_v: float,
) -> float:
    """Yu et al. 2002 laminar-liquid/turbulent-vapor Martinelli parameter.

    This is Eq. (9) in Yu, France, Wambsganss, and Hull (2002) for the
    low-mass-flux 2.98 mm horizontal water tube data. It is intentionally kept
    separate from ``martinelli_parameter_tt`` because the phase-only Reynolds
    regimes differ.
    """
    x = min(max(float(quality), 1.0e-9), 1.0 - 1.0e-9)
    Re_l = max(float(Re_l), 1.0e-12)
    Re_v = max(float(Re_v), 1.0e-12)
    return 18.65 * (rho_v / rho_l) ** 0.5 * ((1.0 - x) / x) * Re_v**0.1 / Re_l**0.5


def chisholm_two_phase_multiplier(X: float, C: float = 12.0) -> float:
    """Chisholm pressure-gradient multiplier ``phi_l^2`` from Martinelli ``X``."""
    X = max(float(X), 1.0e-12)
    return 1.0 + C / X + 1.0 / X**2


def yu2002_small_channel_pressure_multiplier(X: float) -> float:
    """Yu et al. 2002 fitted small-channel water multiplier ``phi_l^2``.

    Yu et al. report that the standard Chisholm form over-predicted their
    small horizontal water-tube data with an RMS error of 33%, while this
    fitted form gave about 7% RMS error over their plotted dataset.
    """
    X = max(float(X), 1.0e-12)
    return X**-1.9


def yu2002_modified_anl_boiling_htc(
    *,
    p_Pa: float,
    mass_flux_kg_m2_s: float,
    diameter_m: float,
    heat_flux_W_m2: float,
    fluid: str = "Water",
) -> float:
    """Yu et al. 2002 modified ANL small-channel boiling HTC correlation.

    Implements Eq. (16)-(18) from Yu et al. for water in a 2.98 mm horizontal
    tube. The correlation is intended for nucleate-type boiling before
    transition boiling and is independent of local quality in the paper's
    fitted form.
    """
    if diameter_m <= 0.0 or mass_flux_kg_m2_s <= 0.0 or heat_flux_W_m2 < 0.0:
        raise ValueError("diameter, mass flux, and heat flux must be non-negative with D,G > 0")
    sat = saturation_state(fluid, p_Pa)
    Bo = heat_flux_W_m2 / (mass_flux_kg_m2_s * sat.h_fg_J_kg) if heat_flux_W_m2 > 0.0 else 0.0
    We_l = mass_flux_kg_m2_s**2 * diameter_m / (sat.rho_l_kg_m3 * sat.sigma_N_m)
    density_ratio = sat.rho_l_kg_m3 / sat.rho_v_kg_m3
    return 6.4e6 * (Bo**2 * We_l) ** 0.27 * density_ratio**-0.2


def gungor_winterton_boiling_htc(
    *,
    p_Pa: float,
    mass_flux_kg_m2_s: float,
    diameter_m: float,
    quality: float,
    heat_flux_W_m2: float,
    fluid: str = "Water",
) -> float:
    """Gungor-Winterton saturated flow-boiling HTC for tubes/annuli.

    Reference: Gungor and Winterton, Int. J. Heat Mass Transfer, 1986.
    Horizontal low-Froude correction is intentionally omitted in this first
    straight-pipe PoC; callers should treat this as a vertical/high-Fr closure.
    """
    if diameter_m <= 0.0 or mass_flux_kg_m2_s <= 0.0 or heat_flux_W_m2 < 0.0:
        raise ValueError("diameter, mass flux, and heat flux must be non-negative with D,G > 0")
    sat = saturation_state(fluid, p_Pa)
    x = min(max(float(quality), 1.0e-9), 1.0 - 1.0e-9)
    Re_l = mass_flux_kg_m2_s * (1.0 - x) * diameter_m / sat.mu_l_Pa_s
    h_l = 0.023 * Re_l**0.8 * sat.pr_l**0.4 * sat.k_l_W_m_K / diameter_m
    x_tt = martinelli_parameter_tt(
        x, sat.rho_l_kg_m3, sat.rho_v_kg_m3, sat.mu_l_Pa_s, sat.mu_v_Pa_s
    )
    Bo = heat_flux_W_m2 / (mass_flux_kg_m2_s * sat.h_fg_J_kg) if heat_flux_W_m2 > 0.0 else 0.0
    enhancement = 1.0 + 24000.0 * Bo**1.16 + 1.37 * (1.0 / x_tt) ** 0.86
    suppression = 1.0 / (1.0 + 1.15e-6 * enhancement**2 * Re_l**1.17)
    if heat_flux_W_m2 <= 0.0:
        return enhancement * h_l
    p_r = p_Pa / sat.p_crit_Pa
    molar_mass_kg_kmol = sat.molar_mass_kg_mol * 1000.0
    h_pool = (
        55.0
        * p_r**0.12
        * (-math.log10(p_r)) ** -0.55
        * molar_mass_kg_kmol**-0.5
        * heat_flux_W_m2**0.67
    )
    return enhancement * h_l + suppression * h_pool


def muller_steinhagen_heck_friction_gradient(
    *,
    p_Pa: float,
    mass_flux_kg_m2_s: float,
    diameter_m: float,
    quality: float,
    fluid: str = "Water",
) -> float:
    """Muller-Steinhagen and Heck two-phase frictional pressure gradient.

    Returns positive ``-dp/dz`` in Pa/m for a smooth circular tube.
    """
    if diameter_m <= 0.0 or mass_flux_kg_m2_s <= 0.0:
        raise ValueError("diameter and mass flux must be positive")
    sat = saturation_state(fluid, p_Pa)
    x = min(max(float(quality), 0.0), 1.0)
    Re_l0 = mass_flux_kg_m2_s * diameter_m / sat.mu_l_Pa_s
    Re_v0 = mass_flux_kg_m2_s * diameter_m / sat.mu_v_Pa_s
    f_l0 = darcy_friction_smooth_pipe(Re_l0)
    f_v0 = darcy_friction_smooth_pipe(Re_v0)
    A = f_l0 * mass_flux_kg_m2_s**2 / (2.0 * diameter_m * sat.rho_l_kg_m3)
    B = f_v0 * mass_flux_kg_m2_s**2 / (2.0 * diameter_m * sat.rho_v_kg_m3)
    interpolation = A + 2.0 * (B - A) * x
    return interpolation * (1.0 - x) ** (1.0 / 3.0) + B * x**3


def homogeneous_acceleration_pressure_gradient(
    *,
    mass_flux_kg_m2_s: float,
    p_Pa: float,
    quality_gradient_1_m: float,
    fluid: str = "Water",
) -> float:
    """HEM acceleration pressure gradient contribution, positive ``-dp/dz``."""
    sat = saturation_state(fluid, p_Pa)
    dv_dx = (1.0 / sat.rho_v_kg_m3 - 1.0 / sat.rho_l_kg_m3) * quality_gradient_1_m
    return mass_flux_kg_m2_s**2 * dv_dx


# ---------------------------------------------------------------------------
# Grant / Chisholm shell-side two-phase pressure-drop multiplier
# ---------------------------------------------------------------------------
# NOTE the name: this is DISTINCT from `chisholm_two_phase_multiplier(X, C)`
# above, which is the Lockhart-Martinelli C-form for flow INSIDE tubes. The
# function below is the shell-side crossflow form, keyed on bundle geometry.

# Grant's B coefficients for shell-side two-phase crossflow, by flow path.
# Source: Doo (2005), "A Modelling and Experimental Study of Evaporating
# Two-Phase Flow on the Shellside of Shell-and-Tube Heat Exchangers", Univ. of
# Strathclyde PhD, sec. 2.4, reporting Grant (1977) and Grant et al. (1986).
# Held at docs/reference/Doo2005.pdf (+ .md text extraction).
GRANT_SHELLSIDE_B = {
    "segmental_baffle_vertical": 1.0,      # vertical up-and-down between segmental baffles
    "ideal_bank_vertical_up": 3.0,
    "ideal_bank_horizontal_rotated_square": 0.6,
    "ideal_bank_horizontal_rotated_triangular": 0.35,
    "ideal_bank_horizontal_square": 0.28,
    "horizontal_side_to_side_spray": 0.75,      # spray / bubbly
    "horizontal_side_to_side_stratified": 0.35,  # stratified / stratified-spray
}


def chisholm_gamma(*, p_Pa: float, fluid: str) -> float:
    """Chisholm's physical-property coefficient Gamma.

        Gamma^2 = (dp/dz)_all-vapour / (dp/dz)_all-liquid
        Gamma   = sqrt(rho_l / rho_v) * (mu_v / mu_l)^0.1

    Doo (2005) eq. 2.54.
    """
    sat = saturation_state(fluid, p_Pa)
    return math.sqrt(sat.rho_l_kg_m3 / sat.rho_v_kg_m3) * (
        sat.mu_v_Pa_s / sat.mu_l_Pa_s
    ) ** 0.1


def grant_shellside_B(*, flow_path: str = "segmental_baffle_vertical",
                      gamma: float | None = None, n: float = 0.2) -> float:
    """Grant's B coefficient for the shell-side two-phase multiplier.

    ``flow_path`` selects from :data:`GRANT_SHELLSIDE_B`. The default,
    ``"segmental_baffle_vertical"`` (B = 1.0), is the case matching this
    codebase's segmentally-baffled geometry, where the coolant is driven
    vertically up and down across the bundle between successive baffles.

    ``flow_path=None`` falls back to Grant's general expression
    (Doo 2005 eq. 2.59), for geometries with no tabulated value:

        B = (2^(2-n) - 2) / (Gamma + 1)

    B here is GEOMETRIC. Chisholm's in-tube B table keyed on (Gamma, mass flux)
    is a different correlation family and must not be substituted.
    """
    if flow_path is not None:
        try:
            return GRANT_SHELLSIDE_B[flow_path]
        except KeyError:
            raise ValueError(
                "unknown flow_path %r; expected one of %s"
                % (flow_path, sorted(GRANT_SHELLSIDE_B))
            ) from None
    if gamma is None:
        raise ValueError("gamma is required when flow_path is None")
    return (2.0 ** (2.0 - n) - 2.0) / (gamma + 1.0)


def grant_chisholm_shellside_multiplier(
    *,
    p_Pa: float,
    quality: float,
    fluid: str = "Water",
    flow_path: str = "segmental_baffle_vertical",
    n: float = 0.2,
) -> float:
    """Shell-side two-phase multiplier phi^2, referenced to the ALL-LIQUID drop.

        phi^2 = 1 + (Gamma^2 - 1) * [ B x^((2-n)/2) (1-x)^((2-n)/2) + x^(2-n) ]

    Verified against the primary source: Chisholm, "Pressure gradients due to
    friction during the flow of evaporating two-phase mixtures in smooth tubes
    and channels", Int. J. Heat Mass Transfer 16 (1973) 347-358, eq. (26)
    (docs/reference/chisholm1973_2.pdf) -- identical, with q the dryness
    fraction. Also given as Doo (2005) eq. 2.55. Multiply an all-liquid
    shell-side pressure drop by this to get the two-phase drop. Recovers 1.0 at
    x=0 and Gamma^2 at x=1 exactly.

    ``n`` is the Reynolds exponent of the bundle's single-phase friction factor
    (Grant fitted 0.462 on his test geometry; 0.2, the smooth-turbulent value,
    is kept as the default here). ``flow_path`` selects B -- see
    :func:`grant_shellside_B`.
    """
    x = min(max(float(quality), 0.0), 1.0)
    gamma = chisholm_gamma(p_Pa=p_Pa, fluid=fluid)
    B = grant_shellside_B(flow_path=flow_path, gamma=gamma, n=n)
    e = (2.0 - n) / 2.0
    return 1.0 + (gamma ** 2 - 1.0) * (
        B * x ** e * (1.0 - x) ** e + x ** (2.0 - n)
    )


def chisholm_intube_B(*, gamma: float, mass_flux_kg_m2_s: float) -> float:
    """Chisholm's B for flow INSIDE TUBES, from the Baroczy correlation.

    Chisholm (1973), Int. J. Heat Mass Transfer 16, 347-358, eqs. (31)-(33)
    (docs/reference/chisholm1973_2.pdf). Mass velocity G must be in kg/m2s::

        Gamma < 9.5        B = 55 / sqrt(G)                 (31)
        9.5 < Gamma < 28   B = 520 / (Gamma * sqrt(G))      (32)
        28 < Gamma         B = 15000 / (Gamma^2 * sqrt(G))  (33)

    **Do not use this for a baffled shell.** It is the in-tube family; the
    shell-side multiplier takes a GEOMETRIC B instead -- see
    :func:`grant_shellside_B`. Provided for tube-side work and for comparison.

    The paper notes that for G > 1900 kg/m2s with Gamma < 9.5, eq. (31) returns
    values below Chisholm's own earlier correlation.
    """
    G = max(float(mass_flux_kg_m2_s), 1.0e-12)
    if gamma < 9.5:
        return 55.0 / math.sqrt(G)
    if gamma < 28.0:
        return 520.0 / (gamma * math.sqrt(G))
    return 15000.0 / (gamma ** 2 * math.sqrt(G))

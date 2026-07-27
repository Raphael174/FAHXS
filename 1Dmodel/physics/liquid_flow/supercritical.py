"""Supercritical-pressure forced-convection heat-transfer closures.

First client of the closure registry (registry.py). Every correlation here is a
property-ratio-corrected Nusselt form: a base turbulent Nu times wall-to-bulk
property ratios that capture the steep property variation through the pseudo-
critical region. All are registered at import.

Sourcing discipline (same as the rest of this package): coefficients are only
hardcoded after visual verification against the rendered source page, because
plain-text PDF extraction silently drops minus signs in exponents. The verified
forms and their provenance:

- ``mccarthy_wolf_nu`` / ``taylor_nu`` -- the rocket-regenerative-cooling
  property-ratio family (McCarthy-Wolf 1960, Taylor 1968/NERVA), as tabulated
  and compared against 2992 supercritical-H2 heated-tube measurements by Locke &
  Landrum (2008), J. Propulsion & Power 24(1):94-103
  (docs/reference/2008_Locke_Landrum_...md, "VERIFIED CORRELATIONS"). These are
  the PRIMARY closures for this project's LN2 case: rocket-cooling channels run
  at Re up to >1e6, so the LN2 case Re of 1e5-6e5 is squarely inside their
  validated envelope (unlike Cheng2020's 7000-27000). Fit on H2 but explicitly
  fluid-transferable via the (T_s/T_b) property-ratio mechanism, hence
  FLUID_ANY. Locke & Landrum's key caveat is honored downstream: these
  bulk-reference forms OVERPREDICT (spike) around the pseudo-critical line, so
  the pseudo_critical regime is flagged low-confidence regardless of closure.
- ``cheng2020_supercritical_nu`` -- Cheng et al. (2020), N2-specific, vertical
  20 mm tube, 7.5-9.0 MPa (P/Pc 2.21-2.65). Now a NARROW-niche closure: it wins
  selection only inside its own 7000-27000 Re window (low-flow N2), where being
  N2-specific beats the generic rocket-cooling forms. Above that Re (our helical
  case) the registry prefers the in-range McCarthy-Wolf/Taylor.
- ``krasnoshchekov_protopopov_nu`` -- K-P (Ting Part II Eq. 9-10, Nu_0 =
  Petukhov-Kirillov Eq. 1). Generic property-ratio fallback; CO2/water fit
  envelope, tier structural_extrapolation.
- conservative Gnielinski-at-bulk bound (``_conservative_bound_htc``) -- always
  in range, the guaranteed last-resort fallback.

- ``wang2023_eckert_split_nu`` -- Wang et al. (2023), N2-specific, 4.57mm
  vertical tube, 3.5 MPa (P/Pc~1.03-1.1). Verified from the rendered page
  (docs/reference/Wang2023.pdf, Eq. 8-10) -- signs matched the original text
  extraction exactly this time, no drops. Registered but validity is
  ``p_reduced in (1.0, 1.1)``, so it NEVER wins at our 80 bar target
  (P/Pc~2.36); kept as a niche near-critical cross-check only, same pattern as
  Cheng2020's own Re niche.

Cook's methane correlation (Urbano-Nasuti 2013 Eq. 15,
Nu=0.0022 Re^0.8 Pr^0.4 (T_b/T_w)^0.45) is the same property-ratio family as
McCarthy-Wolf and adds nothing beyond it for N2; noted for reference only, not
registered.
"""

from __future__ import annotations

import math

from hps_combustor.physics.liquid_flow.coolprop_state_cache import get_cached_state
from hps_combustor.physics.liquid_flow.correlations import (
    darcy_friction_smooth_pipe,
    liquid_single_phase_nusselt,
)
from hps_combustor.physics.liquid_flow.registry import (
    FLUID_ANY,
    TIER_CONSERVATIVE_BOUND,
    TIER_STRUCTURAL_EXTRAPOLATION,
    TIER_VALIDATED_IN_RANGE,
    ClosureContext,
    ClosureRecord,
    register,
)
from hps_combustor.physics.liquid_flow.regime import (
    PSEUDO_CRITICAL,
    SUPERCRITICAL_GAS_LIKE,
    SUPERCRITICAL_LIQUID_LIKE,
)

_ALL_SUPERCRITICAL = frozenset(
    {SUPERCRITICAL_LIQUID_LIKE, PSEUDO_CRITICAL, SUPERCRITICAL_GAS_LIKE}
)


def _bulk_reynolds(ctx: ClosureContext) -> float:
    return ctx.mass_flux_kg_m2_s * ctx.diameter_m / ctx.mu_b


def _wall_properties(ctx: ClosureContext):
    """(mu_w, k_w, h_w) at the lagged wall temperature, or None if unavailable
    (first node/sweep) -- callers fall back to unit property ratios then."""
    if ctx.wall_temp_K is None:
        return None
    st = get_cached_state(ctx.fluid).wall_state_tp(ctx.wall_temp_K, ctx.p_Pa)
    return float(st.viscosity()), float(st.conductivity()), float(st.hmass())


def petukhov_kirillov_nu0(Re_b: float, Pr_b: float) -> float:
    """Base turbulent Nu_0 (Petukhov-Kirillov, Ting Part II Eq. 1) using the
    Filonenko friction factor. darcy_friction_smooth_pipe returns exactly the
    Filonenko form 1/(0.79 ln Re - 1.64)^2 = (1.82 log10 Re - 1.64)^-2 for
    turbulent Re, so it is reused here for xi."""
    xi = darcy_friction_smooth_pipe(Re_b)
    num = (xi / 8.0) * Re_b * Pr_b
    den = 12.7 * math.sqrt(xi / 8.0) * (Pr_b ** (2.0 / 3.0) - 1.0) + 1.07
    return num / den


def mccarthy_wolf_nu(ctx: ClosureContext) -> float:
    """McCarthy-Wolf (1960) supercritical-cryogen HTC (Locke-Landrum Eq. 2):

        Nu_b = 0.025 Re_b^0.8 Pr_b^0.4 (T_s/T_b)^(-0.55)

    bulk properties, T_s = wall (surface) temperature. The (T_s/T_b) ratio
    carries the property-variation physics, making the H2-fit form transferable
    to other supercritical cryogens. Falls back to the ratio=1 (Dittus-Boelter-
    like) value when the wall temperature is unavailable (first node/sweep).
    Returns HTC = Nu_b k_b / d.
    """
    Re_b = _bulk_reynolds(ctx)
    ratio = 1.0 if ctx.wall_temp_K is None else (ctx.wall_temp_K / ctx.T_bulk_K) ** (-0.55)
    Nu_b = 0.025 * Re_b**0.8 * ctx.Pr_b**0.4 * ratio
    return Nu_b * ctx.k_b / ctx.diameter_m


def taylor_nu(ctx: ClosureContext) -> float:
    """Taylor (1968, NERVA) supercritical-cryogen HTC (Locke-Landrum Eq. 8):

        Nu_b = 0.023 Re_b^0.8 Pr_b^0.4 (T_s/T_b)^(-(0.57 - 1.59/(x/D)))

    An axially-varying property-ratio exponent that adds an entrance correction
    (x/D = flow length from heating start over diameter). For a long channel the
    1.59/(x/D) term vanishes and the exponent -> -0.57 (near McCarthy-Wolf's
    -0.55). ``ctx.x_over_D`` supplies the local x/D; when absent, the fully
    developed limit (-0.57) is used. Better than McCarthy-Wolf for low x/D and
    reduced pressure > 4.
    """
    Re_b = _bulk_reynolds(ctx)
    if ctx.wall_temp_K is None:
        ratio = 1.0
    else:
        if ctx.x_over_D is not None and ctx.x_over_D > 0.0:
            exponent = -(0.57 - 1.59 / ctx.x_over_D)
        else:
            exponent = -0.57
        ratio = (ctx.wall_temp_K / ctx.T_bulk_K) ** exponent
    Nu_b = 0.023 * Re_b**0.8 * ctx.Pr_b**0.4 * ratio
    return Nu_b * ctx.k_b / ctx.diameter_m


def wang2023_eckert_split_nu(ctx: ClosureContext) -> float:
    """Wang et al. (2023) supercritical-N2 HTC, verified Eq. 8-10:

        Pr_b_bar = (h_w - h_b)/(T_w - T_b) * mu_b/k_b     (effective Prandtl)
        E = (T_crit - T_b) / (T_w - T_b)                  (Eckert regime split,
                                                             T_crit = TRUE
                                                             critical temperature,
                                                             not the
                                                             pseudo-critical T_pc)

        0<=E<=1: Nu_b = 104.85 Re_b^0.26 Pr_b_bar^-0.083 (rho_w/rho_b)^-0.013
                        (mu_w/mu_b)^1.02 (k_w/k_b)^1.39
        E<0:     Nu_b = 124.34 Re_b^0.02  Pr_b_bar^-0.16  (rho_w/rho_b)^0.63
                        (mu_w/mu_b)^-1.05 (k_w/k_b)^0.75

    E>1 (cold liquid phase) is outside the paper's studied range; this
    implementation defaults to the 0<=E<=1 branch there rather than
    extrapolating into a regime the paper explicitly excludes. Falls back to
    the paper's own pseudo-critical-adjacent default (E=0.5, i.e. the first
    branch) when no wall temperature is available yet.
    """
    Re_b = _bulk_reynolds(ctx)
    wall = _wall_properties(ctx)
    if wall is None or ctx.wall_temp_K == ctx.T_bulk_K:
        rho_ratio = mu_ratio = k_ratio = 1.0
        Pr_bar = ctx.Pr_b
        E = 0.5
    else:
        mu_w, k_w, h_w = wall
        rho_w = get_cached_state(ctx.fluid).wall_state_tp(ctx.wall_temp_K, ctx.p_Pa).rhomass()
        Pr_bar = (h_w - ctx.h_J_kg) / (ctx.wall_temp_K - ctx.T_bulk_K) * ctx.mu_b / ctx.k_b
        Pr_bar = max(Pr_bar, 1.0e-6)
        rho_ratio = rho_w / ctx.rho_b
        mu_ratio = mu_w / ctx.mu_b
        k_ratio = k_w / ctx.k_b
        T_crit = get_cached_state(ctx.fluid).T_crit_K
        E = (T_crit - ctx.T_bulk_K) / (ctx.wall_temp_K - ctx.T_bulk_K)

    if E < 0.0:
        Nu_b = (
            124.34 * Re_b**0.02 * Pr_bar**(-0.16)
            * rho_ratio**0.63 * mu_ratio**(-1.05) * k_ratio**0.75
        )
    else:  # 0 <= E <= 1 branch (also used as the E>1/no-wall-temp default)
        Nu_b = (
            104.85 * Re_b**0.26 * Pr_bar**(-0.083)
            * rho_ratio**(-0.013) * mu_ratio**1.02 * k_ratio**1.39
        )
    return Nu_b * ctx.k_b / ctx.diameter_m


def cheng2020_supercritical_nu(ctx: ClosureContext) -> float:
    """Cheng et al. (2020) supercritical-N2 HTC. Verified Eq. 18-19/16:

        Nu_b = 0.023 Re_b^0.9306 Pr_b^1.3873 (Gr*)^(-0.0013) c_b
        c_b  = (mu_b / mu_w)^(-0.8709)
        Gr*  = g alpha_b q d^4 / (k_b nu_b^2),  nu_b = mu_b / rho_b

    All bulk except mu_w. The Gr* exponent is tiny (-0.0013) so the buoyancy
    term is near-unity and q=0 is harmless once clamped; alpha_b is the isobaric
    expansion coefficient at bulk. Returns HTC = Nu_b k_b / d.
    """
    Re_b = _bulk_reynolds(ctx)
    nu_b = ctx.mu_b / ctx.rho_b
    alpha_b = abs(get_cached_state(ctx.fluid).flash_ph(
        ctx.p_Pa, ctx.h_J_kg).isobaric_expansion_coefficient())
    q = max(ctx.heat_flux_W_m2, 1.0e-6)
    gr_star = 9.80665 * alpha_b * q * ctx.diameter_m**4 / (ctx.k_b * nu_b**2)
    gr_star = max(gr_star, 1.0e-12)

    wall = _wall_properties(ctx)
    c_b = 1.0 if wall is None else (ctx.mu_b / wall[0]) ** (-0.8709)

    Nu_b = (
        0.023
        * Re_b**0.9306
        * ctx.Pr_b**1.3873
        * gr_star ** (-0.0013)
        * c_b
    )
    return Nu_b * ctx.k_b / ctx.diameter_m


def krasnoshchekov_protopopov_nu(ctx: ClosureContext) -> float:
    """Krasnoshchekov-Protopopov generic supercritical HTC (Ting Part II Eq. 9):

        Nu = Nu_0 (mu_w/mu_b)^(-0.11) (k_w/k_b)^0.33 (cbar_p/cp_b)^0.35
        cbar_p = (h_w - h_b) / (T_w - T_b)

    with Nu_0 the Petukhov-Kirillov base. Falls back to unit property ratios and
    cbar_p = cp_b when the wall state is unavailable (first node/sweep).
    """
    Re_b = _bulk_reynolds(ctx)
    Nu_0 = petukhov_kirillov_nu0(Re_b, ctx.Pr_b)
    wall = _wall_properties(ctx)
    if wall is None or ctx.wall_temp_K == ctx.T_bulk_K:
        factor = 1.0
    else:
        mu_w, k_w, h_w = wall
        cbar_p = (h_w - ctx.h_J_kg) / (ctx.wall_temp_K - ctx.T_bulk_K)
        # cbar_p can go negative if the (h,T) pair straddles a nonmonotonic
        # region; clamp to a small positive so the fractional power is defined.
        cp_ratio = max(cbar_p, 1.0e-6) / ctx.cp_b
        factor = (
            (mu_w / ctx.mu_b) ** (-0.11)
            * (k_w / ctx.k_b) ** 0.33
            * cp_ratio**0.35
        )
    return Nu_0 * factor * ctx.k_b / ctx.diameter_m


def _conservative_bound_htc(ctx: ClosureContext) -> float:
    """Plain Gnielinski at bulk properties -- no property-ratio correction. The
    guaranteed fallback: any fluid, any geometry, always in range. Under-
    predicts near the pseudo-critical HTC peak (a conservative bound, not an
    accurate closure), which is why it is the lowest-tier last resort."""
    Re_b = _bulk_reynolds(ctx)
    Nu = liquid_single_phase_nusselt(Re_b, ctx.Pr_b)
    return Nu * ctx.k_b / ctx.diameter_m


# --- registration -----------------------------------------------------------
# P/Pc ranges use "p_reduced"; Re uses "Re_b"; these keys match the operating_
# point dict dispatch.py assembles. Geometry: all are straight-tube fits;
# helical_coil is tagged as an accepted mild extrapolation (locally tube-like,
# secondary-flow enhancement modest and unquantified for a plain circular coil);
# shell_crossflow is deliberately NOT tagged (no validated supercritical shell-
# side closure exists -- an honest gap, flagged at dispatch if requested).
# priority breaks ties between equally-ranked closures (higher wins).

register(ClosureRecord(
    name="mccarthy_wolf",
    regime_tags=_ALL_SUPERCRITICAL,
    geometry_tags=frozenset({"straight_tube", "helical_coil"}),
    orientation_tags=frozenset({"vertical", "horizontal", "any"}),
    fluid_scope=frozenset({FLUID_ANY}),
    # Rocket-cooling-channel Re envelope (Locke-Landrum data span; up to >1e6).
    validity={"Re_b": (1.0e4, 1.0e7)},
    provenance="Locke-Landrum 2008 Eq. 2 (McCarthy-Wolf 1960), verified",
    tier=TIER_VALIDATED_IN_RANGE,
    callable=mccarthy_wolf_nu,
    priority=2,  # canonical baseline: default pick among the generic in-range forms
))

register(ClosureRecord(
    name="taylor_nerva",
    regime_tags=_ALL_SUPERCRITICAL,
    geometry_tags=frozenset({"straight_tube", "helical_coil"}),
    orientation_tags=frozenset({"vertical", "horizontal", "any"}),
    fluid_scope=frozenset({FLUID_ANY}),
    validity={"Re_b": (1.0e4, 1.0e7)},
    provenance="Locke-Landrum 2008 Eq. 8 (Taylor 1968/NERVA), verified; best for P_r>4, low x/D",
    tier=TIER_VALIDATED_IN_RANGE,
    callable=taylor_nu,
    priority=1,  # entrance-corrected alternative to McCarthy-Wolf
))

register(ClosureRecord(
    name="wang2023_eckert_split",
    regime_tags=_ALL_SUPERCRITICAL,
    geometry_tags=frozenset({"straight_tube", "helical_coil"}),
    orientation_tags=frozenset({"vertical"}),
    fluid_scope=frozenset({"Nitrogen"}),
    validity={"p_reduced": (1.0, 1.1)},
    provenance="Wang2023 (Int. J. Thermal Sci. 184, 108001), verified Eq. 8-10",
    tier=TIER_VALIDATED_IN_RANGE,
    callable=wang2023_eckert_split_nu,
    priority=3,  # matches cheng2020: wins its own narrow (near-critical) niche
))

register(ClosureRecord(
    name="cheng2020_supercritical",
    regime_tags=_ALL_SUPERCRITICAL,
    geometry_tags=frozenset({"straight_tube", "helical_coil"}),
    orientation_tags=frozenset({"vertical"}),
    fluid_scope=frozenset({"Nitrogen"}),
    validity={"p_reduced": (2.21, 2.65), "Re_b": (7000.0, 27000.0)},
    provenance="cheng2020 (Int. J. Thermal Sci. 152, 106327), verified Eq. 18-19",
    tier=TIER_VALIDATED_IN_RANGE,
    callable=cheng2020_supercritical_nu,
    priority=3,  # highest so it wins its own niche (in-range fluid-specific N2)
))

register(ClosureRecord(
    name="krasnoshchekov_protopopov",
    regime_tags=_ALL_SUPERCRITICAL,
    geometry_tags=frozenset({"straight_tube", "helical_coil"}),
    orientation_tags=frozenset({"vertical", "horizontal", "any"}),
    fluid_scope=frozenset({FLUID_ANY}),
    validity={"Re_b": (20000.0, 860000.0), "Pr_b": (0.85, 65.0)},
    provenance="Ting Part II Eq. 9-10 (K-P), Nu_0 = Petukhov-Kirillov Eq. 1",
    tier=TIER_STRUCTURAL_EXTRAPOLATION,
    callable=krasnoshchekov_protopopov_nu,
))

register(ClosureRecord(
    name="gnielinski_bulk_bound",
    regime_tags=_ALL_SUPERCRITICAL,
    geometry_tags=frozenset({"straight_tube", "helical_coil", "shell_crossflow"}),
    orientation_tags=frozenset({"vertical", "horizontal", "any"}),
    fluid_scope=frozenset({FLUID_ANY}),
    validity={},  # no property-ratio correction claimed -> nothing to violate
    provenance="Gnielinski single-phase at bulk properties (conservative bound)",
    tier=TIER_CONSERVATIVE_BOUND,
    callable=_conservative_bound_htc,
))

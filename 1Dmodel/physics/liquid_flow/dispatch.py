"""Configuration-agnostic coolant state and correlation dispatch.

This module is the integration layer between HX solvers and coolant physics.
The maintained solvers still use their legacy helium-specific calls by default;
new liquid-coolant work should enter through this dispatcher so the same state
closure and correlations can be shared by straight-pipe, helical, shell-side,
steady, and transient validation paths.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import CoolProp.CoolProp as CP

from hps_combustor.core.thermo import (
    ThermoState as CoolantState,
    coolant_state_from_Tp,
    coolant_state_from_ph,
    coolant_inlet_state,
)
from hps_combustor.physics.liquid_flow.chf import groeneveld_2006_chf
from hps_combustor.physics.liquid_flow.coolprop_state_cache import (
    coolprop_fluid_string,
    get_cached_state,
)
from hps_combustor.physics.liquid_flow.correlations import (
    bergles_rohsenow_onb_wall_superheat,
    darcy_friction_smooth_pipe,
    equilibrium_state_ph,
    gungor_winterton_boiling_htc,
    liquid_single_phase_nusselt,
    muller_steinhagen_heck_friction_gradient,
    post_chf_dispersed_flow_htc,
    saturation_state,
    two_phase_sound_speed,
)
from hps_combustor.physics.liquid_flow import regime as _regime
from hps_combustor.physics.liquid_flow import supercritical as _supercritical  # noqa: F401 (registers closures)
from hps_combustor.physics.liquid_flow.registry import (
    ClosureContext,
    ExtrapolationReport,
    select_supercritical,
)

# CoolantState, coolant_state_from_Tp, coolant_state_from_ph, and
# coolant_inlet_state now live in hps_combustor.core.thermo (Stage A of
# docs/solver_design/FV_CORE_REWORK_PLAN.md) and are re-exported here
# unchanged so every existing importer of this module keeps working. This is
# a pure relocation: the CoolProp call sequences are identical to what used
# to be defined inline in this file.

# Quality window, centered on quality=0 (saturated-liquid boiling onset), over
# which the subcooled single-phase-liquid closure and the saturated two-phase
# closure are blended instead of hard-switched. Purely a numerical smoothing
# choice (removes a spurious ~5x HTC jump in one march node) - not a physical
# onset-of-nucleate-boiling (ONB) criterion (see
# ``bergles_rohsenow_onb_wall_superheat`` for the actual physical gate), and
# not applied at the quality=1 (complete-vaporization) boundary.
BOILING_ONSET_BLEND_HALF_WIDTH = 0.02


@dataclass(frozen=True)
class CoolantClosureResult:
    state: CoolantState
    htc_W_m2_K: float
    dpdz_friction_Pa_m: float
    chf_W_m2: float | None
    chf_margin: float | None
    sound_speed_m_s: float | None = None
    # Bergles-Rohsenow ONB wall-superheat margin [K]: (estimated wall
    # superheat) - (ONB threshold superheat at this heat flux/pressure).
    # Positive => nucleate boiling is physically expected even though bulk
    # quality is still subcooled (state.quality < 0). None when heat flux is
    # non-positive or outside the subcooled branch (not evaluated).
    onb_wall_superheat_margin_K: float | None = None
    # Regime/closure provenance (populated in the supercritical branch; None on
    # the subcritical path, which keeps its inline closure selection).
    regime: str | None = None
    closure_name: str | None = None
    extrapolation_report: ExtrapolationReport | None = None
    # Supercritical heat-transfer-deterioration (HTD) onset flag + the buoyancy
    # parameter it was judged on (Hall-Jackson). Detection only -- no degraded-
    # HTC magnitude model (see regime.buoyancy_parameter docstring).
    htd_risk: bool = False
    buoyancy_parameter: float | None = None
    # Diagnostic-only cross-check HTC from an alternate registered closure
    # (currently: RPE Eq. 8-24 Dittus-Boelter/Colburn vs. the active
    # Gnielinski subcooled-liquid HTC) -- reported alongside htc_W_m2_K,
    # NEVER used in its place. Populated only in the subcooled-liquid branch;
    # None everywhere else. See correlations.py:dittus_boelter_colburn_liquid_
    # nusselt's docstring for measured divergence from the active closure.
    cross_check_closure_name: str | None = None
    cross_check_htc_W_m2_K: float | None = None


def _subcooled_liquid_cross_check(
    *,
    state: CoolantState,
    mass_flux_kg_m2_s: float,
    hydraulic_diameter_m: float,
    heat_flux_W_m2: float,
    wall_temp_K: float | None,
    geometry: str,
    orientation: str,
    x_over_D: float | None,
) -> tuple[str | None, float | None]:
    """Diagnostic-only alternate HTC for the subcooled-liquid branch (never
    used to compute the returned htc_W_m2_K). Returns (closure_name, htc) or
    (None, None) if nothing is registered / applicable -- must never raise,
    this is a reporting sidecar, not part of the solved physics."""
    try:
        ctx = ClosureContext(
            fluid=state.fluid,
            p_Pa=state.p_Pa,
            h_J_kg=state.h_J_kg,
            T_bulk_K=state.T_K,
            rho_b=state.rho_kg_m3,
            mu_b=state.mu_Pa_s,
            k_b=state.k_W_m_K,
            cp_b=state.cp_J_kg_K,
            Pr_b=state.Pr,
            mass_flux_kg_m2_s=mass_flux_kg_m2_s,
            diameter_m=hydraulic_diameter_m,
            heat_flux_W_m2=heat_flux_W_m2,
            wall_temp_K=wall_temp_K,
            x_over_D=x_over_D,
        )
        record, _report = select_supercritical(
            regime="subcooled_liquid",
            geometry=geometry,
            orientation=orientation,
            fluid=state.fluid,
            operating_point={
                "Re_b": mass_flux_kg_m2_s * hydraulic_diameter_m / state.mu_Pa_s,
                "Pr_b": state.Pr,
            },
        )
        return record.name, float(record.callable(ctx))
    except (LookupError, ValueError, ZeroDivisionError):
        return None, None


def evaluate_coolant_closure(
    *,
    coolant_prop,
    p_Pa: float,
    h_J_kg: float,
    mass_flux_kg_m2_s: float,
    hydraulic_diameter_m: float,
    heat_flux_W_m2: float,
    lut_path: str | Path | None = None,
    wall_temp_K: float | None = None,
    geometry: str = "straight_tube",
    orientation: str = "vertical",
    x_over_D: float | None = None,
) -> CoolantClosureResult:
    """Evaluate state, cold-side HTC, pressure drop, and optional CHF margin.

    ``heat_flux_W_m2`` is the wall heat flux into the coolant. For the
    equilibrium_liquid model this covers subcritical liquid/boiling (the
    validated water closures) AND supercritical-pressure forced convection
    (``supercritical.py`` closures via the registry) -- the branch is chosen
    from the state's pressure relative to p_crit, transparently. ``wall_temp_K``
    (lagged coolant-side wall temperature), ``geometry``, ``orientation`` and
    ``x_over_D`` (flow length / diameter from heating start, for Taylor's
    entrance correction) feed the supercritical property-ratio corrections and
    closure selection; they are ignored on the subcritical and single-phase
    paths.
    """
    model = getattr(coolant_prop, "coolant_model", "single_phase_coolprop")
    fluid = getattr(coolant_prop, "coolant", "Helium")
    backend = getattr(coolant_prop, "liquid_property_backend", "HEOS")
    state = coolant_state_from_ph(fluid, p_Pa, h_J_kg, model, backend=backend)

    if model == "equilibrium_liquid" and state.is_supercritical:
        return _supercritical_closure(
            coolant_prop=coolant_prop,
            state=state,
            fluid=fluid,
            backend=backend,
            p_Pa=p_Pa,
            h_J_kg=h_J_kg,
            mass_flux_kg_m2_s=mass_flux_kg_m2_s,
            hydraulic_diameter_m=hydraulic_diameter_m,
            heat_flux_W_m2=heat_flux_W_m2,
            wall_temp_K=wall_temp_K,
            geometry=geometry,
            orientation=orientation,
            x_over_D=x_over_D,
        )

    if model == "single_phase_coolprop":
        if hydraulic_diameter_m <= 0.0 or mass_flux_kg_m2_s <= 0.0:
            raise ValueError("hydraulic diameter and mass flux must be positive")
        Re = mass_flux_kg_m2_s * hydraulic_diameter_m / state.mu_Pa_s
        f = darcy_friction_smooth_pipe(Re)
        Nu = liquid_single_phase_nusselt(Re, state.Pr, f)
        htc = Nu * state.k_W_m_K / hydraulic_diameter_m
        dpdz = f * mass_flux_kg_m2_s**2 / (2.0 * hydraulic_diameter_m * state.rho_kg_m3)
        c = get_cached_state(fluid).flash_ph(p_Pa, h_J_kg).speed_sound()
        return CoolantClosureResult(
            state=state,
            htc_W_m2_K=float(htc),
            dpdz_friction_Pa_m=float(dpdz),
            chf_W_m2=None,
            chf_margin=None,
            sound_speed_m_s=float(c),
        )

    if model != "equilibrium_liquid":
        raise ValueError(f"unknown coolant model: {model!r}")

    # Backend-tagged fluid string (e.g. "BICUBIC::Water") for every
    # CoolProp/correlations call below - opt-in via coolant_prop.
    # liquid_property_backend, default "HEOS" (exact, unchanged behavior).
    fluid = coolprop_fluid_string(fluid, backend)

    def _single_phase_htc_dpdz():
        Re = mass_flux_kg_m2_s * hydraulic_diameter_m / state.mu_Pa_s
        f = darcy_friction_smooth_pipe(Re)
        Nu = liquid_single_phase_nusselt(Re, state.Pr, f)
        htc = Nu * state.k_W_m_K / hydraulic_diameter_m
        dpdz = f * mass_flux_kg_m2_s**2 / (2.0 * hydraulic_diameter_m * state.rho_kg_m3)
        return htc, dpdz

    if state.quality > 1.0:
        htc, dpdz = _single_phase_htc_dpdz()
        c = get_cached_state(fluid).flash_ph(p_Pa, h_J_kg).speed_sound()
        return CoolantClosureResult(
            state=state,
            htc_W_m2_K=float(htc),
            dpdz_friction_Pa_m=float(dpdz),
            chf_W_m2=None,
            chf_margin=None,
            sound_speed_m_s=float(c),
        )

    if state.quality < -BOILING_ONSET_BLEND_HALF_WIDTH:
        htc, dpdz = _single_phase_htc_dpdz()
        c = get_cached_state(fluid).flash_ph(p_Pa, h_J_kg).speed_sound()
        chf = None
        margin = None
        chf_model = getattr(coolant_prop, "liquid_chf_model", "groeneveld_2006")
        if chf_model == "groeneveld_2006" and lut_path is not None:
            try:
                chf = groeneveld_2006_chf(
                    p_Pa=p_Pa,
                    mass_flux_kg_m2_s=mass_flux_kg_m2_s,
                    quality=max(state.quality, -0.5),
                    diameter_m=hydraulic_diameter_m,
                    lut_path=lut_path,
                )
            except ValueError:
                chf = None
            margin = chf / heat_flux_W_m2 if (chf is not None and heat_flux_W_m2 > 0.0) else None
        elif chf_model not in ("none", "groeneveld_2006"):
            raise ValueError(f"unsupported liquid CHF model: {chf_model!r}")

        onb_margin = None
        if heat_flux_W_m2 > 0.0 and htc > 0.0:
            sat_local = saturation_state(fluid, p_Pa)
            T_wall_est = state.T_K + heat_flux_W_m2 / htc
            dT_wall_superheat = T_wall_est - sat_local.T_sat_K
            dT_onb = bergles_rohsenow_onb_wall_superheat(p_Pa=p_Pa, heat_flux_W_m2=heat_flux_W_m2)
            onb_margin = dT_wall_superheat - dT_onb

        cross_check_name, cross_check_htc = _subcooled_liquid_cross_check(
            state=state,
            mass_flux_kg_m2_s=mass_flux_kg_m2_s,
            hydraulic_diameter_m=hydraulic_diameter_m,
            heat_flux_W_m2=heat_flux_W_m2,
            wall_temp_K=wall_temp_K,
            geometry=geometry,
            orientation=orientation,
            x_over_D=x_over_D,
        )

        return CoolantClosureResult(
            state=state,
            htc_W_m2_K=float(htc),
            dpdz_friction_Pa_m=float(dpdz),
            chf_W_m2=None if chf is None else float(chf),
            chf_margin=None if margin is None else float(margin),
            sound_speed_m_s=float(c),
            onb_wall_superheat_margin_K=onb_margin,
            cross_check_closure_name=cross_check_name,
            cross_check_htc_W_m2_K=cross_check_htc,
        )

    # -BOILING_ONSET_BLEND_HALF_WIDTH <= quality <= 1.0: saturated two-phase
    # region (possibly blended against the subcooled closure near quality=0).
    dp_model = getattr(coolant_prop, "liquid_pressure_drop_model", "muller_steinhagen_heck")
    if dp_model != "muller_steinhagen_heck":
        raise ValueError(f"unsupported liquid pressure-drop model: {dp_model!r}")
    dpdz_two_phase = muller_steinhagen_heck_friction_gradient(
        p_Pa=p_Pa,
        mass_flux_kg_m2_s=mass_flux_kg_m2_s,
        diameter_m=hydraulic_diameter_m,
        quality=state.quality,
        fluid=fluid,
    )

    # CHF margin must be known BEFORE picking the HTC correlation: if CHF is
    # exceeded (margin < 1), nucleate/convective boiling has broken down (DNB
    # or dryout) and Gungor-Winterton's optimistic HTC no longer applies.
    chf = None
    margin = None
    chf_model = getattr(coolant_prop, "liquid_chf_model", "groeneveld_2006")
    if chf_model == "groeneveld_2006" and lut_path is not None:
        try:
            chf = groeneveld_2006_chf(
                p_Pa=p_Pa,
                mass_flux_kg_m2_s=mass_flux_kg_m2_s,
                quality=state.quality,
                diameter_m=hydraulic_diameter_m,
                lut_path=lut_path,
            )
        except ValueError:
            # Lookup point outside the LUT's tabulated range (pressure,
            # mass flux, or quality) - CHF margin unavailable at this node,
            # not a fatal condition. Falls back to the boiling HTC below.
            chf = None
        margin = chf / heat_flux_W_m2 if (chf is not None and heat_flux_W_m2 > 0.0) else None
    elif chf_model not in ("none", "groeneveld_2006"):
        raise ValueError(f"unsupported liquid CHF model: {chf_model!r}")

    if margin is not None and margin < 1.0:
        # Post-CHF: see post_chf_dispersed_flow_htc's docstring for the
        # conservative-simplification caveat (treats the coolant as 100%
        # vapor at the full mass flux; not a validated dispersed-flow
        # correlation).
        htc_two_phase = post_chf_dispersed_flow_htc(
            p_Pa=p_Pa,
            mass_flux_kg_m2_s=mass_flux_kg_m2_s,
            diameter_m=hydraulic_diameter_m,
            fluid=fluid,
        )
    else:
        htc_model = getattr(coolant_prop, "liquid_heat_transfer_model", "gungor_winterton")
        if htc_model != "gungor_winterton":
            raise ValueError(f"unsupported liquid heat-transfer model: {htc_model!r}")
        # gungor_winterton_boiling_htc internally clamps quality to
        # [1e-9, 1-1e-9], so calling it with slightly negative quality
        # (inside the blend window) is safe - it degrades to the x->0 limit.
        htc_two_phase = gungor_winterton_boiling_htc(
            p_Pa=p_Pa,
            mass_flux_kg_m2_s=mass_flux_kg_m2_s,
            diameter_m=hydraulic_diameter_m,
            quality=state.quality,
            heat_flux_W_m2=max(heat_flux_W_m2, 0.0),
            fluid=fluid,
        )

    c_two_phase = two_phase_sound_speed(
        p_Pa=p_Pa,
        void_fraction=state.void_fraction,
        rho_mix_kg_m3=state.rho_kg_m3,
        fluid=fluid,
    )

    if state.quality < BOILING_ONSET_BLEND_HALF_WIDTH:
        htc_liquid, dpdz_liquid = _single_phase_htc_dpdz()
        # NOT CP.PropsSI(..., 'P', p_Pa, 'H', h_J_kg, ...): the actual (p, h)
        # state here can be genuinely two-phase (0 <= quality < blend width),
        # and CoolProp correctly refuses a (P,H) SPEED_OF_SOUND query inside
        # the dome (sound speed is not uniquely defined there via the real
        # EOS - that is exactly why Wood's equation exists). Use the
        # saturated-liquid limit (Q=0, always single-phase and well-posed)
        # as the blend's liquid-side anchor instead.
        c_liquid = get_cached_state(fluid).saturated_liquid(p_Pa).speed_sound()
        t = (state.quality + BOILING_ONSET_BLEND_HALF_WIDTH) / (2.0 * BOILING_ONSET_BLEND_HALF_WIDTH)
        t = min(max(t, 0.0), 1.0)
        weight = t * t * (3.0 - 2.0 * t)  # smoothstep: zero slope at both ends
        htc = (1.0 - weight) * htc_liquid + weight * htc_two_phase
        dpdz = (1.0 - weight) * dpdz_liquid + weight * dpdz_two_phase
        c = (1.0 - weight) * c_liquid + weight * c_two_phase
    else:
        htc = htc_two_phase
        dpdz = dpdz_two_phase
        c = c_two_phase

    return CoolantClosureResult(
        state=state,
        htc_W_m2_K=float(htc),
        dpdz_friction_Pa_m=float(dpdz),
        chf_W_m2=None if chf is None else float(chf),
        chf_margin=None if margin is None else float(margin),
        sound_speed_m_s=float(c),
    )


def _supercritical_closure(
    *,
    coolant_prop,
    state: CoolantState,
    fluid: str,
    backend: str,
    p_Pa: float,
    h_J_kg: float,
    mass_flux_kg_m2_s: float,
    hydraulic_diameter_m: float,
    heat_flux_W_m2: float,
    wall_temp_K: float | None,
    geometry: str,
    orientation: str,
    x_over_D: float | None = None,
) -> CoolantClosureResult:
    """Supercritical-pressure forced-convection closure (no dome, no boiling).

    Selects a property-ratio Nusselt correlation from the registry by regime,
    geometry, fluid and operating point; friction is single-phase Darcy at bulk;
    sound speed is the real-EOS value (Wood's equation never applies -- there is
    no two-phase mixture); CHF/quality are undefined and returned as None/NaN.
    HTD (heat-transfer-deterioration) risk is flagged via the McEligot-Jackson
    buoyancy AND flow-acceleration parameters (Urbano & Nasuti 2013) -- both are
    assumption-validity checks on the forced-convection property-ratio closures,
    not a deteriorated-HTC magnitude model.
    """
    if hydraulic_diameter_m <= 0.0 or mass_flux_kg_m2_s <= 0.0:
        raise ValueError("hydraulic diameter and mass flux must be positive")
    fluid_cp = coolprop_fluid_string(fluid, backend)

    Re_b = mass_flux_kg_m2_s * hydraulic_diameter_m / state.mu_Pa_s
    ctx = ClosureContext(
        fluid=fluid_cp,
        p_Pa=p_Pa,
        h_J_kg=h_J_kg,
        T_bulk_K=state.T_K,
        rho_b=state.rho_kg_m3,
        mu_b=state.mu_Pa_s,
        k_b=state.k_W_m_K,
        cp_b=state.cp_J_kg_K,
        Pr_b=state.Pr,
        mass_flux_kg_m2_s=mass_flux_kg_m2_s,
        diameter_m=hydraulic_diameter_m,
        heat_flux_W_m2=heat_flux_W_m2,
        wall_temp_K=wall_temp_K,
        x_over_D=x_over_D,
    )
    operating_point = {
        "p_reduced": state.p_reduced,
        "Re_b": Re_b,
        "Pr_b": state.Pr,
    }
    record, report = select_supercritical(
        regime=state.phase,
        geometry=geometry,
        orientation=orientation,
        fluid=fluid,
        operating_point=operating_point,
    )
    htc = record.callable(ctx)

    f = darcy_friction_smooth_pipe(Re_b)
    dpdz = f * mass_flux_kg_m2_s**2 / (2.0 * hydraulic_diameter_m * state.rho_kg_m3)
    c = get_cached_state(fluid_cp).flash_ph(p_Pa, h_J_kg).speed_sound()

    # HTD-precursor flags (buoyancy needs only bulk properties -- no wall state
    # required, unlike the old density-defect formula this replaces).
    beta_b = get_cached_state(fluid_cp).flash_ph(p_Pa, h_J_kg).isobaric_expansion_coefficient()
    buoyancy = _regime.buoyancy_parameter(
        beta_bulk_1_K=beta_b,
        heat_flux_W_m2=heat_flux_W_m2,
        k_bulk=state.k_W_m_K,
        rho_bulk=state.rho_kg_m3,
        mu_bulk=state.mu_Pa_s,
        diameter_m=hydraulic_diameter_m,
        Re_bulk=Re_b,
        Pr_bulk=state.Pr,
    )
    acceleration = _regime.acceleration_parameter(
        beta_bulk_1_K=beta_b,
        heat_flux_W_m2=heat_flux_W_m2,
        mass_flux_kg_m2_s=mass_flux_kg_m2_s,
        cp_bulk=state.cp_J_kg_K,
        Re_bulk=Re_b,
        Pr_bulk=state.Pr,
    )
    htd = _regime.htd_risk(buoyancy, acceleration)

    return CoolantClosureResult(
        state=state,
        htc_W_m2_K=float(htc),
        dpdz_friction_Pa_m=float(dpdz),
        chf_W_m2=None,
        chf_margin=None,
        sound_speed_m_s=float(c),
        regime=state.phase,
        closure_name=record.name,
        extrapolation_report=report,
        htd_risk=bool(htd),
        buoyancy_parameter=None if buoyancy is None else float(buoyancy),
    )

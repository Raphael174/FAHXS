"""Shell-side closure selection and the Sieder-Tate property correction.

Two fixes landed together (2026-08-20), both in `main_solve_shellntube.py`:

1. Bell-Delaware's `(mu_b/mu_w)^0.14` term was never being passed a viscosity
   ratio, so it sat at its neutral default of 1.0 — the correlation's own
   property-variation correction was switched off.

2. Which shell-side closure a node used was decided purely by the
   `coolantProp.coolant_model` STRING, not by the coolant's thermodynamic
   state. Any supercritical node in `equilibrium_liquid` mode went to a
   property-ratio closure regardless of how far it sat from the
   pseudo-critical region, while the identical state in
   `single_phase_coolprop` mode always got Bell-Delaware. Selection is now
   made on the state: the property-ratio closure is used only where the
   bulk-to-wall interval actually reaches the pseudo-critical band.

These tests pin the selection logic itself rather than solver output, so they
stay meaningful if the correlations behind them are retuned.
"""
from __future__ import annotations

import numpy as np
import pytest

from hps_combustor.input_data import (
    CorrelationCoefficients,
    coolantProp,
    hotgasProp,
    numericalProp,
    shellTubeProp,
    system_requirements,
)
from hps_combustor.main_solve_shellntube import shellntube_solver
from hps_combustor.physics.liquid_flow.coolprop_state_cache import coolprop_fluid_string
from hps_combustor.physics.liquid_flow.regime import pseudo_critical_temperature


def _solver(coolant="Nitrogen", model="equilibrium_liquid", T_in=100.0, p_in=88e5,
            mdot_c=17.5, N=40):
    return shellntube_solver(
        coolantProp=coolantProp(coolant=coolant, coolant_model=model,
                                mass_flow_c=mdot_c, T_in=T_in, p_in=p_in, T_out=T_in + 30.0),
        hotgasProp=hotgasProp(mass_flow_g=0.115),
        shellTubeProp=shellTubeProp(),
        numericalProp=numericalProp(chemistry_model="frozen"),
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
        N_axial=N, flow_config="co")


# ---------------------------------------------------------------- mu ratio
def test_mu_ratio_neutral_without_a_usable_wall_state():
    s = _solver()
    cool_cp = coolprop_fluid_string("Nitrogen", s._liquid_backend)
    for bad in (None, float("nan"), 0.0, -5.0):
        assert s._shell_mu_ratio("Nitrogen", cool_cp, 1e-4, bad, 88e5) == 1.0


def test_mu_ratio_is_bulk_over_wall_and_clamped():
    s = _solver()
    cool_cp = coolprop_fluid_string("Nitrogen", s._liquid_backend)
    import CoolProp.CoolProp as CP

    mu_b = CP.PropsSI("V", "T", 110.0, "P", 88e5, "Nitrogen")
    mu_w = CP.PropsSI("V", "T", 160.0, "P", 88e5, "Nitrogen")
    got = s._shell_mu_ratio("Nitrogen", cool_cp, mu_b, 160.0, 88e5)
    assert got == pytest.approx(mu_b / mu_w, rel=1e-9)

    lo, hi = s.MU_RATIO_LIMITS
    assert s._shell_mu_ratio("Nitrogen", cool_cp, 1e9, 160.0, 88e5) == hi
    assert s._shell_mu_ratio("Nitrogen", cool_cp, 1e-30, 160.0, 88e5) == lo


def test_mu_ratio_reaches_bell_delaware():
    """h_shell must actually respond to the correction, in the Sieder-Tate
    direction: a hotter wall (higher mu_w for a gas) lowers h."""
    s = _solver(coolant="Helium", model="single_phase_coolprop", T_in=303.15, mdot_c=0.075)
    h_neutral, _ = s._shell_h_at(400.0, wall_temp_K=400.0)          # mu_b == mu_w
    h_hot_wall, _ = s._shell_h_at(400.0, wall_temp_K=900.0)         # mu_w > mu_b
    assert h_hot_wall < h_neutral
    assert h_hot_wall / h_neutral == pytest.approx(1.0, abs=0.15)   # mild, ^0.14


# ------------------------------------------------------- regime dispatch
def test_helium_far_above_T_pc_does_not_need_property_ratio_closure():
    """Helium at 80 bar is supercritical by pressure but T_pc ~ 11 K against a
    300-1400 K march — no pseudo-critical anomaly is reachable."""
    s = _solver(coolant="Helium", T_in=303.15, p_in=80e5, mdot_c=0.075)
    cool_cp = coolprop_fluid_string("Helium", s._liquid_backend)
    T_pc = pseudo_critical_temperature(cool_cp, 80e5)
    assert T_pc < 50.0
    assert not s._needs_property_ratio_closure(cool_cp, 80e5, 303.15, 400.0, None)
    assert not s._needs_property_ratio_closure(cool_cp, 80e5, 1000.0, 1400.0, None)


def test_nitrogen_straddling_T_pc_needs_property_ratio_closure():
    """N2 at 88 bar has T_pc ~ 148 K with bulk ~100-124 K and wall ~164 K, so
    the pseudo-critical transition sits inside the thermal boundary layer."""
    s = _solver()
    cool_cp = coolprop_fluid_string("Nitrogen", s._liquid_backend)
    T_pc = pseudo_critical_temperature(cool_cp, 88e5)
    assert 130.0 < T_pc < 165.0
    assert s._needs_property_ratio_closure(cool_cp, 88e5, 110.0, 170.0, None)
    # bulk and wall both far below T_pc -> ordinary dense single-phase fluid
    assert not s._needs_property_ratio_closure(cool_cp, 88e5, 90.0, 100.0, None)
    # both far above -> ordinary gas-like fluid
    assert not s._needs_property_ratio_closure(cool_cp, 88e5, 400.0, 450.0, None)


def test_selection_latches_one_way_per_node():
    """The choice must not chatter while the lagged wall temperature is still
    climbing off its cold seed."""
    s = _solver()
    cool_cp = coolprop_fluid_string("Nitrogen", s._liquid_backend)
    assert not s._sc_latch[3]
    assert not s._needs_property_ratio_closure(cool_cp, 88e5, 90.0, 100.0, 3)
    assert not s._sc_latch[3]
    assert s._needs_property_ratio_closure(cool_cp, 88e5, 110.0, 170.0, 3)
    assert s._sc_latch[3]
    # now latched: a state that would otherwise fall back keeps the closure
    assert s._needs_property_ratio_closure(cool_cp, 88e5, 90.0, 100.0, 3)
    # and the latch is per node, not global
    assert not s._sc_latch[4]


def test_latch_resets_between_solves():
    s = _solver(N=20)
    s._sc_latch[:] = True
    s.solve(verbose=False, max_sweeps=1)
    assert s._sc_latch.sum() < s.N          # rebuilt from the state, not inherited


# ------------------------------------------------- pressure-march consistency
def test_pressure_march_always_uses_bell_delaware():
    """The pressure march uses Bell-Delaware at EVERY liquid-mode node, never the
    closure's own friction gradient — regardless of which closure supplies h.

    Gungor-Winterton/MSH and the supercritical registry are straight-TUBE
    correlations: they model axial flow along one L_tube-long channel with wall
    skin friction. The real shell-side path crosses the bundle N_baffles+1 times
    through N_tcc rows each — roughly 7.5x the path length on this geometry,
    with form drag rather than skin friction — so they under-predict by about
    25x, far too much to accept as an extrapolation."""
    import CoolProp.CoolProp as CP
    s = _solver()
    h_bulk = CP.PropsSI("H", "T", 105.0, "P", 88e5, "Nitrogen")

    # (a) supercritical node far from T_pc -> Bell-Delaware supplies h as well
    _, dp_fb = s._shell_h_at(105.0, h_c_enthalpy=h_bulk, quality_local=float("nan"),
                             p_local=88e5, wall_temp_K=108.0, node_index=0)
    assert s._sc_bell_fallback_nodes == 1
    assert dp_fb == pytest.approx(s._last_bell["dp_shell"] / s.N, rel=1e-12)

    # (b) supercritical node straddling T_pc -> h comes from the property-ratio
    #     closure, but dp must STILL be Bell-Delaware's
    _, dp_pr = s._shell_h_at(140.0, h_c_enthalpy=CP.PropsSI("H", "T", 140.0, "P", 88e5, "Nitrogen"),
                             quality_local=float("nan"), p_local=88e5,
                             wall_temp_K=170.0, node_index=1)
    assert s._sc_bell_fallback_nodes == 1          # h did NOT come from Bell-Delaware
    assert dp_pr == pytest.approx(s._last_bell["dp_shell"] / s.N, rel=1e-12)
    assert dp_pr > 0.0


def test_bell_delaware_survives_the_two_phase_dome():
    """Inside the dome CoolProp's c_p and transport properties are undefined or
    pathological — a negative c_p makes Pr negative and Pr^(-2/3) COMPLEX, which
    used to crash the Bell-Delaware call outright. Saturated-liquid transport
    properties are substituted while the homogeneous two-phase density (which is
    what actually carries the pressure drop through vaporization) is kept."""
    import CoolProp.CoolProp as CP
    s = _solver(coolant="Water", T_in=300.0, p_in=81e5, mdot_c=0.86)
    p = 78e5
    for x in (0.05, 0.5, 0.95):
        h = CP.PropsSI("H", "P", p, "Q", x, "Water")
        h_c, dp = s._shell_h_at(CP.PropsSI("T", "P", p, "Q", x, "Water"),
                                h_c_enthalpy=h, quality_local=x, p_local=p,
                                wall_temp_K=600.0, node_index=0)
        assert np.isfinite(h_c) and h_c > 0.0
        assert np.isfinite(dp) and dp > 0.0
    # density falls through vaporization, so dp must rise with quality
    def _dp(x):
        h = CP.PropsSI("H", "P", p, "Q", x, "Water")
        return s._shell_h_at(CP.PropsSI("T", "P", p, "Q", x, "Water"),
                             h_c_enthalpy=h, quality_local=x, p_local=p,
                             wall_temp_K=600.0, node_index=0)[1]
    assert _dp(0.9) > _dp(0.1)


def test_bell_delaware_leakage_ratio_validity_is_reported():
    """r_lm = (S_sb+S_tb)/S_m is ~6.5 on this geometry against a fitted range of
    roughly r_lm <= 1, so the leakage corrections are extrapolated. Pinned so
    the diagnostic cannot be dropped silently."""
    s = _solver()
    assert s._bell_r_lm > 1.0
    expected = (s.geom["S_sb"] + s.geom["S_tb"]) / s.geom["S_m"]
    assert s._bell_r_lm == pytest.approx(expected, rel=1e-12)

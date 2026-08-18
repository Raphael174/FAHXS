"""Stage A of docs/solver_design/FV_CORE_REWORK_PLAN.md.

Proves the new ``hps_combustor.core.thermo`` backends reproduce the exact
CoolProp call patterns the legacy code paths use, and that
``physics/liquid_flow/dispatch.py``'s ``CoolantState``/state-construction
re-exports are the SAME objects (not independent reimplementations) after the
Stage A relocation -- the mechanism behind the bit-identical acceptance gate.
"""
from __future__ import annotations

import CoolProp.CoolProp as CP
import pytest

from hps_combustor.core.thermo import (
    IdealGasBackend,
    RealFluidBackend,
    ThermoState,
    coolant_state_from_Tp,
    coolant_state_from_ph,
)
from hps_combustor.physics.liquid_flow import dispatch


def test_dispatch_reexports_are_identical_objects():
    """CoolantState/state-construction functions must be the SAME objects
    core.thermo defines, not parallel copies -- this is what makes the
    Stage A move behaviorally inert for every existing importer."""
    assert dispatch.CoolantState is ThermoState
    assert dispatch.coolant_state_from_Tp is coolant_state_from_Tp
    assert dispatch.coolant_state_from_ph is coolant_state_from_ph


def test_real_fluid_backend_state_pT_matches_raw_propsi():
    backend = RealFluidBackend()
    state = backend.state_pT("Helium", 300.0, 80e5)
    assert state.rho_kg_m3 == pytest.approx(
        CP.PropsSI("D", "T", 300.0, "P", 80e5, "Helium")
    )
    assert state.h_J_kg == pytest.approx(CP.PropsSI("H", "T", 300.0, "P", 80e5, "Helium"))
    assert state.phase == "single_phase"


def test_real_fluid_backend_state_ph_water_two_phase():
    backend = RealFluidBackend()
    p = 5e5
    h_sat_l = CP.PropsSI("H", "P", p, "Q", 0.0, "Water")
    h_sat_v = CP.PropsSI("H", "P", p, "Q", 1.0, "Water")
    h_mid = 0.5 * (h_sat_l + h_sat_v)
    state = backend.state_ph("Water", p, h_mid)
    assert state.phase == "two_phase"
    assert 0.0 < state.quality < 1.0


def test_ideal_gas_backend_matches_exact_propsi_call_granularity():
    """main_solve.py's inline gas-path calls one property at a time (e.g.
    only density at one march point, only enthalpy at another) -- confirm the
    backend's per-property getters return exactly what the raw PropsSI call
    at that call site returns today, so a future call-site swap changes
    nothing numerically."""
    backend = IdealGasBackend()
    fluid, T, p = "Helium", 450.0, 78e5
    assert backend.density(fluid, T, p) == CP.PropsSI("D", "T", T, "P", p, fluid)
    assert backend.enthalpy(fluid, T, p) == CP.PropsSI("H", "T", T, "P", p, fluid)
    assert backend.viscosity(fluid, T, p) == CP.PropsSI("V", "T", T, "P", p, fluid)
    assert backend.conductivity(fluid, T, p) == CP.PropsSI("L", "T", T, "P", p, fluid)
    assert backend.cp(fluid, T, p) == CP.PropsSI("C", "T", T, "P", p, fluid)
    assert backend.cv(fluid, T, p) == CP.PropsSI("CVMASS", "T", T, "P", p, fluid)
    assert backend.compressibility(fluid, T, p) == CP.PropsSI("Z", "T", T, "P", p, fluid)
    assert backend.molar_mass(fluid) == CP.PropsSI("MOLAR_MASS", fluid)


def test_ideal_gas_backend_state_pT_matches_coolant_state_from_Tp():
    import math

    backend = IdealGasBackend()
    a = backend.state_pT("Helium", 350.0, 80e5)
    b = coolant_state_from_Tp("Helium", 350.0, 80e5)
    assert math.isnan(a.quality) and math.isnan(b.quality)
    a = a.__class__(**{**a.__dict__, "quality": 0.0})
    b = b.__class__(**{**b.__dict__, "quality": 0.0})
    assert a == b


def test_p_crit_consistent_between_backends():
    real = RealFluidBackend().p_crit("Nitrogen")
    ideal = IdealGasBackend().p_crit("Nitrogen")
    assert real == pytest.approx(ideal, rel=1e-9)


def test_rpe_dittus_boelter_cross_check_registered_and_reported():
    """The Eq. 8-24 cross-check (added 2026-07-31, user-requested sanity
    check) must be reported alongside the active subcooled-liquid HTC, never
    used in its place."""
    from hps_combustor.input_data import coolantProp

    coolant_prop = coolantProp(
        coolant="Water", coolant_model="equilibrium_liquid", T_in=300.0, p_in=20e5
    )
    h_in = CP.PropsSI("H", "T", 300.0, "P", 20e5, "Water")
    result = dispatch.evaluate_coolant_closure(
        coolant_prop=coolant_prop,
        p_Pa=20e5,
        h_J_kg=h_in,
        mass_flux_kg_m2_s=2000.0,
        hydraulic_diameter_m=6e-3,
        heat_flux_W_m2=5e5,
    )
    assert result.cross_check_closure_name == "rpe_dittus_boelter_8_24"
    assert result.cross_check_htc_W_m2_K is not None
    assert result.cross_check_htc_W_m2_K > 0.0
    # Cross-check must never overwrite the active (Gnielinski) HTC.
    from hps_combustor.physics.liquid_flow.correlations import (
        dittus_boelter_colburn_liquid_nusselt,
        liquid_single_phase_nusselt,
    )

    Re = 2000.0 * 6e-3 / result.state.mu_Pa_s
    active_expected = (
        liquid_single_phase_nusselt(Re, result.state.Pr) * result.state.k_W_m_K / 6e-3
    )
    assert result.htc_W_m2_K == pytest.approx(active_expected, rel=1e-9)
    assert result.htc_W_m2_K != pytest.approx(result.cross_check_htc_W_m2_K, rel=1e-2)


def test_rpe_dittus_boelter_closure_matches_hand_formula():
    from hps_combustor.physics.liquid_flow.correlations import (
        dittus_boelter_colburn_liquid_nusselt,
    )

    Re, Pr = 1.0e5, 3.0
    Nu = dittus_boelter_colburn_liquid_nusselt(Re, Pr)
    assert Nu == pytest.approx(0.023 * Re**0.8 * Pr ** (1.0 / 3.0), rel=1e-12)

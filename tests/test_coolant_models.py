import CoolProp.CoolProp as CP
import numpy as np
import pytest

from hps_combustor.input_data import coolantProp
from hps_combustor.physics.liquid_flow.dispatch import (
    coolant_inlet_state,
    coolant_state_from_ph,
    evaluate_coolant_closure,
)
from hps_combustor.physics.liquid_flow.governing_equations import (
    HeatedChannelCase,
    HeatedChannelProfileCase,
    heated_channel_cell_fields,
    solve_steady_heated_channel,
    solve_steady_heated_channel_on_hx_grid,
    solve_steady_heated_channel_profile,
    summarize_heated_channel_result,
)
from hps_combustor.physics.liquid_flow.correlations import saturation_state


def test_default_coolant_prop_keeps_legacy_single_phase_helium_model():
    props = coolantProp()

    assert props.coolant == "Helium"
    assert props.coolant_model == "single_phase_coolprop"

    state = coolant_inlet_state(props)

    assert state.fluid == "Helium"
    assert state.model == "single_phase_coolprop"
    assert state.phase == "single_phase"
    assert state.T_K == pytest.approx(props.T_in)
    assert state.p_Pa == pytest.approx(props.p_in)
    assert state.rho_kg_m3 > 0.0
    assert state.cp_J_kg_K > 0.0
    assert state.Pr > 0.0


def test_equilibrium_liquid_state_from_ph_tracks_saturated_quality():
    p = 5.0e6
    sat = saturation_state("Water", p)
    h = sat.h_l_J_kg + 0.2 * sat.h_fg_J_kg

    state = coolant_state_from_ph("Water", p, h, "equilibrium_liquid")

    assert state.fluid == "Water"
    assert state.model == "equilibrium_liquid"
    assert state.phase == "two_phase"
    assert state.T_K == pytest.approx(sat.T_sat_K)
    assert state.quality == pytest.approx(0.2)
    assert state.void_fraction > 0.0
    assert state.rho_kg_m3 > 0.0


def test_equilibrium_liquid_closure_returns_valid_htc_dp_and_chf_margin():
    p = 5.0e6
    sat = saturation_state("Water", p)
    h = sat.h_l_J_kg + 0.05 * sat.h_fg_J_kg
    props = coolantProp(coolant="Water", coolant_model="equilibrium_liquid")

    closure = evaluate_coolant_closure(
        coolant_prop=props,
        p_Pa=p,
        h_J_kg=h,
        mass_flux_kg_m2_s=636.6197723675814,
        hydraulic_diameter_m=0.02,
        heat_flux_W_m2=1.0e5,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert closure.state.quality == pytest.approx(0.05)
    assert closure.htc_W_m2_K == pytest.approx(21235.036, rel=2.0e-4)
    assert closure.dpdz_friction_Pa_m > 0.0
    assert closure.chf_W_m2 is not None
    assert closure.chf_W_m2 > 1.0e5
    assert closure.chf_margin == pytest.approx(closure.chf_W_m2 / 1.0e5)


def test_equilibrium_liquid_closure_allows_fully_vapor_expansion_state():
    p = 2.0e5
    sat = saturation_state("Water", p)
    props = coolantProp(coolant="Water", coolant_model="equilibrium_liquid")

    closure = evaluate_coolant_closure(
        coolant_prop=props,
        p_Pa=p,
        h_J_kg=sat.h_v_J_kg + 2.0e4,
        mass_flux_kg_m2_s=100.0,
        hydraulic_diameter_m=0.003,
        heat_flux_W_m2=5.0e4,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert closure.state.phase == "vapor"
    assert closure.state.quality > 1.0
    assert closure.htc_W_m2_K > 0.0
    assert closure.dpdz_friction_Pa_m > 0.0


def test_single_phase_coolprop_closure_works_from_pressure_enthalpy():
    # Pinned explicitly: this test's intent is the single_phase_coolprop
    # closure branch, not whatever coolantProp()'s current default happens to
    # be (that default is shared, mutable project state — see
    # test_steady_baseline_regression.py's module docstring for why relying
    # on it silently broke a different test).
    props = coolantProp(coolant="Helium", coolant_model="single_phase_coolprop")
    h = CP.PropsSI("H", "T", props.T_in, "P", props.p_in, props.coolant)

    closure = evaluate_coolant_closure(
        coolant_prop=props,
        p_Pa=props.p_in,
        h_J_kg=h,
        mass_flux_kg_m2_s=200.0,
        hydraulic_diameter_m=0.01,
        heat_flux_W_m2=0.0,
    )

    assert closure.state.phase == "single_phase"
    assert closure.htc_W_m2_K > 0.0
    assert closure.dpdz_friction_Pa_m > 0.0
    assert closure.chf_W_m2 is None
    assert closure.chf_margin is None


def test_heated_channel_solver_matches_saturated_water_regression():
    p = 5.0e6
    sat = saturation_state("Water", p)
    case = HeatedChannelCase(
        fluid="Water",
        length_m=1.0,
        hydraulic_diameter_m=0.02,
        mass_flow_kg_s=0.20,
        p_in_Pa=p,
        h_in_J_kg=sat.h_l_J_kg + 0.05 * sat.h_fg_J_kg,
        heat_flux_W_m2=1.0e5,
        n_cells=80,
        coolant_model="equilibrium_liquid",
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    result = solve_steady_heated_channel(case)

    assert abs(result.energy_residual_J_kg) < 1.0e-6
    assert result.outlet_quality == pytest.approx(0.0691987, rel=2.0e-5)
    assert result.pressure_drop_Pa == pytest.approx(1064.67, rel=2.0e-4)
    assert result.min_chf_margin > 1.0
    assert result.htc_W_m2_K[-1] > 0.0


def test_heated_channel_profile_solver_matches_uniform_channel_api():
    p = 5.0e6
    sat = saturation_state("Water", p)
    uniform = HeatedChannelCase(
        fluid="Water",
        length_m=1.0,
        hydraulic_diameter_m=0.02,
        mass_flow_kg_s=0.20,
        p_in_Pa=p,
        h_in_J_kg=sat.h_l_J_kg + 0.05 * sat.h_fg_J_kg,
        heat_flux_W_m2=1.0e5,
        n_cells=40,
        coolant_model="equilibrium_liquid",
        lut_path="docs/reference/external/2006LUTdata.txt",
    )
    profile = HeatedChannelProfileCase(
        coolant_prop=coolantProp(coolant="Water", coolant_model="equilibrium_liquid"),
        z_edges_m=np.linspace(0.0, uniform.length_m, uniform.n_cells + 1),
        hydraulic_diameter_m=uniform.hydraulic_diameter_m,
        flow_area_m2=np.pi * uniform.hydraulic_diameter_m**2 / 4.0,
        heated_perimeter_m=np.pi * uniform.hydraulic_diameter_m,
        mass_flow_kg_s=uniform.mass_flow_kg_s,
        p_in_Pa=uniform.p_in_Pa,
        h_in_J_kg=uniform.h_in_J_kg,
        heat_flux_W_m2=uniform.heat_flux_W_m2,
        lut_path=uniform.lut_path,
    )

    result_uniform = solve_steady_heated_channel(uniform)
    result_profile = solve_steady_heated_channel_profile(profile)

    np.testing.assert_allclose(result_profile.p_Pa, result_uniform.p_Pa)
    np.testing.assert_allclose(result_profile.h_J_kg, result_uniform.h_J_kg)
    np.testing.assert_allclose(result_profile.quality, result_uniform.quality)
    np.testing.assert_allclose(result_profile.htc_W_m2_K, result_uniform.htc_W_m2_K)


def test_heated_channel_profile_accepts_wall_heat_per_length_input():
    p = 2.0e5
    sat = saturation_state("Water", p)
    z_edges = np.linspace(0.0, 0.5, 21)
    diameter = np.full(20, 0.003)
    perimeter = np.pi * diameter
    heat_flux = 1.0e5 + 5.0e4 * np.linspace(0.0, 1.0, 20)
    heat_per_length = heat_flux * perimeter
    case = HeatedChannelProfileCase(
        coolant_prop=coolantProp(coolant="Water", coolant_model="equilibrium_liquid"),
        z_edges_m=z_edges,
        hydraulic_diameter_m=diameter,
        flow_area_m2=np.pi * diameter**2 / 4.0,
        heated_perimeter_m=perimeter,
        mass_flow_kg_s=8.0e-4,
        p_in_Pa=p,
        h_in_J_kg=sat.h_l_J_kg - 1.0e4,
        heat_per_length_W_m=heat_per_length,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    result = solve_steady_heated_channel_profile(case)

    expected_heat = float(np.sum(heat_per_length * np.diff(z_edges)))
    assert result.heat_rate_W == pytest.approx(expected_heat)
    assert abs(result.energy_residual_J_kg) < 1.0e-6
    assert result.quality[0] < 0.0
    assert result.quality[-1] > result.quality[0]
    assert result.pressure_drop_Pa > 0.0
    np.testing.assert_allclose(result.heat_flux_W_m2[:-1], heat_flux)


def test_heated_channel_profile_accepts_wall_heat_per_segment_input():
    p = 2.0e5
    sat = saturation_state("Water", p)
    z_edges = np.array([0.0, 0.04, 0.10, 0.19, 0.31, 0.50])
    dz = np.diff(z_edges)
    n_cells = dz.size
    diameter = 0.003 + 2.0e-4 * np.linspace(0.0, 1.0, n_cells)
    perimeter = np.pi * diameter
    heat_flux = 8.0e4 + 4.0e4 * np.linspace(0.0, 1.0, n_cells)
    heat_per_segment = heat_flux * perimeter * dz
    common = dict(
        coolant_prop=coolantProp(coolant="Water", coolant_model="equilibrium_liquid"),
        z_edges_m=z_edges,
        hydraulic_diameter_m=diameter,
        flow_area_m2=np.pi * diameter**2 / 4.0,
        heated_perimeter_m=perimeter,
        mass_flow_kg_s=8.0e-4,
        p_in_Pa=p,
        h_in_J_kg=sat.h_l_J_kg - 1.0e4,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    by_segment = solve_steady_heated_channel_profile(
        HeatedChannelProfileCase(heat_per_segment_W=heat_per_segment, **common)
    )
    by_flux = solve_steady_heated_channel_profile(
        HeatedChannelProfileCase(heat_flux_W_m2=heat_flux, **common)
    )

    assert by_segment.heat_rate_W == pytest.approx(float(np.sum(heat_per_segment)))
    np.testing.assert_allclose(by_segment.heat_flux_W_m2, by_flux.heat_flux_W_m2)
    np.testing.assert_allclose(by_segment.h_J_kg, by_flux.h_J_kg)
    np.testing.assert_allclose(by_segment.p_Pa, by_flux.p_Pa)


def test_hx_grid_liquid_adapter_coflow_matches_profile_solver():
    p = 2.0e5
    sat = saturation_state("Water", p)
    z_edges = np.array([0.0, 0.04, 0.10, 0.19, 0.31, 0.50])
    dz = np.diff(z_edges)
    n_cells = dz.size
    diameter = 0.003 + 2.0e-4 * np.linspace(0.0, 1.0, n_cells)
    heat_flux = 8.0e4 + 4.0e4 * np.linspace(0.0, 1.0, n_cells)
    common = dict(
        coolant_prop=coolantProp(coolant="Water", coolant_model="equilibrium_liquid"),
        z_edges_m=z_edges,
        hydraulic_diameter_m=diameter,
        flow_area_m2=np.pi * diameter**2 / 4.0,
        heated_perimeter_m=np.pi * diameter,
        mass_flow_kg_s=8.0e-4,
        p_in_Pa=p,
        h_in_J_kg=sat.h_l_J_kg - 1.0e4,
        heat_per_segment_W=heat_flux * np.pi * diameter * dz,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    profile = solve_steady_heated_channel_profile(HeatedChannelProfileCase(**common))
    hx = solve_steady_heated_channel_on_hx_grid(coolant_enters_at="z_min", **common)

    np.testing.assert_allclose(hx.flow_result.p_Pa, profile.p_Pa)
    np.testing.assert_allclose(hx.node_fields_hx_order["p_Pa"], profile.p_Pa)
    np.testing.assert_allclose(hx.cell_fields_hx_order["heat_flux_W_m2"], heat_flux)
    assert hx.diagnostics.pressure_drop_Pa == pytest.approx(profile.pressure_drop_Pa)


def test_hx_grid_liquid_adapter_counterflow_maps_fields_to_hx_order():
    p = 2.0e5
    sat = saturation_state("Water", p)
    z_edges = np.array([0.0, 0.04, 0.10, 0.19, 0.31, 0.50])
    dz = np.diff(z_edges)
    n_cells = dz.size
    diameter = 0.003 + 2.0e-4 * np.linspace(0.0, 1.0, n_cells)
    heat_flux = 8.0e4 + 4.0e4 * np.linspace(0.0, 1.0, n_cells)
    heat_per_segment = heat_flux * np.pi * diameter * dz

    hx = solve_steady_heated_channel_on_hx_grid(
        coolant_prop=coolantProp(coolant="Water", coolant_model="equilibrium_liquid"),
        z_edges_m=z_edges,
        hydraulic_diameter_m=diameter,
        flow_area_m2=np.pi * diameter**2 / 4.0,
        heated_perimeter_m=np.pi * diameter,
        mass_flow_kg_s=8.0e-4,
        p_in_Pa=p,
        h_in_J_kg=sat.h_l_J_kg - 1.0e4,
        heat_per_segment_W=heat_per_segment,
        coolant_enters_at="z_max",
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert hx.coolant_enters_at == "z_max"
    assert hx.node_fields_hx_order["z_m"][0] == pytest.approx(z_edges[0])
    assert hx.node_fields_hx_order["z_m"][-1] == pytest.approx(z_edges[-1])
    assert hx.node_fields_hx_order["p_Pa"][-1] == pytest.approx(p)
    assert hx.node_fields_hx_order["h_J_kg"][-1] == pytest.approx(sat.h_l_J_kg - 1.0e4)
    assert hx.node_fields_hx_order["h_J_kg"][0] > hx.node_fields_hx_order["h_J_kg"][-1]
    np.testing.assert_allclose(hx.cell_fields_hx_order["heat_flux_W_m2"], heat_flux)
    assert hx.diagnostics.heat_rate_W == pytest.approx(float(np.sum(heat_per_segment)))
    assert hx.diagnostics.energy_residual_ok is True


def test_heated_channel_cell_fields_and_diagnostics_for_boiling_case():
    p = 5.0e6
    sat = saturation_state("Water", p)
    case = HeatedChannelCase(
        fluid="Water",
        length_m=1.0,
        hydraulic_diameter_m=0.02,
        mass_flow_kg_s=0.20,
        p_in_Pa=p,
        h_in_J_kg=sat.h_l_J_kg + 0.05 * sat.h_fg_J_kg,
        heat_flux_W_m2=1.0e5,
        n_cells=16,
        coolant_model="equilibrium_liquid",
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    result = solve_steady_heated_channel(case)
    fields = heated_channel_cell_fields(result)
    diag = summarize_heated_channel_result(result, min_pressure_Pa=case.min_pressure_Pa)

    expected_keys = {
        "z_m",
        "p_Pa",
        "h_J_kg",
        "T_K",
        "quality",
        "void_fraction",
        "rho_kg_m3",
        "htc_W_m2_K",
        "dpdz_friction_Pa_m",
        "dpdz_acceleration_Pa_m",
        "chf_W_m2",
        "chf_margin",
        "heat_flux_W_m2",
    }
    assert set(fields) == expected_keys
    assert all(value.shape == (case.n_cells,) for value in fields.values())
    assert diag.boiling_reached is True
    assert diag.dryout_or_vapor_reached is False
    assert diag.chf_margin_below_limit is False
    assert diag.energy_residual_ok is True
    assert diag.pressure_drop_Pa == pytest.approx(result.pressure_drop_Pa)
    assert diag.min_chf_margin == pytest.approx(result.min_chf_margin)


def test_heated_channel_diagnostics_flag_vapor_and_absent_chf_margin():
    p = 2.0e5
    sat = saturation_state("Water", p)
    case = HeatedChannelCase(
        fluid="Water",
        length_m=0.2,
        hydraulic_diameter_m=0.003,
        mass_flow_kg_s=8.0e-4,
        p_in_Pa=p,
        h_in_J_kg=sat.h_v_J_kg + 1.0e4,
        heat_flux_W_m2=2.0e5,
        n_cells=8,
        coolant_model="equilibrium_liquid",
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    result = solve_steady_heated_channel(case)
    diag = summarize_heated_channel_result(
        result,
        min_pressure_Pa=case.min_pressure_Pa,
        chf_margin_limit=1.0e9,
    )

    assert diag.boiling_reached is False
    assert diag.dryout_or_vapor_reached is True
    assert diag.chf_margin_below_limit is False
    assert np.isnan(diag.min_chf_margin)
    assert diag.max_quality > 1.0


def test_heated_channel_profile_rejects_ambiguous_heat_inputs():
    props = coolantProp(coolant="Water", coolant_model="equilibrium_liquid")
    base = dict(
        coolant_prop=props,
        z_edges_m=np.array([0.0, 1.0]),
        hydraulic_diameter_m=0.01,
        flow_area_m2=np.pi * 0.01**2 / 4.0,
        heated_perimeter_m=np.pi * 0.01,
        mass_flow_kg_s=0.01,
        p_in_Pa=1.0e5,
        T_in_K=370.0,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    with pytest.raises(ValueError, match="exactly one"):
        solve_steady_heated_channel_profile(HeatedChannelProfileCase(**base))
    with pytest.raises(ValueError, match="exactly one"):
        solve_steady_heated_channel_profile(
            HeatedChannelProfileCase(
                heat_flux_W_m2=1.0e5,
                heat_per_segment_W=10.0,
                **base,
            )
        )


def test_heated_channel_solver_accepts_profile_heat_flux_and_subcooled_inlet():
    p = 2.0e5
    sat = saturation_state("Water", p)
    q_profile = 1.0e5 + 5.0e4 * np.linspace(0.0, 1.0, 20)
    case = HeatedChannelCase(
        fluid="Water",
        length_m=0.5,
        hydraulic_diameter_m=0.003,
        mass_flow_kg_s=8.0e-4,
        p_in_Pa=p,
        h_in_J_kg=sat.h_l_J_kg - 1.0e4,
        heat_flux_W_m2=q_profile,
        n_cells=20,
        coolant_model="equilibrium_liquid",
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    result = solve_steady_heated_channel(case)

    expected_heat = float(q_profile.mean() * 3.141592653589793 * case.hydraulic_diameter_m * case.length_m)
    assert result.heat_rate_W == pytest.approx(expected_heat)
    assert abs(result.energy_residual_J_kg) < 1.0e-6
    assert result.quality[0] < 0.0
    assert result.quality[-1] > result.quality[0]
    assert result.pressure_drop_Pa > 0.0

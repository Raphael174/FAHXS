import math

import numpy as np
import pytest

from hps_combustor.physics.liquid_flow.chf import (
    interpolate_chf_table,
    local_chf_diameter_correction,
    load_groeneveld_2006_lut,
)
from hps_combustor.physics.liquid_flow.correlations import (
    chisholm_two_phase_multiplier,
    darcy_friction_smooth_pipe,
    equilibrium_state_ph,
    gungor_winterton_boiling_htc,
    homogeneous_void_fraction,
    martinelli_parameter_laminar_liquid_turbulent_vapor,
    muller_steinhagen_heck_friction_gradient,
    saturation_state,
    yu2002_modified_anl_boiling_htc,
    yu2002_small_channel_pressure_multiplier,
)
from hps_combustor.validation.liquid_boiling_straight_pipe import (
    generate_yu2002_validation_report,
    saturated_water_reference_case,
    solve_steady_straight_pipe,
)
from hps_combustor.validation.liquid_hx_imposed_duty import (
    generate_imposed_duty_report,
    imposed_duty_reference_case,
)
from hps_combustor.validation.liquid_validation_matrix import run_validation_matrix


def test_equilibrium_state_quality_and_void_fraction():
    sat = saturation_state("Water", 1.0e5)
    h_mid = sat.h_l_J_kg + 0.5 * sat.h_fg_J_kg

    state = equilibrium_state_ph(1.0e5, h_mid, "Water")
    expected_alpha = homogeneous_void_fraction(0.5, sat.rho_l_kg_m3, sat.rho_v_kg_m3)

    assert state.phase == "two_phase"
    assert state.quality == pytest.approx(0.5)
    assert state.void_fraction == pytest.approx(expected_alpha)
    assert state.rho_kg_m3 == pytest.approx(
        1.0 / (0.5 / sat.rho_v_kg_m3 + 0.5 / sat.rho_l_kg_m3)
    )
    assert state.void_fraction > 0.99


def test_gungor_winterton_regression_for_water():
    htc = gungor_winterton_boiling_htc(
        p_Pa=5.0e6,
        mass_flux_kg_m2_s=636.6197723675814,
        diameter_m=0.02,
        quality=0.05,
        heat_flux_W_m2=1.0e5,
        fluid="Water",
    )

    assert htc == pytest.approx(21235.036, rel=2.0e-4)


def test_muller_steinhagen_heck_limits_are_single_phase_endpoints():
    p = 5.0e6
    G = 500.0
    D = 0.02
    sat = saturation_state("Water", p)

    liquid_limit = muller_steinhagen_heck_friction_gradient(
        p_Pa=p, mass_flux_kg_m2_s=G, diameter_m=D, quality=0.0, fluid="Water"
    )
    vapor_limit = muller_steinhagen_heck_friction_gradient(
        p_Pa=p, mass_flux_kg_m2_s=G, diameter_m=D, quality=1.0, fluid="Water"
    )

    assert liquid_limit > 0.0
    assert vapor_limit > liquid_limit
    f_l0 = darcy_friction_smooth_pipe(G * D / sat.mu_l_Pa_s)
    assert liquid_limit == pytest.approx(f_l0 * G**2 / (2.0 * D * sat.rho_l_kg_m3))


def test_yu2002_small_channel_correlations_regressions():
    sat = saturation_state("Water", 2.0e5)
    x = 0.6
    Re_l = 900.0
    Re_v = 12000.0
    X = martinelli_parameter_laminar_liquid_turbulent_vapor(
        quality=x,
        rho_l=sat.rho_l_kg_m3,
        rho_v=sat.rho_v_kg_m3,
        Re_l=Re_l,
        Re_v=Re_v,
    )

    assert X == pytest.approx(0.0366863, rel=2.0e-4)
    assert chisholm_two_phase_multiplier(0.1, C=12.0) == pytest.approx(221.0)
    assert yu2002_small_channel_pressure_multiplier(0.1) == pytest.approx(79.4328, rel=2.0e-5)

    htc = yu2002_modified_anl_boiling_htc(
        p_Pa=2.0e5,
        mass_flux_kg_m2_s=103.0,
        diameter_m=0.00298,
        heat_flux_W_m2=1.0e5,
        fluid="Water",
    )
    assert htc == pytest.approx(22491.3, rel=2.0e-4)


def test_groeneveld_chf_diameter_correction_and_table_interpolation():
    assert local_chf_diameter_correction(1.0e6, 0.008) == pytest.approx(1.0e6)
    assert local_chf_diameter_correction(1.0e6, 0.032) == pytest.approx(0.5e6)

    p_axis = np.array([1.0, 2.0])
    g_axis = np.array([100.0, 200.0])
    x_axis = np.array([0.0, 1.0])
    table = np.empty((2, 2, 2))
    for i, p in enumerate(p_axis):
        for j, g in enumerate(g_axis):
            for k, x in enumerate(x_axis):
                table[i, j, k] = p + 0.01 * g + 10.0 * x

    assert interpolate_chf_table(
        p_MPa=1.5,
        mass_flux_kg_m2_s=150.0,
        quality=0.5,
        pressures_MPa=p_axis,
        mass_fluxes_kg_m2_s=g_axis,
        qualities=x_axis,
        chf_kW_m2=table,
    ) == pytest.approx((1.5 + 1.5 + 5.0) * 1000.0)


def test_groeneveld_2006_lut_file_matches_paper_page9_values():
    p_axis, g_axis, x_axis, table = load_groeneveld_2006_lut(
        "docs/reference/external/2006LUTdata.txt"
    )

    i_p = int(np.where(np.isclose(p_axis, 0.10))[0][0])
    i_g0 = int(np.where(np.isclose(g_axis, 0.0))[0][0])
    i_g50 = int(np.where(np.isclose(g_axis, 50.0))[0][0])
    i_x0 = int(np.where(np.isclose(x_axis, 0.00))[0][0])
    i_x50 = int(np.where(np.isclose(x_axis, 0.50))[0][0])

    assert table.shape == (15, 21, 23)
    assert table[i_p, i_g0, i_x0] == pytest.approx(1142.0)
    assert table[i_p, i_g0, i_x50] == pytest.approx(123.0)
    assert table[i_p, i_g50, i_x0] == pytest.approx(1570.0)


def test_straight_pipe_reference_case_energy_and_ranges():
    case = saturated_water_reference_case()
    result = solve_steady_straight_pipe(case)

    expected_heat = math.pi * case.diameter_m * case.length_m * case.heat_flux_W_m2
    expected_dh = expected_heat / case.mass_flow_kg_s

    assert result.heat_rate_W == pytest.approx(expected_heat)
    assert result.h_J_kg[-1] - result.h_J_kg[0] == pytest.approx(expected_dh)
    assert abs(result.energy_residual_J_kg) < 1.0e-6
    assert result.outlet_quality == pytest.approx(0.0691987, rel=2.0e-5)
    assert result.pressure_drop_Pa == pytest.approx(1064.67, rel=2.0e-4)
    assert np.all(result.htc_W_m2_K > 0.0)
    assert np.all(np.diff(result.h_J_kg) > 0.0)


def test_yu2002_validation_report_generates_metrics(tmp_path):
    summary = generate_yu2002_validation_report(tmp_path)

    assert summary["pressure_multiplier_mean_abs_rel_error_yu2002_fit"] < 0.10
    assert summary["pressure_multiplier_mean_abs_rel_error_chisholm"] > 0.25
    assert summary["htc_fig10_digitized_mean_abs_rel_error"] < 0.10
    assert summary["chf_digitized_monotonic_decrease_with_quality"] is True
    for name in summary["outputs"]:
        assert (tmp_path / name).exists()


def test_hx_imposed_duty_reference_case_has_nonuniform_integration_inputs():
    case = imposed_duty_reference_case()
    dz = np.diff(case.z_edges_m)
    diameter = np.asarray(case.hydraulic_diameter_m, dtype=float)
    heat_per_segment = np.asarray(case.heat_per_segment_W, dtype=float)

    assert dz.size == 48
    assert np.ptp(dz) > 0.0
    assert np.ptp(diameter) > 0.0
    assert np.all(heat_per_segment > 0.0)


def test_hx_imposed_duty_report_generates_adapter_artifacts(tmp_path):
    summary = generate_imposed_duty_report(tmp_path)

    assert summary["n_cells"] == 48
    assert summary["heat_rate_W"] > 0.0
    assert summary["pressure_drop_Pa"] > 0.0
    assert summary["outlet_T_K"] > summary["inlet_T_K"]
    assert summary["energy_residual_ok"] is True
    assert summary["chf_margin_below_limit"] is False
    for name in summary["outputs"]:
        assert (tmp_path / name).exists()


def test_liquid_validation_matrix_generates_readiness_gate(tmp_path):
    report = run_validation_matrix(
        output=tmp_path / "liquid_validation_matrix.json",
        artifact_root=tmp_path,
    )

    assert report["checks"]["all_passed"] is True
    assert report["scope"]["steady_governing_state"].startswith("p,h")
    assert "transient boiling/liquid finite-volume model" in report["scope"]["not_yet_validated"]
    assert (tmp_path / "liquid_validation_matrix.json").exists()
    assert (tmp_path / "liquid_solver_postprocess_audit.json").exists()

import numpy as np

from hps_combustor.validation import coupled_bangbang_hx as coupled
from hps_combustor.validation.pressurant_bangbang_sizing import (
    FeedDesign,
    PressurantSystemConfig,
    run_design,
)


def test_build_hx_inputs_from_system_history_maps_boundaries():
    system_config = PressurantSystemConfig(t_end_s=0.02, dt_s=0.01)
    design = FeedDesign(
        n_branches=2,
        orifice_diameter_m=3.0e-3,
        valve_equivalent_diameter_m=3.0e-3,
        control_frequency_Hz=20.0,
    )
    _summary, history = run_design(system_config, design)
    cfg = coupled.CoupledBangBangHxConfig(t_end_s=0.02, hx_nodes=4, hx_save_points=3)

    inputs = coupled.build_hx_inputs_from_system_history(cfg, history)
    transient = inputs["transient"]

    assert inputs["combustor"].HX_config == "shellntube"
    assert transient.coolant_momentum_model == "quasi_steady"
    assert transient.schedule_mass_flow_c[0] == (
        float(history["time_s"][0]),
        float(history["helium_mdot_kg_s"][0]),
    )
    assert transient.schedule_p_c_in[-1] == (
        float(history["time_s"][-1]),
        float(history["line_pressure_before_hx_bar"][-1] * 1.0e5),
    )
    assert transient.schedule_p_c_out[-1] == (
        float(history["time_s"][-1]),
        float(history["water_tank_pressure_bar"][-1] * 1.0e5),
    )
    assert transient.schedule_T_c_in[-1] == (
        float(history["time_s"][-1]),
        float(history["supply_temperature_K"][-1]),
    )


def test_build_hx_inputs_can_select_low_mach_momentum():
    system_config = PressurantSystemConfig(t_end_s=0.02, dt_s=0.01)
    design = FeedDesign(
        n_branches=2,
        orifice_diameter_m=3.0e-3,
        valve_equivalent_diameter_m=3.0e-3,
        control_frequency_Hz=20.0,
    )
    _summary, history = run_design(system_config, design)
    cfg = coupled.CoupledBangBangHxConfig(
        t_end_s=0.02,
        hx_nodes=4,
        hx_save_points=3,
        coolant_momentum_model="low_mach",
        low_mach_mdot_cap_kg_s=0.42,
    )

    inputs = coupled.build_hx_inputs_from_system_history(cfg, history)

    assert inputs["transient"].coolant_momentum_model == "low_mach"
    assert inputs["transient"].transient_coolant_outlet_pressure == history["water_tank_pressure_bar"][0] * 1.0e5
    expected_cap = np.minimum(np.maximum(history["helium_mdot_kg_s"], 0.0), 0.42)
    assert inputs["coolant"].mass_flow_c == float(np.nanmax(expected_cap))
    assert inputs["transient"].schedule_mass_flow_c[0] == (float(history["time_s"][0]), float(expected_cap[0]))
    assert inputs["transient"].schedule_mass_flow_c[-1] == (float(history["time_s"][-1]), float(expected_cap[-1]))


def test_build_hx_inputs_can_select_helical_low_mach_geometry():
    system_config = PressurantSystemConfig(t_end_s=0.02, dt_s=0.01)
    design = FeedDesign(
        n_branches=2,
        orifice_diameter_m=3.0e-3,
        valve_equivalent_diameter_m=3.0e-3,
        control_frequency_Hz=20.0,
    )
    _summary, history = run_design(system_config, design)
    cfg = coupled.CoupledBangBangHxConfig(
        t_end_s=0.02,
        hx_config="shellnHelicalTube",
        hx_nodes=7,
        coolant_momentum_model="low_mach",
        helical_inner_diameter_m=14.0e-3,
        helical_wall_thickness_m=1.25e-3,
        helical_pipe_length_m=12.0,
        helical_centerline_diameter_m=65.0e-3,
        helical_coil_gap_m=4.0e-3,
    )

    inputs = coupled.build_hx_inputs_from_system_history(cfg, history)

    assert inputs["combustor"].HX_config == "shellnHelicalTube"
    assert inputs["combustor"].material_HX == "ST316L"
    assert inputs["combustor"].Dh_coil == 14.0e-3
    assert inputs["combustor"].thickness_coil_wall == 1.25e-3
    assert inputs["combustor"].Nusselt_correction == 0.28
    assert inputs["transient"].n_axial == 7
    assert inputs["transient"].skip_steady_reference_probe is True
    assert inputs["transient"].coolant_momentum_model == "low_mach"
    tube_od = 14.0e-3 + 2.0 * 1.25e-3
    expected_gap = 0.5 * (inputs["combustor"].inner_diameter - 65.0e-3 - tube_od)
    assert inputs["combustor"].gap_shell2coil == expected_gap
    assert inputs["numerical"].L_HX_max > 0.0


def test_hot_gas_schedule_tracks_helium_flow_and_preserves_of():
    cfg = coupled.CoupledBangBangHxConfig(
        lox_mdot_kg_s=0.06,
        diesel_mdot_kg_s=0.03,
        hot_gas_nominal_helium_mdot_kg_s=0.15,
        hot_gas_max_mdot_kg_s=0.09,
    )
    schedule = coupled.hot_gas_schedule_from_helium(
        cfg,
        np.array([0.0, 1.0, 2.0]),
        np.array([0.0, 0.075, 0.15]),
    )

    np.testing.assert_allclose(schedule["hot_gas_mdot_kg_s"], [0.0, 0.045, 0.09])
    np.testing.assert_allclose(schedule["lox_mdot_kg_s"], 2.0 * schedule["diesel_mdot_kg_s"])


def test_run_coupled_case_writes_outputs_with_mocked_hx(monkeypatch, tmp_path):
    t = np.array([0.0, 0.01, 0.02])
    x = np.array([0.0, 0.1])
    time_series = {
        "t": t,
        "x": x,
        "fields": {
            "T_g": np.array([[1000.0, 900.0], [990.0, 890.0], [980.0, 880.0]]),
            "T_c": np.array([[100.0, 120.0], [101.0, 121.0], [102.0, 122.0]]),
            "T_wg": np.array([[300.0, 310.0], [305.0, 315.0], [310.0, 320.0]]),
            "Tbar": np.array([[298.0, 306.0], [303.0, 311.0], [308.0, 316.0]]),
            "T_wc": np.array([[295.0, 300.0], [300.0, 305.0], [305.0, 310.0]]),
            "p_c": np.array([[80.0e5, 75.0e5], [79.0e5, 74.0e5], [78.0e5, 73.0e5]]),
        },
        "scalars": {
            "T_c_out": np.array([120.0, 121.0, 122.0]),
            "T_g_out": np.array([900.0, 890.0, 880.0]),
            "T_wall_max": np.array([310.0, 315.0, 320.0]),
            "T_wall_min": np.array([295.0, 300.0, 305.0]),
            "Q_hot_kW": np.array([10.0, 11.0, 12.0]),
            "dp_g_total_Pa": np.array([1000.0, 1200.0, 1400.0]),
            "dp_shell_total_Pa": np.array([5000.0, 5500.0, 6000.0]),
            "mdot_c_inlet_face": np.array([0.1, 0.1, 0.1]),
            "mdot_c_outlet_face": np.array([0.09, 0.095, 0.1]),
        },
    }

    class DummySolver:
        pass

    solver = DummySolver()
    solver.time_series = time_series

    def fake_run_transient(_inputs):
        return solver, {"T_c_out_final": 122.0}

    monkeypatch.setattr(coupled, "run_transient", fake_run_transient)

    cfg = coupled.CoupledBangBangHxConfig(
        output_dir=str(tmp_path),
        t_end_s=0.02,
        system_dt_s=0.01,
        hx_nodes=2,
        hx_save_points=3,
    )
    payload = coupled.run_coupled_case(cfg)

    assert payload["coupled_diagnostics"]["hx_final_T_c_out_K"] == 122.0
    assert payload["coupled_diagnostics"]["hx_max_flowing_T_c_out_K"] == 122.0
    assert payload["coupled_diagnostics"]["hx_peak_T_wall_hot_face_K"] == 320.0
    assert payload["coupled_diagnostics"]["hx_peak_T_wall_hot_face_time_s"] == 0.02
    assert payload["coupled_diagnostics"]["hx_max_hot_gas_pressure_drop_bar"] == 0.014
    assert payload["coupled_diagnostics"]["hx_max_shell_pressure_drop_estimate_bar"] == 0.06
    assert "hx_max_coolant_pressure_drop_bar" not in payload["coupled_diagnostics"]
    assert "hx_max_coolant_thermodynamic_pressure_span_bar" in payload["coupled_diagnostics"]
    assert (tmp_path / "system_timeseries.csv").exists()
    assert (tmp_path / "hx_boundary_schedule.csv").exists()
    assert (tmp_path / "water_outlet_orifice_sweep.json").exists()
    assert (tmp_path / "coupled_timeseries.csv").exists()
    assert (tmp_path / "hx_transient_timeseries.npz").exists()
    assert (tmp_path / "coupled_dashboard.html").exists()
    assert "0D pressurant/feed/tank surrogate" in (tmp_path / "README.md").read_text(encoding="utf-8")


def test_run_coupled_case_system_only_skips_hx(tmp_path):
    cfg = coupled.CoupledBangBangHxConfig(
        output_dir=str(tmp_path),
        t_end_s=0.02,
        system_dt_s=0.01,
        run_hx=False,
    )

    payload = coupled.run_coupled_case(cfg)

    assert payload["hx_runtime_s"] is None
    assert (tmp_path / "system_timeseries.csv").exists()
    assert (tmp_path / "water_outlet_orifice_sweep.json").exists()
    assert not (tmp_path / "coupled_timeseries.csv").exists()
    assert (tmp_path / "coupled_dashboard.html").exists()


def test_settings_round_trip(tmp_path):
    path = tmp_path / "settings.json"
    cfg = coupled.CoupledBangBangHxConfig(
        output_dir=str(tmp_path / "out"),
        t_end_s=1.5,
        hx_nodes=6,
        run_hx=False,
    )
    system = PressurantSystemConfig(helium_tank_volume_m3=0.3, t_end_s=1.5)
    design = FeedDesign(
        n_branches=3,
        orifice_diameter_m=3.0e-3,
        valve_equivalent_diameter_m=3.5e-3,
        control_frequency_Hz=40.0,
    )

    coupled.write_settings(path, cfg, system, design)
    loaded_cfg, loaded_system, loaded_design = coupled.load_settings(path)

    assert loaded_cfg == cfg
    assert loaded_system == system
    assert loaded_design == design

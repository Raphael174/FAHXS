import numpy as np

from hps_combustor.validation.pressurant_bangbang_sizing import (
    FeedDesign,
    PressurantSystemConfig,
    branch_effective_cda,
    commanded_open_branches,
    helium_orifice_mdot,
    merged_feed_mdot,
    run_design,
    water_exit_flow,
    water_exit_orifice_area,
)


def test_water_exit_orifice_is_sized_to_target_flow():
    config = PressurantSystemConfig()
    area = water_exit_orifice_area(config)
    q = water_exit_flow(config, config.target_water_tank_pressure_Pa, area)

    np.testing.assert_allclose(q, config.target_water_flow_m3_s, rtol=1e-12)


def test_water_exit_orifice_can_be_fixed_by_diameter():
    config = PressurantSystemConfig(water_exit_orifice_diameter_m=20.0e-3)
    area = water_exit_orifice_area(config)

    np.testing.assert_allclose(area, np.pi * (20.0e-3) ** 2 / 4.0)


def test_series_branch_cda_is_below_each_restriction():
    config = PressurantSystemConfig(feed_Cd=0.8)
    design = FeedDesign(
        n_branches=2,
        orifice_diameter_m=3.0e-3,
        valve_equivalent_diameter_m=4.0e-3,
        control_frequency_Hz=10.0,
    )

    cda = branch_effective_cda(config, design)
    cda_orifice = config.feed_Cd * np.pi * design.orifice_diameter_m**2 / 4.0
    cda_valve = config.feed_Cd * np.pi * design.valve_equivalent_diameter_m**2 / 4.0

    assert cda < cda_orifice
    assert cda < cda_valve


def test_helium_orifice_flow_requires_positive_pressure_ratio():
    cda = 1.0e-6
    assert helium_orifice_mdot(
        upstream_pressure=70.0e5,
        downstream_pressure=80.0e5,
        upstream_temperature=100.0,
        cda=cda,
    ) == 0.0
    assert helium_orifice_mdot(
        upstream_pressure=400.0e5,
        downstream_pressure=70.0e5,
        upstream_temperature=100.0,
        cda=cda,
    ) > 0.0


def test_staged_bangbang_command_uses_available_branches():
    assert commanded_open_branches(82.0e5, target_pressure=80.0e5, hysteresis=0.25e5, n_branches=3) == 0
    assert commanded_open_branches(80.0e5, target_pressure=80.0e5, hysteresis=0.25e5, n_branches=3) == 2
    assert commanded_open_branches(79.0e5, target_pressure=80.0e5, hysteresis=0.25e5, n_branches=3) == 3


def test_common_hx_pressure_loss_reduces_merged_feed_flow():
    base = PressurantSystemConfig(hx_pressure_loss_nominal_Pa=0.0)
    lossy = PressurantSystemConfig(
        hx_pressure_loss_nominal_Pa=5.0e5,
        hx_nominal_helium_mdot_kg_s=0.2,
    )
    design = FeedDesign(
        n_branches=2,
        orifice_diameter_m=3.5e-3,
        valve_equivalent_diameter_m=4.0e-3,
        control_frequency_Hz=20.0,
    )
    cda = branch_effective_cda(base, design)

    mdot_base = merged_feed_mdot(
        base,
        open_branches=2,
        supply_pressure=400.0e5,
        supply_temperature=100.0,
        tank_pressure=70.0e5,
        branch_cda=cda,
    )
    mdot_lossy = merged_feed_mdot(
        lossy,
        open_branches=2,
        supply_pressure=400.0e5,
        supply_temperature=100.0,
        tank_pressure=70.0e5,
        branch_cda=cda,
    )

    assert mdot_lossy < mdot_base


def test_short_bangbang_design_tracks_pressure_and_water_flow():
    config = PressurantSystemConfig(t_end_s=5.0, dt_s=0.005)
    design = FeedDesign(
        n_branches=2,
        orifice_diameter_m=3.5e-3,
        valve_equivalent_diameter_m=4.0e-3,
        control_frequency_Hz=20.0,
    )

    summary, _history = run_design(config, design)

    assert abs(summary.mean_pressure_bar - 70.0) < 1.0
    assert abs(summary.mean_water_flow_L_s - 30.0) < 0.3
    assert summary.helium_used_kg > 0.0
    assert summary.final_supply_pressure_bar < 400.0
    assert summary.final_supply_temperature_K < 100.0

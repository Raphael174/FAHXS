import numpy as np

from hps_combustor.transient_core import (
    AxialGrid,
    TransientStateLayout,
    build_helical_core_geometry,
    build_shelltube_core_geometry,
    collect_transient_schedule_times,
    conductance_from_h,
    conservative_mass_energy_step,
    coolant_heat_capacity_J_K,
    coolprop_state_from_mass_energy,
    ShellTubeFluidProperties,
    ShellTubeHotGasMarch,
    ShellTubeShellFilm,
    ShellTubeWallFlux,
    energy_audit,
    equilibrium_gas_state_provider,
    fixed_time_grid,
    flow_direction_from_config,
    fpv_gas_state_provider,
    helical_coolant_film,
    helical_wall_flux,
    interp_schedule,
    implicit_upwind_step,
    implicit_wall_coolant_step,
    initial_mass_energy_from_TP,
    integrate_wall_coolant_fixed_step,
    oxygen_gas_state_provider,
    quasi_steady_face_mdot,
    residence_time_s,
    run_shelltube_transient_core,
    schedule_times,
    semi_implicit_wall_compressible_coolant_step,
    shelltube_conductance_from_h,
    shelltube_coolant_heat_capacity_J_K,
    shelltube_flow_direction,
    shelltube_hot_gas_march,
    shelltube_shell_film,
    shelltube_step_inputs,
    shelltube_tube_gas_film,
    shelltube_wall_flux,
    shelltube_wall_heat_capacity_J_K,
    timescale_audit,
    WallCoolantStepInputs,
    wall_heat_capacity_J_K,
    wall_time_constant_s,
)


def test_transient_state_layout_round_trip():
    layout = TransientStateLayout(3)
    wall = np.array([300.0, 310.0, 320.0])
    coolant = np.array([100.0, 110.0, 120.0])

    y = layout.pack(wall, coolant)
    state = layout.unpack(y)

    np.testing.assert_allclose(state.Tbar_wall, wall)
    np.testing.assert_allclose(state.T_coolant, coolant)
    assert layout.size == 6


def test_axial_grid_uniform_geometry_and_indices():
    grid = AxialGrid.uniform(
        length=2.0,
        n_cells=4,
        coolant_area=0.01,
        wall_area=0.002,
        hot_perimeter=0.2,
        coolant_perimeter=0.1,
        flow_direction=-1,
    )

    np.testing.assert_allclose(grid.dx, np.full(4, 0.5))
    np.testing.assert_allclose(grid.x_centers, [0.25, 0.75, 1.25, 1.75])
    np.testing.assert_allclose(grid.coolant_volume, np.full(4, 0.005))
    np.testing.assert_allclose(grid.wall_volume, np.full(4, 0.001))
    assert grid.length == 2.0
    assert grid.inlet_index == 3
    assert grid.outlet_index == 0


def test_axial_grid_converts_linear_heat_to_cell_heat():
    grid = AxialGrid(
        x_edges=np.array([0.0, 0.1, 0.4]),
        coolant_area=np.array([1.0, 1.0]),
        wall_area=np.array([1.0, 1.0]),
        hot_perimeter=np.array([1.0, 1.0]),
        coolant_perimeter=np.array([1.0, 1.0]),
    )

    np.testing.assert_allclose(
        grid.heat_rate_from_linear(np.array([10.0, 20.0])),
        np.array([1.0, 6.0]),
    )


def test_axial_grid_rejects_nonphysical_geometry():
    with np.testing.assert_raises(ValueError):
        AxialGrid.uniform(
            length=1.0,
            n_cells=2,
            coolant_area=0.0,
            wall_area=1.0,
            hot_perimeter=1.0,
            coolant_perimeter=1.0,
        )


def test_diagnostics_residence_time_and_zero_flow():
    rho = np.array([2.0, 4.0])
    volume = np.array([0.5, 0.25])

    assert residence_time_s(rho, volume, mdot=0.5) == 4.0
    assert np.isinf(residence_time_s(rho, volume, mdot=0.0))


def test_diagnostics_wall_time_constant_handles_insulated_cells():
    tau = wall_time_constant_s(
        wall_heat_capacity=np.array([100.0, 50.0]),
        hot_conductance_W_K=np.array([5.0, 0.0]),
        coolant_conductance_W_K=np.array([15.0, 0.0]),
    )

    np.testing.assert_allclose(tau[0], 5.0)
    assert np.isinf(tau[1])


def test_diagnostics_energy_audit_scales_residual():
    good = energy_audit(
        residual_J=1.0e-3,
        heat_added_J=1000.0,
        advective_energy_in_J=500.0,
        advective_energy_out_J=400.0,
        relative_tol=1.0e-5,
    )
    bad = energy_audit(
        residual_J=1.0,
        heat_added_J=1000.0,
        relative_tol=1.0e-5,
    )

    assert good.passes
    assert good.scale_J == 1000.0
    assert not bad.passes


def test_diagnostics_timescale_audit_flags_slow_coolant_not_fast_hot_gas():
    audit = timescale_audit(
        coolant_residence_s=0.2,
        wall_tau_s=np.array([1.0, 2.0]),
        hot_residence_s=0.0006,
        boundary_tau_s=0.5,
        warning_ratio=0.05,
    )

    assert audit.coolant_to_wall_ratio == 0.2
    assert audit.coolant_to_boundary_ratio == 0.4
    assert not audit.coolant_quasi_steady_ok
    assert audit.hot_quasi_steady_ok


def test_helical_adapter_builds_total_parallel_geometry():
    geom = build_helical_core_geometry(
        pipe_length=6.0,
        n_cells=3,
        tube_inner_diameter=0.004,
        wall_thickness=0.001,
        n_parallel=2,
        flow_config="counter",
    )
    grid = geom.grid

    assert grid.flow_direction == -1
    assert grid.inlet_index == 2
    assert grid.outlet_index == 0
    assert geom.tube_outer_diameter == 0.006
    np.testing.assert_allclose(grid.dx, np.full(3, 2.0))
    np.testing.assert_allclose(
        grid.coolant_area,
        np.full(3, 2.0 * np.pi * 0.004**2 / 4.0),
    )
    np.testing.assert_allclose(
        grid.wall_area,
        np.full(3, 2.0 * np.pi * (0.006**2 - 0.004**2) / 4.0),
    )
    np.testing.assert_allclose(grid.hot_perimeter, np.full(3, 2.0 * np.pi * 0.006))
    np.testing.assert_allclose(grid.coolant_perimeter, np.full(3, 2.0 * np.pi * 0.004))


def test_helical_adapter_inventory_and_conductance_helpers():
    geom = build_helical_core_geometry(
        pipe_length=1.0,
        n_cells=2,
        tube_inner_diameter=0.01,
        wall_thickness=0.001,
        n_parallel=1,
    )
    grid = geom.grid

    Cw = wall_heat_capacity_J_K(grid, density=8000.0, cp=500.0)
    Cc = coolant_heat_capacity_J_K(grid, density=np.array([2.0, 3.0]), cp=5000.0)
    G = conductance_from_h(grid, h_coolant=np.array([100.0, 200.0]))

    np.testing.assert_allclose(Cw, 8000.0 * 500.0 * grid.wall_volume)
    np.testing.assert_allclose(Cc, np.array([2.0, 3.0]) * 5000.0 * grid.coolant_volume)
    np.testing.assert_allclose(G, np.array([100.0, 200.0]) * grid.coolant_perimeter * grid.dx)


def test_helical_adapter_rejects_unknown_flow_config():
    assert flow_direction_from_config("co") == 1
    assert flow_direction_from_config("counter") == -1
    with np.testing.assert_raises(ValueError):
        flow_direction_from_config("sideways")


def test_helical_coolant_film_uses_total_area_and_existing_dispatchers(monkeypatch):
    from hps_combustor.transient_core import adapters_helical

    friction_calls = []
    nusselt_calls = []

    def fake_friction(selector, Re, Dh, Rc, roughness, x, error_factor, corrCoeffs):
        friction_calls.append((selector, Re, Dh, Rc, roughness, x, error_factor, corrCoeffs))
        return 0.02 * error_factor

    def fake_nusselt(selector, Re, Pr, d, R, f_fd, x, error_factor, corrCoeffs):
        nusselt_calls.append((selector, Re, Pr, d, R, f_fd, x, error_factor, corrCoeffs))
        return 100.0 * error_factor

    monkeypatch.setattr(adapters_helical, "dispatch_friction_coil", fake_friction)
    monkeypatch.setattr(adapters_helical, "dispatch_nu_coil", fake_nusselt)

    geom = build_helical_core_geometry(
        pipe_length=1.0,
        n_cells=2,
        tube_inner_diameter=0.01,
        wall_thickness=0.001,
        n_parallel=2,
    )
    film = helical_coolant_film(
        geom,
        mdot_total=0.2,
        rho=np.array([2.0, 4.0]),
        mu=1.0e-5,
        k=np.array([0.1, 0.2]),
        cp=5000.0,
        friction_selector="fake_f",
        nusselt_selector="fake_nu",
        roughness=1.0e-6,
        coil_radius=0.05,
        corrCoeffs=object(),
        x_for_developing=np.array([0.1, 0.2]),
        friction_error_factor=1.5,
        nusselt_error_factor=2.0,
    )

    expected_u = 0.2 / (np.array([2.0, 4.0]) * geom.grid.coolant_area)
    expected_re = np.array([2.0, 4.0]) * expected_u * 0.01 / 1.0e-5
    expected_pr = 5000.0 * 1.0e-5 / np.array([0.1, 0.2])
    expected_h = np.array([200.0, 200.0]) * np.array([0.1, 0.2]) / 0.01

    np.testing.assert_allclose(film.velocity_m_s, expected_u)
    np.testing.assert_allclose(film.reynolds, expected_re)
    np.testing.assert_allclose(film.prandtl, expected_pr)
    np.testing.assert_allclose(film.friction_factor, np.full(2, 0.03))
    np.testing.assert_allclose(film.nusselt, np.full(2, 200.0))
    np.testing.assert_allclose(film.h_W_m2K, expected_h)
    np.testing.assert_allclose(
        film.conductance_W_K,
        expected_h * geom.grid.coolant_perimeter * geom.grid.dx,
    )
    assert [call[0] for call in friction_calls] == ["fake_f", "fake_f"]
    assert [call[0] for call in nusselt_calls] == ["fake_nu", "fake_nu"]


def test_helical_wall_flux_converts_single_tube_per_length_to_total_cell_heat(monkeypatch):
    from hps_combustor.transient_core import adapters_helical

    calls = []

    class FakeConduction:
        def __init__(self, **kwargs):
            calls.append(kwargs)
            self.dx = kwargs["dx"]

        def fluxes_at_Tbar(self, T_bar, h_g_rad=None):
            return {
                "dq_hot__dx": 10.0 * self.dx,
                "dq_cold__dx": 7.0 * self.dx,
                "T_wg": T_bar + 5.0,
                "T_wc": T_bar - 2.0,
                "h_g_rad": h_g_rad,
                "q_w_rad": 3.0,
                "k_w": 20.0,
            }

    monkeypatch.setattr(
        adapters_helical,
        "OneDimensionalSteadyConduction_ShellnHelicalTube",
        FakeConduction,
    )

    geom = build_helical_core_geometry(
        pipe_length=1.0,
        n_cells=2,
        tube_inner_diameter=0.01,
        wall_thickness=0.001,
        n_parallel=3,
    )
    result = helical_wall_flux(
        geom,
        Tbar_wall=np.array([300.0, 310.0]),
        T_coolant=np.array([100.0, 110.0]),
        T_gas=np.array([900.0, 850.0]),
        h_gas=np.array([50.0, 60.0]),
        h_coolant=np.array([1000.0, 1200.0]),
        wall_conductivity_at_T=lambda T: 20.0,
        h_g_rad=np.array([5.0, 6.0]),
    )

    # Fake dq__dx is 10*dx and 7*dx for a single tube. The adapter converts to
    # total W per cell by multiplying by n_parallel*dx.
    np.testing.assert_allclose(result.dq_hot_per_length_W_m, np.full(2, 5.0))
    np.testing.assert_allclose(result.dq_cold_per_length_W_m, np.full(2, 3.5))
    np.testing.assert_allclose(result.hot_heat_W, np.full(2, 5.0 * 3.0 * 0.5))
    np.testing.assert_allclose(result.cold_heat_W, np.full(2, 3.5 * 3.0 * 0.5))
    np.testing.assert_allclose(result.T_wg, [305.0, 315.0])
    np.testing.assert_allclose(result.T_wc, [298.0, 308.0])
    np.testing.assert_allclose(result.h_g_rad, [5.0, 6.0])
    np.testing.assert_allclose(result.q_w_rad, [3.0, 3.0])
    np.testing.assert_allclose(result.k_wall, [20.0, 20.0])
    assert calls[0]["hot_side"] == "outer"
    assert calls[0]["rad_enabled"] is False


def test_zero_flow_soak_matches_lumped_cell_energy():
    T0 = np.array([100.0, 200.0, 300.0])
    rho = np.full(3, 2.0)
    cp = np.full(3, 10.0)
    volume = np.full(3, 0.5)
    heat = np.array([10.0, -5.0, 0.0])

    result = implicit_upwind_step(
        T0, rho, cp, volume, mdot=0.0, T_inlet=90.0, heat_W=heat, dt=2.0
    )

    expected = T0 + 2.0 * heat / (rho * volume * cp)
    np.testing.assert_allclose(result.T_new, expected)
    assert abs(result.energy_residual_J) < 1e-10


def test_implicit_upwind_positive_flow_preserves_uniform_field_without_heat():
    T0 = np.full(5, 120.0)
    props = np.ones(5)

    result = implicit_upwind_step(
        T0, rho=props, cp=props, volume=props, mdot=0.1,
        T_inlet=120.0, heat_W=np.zeros(5), dt=0.5, flow_direction=1
    )

    np.testing.assert_allclose(result.T_new, T0)
    assert abs(result.energy_residual_J) < 1e-12


def test_negative_direction_uses_last_cell_as_inlet_side():
    T0 = np.full(4, 300.0)
    props = np.ones(4)

    result = implicit_upwind_step(
        T0, rho=props, cp=props, volume=props, mdot=1.0,
        T_inlet=100.0, heat_W=np.zeros(4), dt=1.0, flow_direction=-1
    )

    assert result.T_new[-1] < result.T_new[-2] < result.T_new[-3] < result.T_new[-4]
    assert result.T_outlet == result.T_new[0]


def test_constant_wall_source_converges_toward_exponential_duct_solution():
    n = 200
    length = 1.0
    dx = length / n
    perimeter = 0.2
    area = 0.01
    rho = np.full(n, 1.0)
    cp = np.full(n, 1000.0)
    volume = np.full(n, area * dx)
    mdot = 0.05
    h = 50.0
    T_wall = 400.0
    T_in = 300.0
    T = np.full(n, T_in)

    # March pseudo-time until the steady implicit-upwind duct profile stops moving.
    for _ in range(600):
        heat = h * perimeter * dx * (T_wall - T)
        T_next = implicit_upwind_step(
            T, rho, cp, volume, mdot=mdot, T_inlet=T_in,
            heat_W=heat, dt=0.02, flow_direction=1
        ).T_new
        if np.max(np.abs(T_next - T)) < 1e-8:
            T = T_next
            break
        T = T_next

    x = (np.arange(n) + 0.5) * dx
    expected = T_wall - (T_wall - T_in) * np.exp(-h * perimeter * x / (mdot * cp[0]))
    np.testing.assert_allclose(T, expected, rtol=2e-3, atol=2e-2)


def test_single_step_constant_property_energy_balance_with_heat_and_flow():
    T0 = np.array([100.0, 120.0, 140.0, 160.0])
    rho = np.full(4, 4.0)
    cp = np.full(4, 5.0)
    volume = np.full(4, 0.25)
    heat = np.array([1.0, 2.0, 3.0, 4.0])

    result = implicit_upwind_step(
        T0, rho, cp, volume, mdot=0.2, T_inlet=90.0,
        heat_W=heat, dt=0.25, flow_direction=1
    )

    assert abs(result.energy_residual_J) < 1e-10


def test_coolprop_mass_energy_state_round_trip_from_temperature_pressure():
    T = np.array([90.0, 120.0])
    p = np.array([8.0e6, 7.5e6])
    volume = np.array([1.0e-4, 2.0e-4])

    mass, U = initial_mass_energy_from_TP(T, p, volume, "Helium")
    state = coolprop_state_from_mass_energy(mass, U, volume, "Helium")

    np.testing.assert_allclose(state.temperature, T, rtol=1e-10, atol=1e-7)
    np.testing.assert_allclose(state.pressure, p, rtol=1e-9, atol=1e-1)
    np.testing.assert_allclose(state.density, mass / volume)
    assert np.all(state.specific_enthalpy_J_kg > state.specific_internal_energy_J_kg)
    assert np.all(state.cp > 0.0)
    assert np.all(state.mu > 0.0)
    assert np.all(state.k > 0.0)


def test_quasi_steady_face_mdot_responds_to_pressure_and_closed_inlet():
    pressure = np.array([9.0e6, 8.5e6, 8.0e6])
    density = np.array([10.0, 9.0, 8.0])
    resistance = np.array([2.0e12, 2.0e12])

    face = quasi_steady_face_mdot(
        pressure,
        density,
        resistance,
        inlet_pressure=9.0e6,
        outlet_pressure=7.0e6,
        inlet_resistance=1.0e12,
        outlet_resistance=1.0e12,
        inlet_enabled=False,
    )

    assert face[0] == 0.0
    assert np.all(face[1:] > 0.0)
    np.testing.assert_allclose(
        face[1],
        np.sqrt(0.5 * (10.0 + 9.0) * 0.5e6 / 2.0e12),
    )


def test_compressible_mass_energy_step_allows_residual_outflow_after_inlet_shutoff():
    mass = np.array([0.10, 0.10])
    u = np.array([1000.0, 1200.0])
    U = mass * u
    h = np.array([1100.0, 1300.0])
    face_mdot = np.array([0.0, 0.01, 0.01])

    result = conservative_mass_energy_step(
        mass,
        U,
        h,
        face_mdot,
        heat_W=np.zeros(2),
        dt=1.0,
        inlet_enthalpy_J_kg=900.0,
    )

    np.testing.assert_allclose(result.mass_new, [0.09, 0.10])
    assert result.advective_energy_in_J == 0.0
    assert result.advective_energy_out_J == 13.0
    assert abs(result.mass_residual_kg) < 1e-15
    assert abs(result.energy_residual_J) < 1e-12


def test_compressible_mass_energy_step_conserves_closed_domain_with_heat():
    mass = np.array([1.0, 1.0, 1.0])
    U = np.array([100.0, 200.0, 300.0])
    h = np.array([10.0, 20.0, 30.0])
    face_mdot = np.array([0.0, 0.2, -0.1, 0.0])
    heat = np.array([1.0, 2.0, 3.0])

    result = conservative_mass_energy_step(
        mass,
        U,
        h,
        face_mdot,
        heat_W=heat,
        dt=0.5,
        inlet_enthalpy_J_kg=0.0,
    )

    assert abs(np.sum(result.mass_new) - np.sum(mass)) < 1e-15
    np.testing.assert_allclose(
        np.sum(result.internal_energy_new_J) - np.sum(U),
        np.sum(heat) * 0.5,
    )
    assert abs(result.energy_residual_J) < 1e-12


def test_wall_compressible_coolant_step_conserves_total_energy_closed_domain():
    wall = np.array([400.0, 350.0])
    Cw = np.array([10.0, 20.0])
    mass = np.array([1.0, 1.0])
    U = np.array([1000.0, 1200.0])
    Tc = np.array([300.0, 310.0])
    h = np.array([1100.0, 1300.0])
    face_mdot = np.array([0.0, 0.0, 0.0])
    Qhot = np.array([100.0, 50.0])
    G = np.array([5.0, 4.0])

    result = semi_implicit_wall_compressible_coolant_step(
        wall,
        Cw,
        mass,
        U,
        Tc,
        h,
        face_mdot,
        Qhot,
        G,
        dt=0.2,
        inlet_enthalpy_J_kg=1000.0,
    )

    assert np.all(result.heat_wall_to_coolant_W > 0.0)
    assert np.all(result.coolant.internal_energy_new_J > U)
    assert abs(result.coolant.mass_residual_kg) < 1e-15
    assert abs(result.total_energy_residual_J) < 1e-10


def test_wall_compressible_coolant_step_keeps_residual_outflow_when_inlet_closed():
    wall = np.array([300.0, 300.0])
    Cw = np.array([100.0, 100.0])
    mass = np.array([0.10, 0.10])
    U = np.array([100.0, 100.0])
    Tc = np.array([100.0, 100.0])
    h = np.array([1000.0, 1200.0])
    face_mdot = np.array([0.0, 0.01, 0.01])

    result = semi_implicit_wall_compressible_coolant_step(
        wall,
        Cw,
        mass,
        U,
        Tc,
        h,
        face_mdot,
        hot_heat_W=np.zeros(2),
        wall_to_coolant_conductance_W_per_K=np.zeros(2),
        dt=1.0,
        inlet_enthalpy_J_kg=900.0,
    )

    np.testing.assert_allclose(result.coolant.mass_new, [0.09, 0.10])
    assert result.coolant.advective_energy_in_J == 0.0
    assert result.coolant.advective_energy_out_J == 12.0
    assert abs(result.total_energy_residual_J) < 1e-12


def test_wall_coolant_zero_flow_conserves_wall_coolant_energy():
    wall = np.array([400.0, 300.0])
    coolant = np.array([100.0, 120.0])
    Cw = np.array([20.0, 30.0])
    Cc = np.array([10.0, 15.0])
    cp = np.array([5000.0, 5000.0])
    G = np.array([8.0, 4.0])
    Qhot = np.array([50.0, 25.0])

    result = implicit_wall_coolant_step(
        wall,
        coolant,
        Cw,
        Cc,
        cp,
        mdot_coolant=0.0,
        T_coolant_inlet=90.0,
        hot_heat_W=Qhot,
        wall_to_coolant_conductance_W_per_K=G,
        dt=0.2,
    )

    dU = (
        result.wall_internal_energy_new_J
        + result.coolant_internal_energy_new_J
        - result.wall_internal_energy_old_J
        - result.coolant_internal_energy_old_J
    )
    assert np.all(result.T_coolant_new > coolant)
    assert np.all(result.T_wall_new < wall + Qhot * 0.2 / Cw)
    assert abs(dU - np.sum(Qhot) * 0.2) < 1e-10
    assert abs(result.energy_residual_J) < 1e-10


def test_wall_coolant_flow_direction_changes_outlet_side():
    wall = np.array([300.0, 310.0, 320.0])
    coolant = np.array([100.0, 100.0, 100.0])
    Cw = np.full(3, 1.0e6)
    Cc = np.full(3, 50.0)
    cp = np.full(3, 1000.0)
    G = np.full(3, 0.0)
    Qhot = np.zeros(3)

    forward = implicit_wall_coolant_step(
        wall,
        coolant,
        Cw,
        Cc,
        cp,
        mdot_coolant=0.1,
        T_coolant_inlet=200.0,
        hot_heat_W=Qhot,
        wall_to_coolant_conductance_W_per_K=G,
        dt=1.0,
        flow_direction=1,
    )
    reverse = implicit_wall_coolant_step(
        wall,
        coolant,
        Cw,
        Cc,
        cp,
        mdot_coolant=0.1,
        T_coolant_inlet=200.0,
        hot_heat_W=Qhot,
        wall_to_coolant_conductance_W_per_K=G,
        dt=1.0,
        flow_direction=-1,
    )

    assert forward.T_coolant_new[0] > forward.T_coolant_new[-1]
    assert reverse.T_coolant_new[-1] > reverse.T_coolant_new[0]
    assert forward.T_coolant_outlet == forward.T_coolant_new[-1]
    assert reverse.T_coolant_outlet == reverse.T_coolant_new[0]
    assert abs(forward.energy_residual_J) < 1e-7
    assert abs(reverse.energy_residual_J) < 1e-7


def test_wall_coolant_steady_limit_matches_resistance_partition_no_flow():
    wall = np.array([500.0])
    coolant = np.array([100.0])
    Cw = np.array([100.0])
    Cc = np.array([100.0])
    cp = np.array([5000.0])
    G = np.array([10.0])
    Qhot = np.array([1000.0])

    result = implicit_wall_coolant_step(
        wall,
        coolant,
        Cw,
        Cc,
        cp,
        mdot_coolant=0.0,
        T_coolant_inlet=90.0,
        hot_heat_W=Qhot,
        wall_to_coolant_conductance_W_per_K=G,
        dt=5000.0,
    )

    expected_delta = Qhot[0] * Cc[0] / (G[0] * (Cw[0] + Cc[0]))
    assert abs(
        (result.T_wall_new[0] - result.T_coolant_new[0]) - expected_delta
    ) < 1.0
    assert abs(result.energy_residual_J) < 1e-6


def test_fixed_step_integrator_matches_repeated_wall_coolant_steps():
    wall0 = np.array([300.0, 320.0])
    coolant0 = np.array([100.0, 110.0])
    Cw = np.array([100.0, 100.0])
    Cc = np.array([50.0, 50.0])
    cp = np.array([1000.0, 1000.0])
    Qhot = np.array([20.0, 10.0])
    G = np.array([5.0, 2.0])

    def inputs(_t, _wall, _coolant):
        return WallCoolantStepInputs(
            wall_heat_capacity=Cw,
            coolant_heat_capacity=Cc,
            coolant_cp=cp,
            mdot_coolant=0.01,
            T_coolant_inlet=90.0,
            hot_heat_W=Qhot,
            wall_to_coolant_conductance_W_per_K=G,
            flow_direction=1,
        )

    result = integrate_wall_coolant_fixed_step(
        T_wall_initial=wall0,
        T_coolant_initial=coolant0,
        t_eval=np.array([0.0, 0.1, 0.2]),
        step_inputs=inputs,
    )

    first = implicit_wall_coolant_step(
        wall0, coolant0, Cw, Cc, cp, 0.01, 90.0, Qhot, G, 0.1
    )
    second = implicit_wall_coolant_step(
        first.T_wall_new,
        first.T_coolant_new,
        Cw,
        Cc,
        cp,
        0.01,
        90.0,
        Qhot,
        G,
        0.1,
    )

    np.testing.assert_allclose(result.T_wall[1], first.T_wall_new)
    np.testing.assert_allclose(result.T_coolant[1], first.T_coolant_new)
    np.testing.assert_allclose(result.T_wall[2], second.T_wall_new)
    np.testing.assert_allclose(result.T_coolant[2], second.T_coolant_new)
    assert abs(result.energy_residual_J[1]) < 1e-9
    assert abs(result.energy_residual_J[2]) < 1e-9


def test_fixed_step_integrator_keeps_zero_flow_as_local_soak():
    wall0 = np.array([300.0])
    coolant0 = np.array([100.0])

    def inputs(_t, _wall, _coolant):
        return WallCoolantStepInputs(
            wall_heat_capacity=np.array([100.0]),
            coolant_heat_capacity=np.array([100.0]),
            coolant_cp=np.array([1000.0]),
            mdot_coolant=0.0,
            T_coolant_inlet=50.0,
            hot_heat_W=np.array([0.0]),
            wall_to_coolant_conductance_W_per_K=np.array([10.0]),
        )

    result = integrate_wall_coolant_fixed_step(
        T_wall_initial=wall0,
        T_coolant_initial=coolant0,
        t_eval=np.array([0.0, 0.1, 0.2, 0.3]),
        step_inputs=inputs,
    )

    assert np.all(result.T_wall[1:, 0] < wall0[0])
    assert np.all(result.T_coolant[1:, 0] > coolant0[0])
    assert np.all(result.T_coolant[:, 0] > 50.0)
    assert abs(np.sum(result.energy_residual_J)) < 1e-8


def test_fixed_time_grid_inserts_schedule_breakpoints_and_eval_times():
    grid = fixed_time_grid(
        t_end=1.0,
        max_step=0.4,
        t_eval=np.array([0.25, 0.8]),
        schedules=(
            ((0.1, 1.0), (0.4, 2.0), (1.2, 3.0)),
            None,
            ((-0.1, 9.0), (0.8, 4.0)),
        ),
    )

    np.testing.assert_allclose(
        grid,
        np.array([0.0, 0.1, 0.25, 0.4, 0.8, 1.0]),
    )


def test_fixed_time_grid_rejects_invalid_step():
    with np.testing.assert_raises(ValueError):
        fixed_time_grid(t_end=1.0, max_step=0.0)


def test_interp_schedule_matches_flat_held_linear_behavior():
    schedule = ((0.0, 10.0), (2.0, 20.0), (4.0, 0.0))

    assert interp_schedule(None, 1.0, 3.0) == 3.0
    assert interp_schedule(schedule, -1.0, 3.0) == 10.0
    assert interp_schedule(schedule, 1.0, 3.0) == 15.0
    assert interp_schedule(schedule, 3.0, 3.0) == 10.0
    assert interp_schedule(schedule, 5.0, 3.0) == 0.0


def test_schedule_times_collects_sorted_unique_filtered_times():
    times = schedule_times(
        ((0.0, 1.0), (0.5, 2.0), (2.0, 3.0)),
        None,
        ((0.5, 4.0), (1.0, 5.0)),
        t_min=0.25,
        t_max=1.5,
    )

    np.testing.assert_allclose(times, np.array([0.5, 1.0]))


def test_collect_transient_schedule_times_from_named_attributes():
    class Transient:
        schedule_mass_flow_c = ((0.0, 1.0), (1.0, 2.0))
        schedule_T_c_in = ((0.5, 90.0),)
        missing = None

    times = collect_transient_schedule_times(
        Transient(),
        ("schedule_mass_flow_c", "schedule_T_c_in", "schedule_p_c_in"),
        t_min=0.0,
        t_max=0.75,
    )

    np.testing.assert_allclose(times, np.array([0.0, 0.5]))


def test_shelltube_adapter_builds_shell_hold_up_and_total_tube_geometry():
    geom = build_shelltube_core_geometry(
        tube_length=0.3,
        n_cells=3,
        shell_inner_diameter=0.1,
        tube_outer_diameter=0.01,
        wall_thickness=0.001,
        n_tubes=10,
        flow_config="counter",
    )
    grid = geom.grid
    expected_shell_area = np.pi * 0.1**2 / 4.0
    expected_displaced = 10.0 * np.pi * 0.01**2 / 4.0

    assert grid.flow_direction == -1
    assert grid.inlet_index == 2
    assert grid.outlet_index == 0
    assert geom.tube_inner_diameter == 0.008
    np.testing.assert_allclose(grid.dx, np.full(3, 0.1))
    np.testing.assert_allclose(grid.coolant_area, np.full(3, expected_shell_area - expected_displaced))
    np.testing.assert_allclose(
        grid.wall_area,
        np.full(3, 10.0 * np.pi * (0.01**2 - 0.008**2) / 4.0),
    )
    np.testing.assert_allclose(grid.hot_perimeter, np.full(3, 10.0 * np.pi * 0.008))
    np.testing.assert_allclose(grid.coolant_perimeter, np.full(3, 10.0 * np.pi * 0.01))


def test_shelltube_adapter_inventory_and_conductance_helpers():
    geom = build_shelltube_core_geometry(
        tube_length=0.2,
        n_cells=2,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
    )
    grid = geom.grid

    Cw = shelltube_wall_heat_capacity_J_K(grid, density=8200.0, cp=450.0)
    Cc = shelltube_coolant_heat_capacity_J_K(grid, density=np.array([3.0, 4.0]), cp=5200.0)
    G = shelltube_conductance_from_h(grid, h_shell=np.array([100.0, 150.0]))

    np.testing.assert_allclose(Cw, 8200.0 * 450.0 * grid.wall_volume)
    np.testing.assert_allclose(Cc, np.array([3.0, 4.0]) * 5200.0 * grid.coolant_volume)
    np.testing.assert_allclose(G, np.array([100.0, 150.0]) * grid.coolant_perimeter * grid.dx)


def test_shelltube_adapter_rejects_invalid_flow_and_negative_hold_up():
    assert shelltube_flow_direction("co") == 1
    assert shelltube_flow_direction("counter") == -1
    with np.testing.assert_raises(ValueError):
        shelltube_flow_direction("sideways")
    with np.testing.assert_raises(ValueError):
        build_shelltube_core_geometry(
            tube_length=0.1,
            n_cells=1,
            shell_inner_diameter=0.01,
            tube_outer_diameter=0.01,
            wall_thickness=0.001,
            n_tubes=5,
        )


def test_shelltube_shell_film_uses_bell_delaware_interface(monkeypatch):
    from hps_combustor.transient_core import adapters_shelltube

    calls = []

    def fake_bell(geom, Re_s, Pr_s, k_s, cp_s, mu_s, mdot_s, mu_ratio=1.0, corrCoeffs=None):
        calls.append((geom, Re_s, Pr_s, k_s, cp_s, mu_s, mdot_s, mu_ratio, corrCoeffs))
        return {"h_shell": 100.0 + 0.01 * Re_s, "dp_shell": 12.0}

    monkeypatch.setattr(adapters_shelltube, "bell_delaware_shell", fake_bell)

    geom = build_shelltube_core_geometry(
        tube_length=0.2,
        n_cells=2,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
    )
    bell_geom = {"S_m": 0.002, "layout": "triangular30"}
    film = shelltube_shell_film(
        geom,
        bell_geom,
        mdot_shell=0.1,
        rho=np.array([2.0, 3.0]),
        mu=np.array([1.0e-5, 2.0e-5]),
        k=np.array([0.1, 0.2]),
        cp=5000.0,
        corrCoeffs=object(),
        mu_ratio=np.array([1.0, 1.1]),
    )

    expected_mass_flux = np.full(2, 0.1 / 0.002)
    expected_re = 0.006 * expected_mass_flux / np.array([1.0e-5, 2.0e-5])
    expected_pr = 5000.0 * np.array([1.0e-5, 2.0e-5]) / np.array([0.1, 0.2])
    expected_h = 100.0 + 0.01 * expected_re

    np.testing.assert_allclose(film.mass_flux_kg_m2s, expected_mass_flux)
    np.testing.assert_allclose(film.reynolds, expected_re)
    np.testing.assert_allclose(film.prandtl, expected_pr)
    np.testing.assert_allclose(film.h_W_m2K, expected_h)
    np.testing.assert_allclose(film.conductance_W_K, expected_h * geom.grid.coolant_perimeter * geom.grid.dx)
    np.testing.assert_allclose(film.dp_shell_Pa, np.full(2, 12.0))
    assert calls[0][0]["rho_s"] == 2.0
    assert calls[1][0]["rho_s"] == 3.0
    assert calls[1][7] == 1.1


def test_shelltube_shell_film_zero_flow_uses_stagnant_conduction_scale(monkeypatch):
    from hps_combustor.transient_core import adapters_shelltube

    def fail_bell(*_args, **_kwargs):
        raise AssertionError("Bell-Delaware should not be called at zero shell flow")

    monkeypatch.setattr(adapters_shelltube, "bell_delaware_shell", fail_bell)

    geom = build_shelltube_core_geometry(
        tube_length=0.2,
        n_cells=2,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
    )
    film = shelltube_shell_film(
        geom,
        {"S_m": 0.002, "layout": "triangular30"},
        mdot_shell=0.0,
        rho=np.array([2.0, 3.0]),
        mu=np.array([1.0e-5, 2.0e-5]),
        k=np.array([0.12, 0.18]),
        cp=5000.0,
    )

    np.testing.assert_allclose(film.mass_flux_kg_m2s, [0.0, 0.0])
    np.testing.assert_allclose(film.reynolds, [0.0, 0.0])
    np.testing.assert_allclose(film.h_W_m2K, np.array([0.12, 0.18]) / geom.tube_outer_diameter)
    np.testing.assert_allclose(film.dp_shell_Pa, [0.0, 0.0])


def test_shelltube_tube_gas_film_uses_representative_tube_flow_and_smooth_dispatchers(monkeypatch):
    from hps_combustor.transient_core import adapters_shelltube

    friction_calls = []
    nusselt_calls = []

    def fake_friction(Re, roughness, Dh, x=10e10, Re_lo=2300.0, Re_hi=4000.0, error_factor=1.0):
        friction_calls.append((Re, roughness, Dh, x, Re_lo, Re_hi, error_factor))
        return 0.02

    def fake_nusselt(selector, Re, Pr, d, x, f_fd, T_bulk=None, T_wall=None,
                     error_factor=1.0, corrCoeffs=None):
        nusselt_calls.append((selector, Re, Pr, d, x, f_fd, T_bulk, T_wall, error_factor, corrCoeffs))
        return 80.0

    monkeypatch.setattr(adapters_shelltube, "dispatch_friction_tube_straight", fake_friction)
    monkeypatch.setattr(adapters_shelltube, "dispatch_nu_tube_straight", fake_nusselt)

    geom = build_shelltube_core_geometry(
        tube_length=0.2,
        n_cells=2,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
    )
    film = shelltube_tube_gas_film(
        geom,
        mdot_hot_total=0.1,
        rho=np.array([1.0, 2.0]),
        mu=1.0e-5,
        k=np.array([0.1, 0.2]),
        cp=1200.0,
        inside_tube_choice="smooth",
        nusselt_selector="gnielinski_blended",
        roughness=1.0e-6,
        x_for_developing=np.array([0.01, 0.02]),
        T_bulk=np.array([900.0, 800.0]),
        T_wall=np.array([500.0, 510.0]),
    )

    mdot_tube = 0.1 / 5.0
    expected_u = mdot_tube / (np.array([1.0, 2.0]) * geom.single_tube_inner_area)
    expected_re = np.array([1.0, 2.0]) * expected_u * geom.tube_inner_diameter / 1.0e-5
    expected_pr = 1200.0 * 1.0e-5 / np.array([0.1, 0.2])
    expected_h = 80.0 * np.array([0.1, 0.2]) / geom.tube_inner_diameter
    expected_dpdx = 0.02 * np.array([1.0, 2.0]) * expected_u**2 / (2.0 * geom.tube_inner_diameter)

    np.testing.assert_allclose(film.velocity_m_s, expected_u)
    np.testing.assert_allclose(film.reynolds, expected_re)
    np.testing.assert_allclose(film.prandtl, expected_pr)
    np.testing.assert_allclose(film.friction_factor, np.full(2, 0.02))
    np.testing.assert_allclose(film.nusselt, np.full(2, 80.0))
    np.testing.assert_allclose(film.h_W_m2K, expected_h)
    np.testing.assert_allclose(film.hot_conductance_W_K, expected_h * geom.grid.hot_perimeter * geom.grid.dx)
    np.testing.assert_allclose(film.dp_per_length_Pa_m, expected_dpdx)
    assert friction_calls[0][1] == 1.0e-6
    assert nusselt_calls[0][0] == "gnielinski_blended"
    assert nusselt_calls[1][6] == 800.0
    assert nusselt_calls[1][7] == 510.0


def test_shelltube_tube_gas_film_uses_grooved_path_and_calibration_factors(monkeypatch):
    from types import SimpleNamespace
    from hps_combustor.transient_core import adapters_shelltube

    friction_calls = []
    nusselt_calls = []

    def fake_grooved_friction(Re, phi, Re_lo=2000.0, Re_hi=4000.0, error_factor=1.0):
        friction_calls.append((Re, phi, Re_lo, Re_hi, error_factor))
        return 0.03

    def fake_grooved_nusselt(Re, Pr, phi, D_i, x, D_h=None,
                             Re_lo=2000.0, Re_hi=4000.0, error_factor=1.0):
        nusselt_calls.append((Re, Pr, phi, D_i, x, D_h, Re_lo, Re_hi, error_factor))
        return 90.0

    monkeypatch.setattr(adapters_shelltube, "friction_corrugated_tube_vicente", fake_grooved_friction)
    monkeypatch.setattr(adapters_shelltube, "nu_corrugated_tube_vicente", fake_grooved_nusselt)

    geom = build_shelltube_core_geometry(
        tube_length=0.1,
        n_cells=1,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
    )
    coeffs = SimpleNamespace(
        tube_grooved_Nu_factor=1.5,
        tube_grooved_f_factor=2.0,
        Re_transition_lo=2100.0,
        Re_transition_hi=3900.0,
    )
    film = shelltube_tube_gas_film(
        geom,
        mdot_hot_total=0.1,
        rho=1.0,
        mu=1.0e-5,
        k=0.1,
        cp=1200.0,
        inside_tube_choice="grooved",
        nusselt_selector="ignored_for_grooved",
        roughness=0.0,
        corrCoeffs=coeffs,
        corrugation_thickness=0.0002,
        corrugation_pitch=0.002,
    )

    expected_phi = 0.0002**2 / (0.002 * geom.tube_inner_diameter)
    np.testing.assert_allclose(film.friction_factor, [0.06])
    np.testing.assert_allclose(film.nusselt, [135.0])
    assert friction_calls[0][1] == expected_phi
    assert friction_calls[0][2] == 2100.0
    assert friction_calls[0][3] == 3900.0
    assert nusselt_calls[0][2] == expected_phi


def test_shelltube_wall_flux_converts_representative_tube_heat_to_total_cell_heat(monkeypatch):
    from hps_combustor.transient_core import adapters_shelltube

    calls = []

    class FakeConduction:
        def __init__(self, **kwargs):
            calls.append(kwargs)
            self.dx = kwargs["dx"]

        def fluxes_at_Tbar(self, T_bar, h_g_rad=None):
            return {
                "dq_hot__dx": 20.0 * self.dx,
                "dq_cold__dx": 12.0 * self.dx,
                "T_wg": T_bar + 6.0,
                "T_wc": T_bar - 3.0,
                "h_g_rad": h_g_rad,
                "q_w_rad": 0.0,
                "k_w": 16.0,
            }

    monkeypatch.setattr(
        adapters_shelltube,
        "OneDimensionalSteadyConduction_ShellnHelicalTube",
        FakeConduction,
    )

    geom = build_shelltube_core_geometry(
        tube_length=0.3,
        n_cells=3,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
    )
    result = shelltube_wall_flux(
        geom,
        Tbar_wall=np.array([300.0, 310.0, 320.0]),
        T_coolant=np.array([100.0, 110.0, 120.0]),
        T_gas=np.array([900.0, 850.0, 800.0]),
        h_gas=np.array([50.0, 60.0, 70.0]),
        h_shell=np.array([1000.0, 1200.0, 1400.0]),
        wall_conductivity_at_T=lambda T: 16.0,
        h_g_rad=np.array([1.0, 2.0, 3.0]),
    )

    # Fake dq__dx is 20*dx and 12*dx for one representative tube. The adapter
    # converts to total W per cell by multiplying by N_tubes*dx.
    np.testing.assert_allclose(result.dq_hot_per_length_W_m, np.full(3, 2.0))
    np.testing.assert_allclose(result.dq_cold_per_length_W_m, np.full(3, 1.2))
    np.testing.assert_allclose(result.hot_heat_W, np.full(3, 2.0 * 5.0 * 0.1))
    np.testing.assert_allclose(result.cold_heat_W, np.full(3, 1.2 * 5.0 * 0.1))
    np.testing.assert_allclose(result.T_wg, [306.0, 316.0, 326.0])
    np.testing.assert_allclose(result.T_wc, [297.0, 307.0, 317.0])
    np.testing.assert_allclose(result.h_g_rad, [1.0, 2.0, 3.0])
    np.testing.assert_allclose(result.q_w_rad, [0.0, 0.0, 0.0])
    np.testing.assert_allclose(result.k_wall, [16.0, 16.0, 16.0])
    assert calls[0]["hot_side"] == "inner"
    assert calls[0]["rad_enabled"] is False


def test_shelltube_hot_gas_march_advances_representative_tube_enthalpy(monkeypatch):
    from hps_combustor.transient_core import adapters_shelltube

    film_calls = []
    conduction_calls = []
    state_calls = []

    class FakeFilm:
        velocity_m_s = np.array([2.0])
        reynolds = np.array([1000.0])
        prandtl = np.array([0.7])
        friction_factor = np.array([0.02])
        nusselt = np.array([50.0])
        h_W_m2K = np.array([100.0])
        hot_conductance_W_K = np.array([1.0])
        dp_per_length_Pa_m = np.array([5.0])

    def fake_film(geometry, **kwargs):
        film_calls.append((geometry, kwargs))
        return FakeFilm()

    class FakeConduction:
        def __init__(self, **kwargs):
            conduction_calls.append(kwargs)
            self.dx = kwargs["dx"]

        def fluxes_at_Tbar(self, T_bar, h_g_rad=None):
            return {
                "dq_hot__dx": 10.0,
                "dq_cold__dx": 8.0,
                "T_wg": T_bar + 4.0,
                "T_wc": T_bar - 2.0,
                "q_w_rad": 0.0,
                "k_w": 15.0,
            }

    def gas_state_at(h_removed, progress, i):
        state_calls.append((h_removed, progress, i))
        return {
            "T": 900.0 - 0.01 * h_removed,
            "rho": 1.0,
            "mu": 1.0e-5,
            "k": 0.1,
            "cp": 1200.0,
            "progress_source": 0.2,
        }

    monkeypatch.setattr(adapters_shelltube, "shelltube_tube_gas_film", fake_film)
    monkeypatch.setattr(
        adapters_shelltube,
        "OneDimensionalSteadyConduction_ShellnHelicalTube",
        FakeConduction,
    )

    geom = build_shelltube_core_geometry(
        tube_length=0.3,
        n_cells=3,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
    )
    result = shelltube_hot_gas_march(
        geom,
        Tbar_wall=np.array([300.0, 310.0, 320.0]),
        T_coolant=np.array([100.0, 110.0, 120.0]),
        h_shell=np.array([1000.0, 1100.0, 1200.0]),
        mdot_hot_total=0.5,
        gas_state_at=gas_state_at,
        wall_conductivity_at_T=lambda T: 15.0,
        inside_tube_choice="smooth",
        nusselt_selector="gnielinski_blended",
        roughness=1.0e-6,
        h_g_rad=np.array([0.0, 1.0, 2.0]),
    )

    # mdot_tube = 0.5/5 = 0.1 kg/s. Each 0.1 m cell removes
    # dq_hot__dx*dx/mdot_tube = 10*0.1/0.1 = 10 J/kg.
    np.testing.assert_allclose(result.enthalpy_removed_J_kg, [0.0, 10.0, 20.0])
    np.testing.assert_allclose(result.enthalpy_removed_outlet_J_kg, 30.0)
    np.testing.assert_allclose(result.progress_variable, [0.0, 0.01, 0.02])
    np.testing.assert_allclose(result.progress_outlet, 0.03)
    np.testing.assert_allclose(result.T_gas, [900.0, 899.9, 899.8])
    np.testing.assert_allclose(result.h_gas_W_m2K, [100.0, 100.0, 100.0])
    np.testing.assert_allclose(result.wall_flux.hot_heat_W, np.full(3, 10.0 * 5.0 * 0.1))
    np.testing.assert_allclose(result.wall_flux.cold_heat_W, np.full(3, 8.0 * 5.0 * 0.1))
    np.testing.assert_allclose(result.wall_flux.h_g_rad, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(result.T_gas_outlet, 899.7)
    assert conduction_calls[0]["hot_side"] == "inner"
    assert film_calls[0][1]["mdot_hot_total"] == 0.5
    np.testing.assert_allclose(state_calls[-1][:2], (30.0, 0.03))
    assert state_calls[-1][2] == 3


def test_shelltube_fpv_gas_state_provider_wraps_fpv_manifold():
    class FakeFPV:
        def Yc_inlet(self):
            return 0.4

        def state(self, h_removed, progress):
            return (
                900.0 - 0.01 * h_removed,
                1.0,
                1.0e-5,
                0.1,
                1200.0,
                0.2,
                0.1,
                0.03 + progress,
            )

    provider, initial = fpv_gas_state_provider(FakeFPV())
    state = provider(20.0, 0.5, 3)

    assert initial == 0.4
    assert state.T == 899.8
    assert state.rho == 1.0
    assert state.mu == 1.0e-5
    assert state.k == 0.1
    assert state.cp == 1200.0
    assert state.progress_source == 0.53


def test_shelltube_equilibrium_gas_state_provider_wraps_equilibrium_manifold():
    class FakeEquilibrium:
        def at(self, h_removed):
            return (850.0 - h_removed, 1.1, 2.0e-5, 0.2, 1300.0, 0.3, 0.1)

    provider, initial = equilibrium_gas_state_provider(FakeEquilibrium())
    state = provider(5.0, 123.0, 0)

    assert initial == 0.0
    assert state.T == 845.0
    assert state.rho == 1.1
    assert state.mu == 2.0e-5
    assert state.k == 0.2
    assert state.cp == 1300.0
    assert state.progress_source == 0.0


def test_shelltube_oxygen_gas_state_provider_uses_enthalpy_removed(monkeypatch):
    # oxygen_gas_state_provider was relocated to core.hotgas.combustor in
    # Stage D Slice 3 (2026-08-19) -- it calls PropsSI from ITS OWN module
    # namespace now, not adapters_shelltube's (which only re-exports the
    # function object itself, not the name lookups inside it).
    from hps_combustor.core.hotgas import combustor as hotgas_combustor

    calls = []

    def fake_props(output, name1, value1, name2, value2, fluid):
        calls.append((output, name1, value1, name2, value2, fluid))
        if output == "H":
            return 1000.0
        if output == "T":
            assert name1 == "H"
            return 100.0 + 0.01 * value1
        if output == "D":
            return 2.0
        if output == "V":
            return 3.0e-5
        if output == "L":
            return 0.12
        if output == "C":
            return 920.0
        raise AssertionError(output)

    monkeypatch.setattr(hotgas_combustor, "PropsSI", fake_props)

    provider, initial = oxygen_gas_state_provider(
        T_inlet=110.0,
        pressure=8.0e6,
        fluid="Oxygen",
        T_min=95.0,
        T_max=1200.0,
    )
    state = provider(100.0, 0.0, 0)

    assert initial == 0.0
    # h = 1000 - 100, so fake T = 109 K.
    assert state.T == 109.0
    assert state.rho == 2.0
    assert state.mu == 3.0e-5
    assert state.k == 0.12
    assert state.cp == 920.0
    assert state.progress_source == 0.0
    assert calls[0] == ("H", "T", 110.0, "P", 8.0e6, "Oxygen")
    assert calls[1] == ("T", "H", 900.0, "P", 8.0e6, "Oxygen")


def test_shelltube_step_inputs_assembles_wall_coolant_inputs(monkeypatch):
    from hps_combustor.transient_core import adapters_shelltube

    geom = build_shelltube_core_geometry(
        tube_length=0.2,
        n_cells=2,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
        flow_config="counter",
    )

    shell_calls = []
    hot_calls = []

    def fake_shell_film(geometry, bell_geometry, **kwargs):
        shell_calls.append((geometry, bell_geometry, kwargs))
        return ShellTubeShellFilm(
            mass_flux_kg_m2s=np.array([1.0, 1.0]),
            reynolds=np.array([100.0, 100.0]),
            prandtl=np.array([0.7, 0.7]),
            h_W_m2K=np.array([200.0, 300.0]),
            conductance_W_K=np.array([2.0, 3.0]),
            dp_shell_Pa=np.array([5.0, 6.0]),
        )

    def fake_hot_march(geometry, **kwargs):
        hot_calls.append((geometry, kwargs))
        return ShellTubeHotGasMarch(
            T_gas=np.array([900.0, 850.0]),
            h_gas_W_m2K=np.array([100.0, 110.0]),
            gas_velocity_m_s=np.array([2.0, 2.0]),
            reynolds=np.array([1000.0, 1000.0]),
            prandtl=np.array([0.7, 0.7]),
            friction_factor=np.array([0.02, 0.02]),
            nusselt=np.array([50.0, 50.0]),
            dp_per_length_Pa_m=np.array([5.0, 5.0]),
            wall_flux=ShellTubeWallFlux(
                hot_heat_W=np.array([10.0, 20.0]),
                cold_heat_W=np.array([8.0, 16.0]),
                dq_hot_per_length_W_m=np.array([20.0, 40.0]),
                dq_cold_per_length_W_m=np.array([16.0, 32.0]),
                T_wg=np.array([310.0, 320.0]),
                T_wc=np.array([290.0, 300.0]),
                h_g_rad=np.zeros(2),
                q_w_rad=np.zeros(2),
                k_wall=np.array([15.0, 15.0]),
            ),
            enthalpy_removed_J_kg=np.array([0.0, 10.0]),
            progress_variable=np.array([0.1, 0.2]),
            T_gas_outlet=830.0,
            enthalpy_removed_outlet_J_kg=20.0,
            progress_outlet=0.3,
        )

    def coolant_props(T, pressure):
        np.testing.assert_allclose(T, [100.0, 110.0])
        assert pressure == 8.0e6
        return ShellTubeFluidProperties(
            rho=np.array([2.0, 3.0]),
            mu=np.array([1.0e-5, 2.0e-5]),
            k=np.array([0.1, 0.2]),
            cp=np.array([5000.0, 5200.0]),
        )

    monkeypatch.setattr(adapters_shelltube, "shelltube_shell_film", fake_shell_film)
    monkeypatch.setattr(adapters_shelltube, "shelltube_hot_gas_march", fake_hot_march)

    assembled = shelltube_step_inputs(
        geom,
        {"S_m": 0.002},
        Tbar_wall=np.array([300.0, 310.0]),
        T_coolant=np.array([100.0, 110.0]),
        mdot_coolant=0.15,
        T_coolant_inlet=90.0,
        p_coolant=8.0e6,
        mdot_hot_total=0.1,
        gas_state_at=lambda h, y, i: None,
        coolant_properties_at=coolant_props,
        wall_density=8000.0,
        wall_cp=500.0,
        wall_conductivity_at_T=lambda T: 15.0,
        inside_tube_choice="smooth",
        nusselt_selector="gnielinski_blended",
        tube_roughness=1.0e-6,
        progress_initial=0.4,
    )

    inputs = assembled.wall_coolant_inputs
    np.testing.assert_allclose(inputs.wall_heat_capacity, 8000.0 * 500.0 * geom.grid.wall_volume)
    np.testing.assert_allclose(
        inputs.coolant_heat_capacity,
        np.array([2.0, 3.0]) * np.array([5000.0, 5200.0]) * geom.grid.coolant_volume,
    )
    np.testing.assert_allclose(inputs.coolant_cp, [5000.0, 5200.0])
    assert inputs.mdot_coolant == 0.15
    assert inputs.T_coolant_inlet == 90.0
    np.testing.assert_allclose(inputs.hot_heat_W, [10.0, 20.0])
    np.testing.assert_allclose(inputs.wall_to_coolant_conductance_W_per_K, [2.0, 3.0])
    assert inputs.flow_direction == -1
    assert assembled.hot_gas_march.T_gas_outlet == 830.0
    assert shell_calls[0][2]["mdot_shell"] == 0.15
    assert hot_calls[0][1]["mdot_hot_total"] == 0.1
    assert hot_calls[0][1]["progress_initial"] == 0.4


def test_run_shelltube_transient_core_uses_schedules_and_integrates(monkeypatch):
    from hps_combustor.transient_core import adapters_shelltube

    geom = build_shelltube_core_geometry(
        tube_length=0.2,
        n_cells=2,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
        flow_config="co",
    )
    calls = []

    def fake_step_inputs(geometry, bell_geometry, **kwargs):
        calls.append(kwargs)
        hot = np.array([10.0, 20.0])
        shell_film = ShellTubeShellFilm(
            mass_flux_kg_m2s=np.ones(2),
            reynolds=np.ones(2),
            prandtl=np.ones(2),
            h_W_m2K=np.array([100.0, 100.0]),
            conductance_W_K=np.array([2.0, 2.0]),
            dp_shell_Pa=np.zeros(2),
        )
        march = ShellTubeHotGasMarch(
            T_gas=np.array([900.0, 850.0]),
            h_gas_W_m2K=np.array([100.0, 100.0]),
            gas_velocity_m_s=np.array([2.0, 2.0]),
            reynolds=np.ones(2),
            prandtl=np.ones(2),
            friction_factor=np.ones(2),
            nusselt=np.ones(2),
            dp_per_length_Pa_m=np.zeros(2),
            wall_flux=ShellTubeWallFlux(
                hot_heat_W=hot,
                cold_heat_W=np.array([8.0, 16.0]),
                dq_hot_per_length_W_m=np.array([20.0, 40.0]),
                dq_cold_per_length_W_m=np.array([16.0, 32.0]),
                T_wg=np.array([310.0, 320.0]),
                T_wc=np.array([290.0, 300.0]),
                h_g_rad=np.zeros(2),
                q_w_rad=np.zeros(2),
                k_wall=np.ones(2),
            ),
            enthalpy_removed_J_kg=np.array([0.0, 10.0]),
            progress_variable=np.zeros(2),
            T_gas_outlet=830.0,
            enthalpy_removed_outlet_J_kg=20.0,
            progress_outlet=0.0,
        )
        inputs = WallCoolantStepInputs(
            wall_heat_capacity=np.array([100.0, 100.0]),
            coolant_heat_capacity=np.array([50.0, 50.0]),
            coolant_cp=np.array([1000.0, 1000.0]),
            mdot_coolant=kwargs["mdot_coolant"],
            T_coolant_inlet=kwargs["T_coolant_inlet"],
            hot_heat_W=hot,
            wall_to_coolant_conductance_W_per_K=np.array([2.0, 2.0]),
            flow_direction=kwargs["flow_direction"] or geometry.grid.flow_direction,
        )
        return adapters_shelltube.ShellTubeStepInputDiagnostics(
            wall_coolant_inputs=inputs,
            hot_gas_march=march,
            shell_film=shell_film,
            coolant_properties=ShellTubeFluidProperties(
                rho=np.ones(2),
                mu=np.ones(2),
                k=np.ones(2),
                cp=np.array([1000.0, 1000.0]),
            ),
            wall_heat_capacity_J_K=np.array([100.0, 100.0]),
        )

    monkeypatch.setattr(adapters_shelltube, "shelltube_step_inputs", fake_step_inputs)
    provider_calls = []

    def gas_provider_at_time(t):
        provider_calls.append(t)
        return (lambda h, y, i: None), 0.5 + t

    result = run_shelltube_transient_core(
        geom,
        {"S_m": 0.002},
        T_wall_initial=np.array([300.0, 300.0]),
        T_coolant_initial=np.array([100.0, 100.0]),
        t_end=0.3,
        max_step=0.2,
        gas_provider_at_time=gas_provider_at_time,
        coolant_properties_at=lambda T, p: None,
        wall_density=8000.0,
        wall_cp=500.0,
        wall_conductivity_at_T=lambda T: 15.0,
        inside_tube_choice="smooth",
        nusselt_selector="gnielinski_blended",
        tube_roughness=1.0e-6,
        mdot_coolant_default=0.1,
        T_coolant_inlet_default=90.0,
        p_coolant_default=8.0e6,
        mdot_hot_total_default=0.1,
        mdot_coolant_schedule=((0.1, 0.0), (0.3, 0.2)),
        T_coolant_inlet_schedule=((0.2, 95.0),),
    )

    np.testing.assert_allclose(result.integration.t, [0.0, 0.1, 0.2, 0.3])
    assert len(result.step_diagnostics) == 3
    assert len(calls) == 3
    assert calls[0]["mdot_coolant"] == 0.0
    assert calls[1]["mdot_coolant"] == 0.0
    np.testing.assert_allclose(calls[2]["mdot_coolant"], 0.1)
    np.testing.assert_allclose(provider_calls, [0.0, 0.1, 0.2])
    np.testing.assert_allclose(
        [call["progress_initial"] for call in calls],
        [0.5, 0.6, 0.7],
    )
    assert calls[0]["T_coolant_inlet"] == 95.0
    assert result.integration.T_wall[-1, 0] < 300.0
    assert result.integration.T_coolant[-1, 0] > 100.0


def test_shelltube_mass_energy_mode_keeps_residual_outflow_after_inlet_closure(monkeypatch):
    from hps_combustor.transient_core import adapters_shelltube

    geom = build_shelltube_core_geometry(
        tube_length=0.2,
        n_cells=2,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
        flow_config="co",
    )

    monkeypatch.setattr(
        adapters_shelltube,
        "_shelltube_nominal_pressure_drop",
        lambda *args, **kwargs: 1000.0,
    )

    def fake_step_inputs(geometry, bell_geometry, **kwargs):
        hot = np.array([0.0, 0.0])
        shell_film = ShellTubeShellFilm(
            mass_flux_kg_m2s=np.ones(2),
            reynolds=np.ones(2),
            prandtl=np.ones(2),
            h_W_m2K=np.array([100.0, 100.0]),
            conductance_W_K=np.array([2.0, 2.0]),
            dp_shell_Pa=np.array([1000.0, 1000.0]),
        )
        march = ShellTubeHotGasMarch(
            T_gas=np.array([900.0, 850.0]),
            h_gas_W_m2K=np.array([100.0, 100.0]),
            gas_velocity_m_s=np.array([2.0, 2.0]),
            reynolds=np.ones(2),
            prandtl=np.ones(2),
            friction_factor=np.ones(2),
            nusselt=np.ones(2),
            dp_per_length_Pa_m=np.zeros(2),
            wall_flux=ShellTubeWallFlux(
                hot_heat_W=hot,
                cold_heat_W=np.zeros(2),
                dq_hot_per_length_W_m=np.zeros(2),
                dq_cold_per_length_W_m=np.zeros(2),
                T_wg=np.array([300.0, 300.0]),
                T_wc=np.array([100.0, 100.0]),
                h_g_rad=np.zeros(2),
                q_w_rad=np.zeros(2),
                k_wall=np.ones(2),
            ),
            enthalpy_removed_J_kg=np.zeros(2),
            progress_variable=np.zeros(2),
            T_gas_outlet=850.0,
            enthalpy_removed_outlet_J_kg=0.0,
            progress_outlet=0.0,
        )
        inputs = WallCoolantStepInputs(
            wall_heat_capacity=np.array([100.0, 100.0]),
            coolant_heat_capacity=np.array([50.0, 50.0]),
            coolant_cp=np.array([5200.0, 5200.0]),
            mdot_coolant=kwargs["mdot_coolant"],
            T_coolant_inlet=kwargs["T_coolant_inlet"],
            hot_heat_W=hot,
            wall_to_coolant_conductance_W_per_K=np.array([2.0, 2.0]),
            flow_direction=kwargs["flow_direction"] or geometry.grid.flow_direction,
        )
        props = kwargs["coolant_properties_at"](np.array([90.0, 90.0]), kwargs["p_coolant"])
        return adapters_shelltube.ShellTubeStepInputDiagnostics(
            wall_coolant_inputs=inputs,
            hot_gas_march=march,
            shell_film=shell_film,
            coolant_properties=props,
            wall_heat_capacity_J_K=np.array([100.0, 100.0]),
        )

    monkeypatch.setattr(adapters_shelltube, "shelltube_step_inputs", fake_step_inputs)

    result = run_shelltube_transient_core(
        geom,
        {"S_m": 0.002},
        T_wall_initial=np.array([300.0, 300.0]),
        T_coolant_initial=np.array([90.0, 90.0]),
        t_end=0.002,
        max_step=0.001,
        gas_state_at=lambda h, y, i: None,
        coolant_properties_at=lambda T, p: None,
        wall_density=8000.0,
        wall_cp=500.0,
        wall_conductivity_at_T=lambda T: 15.0,
        inside_tube_choice="smooth",
        nusselt_selector="gnielinski_blended",
        tube_roughness=1.0e-6,
        mdot_coolant_default=1.0e-4,
        T_coolant_inlet_default=90.0,
        p_coolant_default=8.0e6,
        mdot_hot_total_default=0.1,
        mdot_coolant_schedule=((0.0, 1.0e-4), (0.001, 0.0)),
        coolant_state_model="mass_energy",
    )

    assert hasattr(result.integration, "coolant_mass_kg")
    np.testing.assert_allclose(result.integration.face_mdot_kg_s[1, 0], 1.0e-4)
    np.testing.assert_allclose(result.integration.face_mdot_kg_s[2, 0], 0.0)
    np.testing.assert_allclose(result.integration.face_mdot_kg_s[2, 1], 0.0)
    assert result.integration.face_mdot_kg_s[2, -1] > 0.0
    assert result.integration.coolant_mass_kg[2].sum() < result.integration.coolant_mass_kg[1].sum()


def test_shelltube_closed_valve_memory_outflow_requires_pressure_drive(monkeypatch):
    from hps_combustor.transient_core import adapters_shelltube

    geom = build_shelltube_core_geometry(
        tube_length=0.2,
        n_cells=2,
        shell_inner_diameter=0.08,
        tube_outer_diameter=0.006,
        wall_thickness=0.001,
        n_tubes=5,
        flow_config="co",
    )

    monkeypatch.setattr(
        adapters_shelltube,
        "_shelltube_nominal_pressure_drop",
        lambda *args, **kwargs: 1000.0,
    )

    def fake_step_inputs(geometry, bell_geometry, **kwargs):
        shell_film = ShellTubeShellFilm(
            mass_flux_kg_m2s=np.ones(2),
            reynolds=np.ones(2),
            prandtl=np.ones(2),
            h_W_m2K=np.array([1.0, 1.0]),
            conductance_W_K=np.array([0.0, 0.0]),
            dp_shell_Pa=np.zeros(2),
        )
        march = ShellTubeHotGasMarch(
            T_gas=np.full(2, np.nan),
            h_gas_W_m2K=np.zeros(2),
            gas_velocity_m_s=np.zeros(2),
            reynolds=np.zeros(2),
            prandtl=np.zeros(2),
            friction_factor=np.zeros(2),
            nusselt=np.zeros(2),
            dp_per_length_Pa_m=np.zeros(2),
            wall_flux=ShellTubeWallFlux(
                hot_heat_W=np.zeros(2),
                cold_heat_W=np.zeros(2),
                dq_hot_per_length_W_m=np.zeros(2),
                dq_cold_per_length_W_m=np.zeros(2),
                T_wg=np.full(2, 300.0),
                T_wc=np.full(2, 100.0),
                h_g_rad=np.zeros(2),
                q_w_rad=np.zeros(2),
                k_wall=np.ones(2),
            ),
            enthalpy_removed_J_kg=np.zeros(2),
            progress_variable=np.zeros(2),
            T_gas_outlet=np.nan,
            enthalpy_removed_outlet_J_kg=0.0,
            progress_outlet=0.0,
        )
        inputs = WallCoolantStepInputs(
            wall_heat_capacity=np.array([100.0, 100.0]),
            coolant_heat_capacity=np.array([50.0, 50.0]),
            coolant_cp=np.array([5200.0, 5200.0]),
            mdot_coolant=kwargs["mdot_coolant"],
            T_coolant_inlet=kwargs["T_coolant_inlet"],
            hot_heat_W=np.zeros(2),
            wall_to_coolant_conductance_W_per_K=np.zeros(2),
            flow_direction=kwargs["flow_direction"] or geometry.grid.flow_direction,
        )
        props = kwargs["coolant_properties_at"](np.array([90.0, 90.0]), kwargs["p_coolant"])
        return adapters_shelltube.ShellTubeStepInputDiagnostics(
            wall_coolant_inputs=inputs,
            hot_gas_march=march,
            shell_film=shell_film,
            coolant_properties=props,
            wall_heat_capacity_J_K=np.array([100.0, 100.0]),
        )

    monkeypatch.setattr(adapters_shelltube, "shelltube_step_inputs", fake_step_inputs)

    class ProgressConfig:
        schedule_p_c_out = ((0.0, 1.0e8), (0.02, 1.0e8))
        insert_schedule_breakpoints = True
        progress_print = False
        transient_coolant_outlet_pressure = None

    result = run_shelltube_transient_core(
        geom,
        {"S_m": 0.002},
        T_wall_initial=np.array([300.0, 300.0]),
        T_coolant_initial=np.array([90.0, 90.0]),
        t_end=0.02,
        max_step=0.01,
        gas_state_at=lambda h, y, i: None,
        coolant_properties_at=lambda T, p: None,
        wall_density=8000.0,
        wall_cp=500.0,
        wall_conductivity_at_T=lambda T: 15.0,
        inside_tube_choice="smooth",
        nusselt_selector="gnielinski_blended",
        tube_roughness=1.0e-6,
        mdot_coolant_default=1.0e-4,
        T_coolant_inlet_default=90.0,
        p_coolant_default=8.0e6,
        mdot_hot_total_default=0.0,
        mdot_coolant_schedule=((0.0, 1.0e-4), (0.01, 0.0)),
        coolant_state_model="mass_energy",
        progress_config=ProgressConfig(),
    )

    np.testing.assert_allclose(result.integration.face_mdot_kg_s[2], 0.0, atol=1.0e-15)


def test_shelltube_low_mach_momentum_update_satisfies_implicit_balance():
    from hps_combustor.transient_core import adapters_shelltube

    mdot = adapters_shelltube._implicit_quadratic_momentum_update(
        old_mdot=0.02,
        pressure_drive=5000.0,
        inertance=30.0,
        resistance_over_density=2.0e5,
        dt=0.01,
    )

    residual = 30.0 * (mdot - 0.02) / 0.01 - 5000.0 + 2.0e5 * mdot * abs(mdot)
    assert abs(residual) < 1.0e-9


def test_shelltube_low_mach_faces_apply_counterflow_boundaries():
    from hps_combustor.transient_core import adapters_shelltube

    faces = adapters_shelltube._shelltube_low_mach_momentum_faces(
        face_old=np.zeros(4),
        pressure=np.array([7.1e6, 7.2e6, 7.3e6]),
        density=np.array([10.0, 10.0, 10.0]),
        resistance=np.array([1.0e7, 1.0e7, 1.0e7]),
        inertance=np.array([20.0, 40.0, 40.0, 20.0]),
        dt=0.01,
        inlet_pressure=7.5e6,
        outlet_pressure=7.0e6,
        flow_direction=-1,
    )

    assert faces[0] < 0.0
    assert faces[-1] < 0.0

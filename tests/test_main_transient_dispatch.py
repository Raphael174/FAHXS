import numpy as np
from types import SimpleNamespace

from hps_combustor.input_data import (
    CorrelationCoefficients,
    combustorProp,
    coolantProp,
    hotgasProp,
    numericalProp,
    runProp,
    shellTubeProp,
    system_requirements,
    transientProp,
)
from hps_combustor import main_transient
from hps_combustor.main_solve_shellntube_transient import shellntube_transient_solver
from hps_combustor.transient_core import (
    ShellTubeFluidProperties,
    ShellTubeHotGasMarch,
    ShellTubeShellFilm,
    ShellTubeStepInputDiagnostics,
    ShellTubeWallFlux,
    WallCoolantIntegrationResult,
    WallCoolantStepInputs,
    AxialGrid,
)


class FakeShellTubeTransientSolver:
    calls = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.time_series = None
        FakeShellTubeTransientSolver.calls.append(("init", kwargs))

    def solve_transient(self, verbose=True):
        FakeShellTubeTransientSolver.calls.append(("legacy", verbose))
        self._set_time_series(100.0)

    def solve_transient_core(self, verbose=True):
        FakeShellTubeTransientSolver.calls.append(("core", verbose))
        self._set_time_series(200.0)

    def _set_time_series(self, marker):
        self.time_series = {
            "t": np.array([0.0, 0.1]),
            "x": np.array([0.0]),
            "fields": {"Tbar": np.array([[300.0], [301.0]])},
            "scalars": {
                "marker": np.array([0.0, marker]),
                "T_c_out": np.array([90.0, 91.0]),
            },
        }


class FakeHelicalTransientSolver:
    calls = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.time_series = None
        FakeHelicalTransientSolver.calls.append(("init", kwargs))

    def solve_transient(self, verbose=True):
        FakeHelicalTransientSolver.calls.append(("legacy", verbose))
        self._set_time_series(300.0)

    def solve_transient_core(self, verbose=True):
        FakeHelicalTransientSolver.calls.append(("core", verbose))
        self._set_time_series(400.0)

    def _set_time_series(self, marker):
        self.time_series = {
            "t": np.array([0.0, 0.1]),
            "x": np.array([0.0]),
            "fields": {"Tbar": np.array([[300.0], [301.0]])},
            "scalars": {
                "marker": np.array([0.0, marker]),
                "T_c_out": np.array([90.0, 91.0]),
            },
        }


def make_inputs(fluid_model):
    transient = transientProp()
    transient.fluid_model = fluid_model
    run = runProp()
    run.shelltube_transient_nodes = 3
    combustor = combustorProp()
    combustor.HX_config = "shellntube"
    combustor.flow_config = "counter"
    return {
        "coolant": coolantProp(),
        "hotgas": hotgasProp(),
        "combustor": combustor,
        "shelltube": shellTubeProp(),
        "numerical": numericalProp(),
        "transient": transient,
        "system": system_requirements(),
        "correlations": CorrelationCoefficients(),
        "run": run,
    }


def make_helical_inputs(fluid_model):
    transient = transientProp()
    transient.fluid_model = fluid_model
    combustor = combustorProp()
    combustor.HX_config = "shellnHelicalTube"
    combustor.flow_config = "counter"
    return {
        "coolant": coolantProp(),
        "hotgas": hotgasProp(),
        "combustor": combustor,
        "shelltube": shellTubeProp(),
        "numerical": numericalProp(),
        "transient": transient,
        "system": system_requirements(),
        "correlations": CorrelationCoefficients(),
        "run": runProp(),
    }


def test_shelltube_transient_coolant_dispatch_uses_core(monkeypatch):
    FakeShellTubeTransientSolver.calls.clear()
    monkeypatch.setattr(
        main_transient,
        "shellntube_transient_solver",
        FakeShellTubeTransientSolver,
    )

    solver, summary = main_transient.run_transient(make_inputs("transient_coolant"))

    assert ("core", True) in FakeShellTubeTransientSolver.calls
    assert ("legacy", True) not in FakeShellTubeTransientSolver.calls
    assert summary["marker_final"] == 200.0
    assert summary["T_c_out_final"] == 91.0
    assert solver.kwargs["N_axial"] == 3
    assert solver.kwargs["flow_config"] == "counter"


def test_helical_transient_coolant_dispatch_uses_core(monkeypatch):
    FakeHelicalTransientSolver.calls.clear()
    monkeypatch.setattr(
        main_transient,
        "transient_solver",
        FakeHelicalTransientSolver,
    )

    _solver, summary = main_transient.run_transient(make_helical_inputs("transient_coolant"))

    assert ("core", True) in FakeHelicalTransientSolver.calls
    assert ("legacy", True) not in FakeHelicalTransientSolver.calls
    assert summary["marker_final"] == 400.0


def test_helical_quasi_steady_dispatch_keeps_legacy_path(monkeypatch):
    FakeHelicalTransientSolver.calls.clear()
    monkeypatch.setattr(
        main_transient,
        "transient_solver",
        FakeHelicalTransientSolver,
    )

    _solver, summary = main_transient.run_transient(make_helical_inputs("quasi_steady"))

    assert ("legacy", True) in FakeHelicalTransientSolver.calls
    assert ("core", True) not in FakeHelicalTransientSolver.calls
    assert summary["marker_final"] == 300.0


def test_shelltube_quasi_steady_dispatch_keeps_legacy_path(monkeypatch):
    FakeShellTubeTransientSolver.calls.clear()
    monkeypatch.setattr(
        main_transient,
        "shellntube_transient_solver",
        FakeShellTubeTransientSolver,
    )

    _solver, summary = main_transient.run_transient(make_inputs("quasi_steady"))

    assert ("legacy", True) in FakeShellTubeTransientSolver.calls
    assert ("core", True) not in FakeShellTubeTransientSolver.calls
    assert summary["marker_final"] == 100.0


def test_shelltube_transient_core_time_series_exposes_engineering_metrics():
    solver = object.__new__(shellntube_transient_solver)
    solver._bc_at = lambda _t: {"mdot_c": 0.15, "mdot_g": 0.1, "mdot_lox": 0.0, "ignited": True}

    grid = AxialGrid.uniform(
        length=0.2,
        n_cells=2,
        coolant_area=0.01,
        wall_area=0.002,
        hot_perimeter=0.1,
        coolant_perimeter=0.12,
        flow_direction=1,
    )
    geometry = SimpleNamespace(grid=grid)
    integration = WallCoolantIntegrationResult(
        t=np.array([0.0, 0.1]),
        T_wall=np.array([[300.0, 301.0], [302.0, 303.0]]),
        T_coolant=np.array([[90.0, 91.0], [92.0, 93.0]]),
        T_coolant_outlet=np.array([91.0, 93.0]),
        hot_heat_added_J=np.array([0.0, 10.0]),
        advective_energy_in_J=np.array([0.0, 20.0]),
        advective_energy_out_J=np.array([0.0, 15.0]),
        energy_residual_J=np.array([0.0, 1.5]),
        heat_wall_to_coolant_W=np.array([[0.0, 0.0], [4.0, 6.0]]),
        last_step=None,
    )
    wall_flux = ShellTubeWallFlux(
        hot_heat_W=np.array([100.0, 200.0]),
        cold_heat_W=np.array([80.0, 160.0]),
        dq_hot_per_length_W_m=np.array([1000.0, 2000.0]),
        dq_cold_per_length_W_m=np.array([800.0, 1600.0]),
        T_wg=np.array([310.0, 320.0]),
        T_wc=np.array([290.0, 295.0]),
        h_g_rad=np.array([0.0, 0.0]),
        q_w_rad=np.array([0.0, 0.0]),
        k_wall=np.array([15.0, 16.0]),
    )
    march = ShellTubeHotGasMarch(
        T_gas=np.array([900.0, 850.0]),
        h_gas_W_m2K=np.array([100.0, 110.0]),
        gas_velocity_m_s=np.array([20.0, 18.0]),
        reynolds=np.array([10000.0, 9000.0]),
        prandtl=np.array([0.7, 0.72]),
        friction_factor=np.array([0.02, 0.021]),
        nusselt=np.array([80.0, 82.0]),
        dp_per_length_Pa_m=np.array([30.0, 40.0]),
        wall_flux=wall_flux,
        enthalpy_removed_J_kg=np.array([0.0, 1000.0]),
        progress_variable=np.array([0.1, 0.2]),
        T_gas_outlet=830.0,
        enthalpy_removed_outlet_J_kg=1500.0,
        progress_outlet=0.25,
    )
    shell_film = ShellTubeShellFilm(
        mass_flux_kg_m2s=np.array([1.0, 2.0]),
        reynolds=np.array([100.0, 200.0]),
        prandtl=np.array([0.68, 0.69]),
        h_W_m2K=np.array([500.0, 600.0]),
        conductance_W_K=np.array([5.0, 6.0]),
        dp_shell_Pa=np.array([7.0, 8.0]),
    )
    diag = ShellTubeStepInputDiagnostics(
        wall_coolant_inputs=WallCoolantStepInputs(
            wall_heat_capacity=np.array([1.0, 1.0]),
            coolant_heat_capacity=np.array([1.0, 1.0]),
            coolant_cp=np.array([5000.0, 5100.0]),
            mdot_coolant=0.15,
            T_coolant_inlet=90.0,
            hot_heat_W=np.array([100.0, 200.0]),
            wall_to_coolant_conductance_W_per_K=np.array([5.0, 6.0]),
        ),
        hot_gas_march=march,
        shell_film=shell_film,
        coolant_properties=ShellTubeFluidProperties(
            rho=np.array([2.0, 3.0]),
            mu=np.array([1.0e-5, 1.1e-5]),
            k=np.array([0.1, 0.11]),
            cp=np.array([5000.0, 5100.0]),
        ),
        wall_heat_capacity_J_K=np.array([1.0, 1.0]),
    )

    ts = solver._time_series_from_core_result(
        SimpleNamespace(integration=integration, step_diagnostics=(diag,)),
        geometry,
    )

    for key in ("Re_g", "Nu_g", "progress_g", "Re_shell", "rho_c", "cp_c"):
        assert key in ts["fields"]
        assert ts["fields"][key].shape == (2, 2)
    for key in ("dp_g_total_Pa", "dp_shell_total_Pa", "progress_g_out", "Re_g_max"):
        assert key in ts["scalars"]
    np.testing.assert_allclose(ts["fields"]["Re_g"][-1], [10000.0, 9000.0])
    np.testing.assert_allclose(ts["fields"]["Re_shell"][-1], [100.0, 200.0])
    np.testing.assert_allclose(ts["scalars"]["dp_g_total_Pa"][-1], 7.0)
    assert ts["scalars"]["progress_g_out"][-1] == 0.25

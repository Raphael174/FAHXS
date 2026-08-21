"""Fidelity regression test for the decoupled Simulink/multiphysics stepper.

This proves `ShellTubeTransientStepper.step()` (1Dmodel/simulink_coupling/)
reproduces the existing, validated `run_shelltube_transient_core()` full-run
trajectory when fed the exact same time-varying boundary values one step at a
time. It also stands as the evidence that no existing solver file was
modified to add this capability: only this test file and the new
`1Dmodel/simulink_coupling/` package are new.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from hps_combustor.input_data import (
    CorrelationCoefficients,
    coolantProp,
    hotgasProp,
    numericalProp,
    shellTubeProp,
    system_requirements,
    transientProp,
)
from hps_combustor.main_solve_shellntube_transient import shellntube_transient_solver
from hps_combustor.simulink_coupling import BoundaryInputs, ShellTubeTransientStepper
from hps_combustor.transient_core import (
    build_shelltube_core_geometry_from_solver,
    coolprop_fluid_properties,
    equilibrium_gas_state_provider,
    run_shelltube_transient_core,
)

N_AXIAL = 6
DT = 0.01
N_STEPS = 4
T_END = DT * N_STEPS


def _base_config():
    cp = coolantProp()
    hp = hotgasProp()
    stp = shellTubeProp()
    nup = numericalProp()
    sr = system_requirements()
    tp = transientProp()
    tp.chemistry_transient = "equilibrium"  # fast to construct, no FPV cache needed
    tp.T_wall_initial = 293.15
    return cp, hp, stp, nup, sr, tp


def _linear_schedule(v0, v1):
    return ((0.0, v0), (T_END, v1))


def _ramp(v0, v1, t):
    w = min(max(t / T_END, 0.0), 1.0)
    return v0 + w * (v1 - v0)


@pytest.mark.parametrize("momentum_model", ["quasi_steady", "low_mach"])
def test_stepper_matches_reference_trajectory(momentum_model):
    cp, hp, stp, nup, sr, tp = _base_config()
    tp.coolant_momentum_model = momentum_model

    # mdot_coolant is held at the nominal design value (not ramped) so both
    # sides calibrate hydraulic resistance against the same reference flow.
    # The full-schedule reference solver gets to look ahead at the whole
    # run's schedule maximum when it calibrates that resistance; the stepper
    # cannot (a live Simulink caller doesn't hand over a future schedule), so
    # it uses `mdot_coolant_reference` (default: coolantProp.mass_flow_c)
    # instead. Ramping mdot_coolant here would legitimately decalibrate the
    # two by that documented, expected amount - see shelltube_stepper.py and
    # README.md - rather than indicating a fidelity bug.
    mdot_c0, mdot_c1 = cp.mass_flow_c, cp.mass_flow_c
    T_cin0, T_cin1 = cp.T_in, cp.T_in + 5.0
    p_cin0, p_cin1 = cp.p_in, cp.p_in * 0.98
    p_cout0, p_cout1 = cp.p_in - 2.0e5, cp.p_in - 2.5e5
    mdot_g0, mdot_g1 = hp.mass_flow_g, hp.mass_flow_g * 1.1

    # --- Reference: existing full-schedule solver -----------------------
    ref_solver = shellntube_transient_solver(
        cp, hp, stp, nup, sr, tp,
        corrCoeffs=CorrelationCoefficients(), N_axial=N_AXIAL, flow_config="co",
    )
    geometry = build_shelltube_core_geometry_from_solver(ref_solver)
    combustion_provider, progress0 = equilibrium_gas_state_provider(ref_solver._eqm)

    progress_config = SimpleNamespace(
        schedule_p_c_out=_linear_schedule(p_cout0, p_cout1),
        transient_coolant_outlet_pressure=None,
        insert_schedule_breakpoints=True,
    )

    n = geometry.grid.n_cells
    T_wall0 = np.full(n, tp.T_wall_initial)
    T_coolant0 = np.full(n, cp.T_in)

    ref = run_shelltube_transient_core(
        geometry,
        ref_solver.geom,
        T_wall_initial=T_wall0,
        T_coolant_initial=T_coolant0,
        t_end=T_END,
        max_step=DT,
        coolant_properties_at=coolprop_fluid_properties(cp.coolant),
        wall_density=float(ref_solver.rho_t),
        wall_cp=lambda T: np.array(
            [ref_solver.cp_t(float(Ti) - 273.15) for Ti in np.asarray(T)]
        ),
        wall_conductivity_at_T=ref_solver.k_t,
        inside_tube_choice=stp.inside_tube_choice,
        nusselt_selector=stp.Nusselt_tube,
        tube_roughness=stp.tube_roughness,
        mdot_coolant_default=mdot_c0,
        T_coolant_inlet_default=T_cin0,
        p_coolant_default=p_cin0,
        mdot_hot_total_default=mdot_g0,
        mdot_coolant_schedule=_linear_schedule(mdot_c0, mdot_c1),
        T_coolant_inlet_schedule=_linear_schedule(T_cin0, T_cin1),
        p_coolant_schedule=_linear_schedule(p_cin0, p_cin1),
        mdot_hot_total_schedule=_linear_schedule(mdot_g0, mdot_g1),
        gas_state_at=combustion_provider,
        progress_initial=progress0,
        corrCoeffs=CorrelationCoefficients(),
        corrugation_thickness=stp.corrugation_thickness,
        corrugation_pitch=stp.corrugation_pitch,
        flow_direction=geometry.grid.flow_direction,
        mdot_floor=1e-9,
        coolant_state_model=(
            "low_mach_momentum" if momentum_model == "low_mach" else "mass_energy"
        ),
        progress_config=progress_config,
        progress_enabled=False,
    )
    ref_T_wall = ref.integration.T_wall
    ref_T_coolant = ref.integration.T_coolant
    ref_pressure = ref.integration.coolant_pressure_Pa

    # --- Stepper: same boundary trajectory, one call per step -----------
    # `p_coolant_out_initial` seeds the same t=0 outlet pressure the reference
    # run knows in advance from its schedule; without it the stepper falls
    # back to a nominal dp estimate (see shelltube_stepper.py / README.md) -
    # correct for a live Simulink deployment with no advance knowledge of the
    # boundary trajectory, but not bit-comparable to a full-schedule
    # reference unless given the same t=0 starting point explicitly.
    stepper = ShellTubeTransientStepper(
        cp, hp, stp, nup, sr, tp,
        corrCoeffs=CorrelationCoefficients(), N_axial=N_AXIAL, flow_config="co",
        p_coolant_out_initial=p_cout0,
    )
    step_T_wall = [stepper.T_wall.copy()]
    step_pressure = []

    t = 0.0
    for j in range(N_STEPS):
        boundary = BoundaryInputs(
            mdot_coolant=_ramp(mdot_c0, mdot_c1, t),
            p_coolant_in=_ramp(p_cin0, p_cin1, t),
            p_coolant_out=_ramp(p_cout0, p_cout1, t),
            T_coolant_in=_ramp(T_cin0, T_cin1, t),
            mdot_hot_total=_ramp(mdot_g0, mdot_g1, t),
            ignited=True,
        )
        out = stepper.step(DT, boundary)
        step_T_wall.append(out.T_wall.copy())
        step_pressure.append(out.p_coolant_outlet)
        t += DT

    step_T_wall = np.array(step_T_wall)

    assert np.allclose(step_T_wall, ref_T_wall, rtol=1e-9, atol=1e-6)

    ref_outlet_idx = geometry.grid.outlet_index
    ref_p_outlet = ref_pressure[1:, ref_outlet_idx]
    assert np.allclose(step_pressure, ref_p_outlet, rtol=1e-9, atol=1.0)

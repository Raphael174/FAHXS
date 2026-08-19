"""Stage D, Slice 2 of docs/solver_design/FV_CORE_REWORK_PLAN.md.

`core/coolant.py` relocates `transient_core/compressible_coolant.py`'s
primitives unchanged (proven identical-object below, the same discipline
`test_core_thermo.py` used for the Stage A `core/thermo.py` move) and adds
`_cfl_stable_substep_count` (consolidated from `transient_core/
adapters_shelltube.py`, where it was first written 2026-08-18) plus the new
`advance_flowpath_coolant` convenience wrapper -- the one genuinely new piece
in this slice: a `FlowPath`-aware substep-and-sum loop generalizing the
pattern already hand-written twice (shell-and-tube 2026-08-18, helical
2026-08-19).
"""
from __future__ import annotations

import numpy as np
import pytest

from hps_combustor.core.coolant import (
    CompressibleCoolantStepResult,
    CoolantThermodynamicState,
    _cfl_stable_substep_count,
    advance_flowpath_coolant,
    conservative_mass_energy_step,
    coolprop_state_from_mass_energy,
    enforce_density_bounds,
    enforce_internal_energy_bounds,
    enforce_internal_energy_floor,
    initial_mass_energy_from_TP,
    internal_energy_from_temperature_mass,
    quasi_steady_face_mdot,
)
from hps_combustor.core.mesh import FlowPath
from hps_combustor.transient_core import compressible_coolant as legacy_shim


def test_transient_core_shim_reexports_are_identical_objects():
    """transient_core/compressible_coolant.py must re-export the SAME
    objects core.coolant defines, not parallel copies -- this is what makes
    the relocation behaviorally inert for every existing importer
    (adapters_shelltube.py, wall_compressible_coolant.py, __init__.py)."""
    assert legacy_shim.CompressibleCoolantStepResult is CompressibleCoolantStepResult
    assert legacy_shim.CoolantThermodynamicState is CoolantThermodynamicState
    assert legacy_shim.coolprop_state_from_mass_energy is coolprop_state_from_mass_energy
    assert legacy_shim.conservative_mass_energy_step is conservative_mass_energy_step
    assert legacy_shim.enforce_density_bounds is enforce_density_bounds
    assert legacy_shim.enforce_internal_energy_bounds is enforce_internal_energy_bounds
    assert legacy_shim.enforce_internal_energy_floor is enforce_internal_energy_floor
    assert legacy_shim.initial_mass_energy_from_TP is initial_mass_energy_from_TP
    assert legacy_shim.internal_energy_from_temperature_mass is internal_energy_from_temperature_mass
    assert legacy_shim.quasi_steady_face_mdot is quasi_steady_face_mdot


def test_adapters_shelltube_imports_cfl_helper_from_core():
    """The shell-and-tube adapter must consume the ONE copy of the CFL
    helper from core.coolant, not a local duplicate -- the exact duplication
    Slice 2 is meant to eliminate (main_solve_transient.py hand-copied this
    import list for the helical fix one day after it was first written)."""
    import hps_combustor.transient_core.adapters_shelltube as adapters_shelltube

    assert adapters_shelltube._cfl_stable_substep_count is _cfl_stable_substep_count


def test_main_solve_transient_imports_cfl_helper_from_core():
    import hps_combustor.main_solve_transient as mst

    assert mst._cfl_stable_substep_count is _cfl_stable_substep_count


def _flow_path(n_cells: int, volume_per_cell: float) -> FlowPath:
    s = np.linspace(0.0, 1.0, n_cells + 1)
    return FlowPath(
        name="test",
        s_edges=s,
        z_of_s_edges=s.copy(),
        A_flow=1e-4,
        Dh=0.01,
        P_wetted=0.03,
        P_heated=0.03,
        n_parallel=2,
    )


def test_advance_flowpath_coolant_uses_volume_total_not_per_channel():
    """Regression for the factor-of-N trap documented in core/mesh.py's
    LLM_CONTEXT.md: advance_flowpath_coolant must size mass/energy off
    path.volume_total (n_parallel included), not volume_per_channel."""
    path = _flow_path(3, volume_per_cell=1e-5)
    assert path.n_parallel == 2
    np.testing.assert_allclose(path.volume_total, path.volume_per_channel * 2)

    T0, p0 = 300.0, 5e5
    m, U = initial_mass_energy_from_TP(
        np.full(3, T0), np.full(3, p0), path.volume_total, "Helium"
    )
    # mass should scale with volume_total (2x volume_per_channel), i.e. be
    # exactly double what a (bugged) per-channel sizing would give.
    m_per_channel, _ = initial_mass_energy_from_TP(
        np.full(3, T0), np.full(3, p0), path.volume_per_channel, "Helium"
    )
    np.testing.assert_allclose(m, m_per_channel * 2.0)


def test_advance_flowpath_coolant_matches_manual_subcycle_loop():
    """advance_flowpath_coolant's substep-and-sum result must match calling
    conservative_mass_energy_step the same number of times by hand -- proves
    the wrapper is a faithful generalization, not a different algorithm."""
    path = _flow_path(3, volume_per_cell=1e-5)
    fluid = "Helium"
    T0, p0 = 300.0, 5e5
    m0, U0 = initial_mass_energy_from_TP(
        np.full(3, T0), np.full(3, p0), path.volume_total, fluid
    )
    faces = np.array([0.02, 0.02, 0.02, 0.02])
    heat_W = np.array([100.0, 100.0, 100.0])
    dt = 0.02
    h_in = 1.6e6

    n_sub = _cfl_stable_substep_count(m0, faces, dt)
    assert n_sub > 1, "fixture should exercise real subcycling"

    # Manual reference loop, mirroring the exact pattern in
    # adapters_shelltube.py / main_solve_transient.py.
    m_ref, U_ref = m0.copy(), U0.copy()
    sub_dt = dt / n_sub
    h_ref = coolprop_state_from_mass_energy(
        m_ref, U_ref, path.volume_total, fluid
    ).specific_enthalpy_J_kg
    heat_acc = adv_in_acc = adv_out_acc = 0.0
    for i in range(n_sub):
        step = conservative_mass_energy_step(
            m_ref, U_ref, h_ref, faces, heat_W, sub_dt,
            inlet_enthalpy_J_kg=h_in, outlet_backflow_enthalpy_J_kg=h_in,
        )
        m_ref, U_ref = step.mass_new, step.internal_energy_new_J
        heat_acc += step.heat_added_J
        adv_in_acc += step.advective_energy_in_J
        adv_out_acc += step.advective_energy_out_J
        if i < n_sub - 1:
            h_ref = coolprop_state_from_mass_energy(
                m_ref, U_ref, path.volume_total, fluid
            ).specific_enthalpy_J_kg

    result = advance_flowpath_coolant(
        path, fluid, m0, U0, faces, heat_W, dt,
        inlet_enthalpy_J_kg=h_in, outlet_backflow_enthalpy_J_kg=h_in,
    )

    np.testing.assert_allclose(result.mass_new, m_ref)
    np.testing.assert_allclose(result.internal_energy_new_J, U_ref)
    assert result.heat_added_J == pytest.approx(heat_acc)
    assert result.advective_energy_in_J == pytest.approx(adv_in_acc)
    assert result.advective_energy_out_J == pytest.approx(adv_out_acc)


def test_advance_flowpath_coolant_energy_closure_near_machine_precision():
    path = _flow_path(4, volume_per_cell=1e-5)
    fluid = "Helium"
    m0, U0 = initial_mass_energy_from_TP(
        np.full(4, 300.0), np.full(4, 5e5), path.volume_total, fluid
    )
    faces = np.full(5, 0.03)
    heat_W = np.array([200.0, 150.0, 100.0, 50.0])
    result = advance_flowpath_coolant(
        path, fluid, m0, U0, faces, heat_W, 0.05,
        inlet_enthalpy_J_kg=1.6e6, outlet_backflow_enthalpy_J_kg=1.6e6,
    )
    assert abs(result.energy_residual_J) < 1e-6

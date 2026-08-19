"""Stage D, Slice 3 of docs/solver_design/FV_CORE_REWORK_PLAN.md.

`core/wall_coolant_coupling.py` relocates `transient_core/
wall_compressible_coolant.py`'s wall-implicit/coolant-explicit step
unchanged. Proves bit-identical equivalence (same objects) to the shim's
re-exports.
"""
from __future__ import annotations

import numpy as np

from hps_combustor.core.coolant import initial_mass_energy_from_TP
from hps_combustor.core.wall_coolant_coupling import (
    WallCompressibleCoolantStepResult,
    semi_implicit_wall_compressible_coolant_step,
)
from hps_combustor.transient_core import wall_compressible_coolant as legacy_shim


def test_transient_core_shim_reexports_are_identical_objects():
    assert legacy_shim.WallCompressibleCoolantStepResult is WallCompressibleCoolantStepResult
    assert (
        legacy_shim.semi_implicit_wall_compressible_coolant_step
        is semi_implicit_wall_compressible_coolant_step
    )


def test_step_conserves_energy_to_machine_precision():
    n = 3
    volume = np.full(n, 1e-5)
    mass, U = initial_mass_energy_from_TP(np.full(n, 300.0), np.full(n, 5e5), volume, "Helium")
    T_c = np.full(n, 300.0)
    h_c = np.full(n, 1.6e6)
    faces = np.full(n + 1, 0.02)
    step = semi_implicit_wall_compressible_coolant_step(
        np.full(n, 350.0),
        np.full(n, 20.0),
        mass,
        U,
        T_c,
        h_c,
        faces,
        np.array([500.0, 500.0, 500.0]),
        np.full(n, 30.0),
        0.01,
        inlet_enthalpy_J_kg=1.6e6,
        outlet_backflow_enthalpy_J_kg=1.6e6,
    )
    assert abs(step.total_energy_residual_J) < 1e-6


def test_zero_dt_is_a_no_op():
    n = 2
    volume = np.full(n, 1e-5)
    mass, U = initial_mass_energy_from_TP(np.full(n, 300.0), np.full(n, 5e5), volume, "Helium")
    step = semi_implicit_wall_compressible_coolant_step(
        np.full(n, 350.0),
        np.full(n, 20.0),
        mass,
        U,
        np.full(n, 300.0),
        np.full(n, 1.6e6),
        np.full(n + 1, 0.0),
        np.full(n, 500.0),
        np.full(n, 30.0),
        0.0,
        inlet_enthalpy_J_kg=1.6e6,
        outlet_backflow_enthalpy_J_kg=1.6e6,
    )
    np.testing.assert_allclose(step.T_wall_new, 350.0)
    np.testing.assert_allclose(step.heat_wall_to_coolant_W, 0.0)

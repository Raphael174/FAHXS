"""Deprecation shim -- relocated to `hps_combustor.core.coolant` (Stage D,
Slice 2 of docs/solver_design/FV_CORE_REWORK_PLAN.md, 2026-08-19).

Pure re-export, not a reimplementation -- every name below is the SAME
object `core.coolant` defines, proven by
`tests/test_core_coolant.py::test_transient_core_shim_reexports_are_identical_objects`.
Existing imports from this module (this file, `adapters_shelltube.py`,
`wall_compressible_coolant.py`, `main_solve_transient.py`,
`transient_core/__init__.py`) continue to work unchanged.
"""

from __future__ import annotations

from hps_combustor.core.coolant import (
    CompressibleCoolantStepResult,
    CoolantThermodynamicState,
    coolprop_state_from_mass_energy,
    conservative_mass_energy_step,
    enforce_density_bounds,
    enforce_internal_energy_bounds,
    enforce_internal_energy_floor,
    initial_mass_energy_from_TP,
    internal_energy_from_temperature_mass,
    quasi_steady_face_mdot,
)

__all__ = [
    "CompressibleCoolantStepResult",
    "CoolantThermodynamicState",
    "coolprop_state_from_mass_energy",
    "conservative_mass_energy_step",
    "enforce_density_bounds",
    "enforce_internal_energy_bounds",
    "enforce_internal_energy_floor",
    "initial_mass_energy_from_TP",
    "internal_energy_from_temperature_mass",
    "quasi_steady_face_mdot",
]

"""Deprecation shim -- relocated to `hps_combustor.core.wall_coolant_coupling`
(Stage D, Slice 3 of docs/solver_design/FV_CORE_REWORK_PLAN.md, 2026-08-19).

Pure re-export, not a reimplementation -- proven by
`tests/test_core_wall_coolant_coupling.py::test_transient_core_shim_reexports_are_identical_objects`.
Existing imports from this module continue to work unchanged.
"""

from __future__ import annotations

from hps_combustor.core.wall_coolant_coupling import (
    WallCompressibleCoolantStepResult,
    semi_implicit_wall_compressible_coolant_step,
)

__all__ = [
    "WallCompressibleCoolantStepResult",
    "semi_implicit_wall_compressible_coolant_step",
]

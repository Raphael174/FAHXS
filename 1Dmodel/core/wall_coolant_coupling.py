"""Coupled wall and compressible-coolant mass/energy stepping — Stage D,
Slice 3 (docs/solver_design/FV_CORE_REWORK_PLAN.md).

Relocated, unchanged, from `transient_core/wall_compressible_coolant.py`
(same "pure infra, move as-is" pattern as `core/coolant.py` for D2 —
independently classified as movable pure infra when Stage D exploration
first inventoried `transient_core/`). `transient_core/
wall_compressible_coolant.py` now re-exports these same objects as a shim,
proven identical-object by
`tests/test_core_wall_coolant_coupling.py::test_transient_core_shim_reexports_are_identical_objects`.

Distinct from `transient_core/wall_coolant.py` (the OLDER temperature-only
fully-implicit 2×2 model, not touched here or by this move — still a
candidate for direct reuse in a future non-compressible-coolant core path).

This is the wall-implicit/coolant-explicit pattern the design doc's §2/§3.5
"do not replace the 2×2 without measuring" caution refers to (see
`docs/solver_design/FV_CORE_REWORK_PLAN.md`'s 2026-08-18 CFL note, item 5):
the wall is solved implicitly against the coolant temperature reconstructed
at the START of the step, not the fully-coupled 2×2 `wall_coolant.py` uses
for the temperature-only model, because `T(m, U)` needs a CoolProp inversion
that cannot be inlined into a 2×2 the way `T` alone can.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .coolant import CompressibleCoolantStepResult, conservative_mass_energy_step


@dataclass(frozen=True)
class WallCompressibleCoolantStepResult:
    """Result from one wall + compressible-coolant step."""

    T_wall_new: np.ndarray
    coolant: CompressibleCoolantStepResult
    heat_wall_to_coolant_W: np.ndarray
    wall_internal_energy_old_J: float
    wall_internal_energy_new_J: float
    hot_heat_added_J: float
    total_energy_residual_J: float


def semi_implicit_wall_compressible_coolant_step(
    T_wall_old,
    wall_heat_capacity,
    coolant_mass,
    coolant_internal_energy_J,
    coolant_temperature,
    coolant_specific_enthalpy_J_kg,
    face_mdot,
    hot_heat_W,
    wall_to_coolant_conductance_W_per_K,
    dt: float,
    *,
    inlet_enthalpy_J_kg: float | None = None,
    outlet_backflow_enthalpy_J_kg: float | None = None,
    mass_floor: float = 1.0e-12,
) -> WallCompressibleCoolantStepResult:
    """Advance wall energy and coolant mass/energy conservatively.

    The wall is integrated with local implicit cooling against the coolant
    temperature reconstructed at the beginning of the step:

    ```text
    Cw (Tw_new - Tw_old)/dt = Qhot - G (Tw_new - Tc_old)
    ```

    The resulting `Qwall_to_coolant = G (Tw_new - Tc_old)` is then used as the
    heat source in the conservative coolant mass/energy update. This is
    first-order in time, but it preserves total wall+coolant energy exactly up
    to the conservative coolant residual.
    """

    Tw_old = _as_1d("T_wall_old", T_wall_old)
    Cw = _as_1d("wall_heat_capacity", wall_heat_capacity)
    Tc = _as_1d("coolant_temperature", coolant_temperature)
    Qhot = _as_1d("hot_heat_W", hot_heat_W)
    G = _as_1d("wall_to_coolant_conductance_W_per_K", wall_to_coolant_conductance_W_per_K)
    _check_same_len(Tw_old, Cw, Tc, Qhot, G)
    if dt < 0.0:
        raise ValueError("dt must be non-negative")
    if np.any(Cw <= 0.0):
        raise ValueError("wall_heat_capacity must be strictly positive")
    if np.any(G < 0.0):
        raise ValueError("wall_to_coolant_conductance_W_per_K must be non-negative")

    if dt == 0.0:
        Tw_new = Tw_old.copy()
        Qwc = np.zeros_like(Tw_old)
    else:
        cw_dt = Cw / dt
        Tw_new = (cw_dt * Tw_old + Qhot + G * Tc) / (cw_dt + G)
        Qwc = G * (Tw_new - Tc)
        U_old = _as_1d("coolant_internal_energy_J", coolant_internal_energy_J)
        min_coolant_U = 0.05 * U_old
        max_removal_W = np.maximum((U_old - min_coolant_U) / dt, 0.0)
        Qwc = np.maximum(Qwc, -max_removal_W)

    coolant = conservative_mass_energy_step(
        coolant_mass,
        coolant_internal_energy_J,
        coolant_specific_enthalpy_J_kg,
        face_mdot,
        Qwc,
        dt,
        inlet_enthalpy_J_kg=inlet_enthalpy_J_kg,
        outlet_backflow_enthalpy_J_kg=outlet_backflow_enthalpy_J_kg,
        mass_floor=mass_floor,
    )

    wall_U_old = float(np.sum(Cw * Tw_old))
    wall_U_new = float(np.sum(Cw * Tw_new))
    hot_heat = float(np.sum(Qhot) * dt)
    total_residual = (
        (wall_U_new - wall_U_old)
        + np.sum(coolant.internal_energy_new_J - _as_1d("coolant_internal_energy_J", coolant_internal_energy_J))
        - (hot_heat + coolant.advective_energy_in_J - coolant.advective_energy_out_J)
    )

    return WallCompressibleCoolantStepResult(
        T_wall_new=Tw_new,
        coolant=coolant,
        heat_wall_to_coolant_W=Qwc,
        wall_internal_energy_old_J=wall_U_old,
        wall_internal_energy_new_J=wall_U_new,
        hot_heat_added_J=hot_heat,
        total_energy_residual_J=float(total_residual),
    )


def _as_1d(name, value) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _check_same_len(*arrays) -> None:
    n = arrays[0].size
    for arr in arrays[1:]:
        if arr.size != n:
            raise ValueError("all array inputs must have the same length")

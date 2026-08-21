"""Coupled wall and coolant finite-volume stepping.

This module is geometry-independent. It advances a wall thermal inventory and a
coolant thermal inventory in each axial cell while treating local wall-coolant
exchange implicitly and coolant advection by implicit upwind.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class WallCoolantStepResult:
    """Result from one coupled wall-coolant finite-volume step."""

    T_wall_new: np.ndarray
    T_coolant_new: np.ndarray
    T_coolant_outlet: float
    heat_wall_to_coolant_W: np.ndarray
    wall_internal_energy_old_J: float
    wall_internal_energy_new_J: float
    coolant_internal_energy_old_J: float
    coolant_internal_energy_new_J: float
    hot_heat_added_J: float
    advective_energy_in_J: float
    advective_energy_out_J: float
    energy_residual_J: float


def implicit_wall_coolant_step(
    T_wall_old,
    T_coolant_old,
    wall_heat_capacity,
    coolant_heat_capacity,
    coolant_cp,
    mdot_coolant,
    T_coolant_inlet,
    hot_heat_W,
    wall_to_coolant_conductance_W_per_K,
    dt,
    *,
    flow_direction=1,
    mdot_floor=1e-12,
    coolant_cp_inlet=None,
) -> WallCoolantStepResult:
    """Advance wall and coolant temperatures by one implicit linear step.

    The solved cell equations are:

    ```text
    Cw_i (Tw_i^{n+1} - Tw_i^n) / dt =
        Qhot_i - G_i (Tw_i^{n+1} - Tc_i^{n+1})

    Cc_i (Tc_i^{n+1} - Tc_i^n) / dt =
        mdot cp_i (Tup_i^{n+1} - Tc_i^{n+1})
      + G_i (Tw_i^{n+1} - Tc_i^{n+1})
    ```

    where `Tup` is the inlet temperature for the upstream boundary or the
    already-solved upstream cell temperature. `hot_heat_W` is positive from hot
    gas into the wall. `wall_to_coolant_conductance_W_per_K` is the local total
    conductance between wall mean temperature and coolant bulk temperature.

    This is a building block for the full HX transient core: geometry adapters
    should compute `hot_heat_W`, `G`, heat capacities, and flow direction.
    """

    T_wall_old = _as_1d("T_wall_old", T_wall_old)
    T_coolant_old = _as_1d("T_coolant_old", T_coolant_old)
    wall_heat_capacity = _as_1d("wall_heat_capacity", wall_heat_capacity)
    coolant_heat_capacity = _as_1d("coolant_heat_capacity", coolant_heat_capacity)
    coolant_cp = _as_1d("coolant_cp", coolant_cp)
    hot_heat_W = _as_1d("hot_heat_W", hot_heat_W)
    conductance = _as_1d(
        "wall_to_coolant_conductance_W_per_K",
        wall_to_coolant_conductance_W_per_K,
    )
    _check_same_len(
        T_wall_old,
        T_coolant_old,
        wall_heat_capacity,
        coolant_heat_capacity,
        coolant_cp,
        hot_heat_W,
        conductance,
    )

    if dt < 0.0:
        raise ValueError("dt must be non-negative")
    if flow_direction not in (-1, 1):
        raise ValueError("flow_direction must be +1 or -1")
    if np.any(wall_heat_capacity <= 0.0) or np.any(coolant_heat_capacity <= 0.0):
        raise ValueError("heat capacities must be strictly positive")
    if np.any(coolant_cp <= 0.0):
        raise ValueError("coolant_cp must be strictly positive")
    if np.any(conductance < 0.0):
        raise ValueError("wall_to_coolant_conductance_W_per_K must be non-negative")

    n = T_wall_old.size
    T_wall_new = np.empty_like(T_wall_old, dtype=float)
    T_coolant_new = np.empty_like(T_coolant_old, dtype=float)
    m = abs(float(mdot_coolant))

    if dt == 0.0:
        T_wall_new[:] = T_wall_old
        T_coolant_new[:] = T_coolant_old
    elif flow_direction == 1:
        for i in range(n):
            upstream = float(T_coolant_inlet) if i == 0 else T_coolant_new[i - 1]
            T_wall_new[i], T_coolant_new[i] = _solve_cell(
                T_wall_old[i],
                T_coolant_old[i],
                wall_heat_capacity[i],
                coolant_heat_capacity[i],
                coolant_cp[i],
                m,
                upstream,
                hot_heat_W[i],
                conductance[i],
                dt,
                mdot_floor,
            )
    else:
        for i in range(n - 1, -1, -1):
            upstream = float(T_coolant_inlet) if i == n - 1 else T_coolant_new[i + 1]
            T_wall_new[i], T_coolant_new[i] = _solve_cell(
                T_wall_old[i],
                T_coolant_old[i],
                wall_heat_capacity[i],
                coolant_heat_capacity[i],
                coolant_cp[i],
                m,
                upstream,
                hot_heat_W[i],
                conductance[i],
                dt,
                mdot_floor,
            )

    inlet_idx = 0 if flow_direction == 1 else n - 1
    outlet_idx = n - 1 if flow_direction == 1 else 0
    cp_boundary = coolant_cp[inlet_idx] if coolant_cp_inlet is None else float(coolant_cp_inlet)
    adv_in = 0.0 if m <= mdot_floor else m * cp_boundary * float(T_coolant_inlet) * dt
    adv_out = 0.0 if m <= mdot_floor else m * coolant_cp[outlet_idx] * T_coolant_new[outlet_idx] * dt

    wall_U_old = float(np.sum(wall_heat_capacity * T_wall_old))
    wall_U_new = float(np.sum(wall_heat_capacity * T_wall_new))
    coolant_U_old = float(np.sum(coolant_heat_capacity * T_coolant_old))
    coolant_U_new = float(np.sum(coolant_heat_capacity * T_coolant_new))
    hot_heat = float(np.sum(hot_heat_W) * dt)
    residual = (
        (wall_U_new + coolant_U_new)
        - (wall_U_old + coolant_U_old)
        - (hot_heat + adv_in - adv_out)
    )

    return WallCoolantStepResult(
        T_wall_new=T_wall_new,
        T_coolant_new=T_coolant_new,
        T_coolant_outlet=float(T_coolant_new[outlet_idx]),
        heat_wall_to_coolant_W=conductance * (T_wall_new - T_coolant_new),
        wall_internal_energy_old_J=wall_U_old,
        wall_internal_energy_new_J=wall_U_new,
        coolant_internal_energy_old_J=coolant_U_old,
        coolant_internal_energy_new_J=coolant_U_new,
        hot_heat_added_J=hot_heat,
        advective_energy_in_J=float(adv_in),
        advective_energy_out_J=float(adv_out),
        energy_residual_J=float(residual),
    )


def _solve_cell(
    T_wall_old,
    T_coolant_old,
    wall_heat_capacity,
    coolant_heat_capacity,
    coolant_cp,
    mdot,
    T_upstream,
    hot_heat_W,
    conductance,
    dt,
    mdot_floor,
):
    cw_dt = wall_heat_capacity / dt
    cc_dt = coolant_heat_capacity / dt
    adv = 0.0 if mdot <= mdot_floor else mdot * coolant_cp

    # Matrix:
    # (cw/dt + G) Tw_new - G Tc_new = cw/dt Tw_old + Qhot
    # -G Tw_new + (cc/dt + adv + G) Tc_new = cc/dt Tc_old + adv Tup
    a11 = cw_dt + conductance
    a12 = -conductance
    a21 = -conductance
    a22 = cc_dt + adv + conductance
    b1 = cw_dt * T_wall_old + hot_heat_W
    b2 = cc_dt * T_coolant_old + adv * T_upstream

    det = a11 * a22 - a12 * a21
    if det <= 0.0 or not np.isfinite(det):
        raise FloatingPointError("invalid local wall-coolant solve")
    T_wall_new = (b1 * a22 - a12 * b2) / det
    T_coolant_new = (a11 * b2 - b1 * a21) / det
    return T_wall_new, T_coolant_new


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

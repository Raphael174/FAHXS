"""Fixed-step driver for transient wall + coolant states."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .wall_coolant import WallCoolantStepResult, implicit_wall_coolant_step


@dataclass(frozen=True)
class WallCoolantStepInputs:
    """Inputs required by one wall/coolant implicit step."""

    wall_heat_capacity: np.ndarray
    coolant_heat_capacity: np.ndarray
    coolant_cp: np.ndarray
    mdot_coolant: float
    T_coolant_inlet: float
    hot_heat_W: np.ndarray
    wall_to_coolant_conductance_W_per_K: np.ndarray
    flow_direction: int = 1
    coolant_cp_inlet: float | None = None


@dataclass(frozen=True)
class WallCoolantIntegrationResult:
    """Time history from `integrate_wall_coolant_fixed_step`."""

    t: np.ndarray
    T_wall: np.ndarray
    T_coolant: np.ndarray
    T_coolant_outlet: np.ndarray
    hot_heat_added_J: np.ndarray
    advective_energy_in_J: np.ndarray
    advective_energy_out_J: np.ndarray
    energy_residual_J: np.ndarray
    heat_wall_to_coolant_W: np.ndarray
    last_step: WallCoolantStepResult | None


StepInputBuilder = Callable[
    [float, np.ndarray, np.ndarray],
    WallCoolantStepInputs,
]


def fixed_time_grid(
    *,
    t_end: float,
    max_step: float,
    t_eval=None,
    schedules=(),
    include_endpoints: bool = True,
) -> np.ndarray:
    """Build a monotonic fixed-step grid with schedule breakpoints inserted.

    Schedules use the project convention `((time_s, value), ...)`. Only the
    time column is used here. Duplicate and out-of-range points are removed.
    """

    if t_end < 0.0:
        raise ValueError("t_end must be non-negative")
    if max_step <= 0.0:
        raise ValueError("max_step must be positive")

    base = np.arange(0.0, float(t_end) + 0.5 * float(max_step), float(max_step))
    points = [base]
    if include_endpoints:
        points.append(np.array([0.0, float(t_end)]))
    if t_eval is not None:
        points.append(_as_1d("t_eval", t_eval))

    for schedule in schedules:
        if not schedule:
            continue
        times = []
        for row in schedule:
            if row is None or len(row) < 1:
                continue
            times.append(float(row[0]))
        if times:
            points.append(np.asarray(times, dtype=float))

    grid = np.concatenate(points)
    grid = grid[np.isfinite(grid)]
    grid = grid[(grid >= 0.0) & (grid <= float(t_end))]
    return np.unique(np.round(grid, decimals=12))


def integrate_wall_coolant_fixed_step(
    *,
    T_wall_initial,
    T_coolant_initial,
    t_eval,
    step_inputs: StepInputBuilder,
    mdot_floor: float = 1e-12,
) -> WallCoolantIntegrationResult:
    """Integrate wall and coolant states over a supplied time grid.

    `step_inputs(t, T_wall, T_coolant)` is called at the start of each interval.
    Geometry adapters should use that callback to compute properties,
    conductances, hot-side heat input, flow direction, and inlet conditions for
    the current state.
    """

    T_wall0 = _as_1d("T_wall_initial", T_wall_initial)
    T_coolant0 = _as_1d("T_coolant_initial", T_coolant_initial)
    if T_wall0.size != T_coolant0.size:
        raise ValueError("T_wall_initial and T_coolant_initial must have the same length")

    t = _as_1d("t_eval", t_eval)
    if t.size < 1:
        raise ValueError("t_eval must not be empty")
    if np.any(np.diff(t) < 0.0):
        raise ValueError("t_eval must be monotonically nondecreasing")

    n_time = t.size
    n_cells = T_wall0.size
    T_wall = np.zeros((n_time, n_cells), dtype=float)
    T_coolant = np.zeros((n_time, n_cells), dtype=float)
    T_coolant_outlet = np.zeros(n_time, dtype=float)
    hot_heat_added_J = np.zeros(n_time, dtype=float)
    advective_energy_in_J = np.zeros(n_time, dtype=float)
    advective_energy_out_J = np.zeros(n_time, dtype=float)
    energy_residual_J = np.zeros(n_time, dtype=float)
    heat_wall_to_coolant_W = np.zeros((n_time, n_cells), dtype=float)

    T_wall[0] = T_wall0
    T_coolant[0] = T_coolant0
    T_coolant_outlet[0] = T_coolant0[-1]
    last_step = None

    for j in range(n_time - 1):
        dt = float(t[j + 1] - t[j])
        inputs = step_inputs(float(t[j]), T_wall[j].copy(), T_coolant[j].copy())
        step = implicit_wall_coolant_step(
            T_wall[j],
            T_coolant[j],
            inputs.wall_heat_capacity,
            inputs.coolant_heat_capacity,
            inputs.coolant_cp,
            inputs.mdot_coolant,
            inputs.T_coolant_inlet,
            inputs.hot_heat_W,
            inputs.wall_to_coolant_conductance_W_per_K,
            dt,
            flow_direction=inputs.flow_direction,
            mdot_floor=mdot_floor,
            coolant_cp_inlet=inputs.coolant_cp_inlet,
        )
        T_wall[j + 1] = step.T_wall_new
        T_coolant[j + 1] = step.T_coolant_new
        T_coolant_outlet[j + 1] = step.T_coolant_outlet
        hot_heat_added_J[j + 1] = step.hot_heat_added_J
        advective_energy_in_J[j + 1] = step.advective_energy_in_J
        advective_energy_out_J[j + 1] = step.advective_energy_out_J
        energy_residual_J[j + 1] = step.energy_residual_J
        heat_wall_to_coolant_W[j + 1] = step.heat_wall_to_coolant_W
        last_step = step

    return WallCoolantIntegrationResult(
        t=t,
        T_wall=T_wall,
        T_coolant=T_coolant,
        T_coolant_outlet=T_coolant_outlet,
        hot_heat_added_J=hot_heat_added_J,
        advective_energy_in_J=advective_energy_in_J,
        advective_energy_out_J=advective_energy_out_J,
        energy_residual_J=energy_residual_J,
        heat_wall_to_coolant_W=heat_wall_to_coolant_W,
        last_step=last_step,
    )


def _as_1d(name: str, value) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr

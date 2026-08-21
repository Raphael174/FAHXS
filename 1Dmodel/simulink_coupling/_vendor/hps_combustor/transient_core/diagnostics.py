"""Diagnostics for transient finite-volume HX states."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class EnergyAudit:
    """Scaled first-law residual."""

    residual_J: float
    scale_J: float
    relative: float
    passes: bool


@dataclass(frozen=True)
class TimescaleAudit:
    """Residence and wall timescale audit for quasi-steady assumptions."""

    coolant_residence_s: float
    wall_tau_min_s: float
    hot_residence_s: float | None
    boundary_tau_s: float | None
    coolant_to_wall_ratio: float
    coolant_to_boundary_ratio: float | None
    hot_to_wall_ratio: float | None
    hot_to_boundary_ratio: float | None
    coolant_quasi_steady_ok: bool
    hot_quasi_steady_ok: bool


def energy_audit(
    residual_J: float,
    *,
    heat_added_J: float = 0.0,
    advective_energy_in_J: float = 0.0,
    advective_energy_out_J: float = 0.0,
    stored_energy_change_J: float = 0.0,
    relative_tol: float = 1e-8,
    absolute_floor_J: float = 1.0,
) -> EnergyAudit:
    """Scale an energy residual against the largest relevant energy term."""

    terms = np.array(
        [
            abs(float(heat_added_J)),
            abs(float(advective_energy_in_J)),
            abs(float(advective_energy_out_J)),
            abs(float(stored_energy_change_J)),
            abs(float(absolute_floor_J)),
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(terms)):
        raise ValueError("energy terms must be finite")
    if relative_tol < 0.0:
        raise ValueError("relative_tol must be non-negative")

    scale = float(np.max(terms))
    residual = float(residual_J)
    if not np.isfinite(residual):
        raise ValueError("residual_J must be finite")
    relative = abs(residual) / scale
    return EnergyAudit(
        residual_J=residual,
        scale_J=scale,
        relative=float(relative),
        passes=bool(relative <= relative_tol),
    )


def residence_time_s(
    density,
    volume,
    mdot,
    *,
    mdot_floor=1e-12,
) -> float:
    """Return total fluid residence time from inventory mass divided by flow."""

    rho = _as_1d("density", density)
    vol = _as_1d("volume", volume)
    _check_same_len(rho, vol)
    if np.any(rho <= 0.0) or np.any(vol <= 0.0):
        raise ValueError("density and volume must be strictly positive")

    m = abs(float(mdot))
    if m <= mdot_floor:
        return float("inf")
    return float(np.sum(rho * vol) / m)


def wall_time_constant_s(
    wall_heat_capacity,
    hot_conductance_W_K,
    coolant_conductance_W_K,
) -> np.ndarray:
    """Return local wall thermal time constants.

    The denominator is the local total conductance removing energy from the wall
    state after linearization around the current condition.
    """

    Cw = _as_1d("wall_heat_capacity", wall_heat_capacity)
    Gh = _as_1d("hot_conductance_W_K", hot_conductance_W_K)
    Gc = _as_1d("coolant_conductance_W_K", coolant_conductance_W_K)
    _check_same_len(Cw, Gh, Gc)
    if np.any(Cw <= 0.0):
        raise ValueError("wall_heat_capacity must be strictly positive")
    if np.any(Gh < 0.0) or np.any(Gc < 0.0):
        raise ValueError("conductances must be non-negative")

    G = Gh + Gc
    tau = np.full_like(Cw, np.inf, dtype=float)
    active = G > 0.0
    tau[active] = Cw[active] / G[active]
    return tau


def timescale_audit(
    *,
    coolant_residence_s: float,
    wall_tau_s,
    hot_residence_s: float | None = None,
    boundary_tau_s: float | None = None,
    warning_ratio: float = 0.05,
) -> TimescaleAudit:
    """Assess whether fluid quasi-steady assumptions are timescale-consistent."""

    wall_tau = _as_1d("wall_tau_s", wall_tau_s)
    if np.any(wall_tau <= 0.0):
        raise ValueError("wall_tau_s must be strictly positive")
    if warning_ratio < 0.0:
        raise ValueError("warning_ratio must be non-negative")

    tau_wall_min = float(np.min(wall_tau))
    tau_c = float(coolant_residence_s)
    if tau_c < 0.0:
        raise ValueError("coolant_residence_s must be non-negative")
    tau_bc = None if boundary_tau_s is None else float(boundary_tau_s)
    tau_h = None if hot_residence_s is None else float(hot_residence_s)
    if tau_bc is not None and tau_bc <= 0.0:
        raise ValueError("boundary_tau_s must be positive when provided")
    if tau_h is not None and tau_h < 0.0:
        raise ValueError("hot_residence_s must be non-negative when provided")

    coolant_to_wall = tau_c / tau_wall_min
    coolant_to_boundary = None if tau_bc is None else tau_c / tau_bc
    hot_to_wall = None if tau_h is None else tau_h / tau_wall_min
    hot_to_boundary = None if (tau_h is None or tau_bc is None) else tau_h / tau_bc

    coolant_ok = coolant_to_wall <= warning_ratio
    if coolant_to_boundary is not None:
        coolant_ok = coolant_ok and coolant_to_boundary <= warning_ratio

    hot_ok = True
    if hot_to_wall is not None:
        hot_ok = hot_ok and hot_to_wall <= warning_ratio
    if hot_to_boundary is not None:
        hot_ok = hot_ok and hot_to_boundary <= warning_ratio

    return TimescaleAudit(
        coolant_residence_s=tau_c,
        wall_tau_min_s=tau_wall_min,
        hot_residence_s=tau_h,
        boundary_tau_s=tau_bc,
        coolant_to_wall_ratio=float(coolant_to_wall),
        coolant_to_boundary_ratio=(
            None if coolant_to_boundary is None else float(coolant_to_boundary)
        ),
        hot_to_wall_ratio=None if hot_to_wall is None else float(hot_to_wall),
        hot_to_boundary_ratio=(
            None if hot_to_boundary is None else float(hot_to_boundary)
        ),
        coolant_quasi_steady_ok=bool(coolant_ok),
        hot_quasi_steady_ok=bool(hot_ok),
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


def _check_same_len(*arrays) -> None:
    n = arrays[0].size
    for arr in arrays[1:]:
        if arr.size != n:
            raise ValueError("all array inputs must have the same length")

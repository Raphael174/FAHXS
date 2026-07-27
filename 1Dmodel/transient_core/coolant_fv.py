"""Transient coolant finite-volume update.

This module is intentionally geometry-independent. Callers provide cell volumes,
fluid properties, heat input per cell, and a flow direction. The update solves a
one-dimensional implicit-upwind temperature step for:

    rho V cp dT_i/dt = mdot cp (T_upwind - T_i) + Q_i

`Q_i` is the total heat rate into coolant cell i [W]. Positive heat warms the
coolant. At zero flow the advection term is removed, so each cell becomes a
stagnant thermal inventory coupled only through its local heat source.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CoolantStepResult:
    """Result from one coolant finite-volume step."""

    T_new: np.ndarray
    T_outlet: float
    internal_energy_old_J: float
    internal_energy_new_J: float
    heat_added_J: float
    advective_energy_in_J: float
    advective_energy_out_J: float
    energy_residual_J: float


def implicit_upwind_step(
    T_old,
    rho,
    cp,
    volume,
    mdot,
    T_inlet,
    heat_W,
    dt,
    *,
    flow_direction=1,
    mdot_floor=1e-12,
    cp_inlet=None,
) -> CoolantStepResult:
    """Advance coolant cell temperatures by one implicit-upwind step.

    Parameters
    ----------
    T_old, rho, cp, volume, heat_W:
        One-dimensional arrays with one entry per coolant cell. `heat_W` is the
        total heat rate into the corresponding cell [W], not per unit length.
    mdot:
        Coolant mass flow rate [kg/s]. Its magnitude is used; direction is
        controlled by `flow_direction`.
    T_inlet:
        Physical inlet temperature at the upstream boundary for the selected
        direction [K].
    dt:
        Time step [s].
    flow_direction:
        `+1` means inlet at cell 0 and flow toward increasing index. `-1` means
        inlet at the last cell and flow toward decreasing index.
    mdot_floor:
        If `abs(mdot) <= mdot_floor`, advection is disabled.
    cp_inlet:
        Optional inlet specific heat [J/kg/K] for energy diagnostics. If omitted,
        the upstream boundary uses the first/last cell `cp`.

    Returns
    -------
    CoolantStepResult
        Includes new temperatures and a first-law residual:

            dU - (Q*dt + H_in*dt - H_out*dt)

        For constant properties this should be near numerical roundoff. With
        temperature-dependent properties it is a diagnostic based on the supplied
        frozen-in-step `cp` values.
    """

    T_old = _as_1d("T_old", T_old)
    rho = _as_1d("rho", rho)
    cp = _as_1d("cp", cp)
    volume = _as_1d("volume", volume)
    heat_W = _as_1d("heat_W", heat_W)
    _check_same_len(T_old, rho, cp, volume, heat_W)

    if dt < 0.0:
        raise ValueError("dt must be non-negative")
    if flow_direction not in (-1, 1):
        raise ValueError("flow_direction must be +1 or -1")
    if np.any(rho <= 0.0) or np.any(cp <= 0.0) or np.any(volume <= 0.0):
        raise ValueError("rho, cp, and volume must be strictly positive")

    n = T_old.size
    thermal_mass = rho * volume * cp
    T_new = np.empty_like(T_old, dtype=float)
    m = abs(float(mdot))

    if dt == 0.0:
        T_new[:] = T_old
    elif m <= mdot_floor:
        T_new[:] = T_old + dt * heat_W / thermal_mass
    elif flow_direction == 1:
        for i in range(n):
            cp_adv = cp[i]
            a = m * cp_adv * dt
            upstream = float(T_inlet) if i == 0 else T_new[i - 1]
            T_new[i] = (thermal_mass[i] * T_old[i] + a * upstream + dt * heat_W[i]) / (
                thermal_mass[i] + a
            )
    else:
        for i in range(n - 1, -1, -1):
            cp_adv = cp[i]
            a = m * cp_adv * dt
            upstream = float(T_inlet) if i == n - 1 else T_new[i + 1]
            T_new[i] = (thermal_mass[i] * T_old[i] + a * upstream + dt * heat_W[i]) / (
                thermal_mass[i] + a
            )

    inlet_idx = 0 if flow_direction == 1 else n - 1
    outlet_idx = n - 1 if flow_direction == 1 else 0
    cp_boundary = cp[inlet_idx] if cp_inlet is None else float(cp_inlet)
    cp_out = cp[outlet_idx]
    adv_in = 0.0 if m <= mdot_floor else m * cp_boundary * float(T_inlet) * dt
    adv_out = 0.0 if m <= mdot_floor else m * cp_out * T_new[outlet_idx] * dt
    U_old = float(np.sum(thermal_mass * T_old))
    U_new = float(np.sum(thermal_mass * T_new))
    heat = float(np.sum(heat_W) * dt)
    residual = (U_new - U_old) - (heat + adv_in - adv_out)

    return CoolantStepResult(
        T_new=T_new,
        T_outlet=float(T_new[outlet_idx]),
        internal_energy_old_J=U_old,
        internal_energy_new_J=U_new,
        heat_added_J=heat,
        advective_energy_in_J=float(adv_in),
        advective_energy_out_J=float(adv_out),
        energy_residual_J=float(residual),
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

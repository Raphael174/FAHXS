"""Compressible coolant mass/energy finite-volume kernel — Stage D, Slice 2
(docs/solver_design/FV_CORE_REWORK_PLAN.md).

Relocated, unchanged, from `transient_core/compressible_coolant.py` (the
module docstring there already described it as "geometry-independent" --
this is a pure move, not a rewrite; `transient_core/compressible_coolant.py`
now re-exports these same objects as a shim, proven identical by
`tests/test_core_coolant.py::test_transient_core_shim_reexports_are_identical_objects`,
the same pattern Stage A used for `core/thermo.py`).

```text
dm_i/dt = mdot_{i-1/2} - mdot_{i+1/2}

dU_i/dt = mdot_{i-1/2} h_up,i-1/2
        - mdot_{i+1/2} h_up,i+1/2
        + Q_i
```

The update is explicit forward Euler in the conserved variables `(m, U)` --
unconditionally unstable once a step exceeds roughly one cell's residence
time. This is a REAL, previously-hit failure mode: see
`_cfl_stable_substep_count`'s docstring and
`docs/solver_design/FV_CORE_REWORK_PLAN.md`'s 2026-08-18/19 notes for the
full story (both the shell-and-tube and helical transient solvers crashed on
their own documented validation cases before this was found and fixed).
`_cfl_stable_substep_count` lived in `transient_core/adapters_shelltube.py`
when first written (2026-08-18) and was then hand-duplicated into
`main_solve_transient.py`'s import list a day later for the helical fix --
this move consolidates it to one place, and both solvers now import it from
here.

`advance_flowpath_coolant` is the one genuinely NEW piece: a
`core/mesh.py::FlowPath`-aware convenience wrapper that subcycles
`conservative_mass_energy_step` using `path.volume_total` (matching
`FlowPath`'s per-channel/total convention -- see `core/mesh.py`'s
`LLM_CONTEXT.md` on the factor-of-N trap), so a future `core/residual.py`
caller (Stage D3) does not need to know about substep bookkeeping. It is a
thin wrapper around the primitives below, not new physics.

Momentum is quasi-steady here: face mass fluxes are algebraic inputs, or can
be estimated from instantaneous pressure differences and hydraulic
resistances (`quasi_steady_face_mdot`) -- see `core/momentum.py` for the
`FlowPath`-generalized version of the resistance-network closure both legacy
transient adapters use.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from CoolProp.CoolProp import PropsSI


@dataclass(frozen=True)
class CompressibleCoolantStepResult:
    """Result from one conservative mass/energy coolant step."""

    mass_new: np.ndarray
    internal_energy_new_J: np.ndarray
    specific_internal_energy_new_J_kg: np.ndarray
    face_mdot: np.ndarray
    face_enthalpy: np.ndarray
    heat_added_J: float
    advective_energy_in_J: float
    advective_energy_out_J: float
    mass_residual_kg: float
    energy_residual_J: float


@dataclass(frozen=True)
class CoolantThermodynamicState:
    """Cellwise state reconstructed from conserved coolant variables."""

    pressure: np.ndarray
    temperature: np.ndarray
    density: np.ndarray
    specific_internal_energy_J_kg: np.ndarray
    specific_enthalpy_J_kg: np.ndarray
    cp: np.ndarray
    mu: np.ndarray
    k: np.ndarray


def coolprop_state_from_mass_energy(mass, internal_energy_J, volume, fluid: str) -> CoolantThermodynamicState:
    """Reconstruct cell thermodynamic/transport state from `m`, `U`, and `V`.

    CoolProp receives density and specific internal energy:

    ```text
    rho_i = m_i / V_i
    u_i   = U_i / m_i
    T,p   = PropsSI(..., "D", rho_i, "U", u_i, fluid)
    ```
    """

    mass = _as_1d("mass", mass)
    internal_energy_J = _as_1d("internal_energy_J", internal_energy_J)
    volume = _as_1d("volume", volume)
    _check_same_len(mass, internal_energy_J, volume)
    if np.any(mass <= 0.0) or np.any(volume <= 0.0):
        raise ValueError("mass and volume must be strictly positive")

    rho = mass / volume
    u = internal_energy_J / mass
    T = np.empty_like(rho)
    p = np.empty_like(rho)
    h = np.empty_like(rho)
    cp = np.empty_like(rho)
    mu = np.empty_like(rho)
    k = np.empty_like(rho)

    for i, (rho_i, u_i) in enumerate(zip(rho, u)):
        T[i] = PropsSI("T", "D", float(rho_i), "U", float(u_i), fluid)
        p[i] = PropsSI("P", "D", float(rho_i), "U", float(u_i), fluid)
        h[i] = PropsSI("H", "D", float(rho_i), "U", float(u_i), fluid)
        cp[i] = _safe_props("C", rho_i, u_i, T[i], p[i], fluid, 5200.0)
        mu[i] = _safe_props("V", rho_i, u_i, T[i], p[i], fluid, 2.0e-5)
        k[i] = _safe_props("L", rho_i, u_i, T[i], p[i], fluid, 0.1)

    return CoolantThermodynamicState(
        pressure=p,
        temperature=T,
        density=rho,
        specific_internal_energy_J_kg=u,
        specific_enthalpy_J_kg=h,
        cp=cp,
        mu=mu,
        k=k,
    )


def initial_mass_energy_from_TP(temperature, pressure, volume, fluid: str) -> tuple[np.ndarray, np.ndarray]:
    """Build initial `mass` and `internal_energy_J` arrays from `T,p,V`."""

    T = _as_1d("temperature", temperature)
    p = _as_1d("pressure", pressure)
    volume = _as_1d("volume", volume)
    _check_same_len(T, p, volume)
    if np.any(volume <= 0.0) or np.any(p <= 0.0):
        raise ValueError("pressure and volume must be strictly positive")

    rho = np.array([PropsSI("D", "T", float(Ti), "P", float(pi), fluid) for Ti, pi in zip(T, p)])
    u = np.array([PropsSI("U", "T", float(Ti), "P", float(pi), fluid) for Ti, pi in zip(T, p)])
    mass = rho * volume
    return mass, mass * u


def internal_energy_from_temperature_mass(temperature, mass, volume, fluid: str) -> np.ndarray:
    """Build cell internal energy from temperature, mass, and volume."""

    T = _as_1d("temperature", temperature)
    mass = _as_1d("mass", mass)
    volume = _as_1d("volume", volume)
    _check_same_len(T, mass, volume)
    rho = mass / volume
    u = np.array([
        PropsSI("U", "T", float(Ti), "D", float(rho_i), fluid)
        for Ti, rho_i in zip(T, rho)
    ])
    return mass * u


def enforce_internal_energy_bounds(
    mass,
    internal_energy_J,
    volume,
    fluid: str,
    *,
    T_floor: float = 60.0,
    T_ceiling: float = 2500.0,
    clip: bool = True,
) -> np.ndarray:
    """Return internal energy checked against a CoolProp-valid temperature range."""

    mass = _as_1d("mass", mass)
    internal_energy_J = _as_1d("internal_energy_J", internal_energy_J)
    volume = _as_1d("volume", volume)
    _check_same_len(mass, internal_energy_J, volume)
    rho = mass / volume
    floor_u = np.array([
        PropsSI("U", "T", float(T_floor), "D", float(rho_i), fluid)
        for rho_i in rho
    ])
    ceiling_u = np.array([
        PropsSI("U", "T", float(T_ceiling), "D", float(rho_i), fluid)
        for rho_i in rho
    ])
    lower = mass * floor_u
    upper = mass * ceiling_u
    out_of_bounds = (internal_energy_J < lower) | (internal_energy_J > upper)
    if np.any(out_of_bounds) and not clip:
        raise FloatingPointError(
            "coolant internal energy left the configured CoolProp temperature "
            f"range ({T_floor:g}-{T_ceiling:g} K)"
        )
    return np.clip(internal_energy_J, lower, upper)


def enforce_internal_energy_floor(*args, **kwargs) -> np.ndarray:
    """Backward-compatible alias for `enforce_internal_energy_bounds`."""

    return enforce_internal_energy_bounds(*args, **kwargs)


def enforce_density_bounds(
    mass,
    internal_energy_J,
    volume,
    *,
    rho_min: float = 1.0e-3,
    rho_max: float = 100.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Clip cell mass to a practical gas-density range, preserving specific energy."""

    mass = _as_1d("mass", mass)
    internal_energy_J = _as_1d("internal_energy_J", internal_energy_J)
    volume = _as_1d("volume", volume)
    _check_same_len(mass, internal_energy_J, volume)
    u = internal_energy_J / mass
    mass_new = np.clip(mass, float(rho_min) * volume, float(rho_max) * volume)
    return mass_new, mass_new * u


def quasi_steady_face_mdot(
    pressure,
    density,
    resistance,
    *,
    inlet_pressure: float | None = None,
    outlet_pressure: float | None = None,
    inlet_resistance: float | None = None,
    outlet_resistance: float | None = None,
    inlet_enabled: bool = True,
    outlet_enabled: bool = True,
) -> np.ndarray:
    """Compute face mass flows from pressure differences.

    Returns an array of length `n_cells + 1`. Positive flow is in increasing
    cell-index direction: face 0 enters cell 0 from the left boundary, and face
    `n` exits the last cell to the right boundary.

    The hydraulic closure is:

    ```text
    mdot = sign(dp) * sqrt(rho_face * |dp| / R)
    ```

    so `R` has units `Pa / (kg/s)^2`. This is deliberately simple: adapters can
    build `R` from Darcy, minor-loss, valve, or calibrated restrictions.
    """

    p = _as_1d("pressure", pressure)
    rho = _as_1d("density", density)
    resistance = _as_1d("resistance", resistance)
    if p.size != rho.size:
        raise ValueError("pressure and density must have the same length")
    n = p.size
    if resistance.size != max(n - 1, 0):
        raise ValueError("resistance must have length n_cells - 1")
    if np.any(rho <= 0.0):
        raise ValueError("density must be strictly positive")
    if np.any(resistance <= 0.0):
        raise ValueError("resistance must be strictly positive")

    face = np.zeros(n + 1, dtype=float)
    if inlet_pressure is not None and inlet_enabled:
        R = _positive_resistance("inlet_resistance", inlet_resistance)
        face[0] = _orifice_mdot(float(inlet_pressure) - p[0], rho[0], R)

    for i in range(n - 1):
        rho_face = 0.5 * (rho[i] + rho[i + 1])
        face[i + 1] = _orifice_mdot(p[i] - p[i + 1], rho_face, resistance[i])

    if outlet_pressure is not None and outlet_enabled:
        R = _positive_resistance("outlet_resistance", outlet_resistance)
        face[-1] = _orifice_mdot(p[-1] - float(outlet_pressure), rho[-1], R)

    return face


def conservative_mass_energy_step(
    mass,
    internal_energy_J,
    specific_enthalpy_J_kg,
    face_mdot,
    heat_W,
    dt: float,
    *,
    inlet_enthalpy_J_kg: float | None = None,
    outlet_backflow_enthalpy_J_kg: float | None = None,
    mass_floor: float = 1.0e-12,
) -> CompressibleCoolantStepResult:
    """Advance cell mass and internal energy with supplied face mass flows.

    `face_mdot` has length `n_cells + 1`; positive direction is from lower to
    higher cell index. Boundary backflow is supported if the corresponding
    boundary enthalpy is supplied.
    """

    m_old = _as_1d("mass", mass)
    U_old = _as_1d("internal_energy_J", internal_energy_J)
    h = _as_1d("specific_enthalpy_J_kg", specific_enthalpy_J_kg)
    face_mdot = _as_1d("face_mdot", face_mdot)
    heat_W = _as_1d("heat_W", heat_W)
    _check_same_len(m_old, U_old, h, heat_W)
    n = m_old.size
    if face_mdot.size != n + 1:
        raise ValueError("face_mdot must have length n_cells + 1")
    if np.any(m_old <= 0.0):
        raise ValueError("mass must be strictly positive")
    if dt < 0.0:
        raise ValueError("dt must be non-negative")

    face_h = _upwind_face_enthalpy(
        h,
        face_mdot,
        inlet_enthalpy_J_kg=inlet_enthalpy_J_kg,
        outlet_backflow_enthalpy_J_kg=outlet_backflow_enthalpy_J_kg,
    )
    m_new = m_old + dt * (face_mdot[:-1] - face_mdot[1:])
    if np.any(m_new <= mass_floor):
        raise FloatingPointError("coolant cell mass became non-positive")

    enthalpy_flux = face_mdot * face_h
    U_new = U_old + dt * (enthalpy_flux[:-1] - enthalpy_flux[1:] + heat_W)
    u_new = U_new / m_new

    mass_residual = float(np.sum(m_new - m_old) - dt * (face_mdot[0] - face_mdot[-1]))
    adv_in = _positive_boundary_energy(face_mdot[0], face_h[0], dt, incoming_left=True)
    adv_in += _positive_boundary_energy(face_mdot[-1], face_h[-1], dt, incoming_left=False)
    adv_out = _outgoing_boundary_energy(face_mdot[0], face_h[0], dt, left=True)
    adv_out += _outgoing_boundary_energy(face_mdot[-1], face_h[-1], dt, left=False)
    heat = float(np.sum(heat_W) * dt)
    energy_residual = float(
        np.sum(U_new - U_old) - (heat + adv_in - adv_out)
    )

    return CompressibleCoolantStepResult(
        mass_new=m_new,
        internal_energy_new_J=U_new,
        specific_internal_energy_new_J_kg=u_new,
        face_mdot=face_mdot.copy(),
        face_enthalpy=face_h,
        heat_added_J=heat,
        advective_energy_in_J=float(adv_in),
        advective_energy_out_J=float(adv_out),
        mass_residual_kg=mass_residual,
        energy_residual_J=energy_residual,
    )


def _cfl_stable_substep_count(
    mass: np.ndarray,
    face_mdot: np.ndarray,
    dt: float,
    *,
    safety: float = 0.25,
) -> int:
    """Number of equal substeps needed to keep the explicit coolant mass/energy
    advection (`conservative_mass_energy_step`, forward Euler in the conserved
    variables) within its stability limit for this macro step.

    This scheme is unconditionally unstable once a step advances a cell by more
    than roughly its residence time `mass / mdot` -- confirmed empirically
    2026-08-18 on the shell-and-tube bang-bang validation case: stable at
    dt/tau <= ~0.2, a fast-growing single-cell spike at dt/tau ~0.4, and a
    `FloatingPointError` from `enforce_internal_energy_bounds` within a handful
    of steps at dt/tau > ~1. Confirmed 2026-08-19 to be the identical root
    cause of the helical transient core's crash too (`main_solve_transient.py`'s
    own inline mass/energy loop -- separate code, same kernel). Originally
    written in `transient_core/adapters_shelltube.py`; consolidated here
    2026-08-19 so both solvers import one copy instead of duplicating it.

    Uses the worst case across cells (max |face flow|, min cell mass) rather
    than a mean, since a single under-resolved cell is enough to trip the
    energy-bounds guard.
    """
    if dt <= 0.0:
        return 1
    m = np.asarray(mass, dtype=float)
    if m.size == 0:
        return 1
    mdot_max = float(np.max(np.abs(np.asarray(face_mdot, dtype=float))))
    if mdot_max <= 0.0:
        return 1
    tau_min = float(np.min(m)) / mdot_max
    if tau_min <= 0.0:
        return 1
    return max(int(np.ceil(dt / (safety * tau_min))), 1)


def advance_flowpath_coolant(
    path,
    fluid: str,
    mass,
    internal_energy_J,
    face_mdot,
    heat_W,
    dt: float,
    *,
    inlet_enthalpy_J_kg: float | None = None,
    outlet_backflow_enthalpy_J_kg: float | None = None,
    mass_floor: float = 1.0e-12,
    safety: float = 0.25,
) -> CompressibleCoolantStepResult:
    """CFL-safe subcycled `conservative_mass_energy_step` for one
    `core/mesh.py::FlowPath`.

    Thin convenience wrapper, not new physics -- exactly the substep-and-sum
    pattern already applied by hand in `transient_core/adapters_shelltube.py`
    (shell-and-tube, 2026-08-18) and `main_solve_transient.py` (helical,
    2026-08-19), generalized so a future `core/residual.py` caller (Stage D3)
    gets it for free instead of reimplementing the loop a third time.

    `face_mdot`/`heat_W` stay frozen across substeps (the same
    quasi-steady-per-macro-step assumption the legacy adapters already use
    for the hot side); only the coolant `(m, U)` state and its reconstructed
    enthalpy are updated each substep, via this module's own
    `coolprop_state_from_mass_energy` -- `fluid` is passed explicitly, same
    convention as every other function in this module (no `ThermoBackend`
    routing here; that is `core/thermo.py`'s concern for a DIFFERENT state
    pair, `(p, h)`, not this module's `(m, U)` conservative pair).

    Uses `path.volume_total` (NOT `volume_per_channel` -- see `core/mesh.py`'s
    factor-of-N convention) as the per-cell volume for mass/energy accounting.
    Extensive diagnostics (`heat_added_J`, advective energies, residuals) are
    summed over substeps so the caller sees one result for the whole macro
    step; intensive/state fields reflect the final substep.
    """
    volume = path.volume_total
    n_sub = _cfl_stable_substep_count(mass, face_mdot, dt, safety=safety)
    sub_dt = dt / n_sub

    m_cur = np.asarray(mass, dtype=float)
    U_cur = np.asarray(internal_energy_J, dtype=float)
    h_cur = np.asarray(
        coolprop_state_from_mass_energy(m_cur, U_cur, volume, fluid).specific_enthalpy_J_kg,
        dtype=float,
    )

    heat_acc = 0.0
    adv_in_acc = 0.0
    adv_out_acc = 0.0
    mass_residual_acc = 0.0
    energy_residual_acc = 0.0
    step = None
    for sub in range(n_sub):
        step = conservative_mass_energy_step(
            m_cur,
            U_cur,
            h_cur,
            face_mdot,
            heat_W,
            sub_dt,
            inlet_enthalpy_J_kg=inlet_enthalpy_J_kg,
            outlet_backflow_enthalpy_J_kg=outlet_backflow_enthalpy_J_kg,
            mass_floor=mass_floor,
        )
        m_cur = step.mass_new
        U_cur = step.internal_energy_new_J
        heat_acc += step.heat_added_J
        adv_in_acc += step.advective_energy_in_J
        adv_out_acc += step.advective_energy_out_J
        mass_residual_acc += step.mass_residual_kg
        energy_residual_acc += step.energy_residual_J
        if sub < n_sub - 1:
            h_cur = np.asarray(
                coolprop_state_from_mass_energy(m_cur, U_cur, volume, fluid).specific_enthalpy_J_kg,
                dtype=float,
            )

    return CompressibleCoolantStepResult(
        mass_new=step.mass_new,
        internal_energy_new_J=step.internal_energy_new_J,
        specific_internal_energy_new_J_kg=step.specific_internal_energy_new_J_kg,
        face_mdot=step.face_mdot,
        face_enthalpy=step.face_enthalpy,
        heat_added_J=heat_acc,
        advective_energy_in_J=adv_in_acc,
        advective_energy_out_J=adv_out_acc,
        mass_residual_kg=mass_residual_acc,
        energy_residual_J=energy_residual_acc,
    )


def _upwind_face_enthalpy(
    cell_h: np.ndarray,
    face_mdot: np.ndarray,
    *,
    inlet_enthalpy_J_kg: float | None,
    outlet_backflow_enthalpy_J_kg: float | None,
) -> np.ndarray:
    n = cell_h.size
    face_h = np.empty(n + 1, dtype=float)

    if face_mdot[0] > 0.0:
        if inlet_enthalpy_J_kg is None:
            raise ValueError("inlet_enthalpy_J_kg is required for positive inlet flow")
        face_h[0] = float(inlet_enthalpy_J_kg)
    else:
        face_h[0] = cell_h[0]

    for i in range(n - 1):
        face_h[i + 1] = cell_h[i] if face_mdot[i + 1] >= 0.0 else cell_h[i + 1]

    if face_mdot[-1] >= 0.0:
        face_h[-1] = cell_h[-1]
    else:
        if outlet_backflow_enthalpy_J_kg is None:
            raise ValueError("outlet_backflow_enthalpy_J_kg is required for outlet backflow")
        face_h[-1] = float(outlet_backflow_enthalpy_J_kg)
    return face_h


def _orifice_mdot(dp: float, rho: float, resistance: float) -> float:
    if dp == 0.0:
        return 0.0
    return float(np.sign(dp) * np.sqrt(max(rho * abs(dp) / resistance, 0.0)))


def _safe_props(output: str, rho: float, u: float, T: float, p: float, fluid: str, fallback: float) -> float:
    try:
        value = PropsSI(output, "D", float(rho), "U", float(u), fluid)
    except Exception:
        try:
            value = PropsSI(output, "T", float(T), "P", max(float(p), 1.0e3), fluid)
        except Exception:
            value = fallback
    if not np.isfinite(value) or value <= 0.0:
        return float(fallback)
    return float(value)


def _positive_resistance(name: str, value: float | None) -> float:
    if value is None or float(value) <= 0.0:
        raise ValueError(f"{name} must be positive")
    return float(value)


def _positive_boundary_energy(mdot: float, h: float, dt: float, *, incoming_left: bool) -> float:
    if incoming_left:
        return float(max(mdot, 0.0) * h * dt)
    return float(max(-mdot, 0.0) * h * dt)


def _outgoing_boundary_energy(mdot: float, h: float, dt: float, *, left: bool) -> float:
    if left:
        return float(max(-mdot, 0.0) * h * dt)
    return float(max(mdot, 0.0) * h * dt)


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

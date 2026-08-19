"""Hot-gas march — Stage D, Slice 3 (docs/solver_design/FV_CORE_REWORK_PLAN.md
section 3.5).

`hot_gas_march` reproduces `transient_core/adapters_shelltube.py::
shelltube_hot_gas_march`'s sequential per-cell algorithm

```text
h_removed_{i+1} = h_removed_i + dq_hot_i * ds_i / mdot_tube
progress_{i+1}  = progress_i  + progress_source_i * ds_i / velocity_i
```

on the new `core/` abstractions instead of `ShellTubeCoreGeometry`:
`FlowPath` (`core/mesh.py`) for geometry, `core/closures.py`'s registered
tube-side Nu/friction closures (Stage D1) instead of the direct
`dispatch_nu_tube_straight`/`nu_corrugated_tube_vicente` calls, and
`core/wall.py::CylindricalWall.fluxes_at_Tbar` (Stage C) instead of
`OneDimensionalSteadyConduction_ShellnHelicalTube`. The gas-state provider
(`gas_state_at`) comes from `core/hotgas/combustor.py` (relocated this
slice) unchanged.

This is genuinely new code, not a relocation — the march itself is not
duplicated infrastructure anywhere, only the pieces it composes are already
built. Equivalence to the legacy march is proven by
`tests/test_core_residual.py` on a hand-built fixture: same `T_gas`,
`h_gas_W_m2K`, `dq_hot`/`dq_cold`, `T_wg`/`T_wc`, and `h_removed`/`progress`
trajectory, exactly (`==`) for both `inside_tube_choice` values, with
radiation and calibration coefficients exercised. One field is
`nusselt`, which this module DERIVES as `Nu = h*D/k` (the registered
closures from `core/closures.py` return `h` directly, not `Nu` — that's
the whole contract) rather than getting it from the underlying dispatch
call the way the legacy march does (`h = Nu*k/D`); the two directions carry
~1e-15 relative floating-point reassociation noise between them, so
`nusselt` alone is checked at `rtol=1e-12` rather than exact equality —
`h_gas_W_m2K`, the physically consequential quantity that actually drives
the wall flux, matches exactly.

**Scope note**: this module does NOT yet decide the two legacy quirks found
investigating the 2026-08-18/19 CFL fixes (`dq_cold` discarded in favor of
`G·(Tbar_new−Tc)`; `mdot_effective = mean(|faces|)` fed to per-cell
closures) — neither lives in the hot-gas march itself (both are in the full
step assembly, `shelltube_step_inputs`, not built here). `cold_heat_W` is
computed and returned by this module (mirroring the legacy march's own
`ShellTubeWallFlux.cold_heat_W`) for whoever assembles the full step to
decide what to do with, exactly as before.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .closures import tube_friction_closure, tube_htc_closure
from .hotgas.combustor import GasStateProvider, _coerce_gas_state
from .mesh import FlowPath
from .wall import CylindricalWall
from ..physics.liquid_flow.registry import ClosureContext


@dataclass(frozen=True)
class HotGasMarchResult:
    """Sequential representative-tube hot-gas march result."""

    T_gas: np.ndarray
    h_gas_W_m2K: np.ndarray
    gas_velocity_m_s: np.ndarray
    reynolds: np.ndarray
    prandtl: np.ndarray
    friction_factor: np.ndarray
    nusselt: np.ndarray
    dp_per_length_Pa_m: np.ndarray
    T_wg: np.ndarray
    T_wc: np.ndarray
    hot_heat_W: np.ndarray
    cold_heat_W: np.ndarray
    dq_hot_per_length_W_m: np.ndarray
    dq_cold_per_length_W_m: np.ndarray
    k_wall: np.ndarray
    enthalpy_removed_J_kg: np.ndarray
    progress_variable: np.ndarray
    T_gas_outlet: float
    enthalpy_removed_outlet_J_kg: float
    progress_outlet: float


def hot_gas_march(
    path: FlowPath,
    wall: CylindricalWall,
    *,
    Tbar_wall,
    T_coolant,
    h_shell,
    mdot_hot_total: float,
    gas_state_at: GasStateProvider,
    inside_tube_choice: str,
    corrCoeffs=None,
    roughness_m: float | None = None,
    corrugation_thickness_m: float = 0.0,
    corrugation_pitch_m: float = 1.0,
    progress_initial: float = 0.0,
    h_g_rad=None,
) -> HotGasMarchResult:
    """March hot gas through one representative tube against current wall/coolant.

    `path` is the hot-stream `FlowPath` (e.g. `core/geometry/shell_and_tube.py::
    build_shelltube_assembly(...).hot`); `roughness_m` defaults to
    `path.roughness` if not given explicitly.
    """
    n = path.n_cells
    Tbar = _cell_array("Tbar_wall", Tbar_wall, n)
    Tc = _cell_array("T_coolant", T_coolant, n)
    hs = _cell_array("h_shell", h_shell, n)
    hrad = np.zeros(n) if h_g_rad is None else _cell_array("h_g_rad", h_g_rad, n)
    if np.any(hs < 0.0) or np.any(hrad < 0.0):
        raise ValueError("heat-transfer coefficients must be non-negative")
    if mdot_hot_total <= 0.0:
        raise ValueError("mdot_hot_total must be positive for hot-gas marching")

    rough = float(path.roughness) if roughness_m is None else float(roughness_m)
    mass_flux = path.mass_flux(mdot_hot_total)
    mdot_tube = abs(float(mdot_hot_total)) / path.n_parallel
    s_centers = path.s_centers
    ds = path.ds
    Dh = _cell_array("path.Dh", path.Dh, n)

    T_g = np.empty(n)
    rho = np.empty(n)
    mu = np.empty(n)
    k = np.empty(n)
    cp = np.empty(n)
    progress_source = np.empty(n)
    h_g = np.empty(n)
    velocity = np.empty(n)
    reynolds = np.empty(n)
    prandtl = np.empty(n)
    friction = np.empty(n)
    nusselt = np.empty(n)
    dpdx = np.empty(n)
    dq_hot = np.empty(n)
    dq_cold = np.empty(n)
    T_wg = np.empty(n)
    T_wc = np.empty(n)
    k_wall = np.empty(n)
    h_removed_history = np.empty(n)
    progress_history = np.empty(n)

    htc_closure = tube_htc_closure(inside_tube_choice)
    friction_closure = tube_friction_closure(inside_tube_choice)

    h_removed = 0.0
    progress = float(progress_initial)

    for i in range(n):
        h_removed_history[i] = h_removed
        progress_history[i] = progress
        state = _coerce_gas_state(gas_state_at(h_removed, progress, i))
        T_g[i] = state.T
        rho[i] = state.rho
        mu[i] = state.mu
        k[i] = state.k
        cp[i] = state.cp
        progress_source[i] = state.progress_source

        velocity[i] = mass_flux[i] / state.rho
        reynolds[i] = state.rho * velocity[i] * Dh[i] / state.mu
        prandtl[i] = state.cp * state.mu / state.k

        ctx = ClosureContext(
            fluid="combustion_products",
            p_Pa=0.0,
            h_J_kg=0.0,
            T_bulk_K=state.T,
            rho_b=state.rho,
            mu_b=state.mu,
            k_b=state.k,
            cp_b=state.cp,
            Pr_b=prandtl[i],
            mass_flux_kg_m2_s=mass_flux[i],
            diameter_m=Dh[i],
            heat_flux_W_m2=0.0,
            wall_temp_K=float(Tbar[i]),
            corrCoeffs=corrCoeffs,
            extra={
                "x_m": float(s_centers[i]),
                "roughness_m": rough,
                "corrugation_thickness_m": corrugation_thickness_m,
                "corrugation_pitch_m": corrugation_pitch_m,
            },
        )
        friction[i] = friction_closure.callable(ctx)
        h_g[i] = htc_closure.callable(ctx)
        nusselt[i] = h_g[i] * Dh[i] / state.k
        dpdx[i] = friction[i] * state.rho * velocity[i] ** 2 / (2.0 * Dh[i])

        faces = wall.fluxes_at_Tbar(
            T_bar=float(Tbar[i]),
            T_g=float(T_g[i]),
            T_c=float(Tc[i]),
            h_g=float(h_g[i]),
            h_c=float(hs[i]),
            h_g_rad=float(hrad[i]),
        )
        dq_hot[i] = faces.q_hot_W_m
        dq_cold[i] = faces.q_cold_W_m
        T_wg[i] = faces.T_wg
        T_wc[i] = faces.T_wc
        k_wall[i] = faces.k_w

        h_removed += dq_hot[i] * ds[i] / mdot_tube
        if velocity[i] > 0.0:
            progress += progress_source[i] * ds[i] / velocity[i]

    outlet_state = _coerce_gas_state(gas_state_at(h_removed, progress, n))

    return HotGasMarchResult(
        T_gas=T_g,
        h_gas_W_m2K=h_g,
        gas_velocity_m_s=velocity,
        reynolds=reynolds,
        prandtl=prandtl,
        friction_factor=friction,
        nusselt=nusselt,
        dp_per_length_Pa_m=dpdx,
        T_wg=T_wg,
        T_wc=T_wc,
        hot_heat_W=dq_hot * path.n_parallel * ds,
        cold_heat_W=dq_cold * path.n_parallel * ds,
        dq_hot_per_length_W_m=dq_hot,
        dq_cold_per_length_W_m=dq_cold,
        k_wall=k_wall,
        enthalpy_removed_J_kg=h_removed_history,
        progress_variable=progress_history,
        T_gas_outlet=float(outlet_state.T),
        enthalpy_removed_outlet_J_kg=float(h_removed),
        progress_outlet=float(progress),
    )


def _cell_array(name: str, value, n_cells: int) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = np.full(n_cells, float(arr))
    if arr.shape != (n_cells,):
        raise ValueError(f"{name} must have shape ({n_cells},)")
    return arr

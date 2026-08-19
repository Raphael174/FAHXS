"""Combustion-chamber hot-gas state providers — Stage D, Slice 3
(docs/solver_design/FV_CORE_REWORK_PLAN.md).

Relocated, unchanged, from `transient_core/adapters_shelltube.py`
(`ShellTubeGasState`, `GasStateProvider`, `fpv_gas_state_provider`,
`equilibrium_gas_state_provider`, `oxygen_gas_state_provider`,
`_coerce_gas_state`) — confirmed geometry-independent before moving (no
shell-and-tube-specific input anywhere in these functions; they wrap the
FPV manifold, the equilibrium/frozen manifold, and CoolProp Oxygen sensible
cooling, all of which are already fluid/geometry-agnostic). Same "pure
infra, move as-is" pattern D2 used for `compressible_coolant.py`;
`transient_core/adapters_shelltube.py` now re-exports these as a shim,
proven identical-object by
`tests/test_core_hotgas_combustor.py::test_adapters_shelltube_shim_reexports_are_identical_objects`.

This is the first real content for `core/hotgas/combustor.py`, the file the
original design doc §3 sketch named for wrapping "chamber Cantera/FPV
provider" — alongside `core/hotgas/nozzle_gas.py` (Stage F groundwork,
already there) and `core/hotgas/prescribed.py` (still not started).

Despite the `ShellTube*` naming (inherited unchanged from the relocation —
renaming is a separate decision, not made here), nothing in this module is
shell-and-tube-specific: `core/residual.py`'s hot-gas march consumes these
providers through the same `gas_state_at(h_removed, progress, i)` interface
regardless of which `FlowPath`-described tube geometry the caller supplies.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from CoolProp.CoolProp import PropsSI


@dataclass(frozen=True)
class ShellTubeGasState:
    """Thermophysical state returned by a hot-gas provider."""

    T: float
    rho: float
    mu: float
    k: float
    cp: float
    progress_source: float = 0.0


GasStateProvider = Callable[[float, float, int], ShellTubeGasState | dict]


def fpv_gas_state_provider(fpv) -> tuple[GasStateProvider, float]:
    """Return a gas-state provider backed by an `FPVManifold`.

    The returned initial progress is `fpv.Yc_inlet()`, matching the maintained
    finite-rate transient convention.
    """

    progress_initial = float(fpv.Yc_inlet())

    def provider(h_removed: float, progress: float, _i: int) -> ShellTubeGasState:
        T, rho, mu, k, cp, _xH2O, _xCO2, omega = fpv.state(h_removed, progress)
        return ShellTubeGasState(T=T, rho=rho, mu=mu, k=k, cp=cp, progress_source=omega)

    return provider, progress_initial


def equilibrium_gas_state_provider(manifold) -> tuple[GasStateProvider, float]:
    """Return a gas-state provider backed by an equilibrium/frozen manifold."""

    def provider(h_removed: float, _progress: float, _i: int) -> ShellTubeGasState:
        T, rho, mu, k, cp, _xH2O, _xCO2 = manifold.at(h_removed)
        return ShellTubeGasState(T=T, rho=rho, mu=mu, k=k, cp=cp, progress_source=0.0)

    return provider, 0.0


def oxygen_gas_state_provider(
    *,
    T_inlet: float,
    pressure: float,
    fluid: str = "Oxygen",
    T_min: float = 95.0,
    T_max: float = 1200.0,
) -> tuple[GasStateProvider, float]:
    """Return a pre-ignition oxygen sensible-cooling gas-state provider.

    `h_removed` is interpreted as specific enthalpy removed from the incoming
    oxygen stream. Temperature is recovered from `(H, P)` through CoolProp and
    clipped to the requested bounds before transport properties are evaluated.
    """

    if T_inlet <= 0.0 or pressure <= 0.0:
        raise ValueError("T_inlet and pressure must be positive")
    if T_min <= 0.0 or T_max <= T_min:
        raise ValueError("temperature bounds are invalid")

    T0 = float(np.clip(T_inlet, T_min, T_max))
    h0 = float(PropsSI("H", "T", T0, "P", pressure, fluid))

    def provider(h_removed: float, _progress: float, _i: int) -> ShellTubeGasState:
        h = h0 - max(float(h_removed), 0.0)
        T = float(PropsSI("T", "H", h, "P", pressure, fluid))
        T = float(np.clip(T, T_min, T_max))
        return ShellTubeGasState(
            T=T,
            rho=float(PropsSI("D", "T", T, "P", pressure, fluid)),
            mu=float(PropsSI("V", "T", T, "P", pressure, fluid)),
            k=float(PropsSI("L", "T", T, "P", pressure, fluid)),
            cp=float(PropsSI("C", "T", T, "P", pressure, fluid)),
            progress_source=0.0,
        )

    return provider, 0.0


def _coerce_gas_state(value) -> ShellTubeGasState:
    if isinstance(value, ShellTubeGasState):
        state = value
    elif isinstance(value, dict):
        state = ShellTubeGasState(
            T=value["T"],
            rho=value["rho"],
            mu=value["mu"],
            k=value["k"],
            cp=value["cp"],
            progress_source=value.get("progress_source", 0.0),
        )
    else:
        raise TypeError("gas_state_at must return ShellTubeGasState or dict")
    vals = np.array([state.T, state.rho, state.mu, state.k, state.cp], dtype=float)
    if not np.all(np.isfinite(vals)):
        raise ValueError("gas state contains non-finite values")
    if state.T <= 0.0 or state.rho <= 0.0 or state.mu <= 0.0:
        raise ValueError("gas T, rho, and mu must be positive")
    if state.k <= 0.0 or state.cp <= 0.0:
        raise ValueError("gas k and cp must be positive")
    return state

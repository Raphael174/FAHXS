"""Conservative state containers for the FV core — Stage D, Slice 2
(docs/solver_design/FV_CORE_REWORK_PLAN.md).

Generalizes `transient_core/state.py::TransientState`/`TransientStateLayout`
(read in full before writing this — it only packs `[Tbar_wall, T_coolant]`,
the OLD temperature-only coolant model) to the conservative `(mass,
internal_energy)` pair the project actually runs today
(`core/coolant.py`/`transient_core/compressible_coolant.py`'s "mass_energy"
coolant state model — the one both the shell-and-tube and helical CFL fixes
targeted this session).

`mass_kg`/`internal_energy_J` are TOTALS per cell (n_parallel channels
included), matching `core/mesh.py::FlowPath.volume_total`'s convention, NOT
`volume_per_channel` — see that module's `LLM_CONTEXT.md` for the
factor-of-N trap this convention exists to prevent.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class WallState:
    """Per-cell lumped wall temperature. Single layer only, matching
    `core/wall.py::CylindricalWall`'s current scope (no multi-layer/coating
    support yet)."""

    Tbar_K: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "Tbar_K", _as_1d("Tbar_K", self.Tbar_K))


@dataclass(frozen=True)
class CoolantState:
    """Per-cell conservative coolant state: total mass and total internal
    energy (both across `n_parallel` channels — see module docstring)."""

    mass_kg: np.ndarray
    internal_energy_J: np.ndarray

    def __post_init__(self) -> None:
        mass = _as_1d("mass_kg", self.mass_kg)
        energy = _as_1d("internal_energy_J", self.internal_energy_J)
        if mass.size != energy.size:
            raise ValueError("mass_kg and internal_energy_J must have equal length")
        if np.any(mass <= 0.0):
            raise ValueError("mass_kg must be strictly positive")
        object.__setattr__(self, "mass_kg", mass)
        object.__setattr__(self, "internal_energy_J", energy)

    @property
    def n_cells(self) -> int:
        return self.mass_kg.size

    @property
    def specific_internal_energy_J_kg(self) -> np.ndarray:
        return self.internal_energy_J / self.mass_kg


@dataclass(frozen=True)
class WallCoolantStateLayout:
    """Pack and unpack `[Tbar_K, mass_kg, internal_energy_J]` state vectors
    for one wall + one coolant stream sharing `n_cells` cells.

    Mirrors `transient_core/state.py::TransientStateLayout`'s pack/unpack
    shape, generalized from `[Tbar, T_coolant]` (2 fields) to
    `[Tbar, mass, U]` (3 fields) for the conservative coolant model.
    """

    n_cells: int

    def __post_init__(self) -> None:
        if int(self.n_cells) <= 0:
            raise ValueError("n_cells must be positive")
        object.__setattr__(self, "n_cells", int(self.n_cells))

    @property
    def wall_slice(self) -> slice:
        return slice(0, self.n_cells)

    @property
    def mass_slice(self) -> slice:
        return slice(self.n_cells, 2 * self.n_cells)

    @property
    def energy_slice(self) -> slice:
        return slice(2 * self.n_cells, 3 * self.n_cells)

    @property
    def size(self) -> int:
        return 3 * self.n_cells

    def pack(self, wall: WallState, coolant: CoolantState) -> np.ndarray:
        if wall.Tbar_K.size != self.n_cells:
            raise ValueError(f"wall.Tbar_K must have shape ({self.n_cells},)")
        if coolant.n_cells != self.n_cells:
            raise ValueError(f"coolant state must have shape ({self.n_cells},)")
        y = np.empty(self.size, dtype=float)
        y[self.wall_slice] = wall.Tbar_K
        y[self.mass_slice] = coolant.mass_kg
        y[self.energy_slice] = coolant.internal_energy_J
        return y

    def unpack(self, y: np.ndarray) -> tuple[WallState, CoolantState]:
        arr = np.asarray(y, dtype=float)
        if arr.ndim != 1 or arr.size != self.size:
            raise ValueError(f"state vector must have shape ({self.size},)")
        wall = WallState(Tbar_K=arr[self.wall_slice])
        coolant = CoolantState(
            mass_kg=arr[self.mass_slice],
            internal_energy_J=arr[self.energy_slice],
        )
        return wall, coolant


def _as_1d(name: str, value) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr

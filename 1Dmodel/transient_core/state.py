"""State packing helpers for transient wall + coolant models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TransientState:
    """Named views of the transient state vector."""

    Tbar_wall: np.ndarray
    T_coolant: np.ndarray


@dataclass(frozen=True)
class TransientStateLayout:
    """Pack and unpack ``[Tbar_wall, T_coolant]`` state vectors."""

    n_cells: int

    def __post_init__(self) -> None:
        if int(self.n_cells) <= 0:
            raise ValueError("n_cells must be positive")
        object.__setattr__(self, "n_cells", int(self.n_cells))

    @property
    def wall_slice(self) -> slice:
        return slice(0, self.n_cells)

    @property
    def coolant_slice(self) -> slice:
        return slice(self.n_cells, 2 * self.n_cells)

    @property
    def size(self) -> int:
        return 2 * self.n_cells

    def pack(self, Tbar_wall: np.ndarray, T_coolant: np.ndarray) -> np.ndarray:
        wall = self._as_cells("Tbar_wall", Tbar_wall)
        coolant = self._as_cells("T_coolant", T_coolant)

        y = np.empty(self.size, dtype=float)
        y[self.wall_slice] = wall
        y[self.coolant_slice] = coolant
        return y

    def unpack(self, y: np.ndarray) -> TransientState:
        state = np.asarray(y, dtype=float)
        if state.ndim != 1 or state.size != self.size:
            raise ValueError(f"state vector must have shape ({self.size},)")
        return TransientState(
            Tbar_wall=state[self.wall_slice],
            T_coolant=state[self.coolant_slice],
        )

    def _as_cells(self, name: str, value: np.ndarray) -> np.ndarray:
        arr = np.asarray(value, dtype=float)
        if arr.ndim != 1 or arr.size != self.n_cells:
            raise ValueError(f"{name} must have shape ({self.n_cells},)")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} contains non-finite values")
        return arr

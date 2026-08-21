"""Geometry-neutral axial grid descriptors for transient HX cells."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class AxialGrid:
    """Axial finite-volume grid and per-cell geometric measures.

    The transient core only needs arrays. Geometry-specific adapters are
    responsible for deriving these arrays from a helical coil or a shell-and-tube
    bundle.
    """

    x_edges: np.ndarray
    coolant_area: np.ndarray
    wall_area: np.ndarray
    hot_perimeter: np.ndarray
    coolant_perimeter: np.ndarray
    flow_direction: int = 1

    def __post_init__(self) -> None:
        x_edges = _as_1d("x_edges", self.x_edges)
        if x_edges.size < 2:
            raise ValueError("x_edges must contain at least two entries")
        if np.any(np.diff(x_edges) <= 0.0):
            raise ValueError("x_edges must be strictly increasing")

        n_cells = x_edges.size - 1
        coolant_area = _as_cells("coolant_area", self.coolant_area, n_cells)
        wall_area = _as_cells("wall_area", self.wall_area, n_cells)
        hot_perimeter = _as_cells("hot_perimeter", self.hot_perimeter, n_cells)
        coolant_perimeter = _as_cells("coolant_perimeter", self.coolant_perimeter, n_cells)

        if self.flow_direction not in (-1, 1):
            raise ValueError("flow_direction must be +1 or -1")
        for name, arr in (
            ("coolant_area", coolant_area),
            ("wall_area", wall_area),
            ("hot_perimeter", hot_perimeter),
            ("coolant_perimeter", coolant_perimeter),
        ):
            if np.any(arr <= 0.0):
                raise ValueError(f"{name} must be strictly positive")

        object.__setattr__(self, "x_edges", x_edges)
        object.__setattr__(self, "coolant_area", coolant_area)
        object.__setattr__(self, "wall_area", wall_area)
        object.__setattr__(self, "hot_perimeter", hot_perimeter)
        object.__setattr__(self, "coolant_perimeter", coolant_perimeter)

    @classmethod
    def uniform(
        cls,
        *,
        length: float,
        n_cells: int,
        coolant_area: float,
        wall_area: float,
        hot_perimeter: float,
        coolant_perimeter: float,
        flow_direction: int = 1,
    ) -> "AxialGrid":
        """Build a uniform grid from scalar geometric measures."""

        if int(n_cells) <= 0:
            raise ValueError("n_cells must be positive")
        if length <= 0.0:
            raise ValueError("length must be positive")
        n = int(n_cells)
        return cls(
            x_edges=np.linspace(0.0, float(length), n + 1),
            coolant_area=np.full(n, float(coolant_area)),
            wall_area=np.full(n, float(wall_area)),
            hot_perimeter=np.full(n, float(hot_perimeter)),
            coolant_perimeter=np.full(n, float(coolant_perimeter)),
            flow_direction=flow_direction,
        )

    @property
    def n_cells(self) -> int:
        return self.x_edges.size - 1

    @property
    def dx(self) -> np.ndarray:
        return np.diff(self.x_edges)

    @property
    def x_centers(self) -> np.ndarray:
        return 0.5 * (self.x_edges[:-1] + self.x_edges[1:])

    @property
    def length(self) -> float:
        return float(self.x_edges[-1] - self.x_edges[0])

    @property
    def coolant_volume(self) -> np.ndarray:
        return self.coolant_area * self.dx

    @property
    def wall_volume(self) -> np.ndarray:
        return self.wall_area * self.dx

    @property
    def inlet_index(self) -> int:
        return 0 if self.flow_direction == 1 else self.n_cells - 1

    @property
    def outlet_index(self) -> int:
        return self.n_cells - 1 if self.flow_direction == 1 else 0

    def heat_rate_from_linear(self, heat_per_length_W_m: np.ndarray) -> np.ndarray:
        """Convert per-length heat rates into per-cell heat rates."""

        qprime = _as_cells("heat_per_length_W_m", heat_per_length_W_m, self.n_cells)
        return qprime * self.dx


def _as_1d(name: str, value: np.ndarray) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _as_cells(name: str, value: np.ndarray, n_cells: int) -> np.ndarray:
    arr = _as_1d(name, value)
    if arr.size != n_cells:
        raise ValueError(f"{name} must have shape ({n_cells},)")
    return arr

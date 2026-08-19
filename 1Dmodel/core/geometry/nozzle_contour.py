"""Nozzle contour generation — Stage F groundwork (see
docs/solver_design/FV_CORE_REWORK_PLAN.md section 5.1).

Standalone: builds ``r(z)``/``A(z)`` from a handful of scalar geometry
parameters. Not yet wired into ``core/mesh.py``'s ``FlowPath`` (that
integration is deferred to when Stage F formally starts, after Stages B-E's
existing-config migration passes its gates) -- this module and
``core/hotgas/nozzle_gas.py`` are built and tested now because they don't
depend on the FV residual/mesh at all, only on scalar/array geometry and
gas-dynamics relations.

Only a conical contour is implemented (the user's "keep it simple for now"
scope, 2026-07-31). Rao/bell contours are a later addition to this same
module -- ``build_conical_contour`` is written so that only the divergent-
section point generation changes for a bell profile; the convergent section
and station bookkeeping are shape-agnostic.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class NozzleContour:
    """Axial station geometry for a nozzle hot-gas path.

    ``z=0`` is the throat; upstream (chamber/convergent side) is negative,
    downstream (divergent/exit side) is positive. This keeps the throat
    index fixed at a known location regardless of how much convergent-side
    geometry is resolved.
    """

    z_m: np.ndarray
    r_m: np.ndarray
    throat_index: int

    @property
    def A_m2(self) -> np.ndarray:
        return np.pi * self.r_m**2

    @property
    def D_m2(self) -> np.ndarray:
        """Local diameter [m] (property name kept short/consistent with A_m2)."""
        return 2.0 * self.r_m

    @property
    def A_t_m2(self) -> float:
        return float(self.A_m2[self.throat_index])

    @property
    def D_t_m(self) -> float:
        return float(self.D_m2[self.throat_index])

    @property
    def expansion_ratio(self) -> np.ndarray:
        return self.A_m2 / self.A_t_m2

    def __post_init__(self) -> None:
        z = np.asarray(self.z_m, dtype=float)
        r = np.asarray(self.r_m, dtype=float)
        if z.shape != r.shape or z.ndim != 1:
            raise ValueError("z_m and r_m must be equal-length 1D arrays")
        if np.any(np.diff(z) <= 0.0):
            raise ValueError("z_m must be strictly increasing")
        if np.any(r <= 0.0):
            raise ValueError("r_m must be strictly positive everywhere")
        i_min = int(np.argmin(r))
        if i_min != self.throat_index:
            raise ValueError(
                f"throat_index={self.throat_index} is not the minimum-radius "
                f"station (minimum is at index {i_min}) -- single-throat "
                f"sanity gate failed (design doc section 5.1)"
            )
        dr = np.diff(r)
        if np.any(dr[: self.throat_index] > 1e-12):
            raise ValueError("convergent section is not monotonically converging")
        if np.any(dr[self.throat_index :] < -1e-12):
            raise ValueError("divergent section is not monotonically diverging")


def build_conical_contour(
    *,
    D_throat_m: float,
    expansion_ratio: float,
    half_angle_div_deg: float = 15.0,
    contraction_ratio: float = 6.0,
    half_angle_conv_deg: float = 30.0,
    n_stations_conv: int = 20,
    n_stations_div: int = 60,
) -> NozzleContour:
    """Conical convergent/divergent contour from throat diameter + expansion ratio.

    Only ``D_throat_m`` and ``expansion_ratio`` were specified by the user for
    the first C2H4/O2 sizing pass (2026-07-31); everything else here is a
    documented, commonly-used DEFAULT ASSUMPTION, not a validated design
    choice:

    - ``half_angle_div_deg=15``: the standard baseline conical-nozzle
      divergence half-angle (roughly matches an 80% Rao-bell's effective
      angle, and is the most commonly cited "simple conical nozzle" default
      in rocket propulsion texts).
    - ``contraction_ratio=6`` (``A_chamber/A_throat``) and
      ``half_angle_conv_deg=30``: reasonable regen-engine defaults, NOT
      derived from any chamber sizing (L*, residence time, injector pattern)
      for this specific engine. The convergent section barely matters for
      peak wall heat flux (RPE: "the largest convective heat flux can be
      expected at the throat... in good agreement with experimental data
      which show the maximum heat flux to occur just before the throat") --
      it exists here mainly to give the area-Mach solver a well-posed
      subsonic starting point.

    Revisit all three once real chamber geometry (contraction ratio, chamber
    length/L*, convergent half-angle) is available.
    """
    if D_throat_m <= 0.0:
        raise ValueError("D_throat_m must be positive")
    if expansion_ratio <= 1.0:
        raise ValueError("expansion_ratio must exceed 1")
    if contraction_ratio <= 1.0:
        raise ValueError("contraction_ratio must exceed 1")

    R_t = 0.5 * D_throat_m
    R_c = R_t * np.sqrt(contraction_ratio)
    R_e = R_t * np.sqrt(expansion_ratio)

    z_conv_start = -(R_c - R_t) / np.tan(np.radians(half_angle_conv_deg))
    z_conv = np.linspace(z_conv_start, 0.0, n_stations_conv, endpoint=False)
    r_conv = R_c - (R_c - R_t) * (z_conv - z_conv_start) / (0.0 - z_conv_start)

    z_exit = (R_e - R_t) / np.tan(np.radians(half_angle_div_deg))
    z_div = np.linspace(0.0, z_exit, n_stations_div)
    r_div = R_t + (R_e - R_t) * (z_div - 0.0) / (z_exit - 0.0)

    z = np.concatenate([z_conv, z_div])
    r = np.concatenate([r_conv, r_div])
    throat_index = n_stations_conv  # first divergent station == throat station
    return NozzleContour(z_m=z, r_m=r, throat_index=throat_index)

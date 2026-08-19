"""Pluggable momentum closures for the FV core — Stage D, Slice 2
(docs/solver_design/FV_CORE_REWORK_PLAN.md section 3.6).

The design doc's §3.6 sketch (`MomentumModel` protocol, `dp/ds` per cell) is
written for the STEADY march driver (`core/drivers/march.py`), which does not
exist yet (Stage E). What this slice actually builds, and what is concretely
testable now, is the TRANSIENT quasi-steady closure both legacy transient
adapters already use: a prescribed uniform face mass flow while a coolant
pump/valve is actively commanding flow, falling back to pressure-driven
orifice discharge at the boundary once the commanded flow drops to zero
(valve closed). The full march-`dp/ds` `MomentumModel` protocol is deferred
to whichever stage actually builds the march driver — building it now with
no real consumer risks guessing the interface wrong.

`quasi_steady_face_mdot` here generalizes
`transient_core/adapters_shelltube.py::_shelltube_quasi_steady_faces` and
`main_solve_transient.py::_helical_quasi_steady_faces`, parameterized by
`core/mesh.py::FlowPath` instead of a raw `flow_direction` int. The
resistance *calibration* (Bell-Delaware-based for shell-and-tube,
coil-arc-length-based for helical) stays the caller's responsibility, same
as `core/closures.py` kept geometry-specific correlation physics in
delegating callables rather than the thin adapter.

**A real divergence between the two legacy functions was found while
generalizing them, and is deliberately NOT reproduced twice**: in the
"commanded flow at/below floor" (valve-closed) branch, the shell-and-tube
version only ever sets the downstream boundary face (interior faces stay
zero); the helical version ADDITIONALLY computes non-zero interior-face
flows from cell-to-cell pressure differences in that same branch — an
inconsistency between the two adapters, not a documented design choice (see
`docs/solver_design/FV_CORE_REWORK_PLAN.md`'s 2026-08-19 note). This module
matches the shell-and-tube behavior (simpler, and the one with full
regression coverage from this session's CFL work) and flags the helical
divergence as tech debt for whoever retires the legacy files (Stage E),
rather than silently picking one without saying so.
"""

from __future__ import annotations

import numpy as np

from .coolant import _orifice_mdot
from .mesh import FlowPath


def quasi_steady_face_mdot(
    path: FlowPath,
    pressure: np.ndarray,
    density: np.ndarray,
    resistance: np.ndarray,
    *,
    mdot_inlet: float,
    outlet_pressure: float,
    mdot_floor: float,
) -> np.ndarray:
    """Face mass flows [kg/s], length `path.n_cells + 1`.

    - `mdot_inlet > mdot_floor`: every face carries the commanded flow
      uniformly, signed by `path.flow_direction`.
    - `mdot_inlet <= mdot_floor`: only the downstream boundary face (the one
      `path.flow_direction` points away from) carries pressure-driven orifice
      discharge; interior faces are zero. Matches
      `_shelltube_quasi_steady_faces`'s low-flow branch exactly (NOT
      `_helical_quasi_steady_faces`'s — see module docstring).

    `resistance` must have length `path.n_cells` (matching both legacy
    functions' signature), though only its last element is read in the
    low-flow branch — a quirk of the legacy interface, preserved here for
    faithful reproduction rather than "fixed" as a drive-by.
    """
    p = np.asarray(pressure, dtype=float)
    rho = np.asarray(density, dtype=float)
    n = path.n_cells
    if p.size != n or rho.size != n:
        raise ValueError(f"pressure and density must have shape ({n},)")
    resistance = np.asarray(resistance, dtype=float)
    if resistance.size != n:
        raise ValueError(f"resistance must have shape ({n},)")

    face = np.zeros(n + 1, dtype=float)
    mdot_cmd = max(float(mdot_inlet), 0.0)
    if mdot_cmd > mdot_floor:
        return np.full(n + 1, mdot_cmd if path.flow_direction == 1 else -mdot_cmd, dtype=float)

    if path.flow_direction == 1:
        face[-1] = max(
            _orifice_mdot(p[-1] - float(outlet_pressure), rho[-1], resistance[-1]), 0.0
        )
    else:
        face[0] = min(
            _orifice_mdot(float(outlet_pressure) - p[0], rho[0], resistance[-1]), 0.0
        )
    return face

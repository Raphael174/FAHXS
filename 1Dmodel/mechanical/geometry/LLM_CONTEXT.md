# mechanical/geometry/ LLM Context

## Scope

Geometry helper functions consumed directly by the maintained legacy solvers
(`main_solve.py`, `main_solve_shellntube.py`) inside their `__init__`
geometry blocks. These are the ORIGINAL, still-live geometry derivations —
not to be confused with the newer `core/geometry/` builders (see cross-
reference below), which reproduce the same physics in a fluid-agnostic,
`FlowPath`/`HXAssembly`-shaped form for the in-progress FV core rework and
are not yet wired into any production solver.

## Contents
| File | Role |
|---|---|
| `helix_geometry.py` | `HelixGeometryRadiusCST`: parametric helical-coil centerline (x,y,z), arc-length integration via Simpson's rule, and curve-fit `s -> x`/`s -> theta` mapping functions for converting coil arc length to axial engine position. Also `compute_Dh_shell` (Salimpour correlation for shell-with-coil hydraulic diameter). |
| `shelltube_geometry.py` | Bell-Delaware shell-side geometry: `compute_bell_delaware_geometry(...)` (crossflow/window/leakage/bypass areas, tube-row counts) and `estimate_tube_count` (bundle circle-packing tube count), consumed by `physics/bell_delaware.py::bell_delaware_shell`. Formulas from Serth, "Process Heat Transfer" (2007) ch. 6. Has a `__main__` self-check block reproducing an EchTherm reference case. |

## Cross-reference: `core/geometry/` (new, parallel implementation)

`core/geometry/shell_and_helical_tube.py` and `core/geometry/shell_and_tube.py`
reproduce this folder's physics for the new fluid-agnostic FV core
(`docs/solver_design/FV_CORE_REWORK_PLAN.md` Stage B). Specifically:

- `core/geometry/shell_and_helical_tube.py::build_helical_assembly` calls
  `HelixGeometryRadiusCST` directly and uses its `func_s_to_x` mapping as
  `z_of_s_edges` — this is exactly the piece that fixes the axial-bookkeeping
  bug flagged in `/CLAUDE.md`'s "Known Sharp Edges" (`main_solve.py`'s legacy
  `_advance_state()` approximates coil arc-length-to-axial-position linearly,
  silently wrong for any `HX_config` other than `"shellnHelicalTube"`; the
  real coil is ~1378 arc-length nodes, not ~100).
- `core/geometry/shell_and_tube.py::build_shelltube_assembly` reproduces
  `main_solve_shellntube.py`'s geometry block plus
  `compute_bell_delaware_geometry` from this folder.

Both `core/geometry/` builders are **built and tested but not yet consumed by
any legacy solver** (per `core/geometry/LLM_CONTEXT.md`, Stage E work). This
folder (`mechanical/geometry/`) remains the live production path until that
wiring lands — do not assume the new builders have superseded these functions.

## Key correctness points

- `HelixGeometryRadiusCST` numerically integrates arc length with `scipy.
  integrate.simpson` inside a Python loop (O(N) Simpson calls, not vectorized)
  then curve-fits linear `s->x` and `s->theta` maps with `scipy.optimize.
  curve_fit` — the returned functions are fitted approximations of the true
  parametric relationship, not exact inversions. `L_max_pipe` (final `s[-1]`)
  is the authoritative total coil length used downstream.
- `compute_bell_delaware_geometry` requires baffle-cut angle terms
  (`arccos` arguments) to be clipped to `[-1, 1]` before use — already done
  internally (`arg = float(np.clip(arg, -1.0, 1.0))`) to avoid NaN from
  floating-point overshoot at the domain boundary.
- `estimate_tube_count` is described as "rough" in its own docstring (bundle
  circle-packing approximation, not exact tube-layout enumeration).

## TODOs

None found (no TODO/FIXME markers in either file).

## Change history

No usable git history for this folder (single initial commit; folder is
currently uncommitted/untracked per `git status`). `shelltube_geometry.py`'s
docstring labels itself "shell-and-tube extension, WP1", indicating it was
added later than `helix_geometry.py`, but no date is recorded in-file.

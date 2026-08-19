# core/geometry/ LLM Context

## Scope

Per-configuration geometry builders that turn `input_data.py` dataclasses
into `core/mesh.py`'s `FlowPath`/`HXAssembly` objects — the layer that
replaces the inline geometry derivations currently duplicated across
`main_solve.py.__init__`/`main_solve_shellntube.py.__init__`. See
`docs/solver_design/FV_CORE_REWORK_PLAN.md` §3.1.

## Contents

| File | Stage | Reproduces | Status |
|---|---|---|---|
| `shell_and_tube.py` | B | `main_solve_shellntube.py.__init__`'s geometry block + `mechanical/geometry/shelltube_geometry.py::compute_bell_delaware_geometry` | Built, tested (`tests/test_core_geometry_builders.py`), NOT yet consumed by any legacy solver |
| `shell_and_helical_tube.py` | B | `main_solve.py.__init__`'s geometry block (lines ~110-148 there) + `mechanical/geometry/helix_geometry.py::HelixGeometryRadiusCST` | Built, tested, NOT yet consumed — own docstring explicitly defers wiring to Stage E |
| `nozzle_contour.py` | F | Nothing — genuinely new capability, no legacy equivalent anywhere in the repo | Built early (2026-07-31), ahead of B-E, standalone (zero dependency on `mesh.py`) |
| `__init__.py` | — | — | Empty/minimal package marker |

## Key correctness points

- `build_shelltube_assembly(shell_tube_prop, combustor_prop, *, n_cells=200)`:
  hot=`tube_gas` (`n_parallel=N_tubes`, per-tube area/perimeter — the
  "divide by N_tubes" convention made structural instead of a per-call-site
  habit), cold=`shell_coolant` (`n_parallel=1`, Bell-Delaware crossflow
  area). Both streams share one uniform axial partition today, so their
  `StreamCoupling` is the identity — but still built through the same
  conservative machinery, so a future non-uniform tube-side grid needs no
  special-casing.
- `build_helical_assembly(combustor_prop, numerical_prop, *,
  n_cells_shell=None)`: hot=`shell` (annular gas passage, `s=z`, uniform
  grid), cold=`coil` (`n_parallel=N_coils`, **`s` is the real coil arc
  length, `z_of_s_edges` comes from `HelixGeometryRadiusCST`'s
  `func_s_to_x`** — this is THE field that fixes the axial-bookkeeping bug
  CLAUDE.md flags: `main_solve.py`'s legacy `_advance_state()` approximates
  this mapping linearly and is silently wrong for any `HX_config` other than
  `"shellnHelicalTube"`). Real coil geometry is ~1378 arc-length nodes for
  this combustor, not ~100 — a placeholder/toy geometry mistake here has
  bitten this project before (CLAUDE.md, 2026-07-13).
- `nozzle_contour.py`: conical contour builder (`theta_conv`, `theta_div`),
  parameterized by throat + expansion ratio. Sanity gates: monotonic
  converging section, single throat, `A/A_t >= 1` everywhere. Convergent
  section geometry defaults are ASSUMED, not this project's real chamber
  geometry — do not treat outputs as design numbers without confirming the
  contour against real geometry first (see `docs/solver_design/
  FV_CORE_REWORK_PLAN.md`'s 2026-07-31 nozzle notes and `validation/
  nozzle_c2h4_o2_bartz_example.py`'s own caveats).

## TODOs (from the design doc, not invented)

- Stage F still needs `nozzle_axial_channels.py` and
  `nozzle_helical_channels.py` (regen-cooling channel geometry) — neither
  exists yet.
- Stage E must actually wire `shell_and_tube.py`/`shell_and_helical_tube.py`
  into the legacy solvers, replacing their inline geometry blocks — not yet
  started.

## Change history

No deep git history available (single initial commit + this project's
ongoing uncommitted work). `nozzle_contour.py` dated 2026-07-31 in the design
doc; the two `shell_and_*.py` builders' exact build date isn't recorded
anywhere, only discovered-already-built 2026-08-18 when Stage D work began
(see `core/LLM_CONTEXT.md`).

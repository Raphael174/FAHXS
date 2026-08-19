# core/ LLM Context

## Scope

`hps_combustor.core` is the NEW fluid-agnostic quasi-1D finite-volume core
from `docs/solver_design/FV_CORE_REWORK_PLAN.md`. It is built up one staged
slice at a time, alongside (not yet replacing) the maintained legacy solvers
at the `1Dmodel/` top level (`main_solve.py`, `main_solve_shellntube.py`, and
their transient counterparts). Nothing here is wired into a production solver
path until its stage's acceptance gate passes — check the design doc's status
banner and dated notes at the top before assuming any piece below is "live".

**Read `docs/solver_design/FV_CORE_REWORK_PLAN.md` first.** It is the
executable design spec — file-by-file architecture (§3), what's reused vs.
new (§2), and the staged gate table (§6, now split into slices D1-D4 for
Stage D specifically — see its 2026-08-18 notes). This file is a map of
what's actually *built*, which lags the design doc's aspirational full
architecture.

## Contents (as of 2026-08-19)

| File | Stage | Status | Role |
|---|---|---|---|
| `thermo.py` | A | **Done, wired into all 4 legacy solvers** | `ThermoState`, `ThermoBackend` protocol, `RealFluidBackend` (CoolProp `(p,h)`), `IdealGasBackend` (legacy fast path, bit-identical to raw `CP.PropsSI`), `ReactingGasBackend` (FPV/equilibrium manifold adapter). `physics/liquid_flow/dispatch.py` re-exports `ThermoState`/`coolant_state_from_*` from here as the SAME objects (not copies). |
| `mesh.py` | B | **Built, tested, not yet wired into legacy solvers** | `FlowPath` (one stream's 1D discretization; per-cell measures are PER CHANNEL, not totals — see "Sharp edges" below), `StreamCoupling`/`build_coupling` (conservative overlap operator between two streams' axial partitions), `HXAssembly` (hot+cold `FlowPath` pair). |
| `wall.py` | C | **Built, tested, not yet wired into legacy solvers** | `CylindricalWall.fluxes_at_Tbar`/`.solve_steady` — generalizes `physics/heat_conduction.py`'s validated quadratic quasi-static radial reconstruction; proven `rel=1e-12` equivalent to it. Adds the NEW analytic rib/fin path (`rectangular_fin_efficiency`, `channel_effective_perimeter`) for ribbed/multi-channel geometries — currently dead code for the two existing configs, exercised only by its own unit tests, will matter starting at Stage F (nozzle regen channels). |
| `closures.py` | D1 | **Built, tested, not yet wired into any solver** | `tube_htc_closure`/`tube_friction_closure(inside_tube_choice)` — forced-name lookup into the registry entries `physics/liquid_flow/gas_closures.py` registers. Thin by design; see that module for the actual physics (which it delegates to, not reimplements). |
| `geometry/` | B/F | Mixed — see its own `LLM_CONTEXT.md` | Per-config `FlowPath`/`HXAssembly` builders (shell-and-tube, shell-and-helical-tube) + nozzle contour (Stage F groundwork, built ahead of schedule). |
| `hotgas/` | F | Groundwork only | Nozzle hot-gas quasi-1D area-Mach + Bartz-Cornelisse HTC. Greenfield — no legacy equivalent exists anywhere in the repo. |
| `state.py`, `momentum.py`, `residual.py`, `assembly.py`, `diagnostics.py`, `drivers/` | D2-D4, E, F, G | **Do not exist yet** | Per the design doc's §3 architecture sketch. Do not assume any of these exist without checking — an earlier version of the design doc's status banner was stale on exactly this point (said `mesh.py` "doesn't exist yet" after it had already been built), which is why this file exists: to keep the built-vs-designed line visible and current. |

## Major changes (chronological, from session work — no deep git history exists;
the repo has one initial commit plus this project's ongoing uncommitted work,
so this list is reconstructed from design-doc dates and memory, not `git log`)

- **2026-07-31**: `core/thermo.py` created (Stage A extraction). `core/
  geometry/nozzle_contour.py` and `core/hotgas/nozzle_gas.py` built early,
  standalone (Stage F groundwork, ahead of B-E, justified because they have
  zero dependency on `mesh.py`/`residual.py`).
- **Between 2026-07-31 and 2026-08-18 (exact date not recorded)**: `core/
  mesh.py`, `core/wall.py`, `core/geometry/shell_and_tube.py`, `core/
  geometry/shell_and_helical_tube.py` built — Stages B and C — with 46
  passing tests, but the design doc's status banner was never updated to
  say so. Discovered stale 2026-08-18 when starting Stage D work.
- **2026-08-18**: Stage A closed for real (all 4 legacy solvers rewired onto
  `IdealGasBackend`, replacing inline `CP.PropsSI` calls — see
  `main_solve.py`, `main_solve_shellntube.py`, `main_solve_transient.py`,
  `main_solve_shellntube_transient.py`). Stage D started; split into slices
  D1-D4 (design doc too large to build/review as one PR). D1 done: `core/
  closures.py` + `physics/liquid_flow/gas_closures.py` (new file, sibling
  package) unify shell-and-tube's tube-side gas Nu/friction correlations into
  the `physics/liquid_flow/registry.py` mechanism, delegating to the existing
  validated correlation functions rather than reimplementing them.

## TODOs / open items (from the design doc + this session, not invented)

- **D2**: `core/state.py` (conservative state vector) + `core/momentum.py`
  (protocol + `QuasiSteadyMomentum`) + a standalone coolant mass/energy
  kernel. Must decide whether to carry forward the CFL-subcycling pattern
  just added to the LEGACY `transient_core/adapters_shelltube.py` (see that
  package's `LLM_CONTEXT.md` — a real, fixed regression, 2026-08-18) or make
  coolant advection genuinely implicit, consistent with §3.5's "linearly
  implicit per cell" direction.
- **D3**: `core/residual.py` — assemble wall + coolant + closures + hot-gas
  march for one time step; gate is reproducing one legacy step's fluxes
  within tolerance on a fixed fixture.
- **D4**: `core/drivers/transient.py` + repoint
  `main_solve_shellntube_transient.py` at it, retiring
  `_run_shelltube_transient_core_mass_energy`.
- D3/D4 must explicitly decide whether to reproduce two legacy shell-and-tube
  quirks found while investigating the CFL regression, or fix them: (a)
  `dq_cold`/`T_wc` are computed by the wall solve and discarded — actual
  coolant heating uses `G·(Tbar_new − Tc)`, driven by the mean wall
  temperature, not the reconstructed cold face; (b) `mdot_effective` is a
  mean-of-faces scalar fed to per-cell closures rather than the local face
  flow. Default: reproduce-then-flag (the stage gate says "reproduce ...
  within tolerance", not redesign).
- `core/wall.py`'s docstring says the driver integrates
  `(rho*cp*A_wall) dT_bar/dt`, but neither `CylindricalWall` nor `FlowPath`
  currently stores `A_wall`/wall density/wall cp — `FlowPath` is NOT a strict
  superset of the legacy `transient_core.AxialGrid` it's meant to replace
  (that also stored `wall_area`/`wall_volume`). Needs a home, as a TOTAL
  (×`n_parallel`), before D3 can build the wall thermal-mass term.
- Bell-Delaware (shell-side correlation) is deliberately NOT wrapped as a
  registry closure — it returns `(h, dp)` from one call against a whole
  geometry dict, not a single scalar from a `ClosureContext`. Needs its own
  two-output closure protocol if/when it's registry-unified; stays a direct
  call (`physics/bell_delaware.py`) for now.
- Stages E (retire the 4 legacy `main_solve*.py` files), F (nozzle configs,
  partially started via `hotgas/`/`geometry/nozzle_contour.py`), and G
  (rename `physics/liquid_flow` → `physics/fluid`, delete `IdealGasBackend`
  if unused) are all still fully ahead, per the design doc §6.

## Sharp edges specific to this package

- **Factor-of-N convention flip vs. the legacy `transient_core.AxialGrid`
  this replaces**: `AxialGrid` stored per-cell areas/perimeters as TOTALS
  (already multiplied by `n_parallel`); `FlowPath` deliberately stores
  PER-CHANNEL values, with `n_parallel` carried separately and
  `FlowPath.mass_flux()` dividing by it exactly once. Critically,
  `AxialGrid.coolant_volume` ≡ `FlowPath.volume_total`, **not**
  `volume_per_channel`. Wiring `volume_per_channel` into a mass/energy
  inventory calculation silently makes it wrong by `N_ch` or `N_tubes` — this
  is exactly the class of bug a 2026-08-18 CFL investigation had to rule out
  before finding the real root cause elsewhere (see `transient_core/
  LLM_CONTEXT.md`).
- `FlowPath.z_of_s_edges` is the field that replaces the fragile
  `main_solve.py::_advance_state()` axial bookkeeping — the bug the
  `HX_config == "shellnHelicalTube"` guard in `main_solve.py` was added to
  catch (CLAUDE.md, 2026-07-13). Here it's DATA, computed once by a geometry
  builder, never re-derived mid-march.
- `core/wall.py`'s rib/fin factor-of-2 convention
  (`rectangular_fin_efficiency`) is deliberate for a rib wetted on BOTH
  faces (between two channels) — matches `physics/heat_conduction.py::
  compute_eta_fin_rectangular`, NOT `compute_fin_efficiency_ch` (factor 1,
  wetted on one side only). Getting this backwards under-predicts peak
  hot-gas-wall temperature on ribbed channels, the quantity that actually
  sizes a regen wall.
- `ClosureContext` (consumed by `closures.py` via the registry) gained two
  new optional fields 2026-08-18: `corrCoeffs` and `extra` (a dict escape
  hatch). Both default to `None`/`{}` and every pre-existing call site is
  unaffected — but any NEW closure needing calibration knobs or a
  closure-specific scalar (raw axial position, roughness, corrugation
  geometry) should use these rather than growing `ClosureContext` with
  narrow one-consumer fields. See `physics/liquid_flow/registry.py`'s
  `LLM_CONTEXT.md` (in `physics/liquid_flow/`) for the full registry
  contract.

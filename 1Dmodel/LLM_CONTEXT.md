# 1Dmodel/ LLM Context

## Scope

The top-level `hps_combustor` package (mapped from `1Dmodel/` in `pyproject.toml`
since `1Dmodel` is not a valid Python identifier). This file covers only the
files directly in this folder: the two user-facing entry points, the four
backend solver classes, the input dataclasses, result packaging, schedule
loading, and the dashboard CLI. Subpackages (`physics/`, `mechanical/`,
`core/`, `transient_core/`, `model_data_process/`, `simulink_coupling/`, …)
have their own context, if present, and are out of scope here.

## Contents

| File | Role |
|---|---|
| `__init__.py` | Empty; marks the package. |
| `input_data.py` | Every input dataclass: `coolantProp`, `hotgasProp`, `combustorProp`, `shellTubeProp`, `numericalProp`, `transientProp`, `runProp`, `system_requirements`, `CorrelationCoefficients` (plus `ToasterProp`, less used). Source of truth for field names used throughout the rest of the repo. |
| `main_steady.py` | User entry point for steady runs (`python -m hps_combustor.main_steady`). Builds inputs, dispatches on `combustorProp.HX_config`, calls `package_steady_run`. |
| `main_transient.py` | User entry point for transient runs (`python -m hps_combustor.main_transient`). Same dispatch shape, plus applies `runProp.schedule_file` via `schedule_inputs.py`. |
| `main_solve.py` | Backend steady solver, class `main_solver` — shell-and-**helical**-tube only (`combustorProp.HX_config == "shellnHelicalTube"`, enforced by a `ValueError` in `__init__`). Also hosts `solve_counterflow_liquid_reference` / `solve_counterflow_physical_reference`, the shooting-method counter-flow references. |
| `main_solve_shellntube.py` | Backend steady solver, class `shellntube_solver` — baffled shell-and-tube (WP2). Hot gas inside straight tubes, coolant on shell side via Bell-Delaware; both inlets prescribed, solved by predictive under-relaxed sweep iteration (a genuine two-point BVP, unlike the helical solver's prescribed-outlet shortcut). |
| `main_solve_transient.py` | Backend transient solver, class `transient_solver(main_solver)` — helical coil. `fluid_model="quasi_steady"` (legacy, lumped wall ODE only) vs `"transient_coolant"` (production, wall + 1D-FV helium mass/energy). Hot gas stays quasi-steady (sub-ms residence, Mach_g < 0.3 design). |
| `main_solve_shellntube_transient.py` | Backend transient solver, class `shellntube_transient_solver(shellntube_solver)`. `solve_transient_core` delegates entirely to `transient_core.adapters_shelltube.run_shelltube_transient_core` — none of this file's own raw call sites are reachable from that path (asymmetric with the helical file, where `solve_transient_core` still calls back into several of its own methods). |
| `result_package.py` | `package_steady_run` / `package_transient_run`: writes JSON metadata + input presets + numeric arrays (+ optional CSV) under `runProp.output_root`, optionally zipped. `load_transient_time_series` reads it back from a folder/`.npz`/`.zip`. Every user-facing run should go through this. |
| `schedule_inputs.py` | `apply_schedule_file`: loads transient BC schedules from CSV (dependency-free) or XLSX (needs pandas/openpyxl), two-sheet (`helium`/`propellants`) or single-sheet with descriptive column names, decimal-comma tolerant. |
| `dashboard_cli.py` | `hps-dashboard` console entry point: builds a self-contained HTML dashboard from a saved `.npz`/run folder/`.zip` via `TransientDashboard`. |

## Recent structural change (Stage A, closed 2026-08-18)

Per `docs/solver_design/FV_CORE_REWORK_PLAN.md`, all four `main_solve*.py`
files now route every gas-mode/GOX raw `CoolProp.CoolProp.PropsSI` call
through `core/thermo.py::IdealGasBackend` — each class instantiates its own
`self._thermo = IdealGasBackend()` in `__init__` (47 call sites converted 1:1,
no expression rewriting, so results are bit-identical before/after). Liquid/
two-phase/supercritical branches were already on the dispatch-routed path and
untouched. Treat `self._thermo` as the correct place to add any new ideal-gas
property lookup in these files — do not add a fresh inline `PropsSI` call.
`IdealGasBackend` exists specifically to make the eventual `core/` rework
(Stages B–E, in progress — see that plan doc) reproduce these files' numbers;
it is a migration scaffold, not meant to be the final architecture.

## Entry points and dispatch

- `main_steady.py` / `main_transient.py` are the ONLY files a normal user run
  should call directly (`python -m hps_combustor.main_steady` /
  `.main_transient`). Both dispatch purely on `combustorProp.HX_config`
  (`"shellnHelicalTube"` → `main_solve.py`/`main_solve_transient.py`,
  `"shellntube"` → `main_solve_shellntube.py`/`main_solve_shellntube_transient.py`).
- `main_steady.py`'s counter-flow branch for the helical config *always*
  shoots on the physical cold-end inlet (`T_in`/`p_in`) — it deliberately does
  not offer the legacy guessed-`T_out`/`p_out` plain march, because that guess
  cannot represent a two-phase state when `coolant_model="equilibrium_liquid"`
  and repeatedly produced silently-wrong ducts or hard crashes.
- Every `main_solve*.py` file keeps a `__main__`-bootstrap block at the top
  (registers the package under the stable alias `_hps` via `runpy`) so it can
  still be run directly (`python main_solve.py`) even though `1Dmodel` is not
  a valid package name on its own. Preserve this block if editing these files.

## `input_data.py` — key fields and gotchas

- `coolantProp.coolant` / `coolant_model`: fluid-agnostic by design — Helium +
  `"single_phase_coolprop"` is the working project baseline (confirmed with
  the user 2026-07-13); `coolant_model="equilibrium_liquid"` switches on the
  `(p,h)` liquid/boiling path, but see "Liquid coolant wiring status" below —
  it does NOT change results for shell-and-tube or either transient solver.
- `combustorProp.HX_config`: `"shellnHelicalTube"` | `"shellntube"` (also
  lists legacy `"coolingjacket"`/`"coolingcoil"` values not implemented by the
  maintained solvers). `main_solver` (`main_solve.py`) raises if this is not
  exactly `"shellnHelicalTube"` — see Known Sharp Edges below.
- `combustorProp.flow_config`: `"co"` | `"counter"`, supported by all four
  maintained solvers.
- `numericalProp.chemistry_model`: `"finite_rate"` (FPV manifold, the required
  production default for the diesel/O2 high-heat-extraction regime) |
  `"equilibrium"` | `"frozen"` (validation-only, used by the fast baseline
  regression fixtures for speed/determinism).
- `transientProp.fluid_model`: `"quasi_steady"` (legacy, wall-only ODE) |
  `"transient_coolant"` (production, wall + compressible coolant FV — routes
  through `solve_transient_core`, which itself hits the CFL fix noted below
  for shell-and-tube).
- `transientProp.solver_method`: `"fixed_step"` is the bounded-cost production
  path for long shell-and-tube counter-flow runs (linearly implicit in the
  local wall-film stiffness); `"BDF"`/`"RK45"` are validation options only —
  their Jacobian probing is expensive when the RHS includes profile
  relaxation.
- `CorrelationCoefficients`: the single source of calibration knobs consumed
  by `optimization/calibrate.py` and design maps elsewhere. **Do not rename
  fields casually** — 21 fields are load-bearing for calibration priors.
- `runProp`: `output_root` (default `"zip_folders"`), `make_archive`,
  `save_csv`, `schedule_file`, and the independent axial-node counts for
  shell-and-tube (`shelltube_steady_nodes=200`, `shelltube_transient_nodes=80`
  — separate from the helical arc-length grid, which is driven by
  `numericalProp.N_arc_steps_per_turn`).

## Liquid coolant wiring status (as of this writing)

- **Wired into the coupled march**: helical steady (`main_solve.py`), when
  `coolantProp.coolant_model == "equilibrium_liquid"`. Co-flow is fully
  self-consistent; counter-flow's prescribed-outlet limitation is resolved by
  `solve_counterflow_liquid_reference()` (adaptive bracket + bisection on the
  hot-end starting enthalpy against the physical `T_in`/`p_in`).
- **Postprocess-only**: shell-and-tube steady (`main_solve_shellntube.py`),
  via the opt-in `shellntube_solver.liquid_coolant_postprocess()` call — not
  in the coupled march itself.
- **Not present at all**: both transient solvers.
- Do not assume flipping `coolant_model` changes solved results anywhere
  other than the helical steady coupled march.

## Known Sharp Edges (apply to files in this folder)

- `main_solver.__init__` (helical) requires `combustorProp.HX_config ==
  "shellnHelicalTube"` and raises `ValueError` otherwise. Its
  `_advance_state()` axial-length bookkeeping silently used a wrong linear-`dx`
  approximation for any other `HX_config` before this guard existed — this bit
  every Phase 0-2 liquid-coolant test. The real helical coil is ~1378
  arc-length nodes (not ~100) for this combustor's geometry; duty is
  ~150-300 kW, not ~20 kW — use realistic scale in any new test/example.
- In shell-and-tube files, tube-side gas quantities are **per tube** — divide
  total hot-gas mass flow by `N_tubes` for per-tube velocity/enthalpy.
- Darcy friction convention is used throughout; do not apply a Fanning-to-
  Darcy factor.
- Shell-and-tube wall conduction must use `hot_side="inner"` (hot gas is
  inside the tubes) — do not regress the hot/cold perimeter mapping.
- Persistent Cantera objects (used inside the combustion setup these files
  call into) must be reset from cached inlet `T/p/Y` before repeated sweeps.
- `data_master` (from `model_data_process/data_processing.py::
  make_solver_data()`) must be fresh per solver instance — do not reuse across
  runs.
- Pre-ignition GOX chilldown (`ignition=0` with LOX/GOX mass flow scheduled)
  uses CoolProp Oxygen properties, a separate path from the He gas backend.
- `main_solve_transient.py`'s own inline mass/energy loop (separate code from
  `transient_core/adapters_shelltube.py`) has the same
  `FloatingPointError: coolant internal energy left the configured CoolProp
  temperature range` signature that was root-caused and fixed for shell-and-
  tube's `solve_transient_core` (an explicit-Euler CFL instability, fixed
  2026-08-18 via substepping in `adapters_shelltube.py`) — **the helical
  version was left alone, out of scope for that fix.** Flag before relying on
  `transient_solver.solve_transient_core` for real helical work.

## TODOs

None found via direct grep of `TODO`/`FIXME` in this folder's own files as of
this writing. Larger, evidenced gaps live in
`docs/solver_design/FV_CORE_REWORK_PLAN.md` (Stage D/E migration of these
four solver files onto `core/residual.py` + `core/drivers/`) rather than as
inline markers.

## Change history

Git history is effectively one initial commit plus a large uncommitted
working tree; do not infer a changelog from `git log`. Evidenced, dated
changes from code comments and `docs/solver_design/FV_CORE_REWORK_PLAN.md`:

- **2026-07-13**: `main_solver` gained the `HX_config == "shellnHelicalTube"`
  guard (previously a silent wrong-length bug for any other value).
- **2026-08-18**: Stage A of the FV core rework — all four `main_solve*.py`
  files rewired off inline `CP.PropsSI` onto `core/thermo.py::IdealGasBackend`
  (`self._thermo`), gated by bit-identical `tests/test_steady_baseline_regression*.py`
  and the newly added `tests/test_transient_baseline_regression.py`.

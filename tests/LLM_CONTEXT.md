# tests/ LLM Context

## Scope

Root-level pytest suite for `hps_combustor` (25 files as of this writing).
Covers steady/transient baseline regression, liquid/boiling coolant physics,
the new `core/` package (FV rework Stages A-D), transient infrastructure, and
the parked Simulink coupling feature. `research/flamelet_kit/tests/` is a
separate, standalone suite — not covered here.

## How to run

```powershell
python -m pip install -e .
python -m pytest tests/ -q --ignore=tests/test_simulink_coupling_stepper.py --ignore=tests/test_simulink_coupling_standalone.py --ignore=tests/test_simulink_coupling_fmu_wrapper.py
```

The three `test_simulink_coupling_*.py` files require `pythonfmu`
(`pip install pythonfmu` or `pip install -e ".[simulink]"`) and are
conventionally excluded from a normal full-suite run via `--ignore` — this is
established practice in this repo's own sessions, not just a suggestion.

## Known pre-existing failures

`test_steady_baseline_regression.py` has **4 pre-existing failures**, unrelated
to any single change in flight: frozen-chemistry mode, obtained values
consistently 0.5-10% higher than expected, not yet root-caused. Do not treat
these as something a given change broke, and do not describe them as fixed —
re-characterize before assuming a new regression is unrelated. See project
memory `baseline-regression-failures-characterized` for the exact
obtained/expected numbers.

## Index by area

### Steady baseline regression (guards against unintended gas-path drift)
| File | Covers |
|---|---|
| `test_steady_baseline_regression.py` | Helical + shell-and-tube, co/counter, `chemistry_model="frozen"` (fast, deterministic). Phase 0 of the liquid integration plan. **4 known pre-existing failures — see above.** |
| `test_steady_baseline_regression_finite_rate.py` | Same 4 cases under `chemistry_model="finite_rate"` (the real production chemistry default) so a finite-rate-only regression doesn't slip through the frozen-mode fixture. |

### Transient baseline regression / infrastructure
| File | Covers |
|---|---|
| `test_transient_baseline_regression.py` | Stage A gate for `main_solve_transient.py` / `main_solve_shellntube_transient.py` rewiring onto `core/thermo.py::IdealGasBackend` — the only pytest coverage that numerically exercises these two files end-to-end (vs. `test_main_transient_dispatch.py`, which mocks them out). |
| `test_transient_core_shelltube.py` | Regression/fix gate for the shell-and-tube `solve_transient_core` CFL instability found and fixed 2026-08-18 (explicit-Euler mass/energy advection unconditionally unstable past ~1 cell residence time per macro step; fixed via CFL-safe substepping in `adapters_shelltube.py`). |
| `test_transient_core_coolant_fv.py` | Unit tests for `transient_core`'s finite-volume building blocks: state layout round-trip, `AxialGrid` geometry/indices/heat conversion, geometry rejection, diagnostics (residence time, zero-flow). |
| `test_transient_progress.py` | `TransientProgressPrinter` terminal-progress output (field selection, enable/disable, custom pressure label). |
| `test_main_transient_dispatch.py` | `main_transient.py`'s dispatch logic only, with `FakeShellTubeTransientSolver`/`FakeHelicalTransientSolver` — proves the right solver class and `solve_transient`/`solve_transient_core` method get called for each `fluid_model`, without running real physics. |
| `test_schedule_inputs.py` | `apply_schedule_file` CSV/Excel parsing, including decimal-comma tolerance. |
| `test_pressurant_bangbang_sizing.py` | Pressurant/orifice sizing helpers (water exit orifice, series-branch CdA, helium orifice flow, staged bang-bang command logic) — supports transient valve-schedule scenarios. |
| `test_coupled_bangbang_hx.py` | Builds HX boundary inputs from a system-level pressurant/bang-bang history; momentum-model and geometry selection; hot-gas schedule tracking; end-to-end coupled-case run with a mocked HX. |

### Liquid / boiling coolant physics
| File | Covers |
|---|---|
| `test_liquid_boiling_poc.py` | Per-correlation unit/regression tests: HEM state, Gungor-Winterton, Müller-Steinhagen-Heck, Yu2002, Groeneveld LUT (incl. exact page-9 match), plus `solve_steady_straight_pipe` energy closure and the validation-matrix runner. |
| `test_coolant_models.py` | Unit tests of the coolant state/closure dispatcher (`physics/liquid_flow/dispatch.py`) — default legacy single-phase model, equilibrium-liquid `(p,h)` state, HTC/dp/CHF closure outputs, vapor-expansion edge case. |
| `test_liquid_coupled_helical.py` | The real coupled `main_solver` march with `coolant_model="equilibrium_liquid"` (helical only): subcooled/boiling co- and counter-flow, grid convergence, cross-check vs. the postprocess bridge. Requires `combustorProp.HX_config="shellnHelicalTube"` and realistic ~1378-node coil scale. |
| `test_liquid_counterflow_reference.py` | `solve_counterflow_liquid_reference()` shooting method (main_solve.py) — HEM-monotonicity check plus a `main_steady.py` dispatch-routing test. **Slow (~20 min)**: adaptive bracket-search-then-bisection needs up to ~40 full solver runs per case at real coil scale. |
| `test_liquid_hx_adapters.py` | Integration tests of the helical/shell-and-tube adapter layer (co/counter duty mapping, energy conservation, bad-input rejection) and `liquid_coolant_postprocess()` on both solvers, using `SimpleNamespace` fixtures rather than a full coupled solve. |
| `test_liquid_flow_shims.py` | Confirms pre-Phase-1 module paths (`physics/liquid_coolant.py` etc.) still re-export correctly and emit `DeprecationWarning` — the only place that deliberately imports the old paths. |

### `core/` package — FV rework Stages A-D (`docs/solver_design/FV_CORE_REWORK_PLAN.md`)
| File | Stage | Covers |
|---|---|---|
| `test_core_thermo.py` | A | `core/thermo.py` backends reproduce legacy CoolProp call patterns; `physics/liquid_flow/dispatch.py`'s re-exports are the SAME objects post-relocation (not reimplementations). |
| `test_core_mesh.py` | B | `FlowPath`/coupling-overlap invariants; the acceptance-gate tests `test_coupling_conserves_energy_helical_to_shell` and `test_coupling_conserves_energy_nonuniform_both_sides` prove exact energy conservation on a non-uniform helical<->shell mapping. |
| `test_core_geometry_builders.py` | B | `core/geometry/` builders reproduce `main_solve.py`/`main_solve_shellntube.py`'s own inline geometry derivations to machine precision (guards against silent divergence before Stage E repoints the legacy solvers). |
| `test_core_wall.py` | C | Quadratic radial reconstruction + `hot_side` orientation match `physics/heat_conduction.py`; independently verifies the new analytic rib/fin path (dead code for existing configs) via hand calculation and limiting cases. |
| `test_core_closures.py` | D1 | The 4 new registered gas-side closures (`physics/liquid_flow/gas_closures.py` via `core/closures.py`) are bit-identical (`==`, not `approx`) to the legacy `shelltube_tube_gas_film` direct calls, across laminar/transitional/turbulent regimes and both `inside_tube_choice` values; includes the caught `Re_lo`/`Re_hi` corrCoeffs-vs-function-default regression. |
| `test_nozzle_hotgas.py` | F groundwork | Area-Mach and adiabatic-wall-temperature relations against closed-form isentropic solutions (tight tolerance — exact analytic checks, not fitted-correlation regressions). |

### Simulink coupling (excluded from normal runs — needs `pythonfmu`)
| File | Covers |
|---|---|
| `test_simulink_coupling_stepper.py` | Fidelity: `ShellTubeTransientStepper.step()` reproduces `run_shelltube_transient_core()`'s full-run trajectory when fed the same time-varying boundary values one step at a time. |
| `test_simulink_coupling_standalone.py` | Copies `simulink_coupling/` to an isolated temp dir, hides this repo's editable `hps_combustor` install from a subprocess, confirms the stepper + FMI wrapper run on the vendored `_vendor/` copy alone. |
| `test_simulink_coupling_fmu_wrapper.py` | `test_fmu_slave_steps_without_error`: wiring check only (FMI variables vs. stepper), can pass via the real install even if the packaged FMU is broken. `test_packaged_fmu_runs_standalone`: builds+unzips the actual `.fmu` and runs it in a subprocess with the real install hidden — the real proof the FMU is self-contained. |

## TODOs

None found via grep of `TODO`/`FIXME` across `tests/*.py`. The one open,
evidenced item affecting this suite is the 4 pre-existing baseline-regression
failures above — not yet root-caused, tracked in project memory rather than
an inline test marker.

## Change history

Git history is one initial commit plus a large uncommitted working tree; not
a real changelog. Evidenced, dated events from test docstrings and
`docs/solver_design/FV_CORE_REWORK_PLAN.md`:

- **2026-08-18**: `test_transient_core_shelltube.py` added as the fix gate for
  the shell-and-tube transient CFL instability. `test_transient_baseline_regression.py`
  and the `core/` Stage-D closure tests (`test_core_closures.py`) also landed
  around this date per the FV core rework plan's Stage A/D notes.
- `test_core_mesh.py`, `test_core_geometry_builders.py`, `test_core_wall.py`
  (Stages B/C) were found already built and passing when Stage D work started
  (documented as a stale-doc correction in the rework plan, same date).

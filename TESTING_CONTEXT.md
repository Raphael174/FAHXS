# Testing And Verification Context

## Dependency Setup

Recommended from repository root:

```powershell
python -m pip install -e .
python -m pip install -r 1Dmodel/requirements.txt
```

Runtime dependencies include `numpy`, `scipy`, `matplotlib`, `CoolProp`, and `cantera`. If Cantera pip wheels fail on the platform, use conda:

```powershell
conda install -c cantera cantera
```

## Existing Automated Tests

The explicit pytest suite is currently for `research/flamelet_kit`:

```powershell
python -m pytest research/flamelet_kit/tests/ -q
```

These tests use Cantera `gri30.yaml` and check grid construction, mixture-fraction stoichiometry, scalar-dissipation profile shape, flamelet boundedness, cache behavior, and cooling-PFR energy balance.

## Steady Helium Baseline Regression (Phase 0 Of The Liquid Integration Plan)

```powershell
python -m pytest tests/test_steady_baseline_regression.py -q
```

Freezes helical + shell-and-tube helium scalars (co/counter,
`chemistry_model="frozen"` for speed/determinism) from before any liquid-
coolant coupling work touched `main_solve.py` / `main_solve_shellntube.py`.
Every phase of `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` must
keep this passing unchanged — it is the guard against unintended gas-path
drift. ~6 s, bit-identical across runs.

Every fluid/geometry field this fixture depends on is pinned explicitly
(`coolant`, `coolant_model`, `mass_flow_c`, `combustorProp.HX_config`, etc.)
rather than left to `input_data.py` dataclass defaults — those defaults are
shared, mutable project state that gets changed for day-to-day
experimentation. An earlier version of this fixture relied on defaults and
broke silently when they changed (see the integration plan's "Configuration-
Drift Incident" section, 2026-07-13) — don't reintroduce that fragility.

## Liquid Coolant / Boiling Tests

Root-level pytest coverage exists for the in-progress liquid-coolant physics,
which lives in `1Dmodel/physics/liquid_flow/` (`correlations.py`, `chf.py`,
`dispatch.py`, `governing_equations.py`, `hx_adapters.py`,
`sanity_checks.py` — the old top-level `physics/liquid_coolant.py` /
`coolant_models.py` / `heated_liquid_channel.py` / `liquid_hx_adapters.py`
paths are now deprecation shims, see
`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` Phase 1). It is wired
into the helical steady coupled march (Phase 2); shell-and-tube remains
postprocess-only (Phase 3, not yet done).

```powershell
python -m pytest tests/test_liquid_boiling_poc.py tests/test_liquid_hx_adapters.py tests/test_coolant_models.py tests/test_liquid_flow_shims.py tests/test_liquid_coupled_helical.py -q
```

- `tests/test_liquid_boiling_poc.py`: unit/regression tests for each
  correlation (HEM state, Gungor-Winterton, Müller-Steinhagen-Heck, Yu2002
  multiplier/HTC, Groeneveld LUT interpolation and exact page-9 match), plus
  `solve_steady_straight_pipe` energy closure and the
  `liquid_validation_matrix` runner.
- `tests/test_liquid_coupled_helical.py`: exercises the real coupled
  `main_solver` march with `coolant_model="equilibrium_liquid"` — subcooled
  co/counter, boiling reached via heating (co-flow), a tractable
  counter-flow case, grid convergence (50 vs 100 arc-steps/turn), and a
  cross-check against the pre-existing postprocess bridge fed the identical
  converged duty profile (co-flow only — see the module docstring for why
  counter-flow is excluded from that specific check). All cases use
  `combustorProp.HX_config="shellnHelicalTube"` explicitly and
  `T_in`/`p_in`/`mass_flow_c` scaled for the real ~1378-node coil (see the
  module docstring's scale note).
- `tests/test_liquid_hx_adapters.py`: integration tests of the helical and
  shell-and-tube adapter layer (co/counter duty mapping, energy conservation,
  bad-input rejection) and of `liquid_coolant_postprocess()` on both solvers
  using `SimpleNamespace` fixtures rather than a full coupled solve.
- `tests/test_coolant_models.py`: unit tests of the coolant state/closure
  dispatcher (`liquid_flow/dispatch.py`).
- `tests/test_liquid_flow_shims.py`: confirms the pre-Phase-1 module paths
  still re-export correctly and emit `DeprecationWarning`. This is the only
  place in the repo that deliberately imports the old paths.

`tests/test_liquid_coupled_helical.py` covers the real `main_solver` end to
end with `coolant_model="equilibrium_liquid"` (helical only). There is still
no equivalent test for `shellntube_solver` (Phase 3, not yet done), and no
transient liquid/boiling test coverage.

`tests/test_liquid_counterflow_reference.py` covers
`solve_counterflow_liquid_reference()` (`main_solve.py`) — the physical
counter-flow shooting reference that converges the march's hot-end starting
enthalpy against the user's `T_in`/`p_in` instead of the legacy `T_out`/
`p_out` guess. Includes a HEM-monotonicity check (the shooting method's
correctness depends on `residual(h_hot_end)` being monotonic) and a
`main_steady.py` dispatch-routing test. This suite is slow (~20 min at the
real ~1378-node coil scale): the adaptive bracket-search-then-bisection
root-find needs up to ~40 full solver runs per case, and each run marches
the whole coil. Cases deliberately use a larger `mass_flow_c` (smaller
relative duty, faster convergence) rather than a boiling-scale case, to keep
this bounded — see the module docstring.

`1Dmodel/validation/water_helical_example.py` is a runnable (not pytest)
confirmed-working water recipe: `run_coflow()` is fast (no shooting needed
in co-flow); `run_counterflow_physical()` uses the shooting reference above
and is correspondingly slow. Use this instead of editing
`coolantProp`/`combustorProp` shared defaults to try water.

Machine-readable validation matrix:

```powershell
python -m hps_combustor.validation.liquid_validation_matrix
```

writes `docs/validation/liquid_validation_matrix.json` (`all_passed: true`
for what is validated); its `scope.not_yet_validated` list explicitly names
the fully coupled production liquid wall march, the transient boiling/liquid
finite-volume model, geometry-specific shell-side/helical-coil boiling
correlations, and dryout/post-CHF model validation.

## Solver Smoke Tests

Use these for changes in `1Dmodel/`:

```powershell
python -m hps_combustor.main_steady
python -m hps_combustor.main_transient
```

Select the backend in `input_data.py` with `combustorProp.HX_config`. These are smoke tests, not deterministic unit tests, and they write a timestamped archive under `runProp.output_root`.

For fast shell-and-tube transient checks, prefer:

- `combustorProp.HX_config = "shellntube"`
- `combustorProp.flow_config in ("co", "counter")`
- `transientProp.solver_method = "fixed_step"`
- `runProp.shelltube_transient_nodes = 16` for smoke, `80` for production-grid timing
- `transientProp.t_end = 5.0`, `max_step = 0.25` for speed checks

Recent reference checks:

- Shell-and-tube steady vs transient wall reconstruction, 16 nodes:
  co/counter x frozen/equilibrium/finite-rate agreed in heat rate within about
  0.4% or better after the `hot_side="inner"` conduction fix.
- Shell-and-tube counter-flow finite-rate, fixed-step, 80 nodes, 5 s simulated:
  about 5.8 s wall time and 26 fluid passes with cached FPV chemistry.

## Focused Verification Ideas

- Correlations: add simple numeric regression tests around dispatcher functions in `1Dmodel/physics/`.
- Geometry: verify diameters, hydraulic diameters, areas, and positive clearances for a known input dataclass.
- Solver changes: compare key scalar outputs from `compute_performance()` for a small, fixed configuration.
- Transient changes: check `max_step=0.25` vs `0.1` on a short case before long runs. `Q_hot` should be much less sensitive than outlet helium during early warm-up.
- Calibration changes: use one synthetic `CalibrationRecord` and mock or minimize solver calls if possible; full MCMC is too expensive for routine checks.
- Plotting changes: instantiate dashboards from a small fake `data_master` rather than running a full solver only to inspect plotting.

## Known Verification Gaps

- No root-level pytest suite covers the helical-coil or shell-and-tube solvers.
- No CI configuration is visible in this workspace.
- Solver smoke tests may be slow because they call Cantera and CoolProp repeatedly.
- First finite-rate run for a new inlet state/grid builds an FPV cache file and is slower; cached runs should be much faster.
- Some output includes non-ASCII encoding artifacts in existing files; avoid treating those as behavioral failures unless the task is documentation cleanup.

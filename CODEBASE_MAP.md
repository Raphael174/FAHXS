# Codebase Map

## Purpose

The repository models a combustor heat exchanger using quasi-1D thermal-fluid solvers. The main use case is sizing and calibrating helium-cooled combustor-HX concepts with combustion gases from liquid fuels or hydrogen in oxygen.

## Top-Level Layout

| Path | Role |
|---|---|
| `pyproject.toml` | Installs `1Dmodel/` as package `hps_combustor`. |
| `run_solver.py` | Compatibility wrapper for `hps_combustor.main_steady`. |
| `1Dmodel/` | Main solver package: inputs, solvers, physics, geometry, materials, plotting. |
| `research/flamelet_kit/` | Standalone flamelet and cooling-PFR research toolkit with tests and docs. |
| `optimization/` | Calibration and optimization wrappers around the main solver. |
| `studies/` | Historical parametric study scripts, outside the package surface. |
| `results/legacy_outputs/` | Historical generated plots and dashboards. New runs write to `zip_folders/` by default. |
| `docs/` | Methodology notes, design artifacts, solver design notes, and references. |
| `archive/` | Legacy solver variants, old input presets, injector/exploratory work, build artifacts. |
| `cache/fpv_manifolds/` | Generated finite-rate FPV manifold cache. Keep for fast repeated runs; safe to regenerate. |

## Core Entry Points

| File | Use |
|---|---|
| `1Dmodel/main_steady.py` | User entry point for steady runs; dispatches from `combustorProp.HX_config`. |
| `1Dmodel/main_transient.py` | User entry point for transient runs; dispatches from `combustorProp.HX_config`. |
| `1Dmodel/main_solve.py` | Backend steady shell-and-helical-tube combustor-HX solver. |
| `1Dmodel/main_solve_shellntube.py` | Backend steady baffled shell-and-tube solver. |
| `1Dmodel/main_solve_transient.py` | Backend transient helical-coil solver. |
| `1Dmodel/main_solve_shellntube_transient.py` | Backend transient shell-and-tube solver. |
| `optimization/calibrate.py` | Bayesian/least-squares calibration of correlation prefactors. |
| `optimization/quick_optimize.py` | Single-objective mass optimization with constraints. |
| `research/flamelet_kit/example_run.py` | Flamelet demo. |
| `research/flamelet_kit/example_cooling.py` | Cooling plug-flow reactor demo. |

## Package Import Model

`pyproject.toml` maps:

```toml
[tool.setuptools.package-dir]
"hps_combustor" = "1Dmodel"
```

After `python -m pip install -e .`, imports should use:

```python
from hps_combustor.main_solve import main_solver
from hps_combustor.input_data import coolantProp, hotgasProp
```

User runs should use `python -m hps_combustor.main_steady` or `python -m hps_combustor.main_transient`. Some backend files still keep direct-run bootstrap code for compatibility while `1Dmodel/` remains the source folder.

## Primary Data Flow

1. Input dataclasses from `1Dmodel/input_data.py` define coolant, hot gas, combustor geometry, numerical controls, system requirements, and correlation coefficients.
2. A solver initializes geometry, material functions, CoolProp helium state, and Cantera combustion state.
3. The solver marches axially or along coil arc length.
4. Each node computes coolant properties, hot-gas properties, friction, Nusselt numbers, radiation, wall conduction, stress, and state derivatives.
5. Node quantities are appended to `data_master`, created by `model_data_process/data_processing.py`.
6. `compute_performance()` or backend summaries derive scalar metrics such as heat duty, pressure drop, mass, stress ratio, and thrust estimates.
7. `result_package.py` writes numeric data, input presets, metadata, and a zip archive under `runProp.output_root`.
8. Dashboard plotting lives in `model_data_process/data_plotting.py` and `data_plotting_transient.py`.

Transient shell-and-tube long runs use `transientProp.solver_method =
"fixed_step"` by default. That path is linearly implicit in local wall-film
stiffness and has bounded cost of roughly one quasi-steady fluid pass per time
step. BDF/Radau remain validation options but can be much slower for
counter-flow because the RHS includes profile relaxation.

Finite-rate chemistry builds an FPV manifold once per inlet state/grid setting
and reuses it from `cache/fpv_manifolds/`. The first build is intentionally more
expensive than a cached run.

## Liquid Coolant / Boiling (In-Progress)

`1Dmodel/physics/liquid_flow/correlations.py`, `chf.py`, `dispatch.py`,
`governing_equations.py`, `hx_adapters.py`, and `sanity_checks.py` implement a
`p,h`-state liquid/boiling coolant path, built against literature under
`docs/reference` with validation artifacts under `docs/validation`. (The old
top-level `physics/liquid_coolant.py` / `coolant_models.py` /
`heated_liquid_channel.py` / `liquid_hx_adapters.py` paths are deprecation
shims — see `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` Phase 1.)

Wiring status (Phase 2 of the integration plan):

- **Helical steady (`main_solve.py`)**: wired into the coupled march itself
  when `coolantProp.coolant_model == "equilibrium_liquid"`. Co-flow is fully
  self-consistent; plain counter-flow's prescribed-outlet limitation is
  resolved by `solve_counterflow_liquid_reference()` (shoots the hot-end
  starting enthalpy against the physical `T_in`/`p_in`, adaptive bracket +
  bisection — see the integration plan's hardening pass).
- **Shell-and-tube (`main_solve_shellntube.py`)**: **corrected 2026-08-19**
  (this section previously said "postprocess-only"; stale) — `self._liquid_mode`
  couples `evaluate_coolant_closure` directly into `_shell_h_at`/
  `_shell_side_march`, the same coupled-march pattern as the helical solver.
  Verified by running Water and LN2/supercritical-N2 cases through it, not
  just reading the code — see `1Dmodel/validation/friday_shelltube_water.py`/
  `friday_shelltube_n2.py`. `shellntube_solver.liquid_coolant_postprocess()`
  is a separate, additional opt-in diagnostic bridge, not the only path.
- Both transient solvers are unaffected (zero liquid/boiling coolant
  presence — confirmed by grep, not assumed, 2026-08-19).
- `1Dmodel/validation/water_helical_example.py`: confirmed-working, standalone
  water/boiling recipe (does not touch shared `coolantProp`/`combustorProp`
  defaults, which remain the Helium/`shellntube` baseline).

See `docs/context/PHYSICS_CONTEXT.md` and
`docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md` for detail.
Relevant tests: `tests/test_liquid_boiling_poc.py`,
`tests/test_liquid_hx_adapters.py`, `tests/test_liquid_flow_shims.py`,
`tests/test_liquid_coupled_helical.py`,
`tests/test_liquid_counterflow_reference.py`.

The ideal-gas compressible governing equations used by the gas/helium march
live in `1Dmodel/physics/gas_flow/governing_equations.py` (moved from the
top-level `physics/governing_equations.py`, now a deprecation shim) — a
completely different equation set from the liquid one above, despite the
shared filename before Phase 1.

## What Is Usually Safe To Ignore

- `archive/` unless recovering old behavior or presets.
- `results/legacy_outputs/` unless comparing historical generated plots.
- `studies/` unless updating a parametric study workflow.
- `research/flamelet_kit/` unless doing standalone flamelet/PFR research.

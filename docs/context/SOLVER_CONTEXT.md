# Solver Context

## Authoritative Solvers

`main_steady.py` and `main_transient.py` are the intended user entry points:

```powershell
python -m hps_combustor.main_steady
python -m hps_combustor.main_transient
```

They dispatch from `combustorProp.HX_config` in `input_data.py`.

`main_solve.py` is the maintained steady shell-and-helical-tube backend. It defines `main_solver`.

`main_solve_shellntube.py` is the maintained steady baffled shell-and-tube backend. It defines `shellntube_solver`.

Older or specialized variants were moved to `archive/legacy_solvers/`. Inspect their headers before reusing logic; do not assume they are as current as `main_solve.py`.

## Inputs

All common configuration lives in `input_data.py`:

| Dataclass | Purpose |
|---|---|
| `coolantProp` | Helium fluid, mass flow, inlet/outlet guesses or boundary conditions. Also carries liquid-coolant selectors (`coolant_model`, `liquid_heat_transfer_model`, `liquid_pressure_drop_model`, `liquid_chf_model`, `liquid_chf_lut_path`) — wired into the helical steady coupled march (see Liquid Coolant Status), still opt-in/unused elsewhere; default `coolant="Helium"` is unchanged. |
| `hotgasProp` | Fuel, oxidizer, O/F, hot-gas pressure, gas mass flow, injection temperatures. |
| `combustorProp` | Helical-coil geometry, flow direction, materials, roughness, correlation selectors. |
| `shellTubeProp` | Baffled shell-and-tube geometry and correlation selectors. |
| `numericalProp` | March spacing, length limits, radiation and chemistry flags, sanity checks. |
| `transientProp` | Transient boundary-condition schedules and integration controls. |
| `system_requirements` | Ambient pressure, thrust target, burn time. |
| `CorrelationCoefficients` | Calibration knobs for heat transfer, friction, radiation, and shell-and-tube factors. |
| `runProp` | User-facing run name, archive output folder, optional schedule file, and shell-tube grid controls. |

Use dataclass instances, not classes:

```python
solver = main_solver(
    coolantProp=coolantProp(),
    hotgasProp=hotgasProp(),
    combustorProp=combustorProp(),
    numericalProp=numericalProp(),
    system_requirements=system_requirements(),
)
```

## Helical-Coil Solver Flow

`main_solver.__init__()`:

1. Stores input dataclasses and correlation coefficients.
2. Selects fuel mechanism through `choose_fuel()`.
3. Computes helical geometry, shell hydraulic diameter, areas, tube wall area, and arc-length mapping.
4. Sets `numericalProp.dx = pi * D_coil / N_arc_steps_per_turn`.
5. Builds WSGGM radiation backend if radiation is enabled.
6. Loads material property interpolators.
7. Initializes coolant state with flow direction:
   - `flow_config == "co"` starts at helium inlet and marches forward.
   - otherwise starts at helium outlet and uses negative coolant march sign.
8. Initializes Cantera hot-gas equilibrium state.
9. Creates a fresh `data_master` via `make_solver_data()`.

`main_solver.solver()` loops while the HX length is below the limit and coolant temperature stays above the CoolProp safety floor. At each node it computes properties, correlations, wall conduction, stress, ODE derivatives, records data, then advances coolant and gas states.

`compute_performance()` should be called after `solver()` before reading scalar outputs. `HX_sizing_brief()` runs global checks, computes performance, prints, and optionally plots.

## Shell-And-Tube Solver Flow

`shellntube_solver` uses predictive sweep iteration:

1. Compute Bell-Delaware geometry once.
2. Cache inlet combustion gas state once.
3. March representative tube-side hot gas against the current shell temperature profile.
4. March shell-side coolant energy against the tube duty profile.
5. Under-relax shell temperature and repeat until convergence.

This is not a simple single-pass co/counter-flow ODE like the helical solver; do not port helical flow assumptions into this solver.

## Transient Solver Flow

Both maintained transient solvers integrate one lumped thickness-mean wall
temperature per axial node. At each time step/RHS evaluation they reconstruct
quasi-steady fluid fields with a `fluid_pass()`.

Important planning update: the next-generation transient solver should not keep
stretching this `fluid_pass()` abstraction. The planned `transient_core` is a
finite-volume transient-coolant architecture:

```text
state(t) = [Tbar_wall(x), T_helium(x)]
hot gas  = quasi-steady initially
```

See `docs/solver_design/TRANSIENT_CORE_IMPLEMENTATION_PLAN.md` for the staged
implementation and validation plan. Until that core is implemented, the notes
below describe the currently maintained quasi-steady-fluid transient solvers.

Current `transient_core` implementation status:

- Implemented: geometry-neutral coolant finite-volume step,
  axial grid descriptor, `[Tbar_wall, T_coolant]` state layout, linear implicit
  wall/coolant step, shared energy/timescale diagnostics, helical
  geometry/inventory adapter, generic fixed-step wall/coolant integrator,
  schedule-breakpoint time-grid builder, shared schedule interpolation helpers,
  shell-and-tube adapter, and shell-and-tube `main_transient.py` dispatch via
  `transientProp.fluid_model = "transient_coolant"`.
- Tested: zero-flow local soak, co/counter-flow advection direction,
  constant-property first-law residuals, helical area/perimeter/inventory
  conversion, coolant-film dispatcher scaling, wall-flux unit conversion,
  shell-and-tube adapter assembly, and short shell-and-tube dispatch smoke runs.
- Pending: helical hot-gas marching wrapper, pressure drop, helical
  transient-core dispatch, shell-side residence refinement, dashboard-field
  review, and broader validation.

Helical transient:

- Uses the same geometry and chemistry setup as `main_solver`.
- Supports co-flow and counter-flow.
- Counter-flow uses a warm-started coolant-profile relaxation.

Shell-and-tube transient:

- Inherits geometry, material, tube-side hydraulics, and Bell-Delaware setup
  from `shellntube_solver`.
- Supports co-flow and counter-flow.
- Production long runs should use `transientProp.solver_method = "fixed_step"`.
  This is a linearly-implicit wall update:
  `Tbar_next = Tbar + dt * R / (1 + dt * lambda)`, where `lambda` is estimated
  from the current hot/cold film conductances.
- BDF/Radau remain useful for validation but can be much slower because they
  repeatedly evaluate the expensive counter-flow RHS to estimate Jacobians.

## Liquid Coolant Status

A liquid/boiling coolant option is under active development. Status per
solver (see `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` for the
phased plan):

- **Helical steady (`main_solve.py`) has a real coupled liquid march** when
  `coolantProp.coolant_model == "equilibrium_liquid"`: coolant state is
  `(p, enthalpy_c)`, evaluated each node through
  `physics/liquid_flow/dispatch.evaluate_coolant_closure()`, integrated via
  `dh/dx = dQ/mdot` and `dp/dx = -friction`. `main_solver.solver()` dispatches
  on `coolant_model` internally — there is no separate liquid entry point.
  Co-flow is fully self-consistent. A sanity-gate report runs automatically
  at the end of a liquid-mode solve and is stored as
  `self.liquid_sanity_report`.
- **Counter-flow physical reference**: the plain liquid counter-flow march
  still starts from the legacy `coolantProp.T_out`/`p_out` guess (a
  single-phase `(T,P)` state — cannot represent a two-phase starting point).
  `solve_counterflow_liquid_reference()` in `main_solve.py` resolves this:
  it shoots the march's hot-end starting **enthalpy** (never temperature) so
  the cold end matches the user's physical `coolantProp.T_in`/`p_in`, using
  an adaptive-bracket-then-bisection root-find (secant was tried first and
  discarded — it overshot into nonphysical states near boiling onset).
  `main_steady.py`'s `run_steady()` dispatches to it automatically when
  `numericalProp.counterflow_physical_steady_reference = True` and
  `coolant_model == "equilibrium_liquid"` (instead of the gas-only
  `solve_counterflow_physical_reference()`, which shoots on temperature and
  is invalid inside the dome). Only enthalpy is shot — the hot-end pressure
  is approximated as `p_in` (friction drop is normally small) — see
  `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`, "Post-Phase-2
  Hardening Pass," for the full derivation and residual limitation.
- **Shell-and-tube (`main_solve_shellntube.py`) is still post-process-only**:
  `shellntube_solver.liquid_coolant_postprocess()` consumes an
  already-converged `dQ` duty profile and stores `self.liquid_coolant`
  diagnostics (`p,h`, quality, void, CHF margin) without changing the wall
  temperatures or coolant march that produced `dQ`; the coupled march still
  computes coolant `T,p` via direct CoolProp `PropsSI` calls regardless of
  `coolant_model`. This is Phase 3 of the integration plan.
- `main_solver.liquid_coolant_postprocess()` still exists on the helical
  solver too (unchanged) as an independent cross-check path — Phase 2's test
  suite (`tests/test_liquid_coupled_helical.py`) uses it to verify the new
  coupled march and the pre-existing bridge agree.
- Neither transient solver (`main_solve_transient.py`,
  `main_solve_shellntube_transient.py`) references liquid/boiling code at all.
- See `docs/context/PHYSICS_CONTEXT.md` (Liquid Coolant / Boiling section) and
  `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md` for correlation
  detail.

## Common Pitfalls

- `data_master` must be a fresh dict per solver instance. Use `make_solver_data()`.
- `CorrelationCoefficients` belongs in `input_data.py`; optimization and calibration depend on those field names.
- `flow_config` changes the coolant state direction and pressure-drop interpretation.
- `numericalProp.equilibrium_dh_gas_ON` is a legacy alias. New logic should prefer `numericalProp.chemistry_model`.
- `main_solve.py` uses Darcy friction convention consistently. Do not multiply Ali/Colebrook friction by four.
- Steady and transient chemistry default to the FPV finite-rate path in `physics/combustion_chemistry/fpv_manifold.py`.
- `equilibrium` and `frozen` remain comparison modes; do not let the legacy `equilibrium_dh_gas_ON` flag override an explicit `finite_rate` setting.
- Solver outputs are often lists in `data_master`; convert to `numpy.array` before vector operations.
- Shell-and-tube hot gas is inside the tubes. `hot_side="inner"` must be passed
  through steady conduction and transient reconstruction; otherwise steady and
  transient heat rates disagree.
- Do not disable `numericalProp.fpv_cache_dir` for normal finite-rate work.
- `main_solver` requires `combustorProp.HX_config == "shellnHelicalTube"` and
  raises otherwise (added 2026-07-13, after this exact mismatch silently
  corrupted axial-length accounting in test fixtures for months). Always pin
  `HX_config` explicitly when constructing `combustorProp` for `main_solver`
  — do not rely on the dataclass default. The real helical coil for this
  combustor's geometry is ~1378 arc-length nodes, not ~100; duty scales
  accordingly (~150-300 kW, not ~20 kW) — do not reuse a short-coil-scale
  test case's `mass_flow_c`/`T_in`/`p_in` assuming a similar duty.
- `coolantProp` has no fluid-specific constant fields (e.g. no `molar_mass`);
  `main_solver` looks up molar mass from CoolProp for whatever `coolant` is
  configured. Do not reintroduce a hardcoded per-fluid property field.

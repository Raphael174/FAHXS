# Transient Core Implementation Plan

Status: planned architecture for the next solver generation. WP1, WP1b, the
geometry-neutral part of WP6, and the geometry/inventory part of WP3 are
implemented; WP2 has started. The geometry-independent coolant finite-volume
kernel, axial grid descriptor, state layout helpers, first linear wall/coolant
implicit step, shared diagnostics, and first helical geometry bridge are
implemented in `1Dmodel/transient_core/` with focused tests in
`tests/test_transient_core_coolant_fv.py`.

Purpose: transform Combustor-HX from a wall-transient solver with
quasi-steady fluids into a near fully transient heat-exchanger solver with:

```text
transient wall/material
transient helium coolant inventory and advection
quasi-steady hot gas, initially
optional transient hot gas as the final extension
```

This document is the implementation and validation plan. It should be read with:

```text
docs/TECHNICAL_REFERENCE.md
docs/context/TRANSIENT_STATUS.md
docs/reference/nellis_klein_heat_transfer_hx_notes.md
```

## 1. Decision Summary

The current transient solver architecture is:

```text
state(t) = wall mean temperature Tbar_w(x)
fluids   = quasi-steady algebraic fluid_pass(Tbar_w, boundary_conditions)
```

That is no longer rigorous for bang-bang helium operation because the helium
fluid inventory has memory. When helium mass flow drops, the residence time can
become comparable to or much larger than the wall time and valve cycle time. At
zero flow, the current quasi-steady coolant march has no physical meaning: the
actual coolant should remain as stagnant thermal mass exchanging heat with the
wall.

The new production architecture should be:

```text
state(t) = [Tbar_w(x), T_c(x)]
hot gas  = quasi-steady march against the current wall/coolant state
```

Keep the existing quasi-steady-fluid transient solver as:

```text
fluid_model = "quasi_steady"
```

and add the new finite-volume core as:

```text
fluid_model = "transient_coolant"
```

Do not start a new repository. The existing codebase contains reusable geometry,
correlations, wall reconstruction, chemistry manifolds, schedules, packaging, and
dashboard infrastructure. The correct move is a new internal transient core with
clean adapters for each HX configuration.

## 2. Scope Boundary

### In Scope For The First Transient Core

- Helium coolant temperature storage and advection.
- Wall/material thermal storage.
- Co-flow and counter-flow through advection direction, not profile relaxation.
- Quasi-steady hot gas with finite-rate/equilibrium/frozen chemistry as today.
- Existing shell-and-helical-tube and shell-and-tube configurations.
- Scheduled helium mass flow, helium inlet temperature, helium inlet pressure,
  hot-side schedule, GOX/ignition schedule.
- Energy accounting including coolant and wall internal energy.
- Result packaging compatible with dashboard reloads.

### Out Of Scope Until Later

- Full compressible momentum and pressure waves.
- Valve/orifice pressure-control dynamics.
- Hot-gas transient storage and transport.
- Transient finite-rate gas chemistry state storage.
- Two-phase LOX chilldown or boiling.
- Radial wall finite-difference model.
- Structural FEA.

The optional final extension is transient hot gas. It is not required to declare
the first transient-core goal successful because hot-side residence is about
0.6 ms in the current regime, much faster than the helium/wall/control dynamics.

## 3. Governing Model

### 3.1 State Variables

For each axial cell `i`:

```text
Tbar_w[i]   wall thickness-mean temperature [K]
T_c[i]      helium bulk temperature in coolant cell [K]
```

Optional future states:

```text
p_c[i]      coolant pressure inventory [Pa]
T_g[i]      hot-gas transient temperature [K]
h_gas[i]    hot-gas enthalpy state [J/kg]
Yc[i]       finite-rate progress variable state [-]
T_w[i,j]    radial wall temperature layers [K]
```

### 3.2 Wall Balance

The wall remains one thermal state per axial node in the first implementation:

```text
rho_w * cp_w * A_wall * dTbar_w/dt = dq_hot__dx - dq_cold__dx
```

The current `fluxes_at_Tbar()` quadratic face reconstruction remains valid for
thin walls if the wall audit passes:

```text
tau_diff_wall << tau_wall
Bi_wall not excessive
```

### 3.3 Coolant Finite-Volume Balance

For cell `i` with coolant volume `V_c[i]`:

```text
rho_c[i] * V_c[i] * cp_c[i] * dT_c[i]/dt =
    mdot_c * cp_up[i] * (T_up[i] - T_c[i])
  + dq_cold_total__dx[i] * dx[i]
```

where:

```text
dq_cold_total__dx = dq_cold__dx              helical, single coil pass
dq_cold_total__dx = dq_cold__dx * N_tubes    shell-and-tube
```

At zero flow:

```text
mdot_c = 0
```

the advection term vanishes and the model becomes local wall/coolant thermal
soak:

```text
rho_c * V_c * cp_c * dT_c/dt = dq_cold_total__dx * dx
```

This is the key physics missing from the current quasi-steady transient solver.

### 3.4 Hot Gas Closure

The hot side is algebraic in the first transient core:

```text
given Tbar_w(x), T_c(x), BC_hot(t):
    march hot gas quasi-steadily
    evaluate chemistry manifold
    evaluate h_g(x)
    evaluate dq_hot__dx(x)
```

This is justified initially by:

```text
tau_hot ~ 0.6 ms
tau_hot << helium/wall/control timescale
```

The model still needs a runtime audit. If hot-side residence becomes too large
relative to the forcing timescale, the solver should warn that transient hot gas
is required.

## 4. Numerical Strategy

### 4.1 Fixed-Step IMEX Core

Use a bounded-cost fixed-step IMEX update:

```text
advection:       implicit upwind
wall exchange:   linearly implicit local source
properties:      evaluated from current or Picard-updated state
hot gas:         quasi-steady pass per Picard iteration
```

Do not use purely explicit coolant advection. In the helical case:

```text
U_c ~ 100 m/s
L ~ 6 m
N ~ 80
dx ~ 0.075 m
explicit CFL dt < 0.75 ms
```

which is smaller than the user's 2.5 ms schedule sampling. Implicit upwind avoids
making the solver impractically slow.

### 4.2 Per-Step Algorithm

For each time step `t_n -> t_{n+1}`:

```text
1. read/interpolate schedules at t_n or t_{n+1/2}
2. unpack state: Tbar_w^n, T_c^n
3. for k = 1..n_picard:
      a. evaluate coolant properties at T_c^k, p_c boundary/estimate
      b. run quasi-steady hot-gas adapter using Tbar_w^k and T_c^k
      c. reconstruct wall faces and fluxes with fluxes_at_Tbar()
      d. assemble coolant implicit-upwind system
      e. update T_c^{k+1}
      f. update Tbar_w^{k+1} with linearly implicit wall source
      g. check Picard convergence
4. store diagnostics and advance
```

Default:

```text
n_picard = 2
```

Single-Picard mode can be allowed for fast screening after validation.

### 4.3 Counter-Flow Handling

In the current quasi-steady transient, counter-flow requires coolant profile
relaxation because both fluid fields are algebraic. In the new transient core,
counter-flow is just the opposite advection direction:

```text
co-flow:      coolant inlet at i = 0
counter-flow: coolant inlet at i = N-1
```

The coolant state carries the memory, so no counter-flow profile relaxation loop
is needed in `transient_coolant` mode.

### 4.4 Low-Flow And Zero-Flow Handling

Define:

```text
mdot_adv_floor = transientProp.min_mdot_for_advection
```

For `abs(mdot_c) <= mdot_adv_floor`:

```text
advection term = 0
inlet boundary not imposed on internal cells
coolant cells exchange heat locally with the wall
```

For small nonzero flow, implicit upwind remains stable and produces a delayed
flush of the stored coolant temperature field.

## 5. Code Organization

Add a new package:

```text
1Dmodel/transient_core/
    __init__.py
    diagnostics.py
    state.py
    grid.py
    coolant_fv.py
    wall_coolant.py
    hot_quasisteady.py
    integrator.py
    schedules.py
    diagnostics.py
    adapters_helical.py
    adapters_shelltube.py
```

Responsibilities:

| Module | Responsibility |
|---|---|
| `state.py` | Pack/unpack state vectors; name slices; initialize wall and coolant states. |
| `grid.py` | Cell centers, lengths, volumes, flow direction, perimeters. |
| `coolant_fv.py` | Implicit upwind helium finite-volume update. |
| `wall_coolant.py` | Geometry-neutral implicit wall/coolant exchange plus coolant advection. |
| `hot_quasisteady.py` | Common interface/protocol for quasi-steady hot-side adapters. |
| `integrator.py` | Generic fixed-step wall/coolant driver; schedule dispatch still pending. |
| `schedules.py` | Shared schedule interpolation and breakpoint-time extraction. |
| `diagnostics.py` | Energy audits, residence time, wall time constants, quasi-steady timescale audits. |
| `adapters_helical.py` | Helical geometry/inventory bridge now; later hot-side/correlation adapter. |
| `adapters_shelltube.py` | Geometry and hot-gas adapter for baffled shell-and-tube. |

Existing files should dispatch into the new core rather than being rewritten
around it:

```text
main_transient.py
main_solve_transient.py
main_solve_shellntube_transient.py
```

The old solvers remain callable for reference:

```text
transientProp.fluid_model = "quasi_steady"
```

The new core is selected with:

```text
transientProp.fluid_model = "transient_coolant"
```

## 6. Input And Output Changes

### 6.1 New `transientProp` Fields

Planned fields:

```python
fluid_model: str = "quasi_steady"  # "quasi_steady" | "transient_coolant" | "transient_both"
coolant_initial_condition: str = "uniform_inlet"  # "uniform_inlet" | "steady"
transient_coolant_picard_iter: int = 2
transient_coolant_picard_tol_K: float = 1e-3
min_mdot_for_advection: float = 1e-8
timescale_warning_ratio: float = 0.05
hot_side_transient_optional: bool = False
```

Keep scheduled mass-flow inputs first. Pressure-controlled valve/orifice physics
is a later boundary-condition module.

### 6.2 Output Additions

Transient packages should add:

```text
T_c_state[t, x]              stored helium-cell temperature
T_c_out[t]                   boundary/outlet reconstructed from coolant state
U_coolant_J[t]               coolant internal energy
U_wall_J[t]                  wall internal energy
energy_residual_J[t]         cumulative energy-balance residual
tau_fluid_s[t]
tau_wall_min_s[t]
tau_hot_s[t]
fluid_model
coolant_advection_direction
picard_iterations[t]
warnings[t]
```

## 7. Geometry Adapter Requirements

Each configuration adapter must provide arrays:

```text
N
dx[i]
x_center[i]
coolant_volume[i]
coolant_flow_area[i]
coolant_heat_perimeter[i]
wall_area_per_length[i]
hot_perimeter[i]
cold_perimeter[i]
flow_direction
```

### 7.1 Helical Adapter

Coolant inventory:

```text
V_c[i] = A_tube_inner * dx_arc[i] * N_parallel_paths
```

Heat perimeter:

```text
P_c = pi * D_inner * N_parallel_paths
```

The current helical grid is along coil arc length. The hot-gas quasi-steady pass
can reuse the current manifold/radiation/correlation logic.

### 7.2 Shell-And-Tube Adapter

Do not use Bell-Delaware `S_m` as the coolant inventory volume. `S_m` is a
maximum crossflow area for heat-transfer and pressure-drop correlations; it is
not shell-side helium hold-up.

Approximate first coolant inventory:

```text
V_shell_fluid =
    shell_internal_volume
  - tube_outer_displaced_volume
  - baffle_solid_volume_approx
  - tubesheet_volume_if_inside_active_domain
```

Then distribute it across the active axial cells. This is sufficient for the
first transient coolant model; later revisions can allocate volume by baffle
compartment.

Shell-side heat perimeter remains:

```text
P_c_total = N_tubes * pi * D_tube_outer
```

For each representative tube:

```text
dq_cold_total__dx = dq_cold_per_tube__dx * N_tubes
```

## 8. Implementation Work Packages

### WP0: Documentation And Requirements Lock

Deliverables:

- This plan.
- Update doc indexes and context notes.
- Add current assumptions to `TECHNICAL_REFERENCE.md` after implementation
  details are confirmed.

Validation gate:

- A future agent can identify target files, equations, validation cases, and
  phase boundaries without reading the full conversation.

### WP1: Pure Finite-Volume Coolant Kernel

Build `transient_core.coolant_fv` independent of HX geometry.

Implementation status:

```text
implemented: implicit_upwind_step()
implemented: CoolantStepResult diagnostics
implemented: zero-flow local soak behavior
implemented: co-flow and counter-flow advection direction
implemented: constant-property first-law residual diagnostic
not yet: variable-property Picard coupling to CoolProp
not yet: integration with wall and HX geometry adapters
```

Inputs:

```text
T_c_old[i]
rho_c[i], cp_c[i]
V_c[i]
mdot_c
T_in
dq_cold_total__dx[i]
dx[i]
flow_direction
dt
```

Outputs:

```text
T_c_new[i]
T_out
advective_energy_in_out
coolant_energy_change
```

Validation:

- Zero-flow soak against analytical lumped cell.
- Pure advection step with no heat transfer preserves a step profile shape with
  expected implicit-upwind numerical diffusion.
- Constant-wall-temperature duct approaches the exponential internal-flow
  solution as grid is refined.
- Energy balance closes for one time step.

### WP1b: Axial Grid Descriptor

Build `transient_core.grid` so geometry adapters pass explicit arrays into the
core rather than relying on solver-specific indexing conventions.

Implementation status:

```text
implemented: AxialGrid
implemented: uniform grid constructor
implemented: cell centers, dx, total length
implemented: coolant and wall volumes
implemented: hot/coolant perimeters
implemented: inlet/outlet index helpers for co-flow and counter-flow
implemented: conversion from W/m heat rates to W per cell
```

Validation:

- Uniform grid geometry gives expected centers, volumes, and length.
- Reverse flow maps inlet to the last cell and outlet to the first cell.
- Per-length heat rates multiply by local `dx`.
- Nonphysical geometry is rejected.

### WP2: Wall-Coolant Coupling Kernel

Build `transient_core.wall_coolant` as the geometry-neutral coupled thermal
kernel, then wrap it with geometry adapters that call the existing
`OneDimensionalSteadyConduction_ShellnHelicalTube.fluxes_at_Tbar()` or the
shell-and-tube equivalent.

Current implementation status:

```text
implemented: implicit_wall_coolant_step()
implemented: WallCoolantStepResult diagnostics
implemented: zero-flow local wall/coolant soak
implemented: co-flow and counter-flow coolant advection direction
implemented: constant-property first-law residual over wall + coolant
implemented: first helical and shell-and-tube wall reconstruction adapters
not yet: nonlinear Picard update of G_i from h_c(T, mdot, p)
not yet: production geometry adapters that compute full Qhot_i and G_i from
         properties, correlations, and hot-side marches at each step
```

The implemented linear cell model is:

```text
Cw_i (Tw_i^{n+1} - Tw_i^n) / dt =
    Qhot_i - G_i (Tw_i^{n+1} - Tc_i^{n+1})

Cc_i (Tc_i^{n+1} - Tc_i^n) / dt =
    mdot_c cp_i (Tup_i^{n+1} - Tc_i^{n+1})
  + G_i (Tw_i^{n+1} - Tc_i^{n+1})
```

The upstream coolant temperature is the inlet boundary at the first upwind cell
or the already solved upstream cell temperature, so the implicit upwind solve
remains a sequence of local 2x2 systems.

Inputs:

```text
Tbar_w_old[i]
T_c[i]
wall heat capacity Cw_i
coolant heat capacity Cc_i
coolant cp_i
mdot_c
T_c_in
hot heat Qhot_i
wall-to-coolant conductance G_i
flow_direction
dt
```

Outputs:

```text
Tbar_w_new[i]
T_c_new[i]
T_c_out
heat_wall_to_coolant_W[i]
U_wall_change
U_coolant_change
first-law residual
```

Validation:

- Zero-flow case conserves combined wall/coolant energy and does not impose an
  inlet boundary inside stagnant coolant.
- Flow direction changes the physical outlet side.
- With no advection, the wall/coolant temperature difference approaches the
  finite-inventory relative equilibrium.
- Combined wall + coolant first-law residual closes near roundoff for
  constant-property test cases.
- Pending adapter validation: constant fluid temperatures settle to the steady
  wall solve from existing geometry-specific conduction objects.

### WP3: Helical Transient-Coolant Adapter

Build `adapters_helical.py`.

Current implementation status:

```text
implemented: flow_direction_from_config()
implemented: build_helical_core_geometry()
implemented: build_helical_core_geometry_from_solver()
implemented: wall_heat_capacity_J_K()
implemented: coolant_heat_capacity_J_K()
implemented: conductance_from_h()
implemented: helical_coolant_film()
implemented: helical_wall_flux()
not yet: helical hot-side quasi-steady march producing Qhot_i
not yet: pressure drop update
not yet: production dispatch from main_transient.py
```

Convention:

```text
Axial coordinate: helical tube arc length.
Grid areas/perimeters: totals across all parallel coils.
Heat rates entering the core: total W per cell, not W/m per single tube.
```

Reuse:

```text
helical geometry setup
coil coolant h/friction correlations
hot-side chemistry manifold
radiation lookup
wall reconstruction
```

Validation:

- Implemented geometry validation: total parallel-coil area/perimeter scaling.
- Implemented inventory validation: wall/coolant heat capacities from grid
  volumes.
- Implemented conductance validation: `G_i = h_c,i * P_coolant,i * dx_i`.
- Implemented coolant-film validation: total-area velocity, Reynolds/Prandtl,
  friction/Nusselt dispatcher calls, `h_c`, and `G_i` scaling.
- Implemented wall-flux validation: existing `fluxes_at_Tbar()` outputs are
  converted from single-tube `W/m` to total `W` per cell using
  `dx_i*N_parallel`.
- High-flow smooth schedule matches current quasi-steady transient within a
  small tolerance.
- Low-flow/bang-bang schedule shows delayed helium outlet packets.
- Counter-flow works without profile relaxation.
- 6 m / 100 m/s residence-time test gives about 0.06 s delay scale.
- Constant-BC transient settles to physical steady reference.

### WP4: Shell-And-Tube Transient-Coolant Adapter

Build `adapters_shelltube.py`.

Current implementation status:

```text
implemented: shelltube_flow_direction()
implemented: build_shelltube_core_geometry()
implemented: build_shelltube_core_geometry_from_solver()
implemented: shelltube_wall_heat_capacity_J_K()
implemented: shelltube_coolant_heat_capacity_J_K()
implemented: shelltube_conductance_from_h()
implemented: shelltube_shell_film()
implemented: shelltube_tube_gas_film()
implemented: shelltube_wall_flux()
implemented: shelltube_hot_gas_march() sequential enthalpy/progress scaffold
implemented: fpv_gas_state_provider()
implemented: equilibrium_gas_state_provider()
implemented: oxygen_gas_state_provider()
implemented: shelltube_step_inputs()
implemented: run_shelltube_transient_core()
implemented: per-step gas_provider_at_time() and mdot_hot_total_at_time() hooks
implemented: no-hot-flow branch for helium-only pre-ignition phases
implemented: shell-and-tube solve_transient_core() bridge
implemented: shell-and-tube production dispatch from main_transient.py
not yet: shell-side residence refinement beyond free-volume baseline
```

Convention:

```text
Axial coordinate: straight tube length.
Coolant hold-up area: shell bore area - total tube outer area.
Wall/perimeters: totals across all tubes.
Heat rates entering the core: total W per cell.
```

Reuse:

```text
EchTherm shell/tube geometry
Bell-Delaware shell-side h
tube-side hot gas correlations
FPV/equilibrium/frozen chemistry
hot_side="inner" wall orientation
```

Validation:

- Implemented hold-up validation: shell-side coolant area is shell bore area
  minus displaced tube outer area and must remain positive.
- Implemented geometry validation: total tube wall area and total hot/cold
  perimeters scale with `N_tubes`.
- Implemented inventory validation: shell coolant and tube wall heat capacities
  use grid volumes.
- Implemented conductance validation: `G_i = h_shell,i * P_outer,total * dx_i`.
- Implemented shell-film validation: Bell-Delaware call receives `G_s`, `Re_s`,
  `Pr_s`, density in the geometry dict, and produces `h_shell`, `G_i`, and
  shell pressure-drop diagnostics.
- Implemented tube-gas-film validation: total hot-side mass flow is divided by
  `N_tubes`, smooth/grooved tube correlation paths are selected correctly, and
  `h_g`, hot conductance, and pressure-drop-per-length diagnostics scale from
  supplied property arrays.
- Implemented wall-flux validation: existing `fluxes_at_Tbar()` is called with
  `hot_side="inner"` and representative-tube `W/m` outputs are converted to
  total `W` per cell with `dx_i*N_tubes`.
- Implemented hot-gas-march validation: an injected gas-state provider receives
  the current enthalpy removed and progress variable, the representative-tube
  heat rate advances `h_removed`, and bundle-total heat rates are returned for
  the transient wall/coolant core.
- Implemented gas-provider validation: FPV wrappers return the finite-rate
  progress source, equilibrium/frozen wrappers return zero progress source, and
  the oxygen provider maps `h_removed` through CoolProp `(H, P)` before
  evaluating sensible-flow properties.
- Implemented step-input validation: coolant properties, shell film, hot-gas
  march, wall heat capacity, coolant heat capacity, heat vector, conductance,
  and flow direction assemble into the `WallCoolantStepInputs` contract used by
  the generic integrator.
- Implemented adapter-runner validation: schedule breakpoints are inserted into
  the fixed-step grid, scheduled values are passed into step assembly at interval
  starts, the generic wall/coolant integrator advances both states, and per-step
  shell/hot-gas diagnostics are retained.
- Implemented dispatch smoke validation: 0.25 s / 4-node finite-rate combustion
  and cold-He no-ignition shell-and-tube runs complete through
  `main_transient.py` with `fluid_model="transient_coolant"`.
- Implemented short-run validation note:
  `docs/validation/TRANSIENT_CORE_SHORT_RUNS_2026-07-09.md`. Two 2 s / 20-node
  shell-and-tube runs completed: finite-rate bang-bang helium and GOX chilldown
  followed by ignition/hot-gas ramp. Runtime was about 1.5 s in-memory and
  about 2.0 s with result packaging.
- Implemented repeatable short-run matrix runner:
  `python -m hps_combustor.validation.transient_core_short_runs`. It writes
  `docs/validation/transient_core_short_run_results.json`.
- Latest matrix result: bang-bang timestep refinement looks reasonable for a
  smoke check (`T_c_out` shift about 3.5 K from `max_step=0.25` to `0.10 s`),
  but GOX/ignition grid refinement from 20 to 40 nodes is not yet converged
  (`T_wall_max` shift about 74 K and large FPV progress-variable shift).
  Production validation still needs a 20/40/80 or 40/80/160 grid study after
  model choices settle.
- Implemented dashboard-style transient-core output fields for shell-and-tube:
  gas velocity/Re/Pr/Nu/friction/dp, shell mass flux/Re/Pr/dp estimate,
  coolant `rho/mu/k/cp`, gas enthalpy removed, and FPV progress variable.
- Implemented zero-shell-flow fallback: when `mdot_shell <= mdot_floor`,
  `shelltube_shell_film()` bypasses Bell-Delaware and uses the stagnant
  conduction scale `h_shell = k_He / D_tube_outer`, `dp_shell = 0`. This keeps
  bang-bang shutoff numerically well-defined but remains a validation item, not
  a final low-flow heat-transfer model.
- 150 g/s, 90 K, 82 bar EchTherm case gives helium residence time about
  15-30 ms by order of magnitude.
- High-flow smooth schedule matches current quasi-steady transient.
- Bang-bang flow captures stagnant shell-side helium heating/cooling and delayed
  outlet response.
- Co-flow/counter-flow differ correctly without relaxation loops.
- Constant-BC transient settles to steady shell-and-tube solver.

### WP5: Unified Integrator And Dispatch

Build `integrator.py` and wire it through `main_transient.py`.

Current implementation status:

```text
implemented: WallCoolantStepInputs
implemented: WallCoolantIntegrationResult
implemented: fixed_time_grid()
implemented: integrate_wall_coolant_fixed_step()
implemented: interp_schedule()
implemented: schedule_times()
implemented: collect_transient_schedule_times()
implemented: history storage for T_wall, T_coolant, outlet T, heat transfer, residuals
implemented: shell-and-tube dispatch from main_transient.py
implemented: result package compatibility through solver.time_series
not yet: helical transient-coolant dispatch from main_transient.py
not yet: dashboard semantic review for new fields
```

Validation:

- Implemented integrator validation: repeated fixed-step integration matches
  repeated direct calls to `implicit_wall_coolant_step()`.
- Implemented zero-flow integration validation: stagnant coolant evolves by
  local wall/coolant heat exchange and is not overwritten by inlet temperature.
- Implemented schedule-grid validation: schedule breakpoints, requested output
  times, and endpoints are inserted; duplicate/out-of-range points are removed.
- Implemented schedule interpolation validation: linear interpolation with
  flat-held endpoints matches the legacy transient solvers.
- `fluid_model="quasi_steady"` still dispatches to current behavior.
- `fluid_model="transient_coolant"` dispatches through `transient_core`.
- Output time vectors and package format remain dashboard-compatible.
- No existing steady solver behavior changes.

### WP6: Diagnostics, Audits, And Packaging

Build/extend diagnostics:

```text
tau_fluid = coolant hold-up / max(abs(mdot_c), floor)
tau_wall = rho_w*cp_w*A_wall/(h_hot*P_hot + h_cold*P_cold)
tau_hot = hot-side residence estimate
tau_BC = local schedule change timescale
energy residual = Q_hot_in_time - Q_cold_out_time - dU_wall - dU_coolant
```

Current implementation status:

```text
implemented: energy_audit()
implemented: residence_time_s()
implemented: wall_time_constant_s()
implemented: timescale_audit()
implemented: EnergyAudit and TimescaleAudit dataclasses
implemented: shell-and-tube transient-core scalar/field histories in result packages
not yet: dashboard exposure review for new transient-coolant histories
not yet: warnings emitted by main_transient.py for transient_core runs
```

Validation:

- Zero flow gives infinite coolant residence time.
- Wall time constants use `Cw/(G_hot + G_cold)` and insulated cells return
  infinity.
- Energy audit scales residuals against the largest relevant energy term.
- Timescale audit flags slow coolant relative to wall/control forcing while
  allowing the current fast hot-gas residence regime.
- Pending integration validation: warnings trigger on quasi-steady mode for
  bang-bang helium.
- Pending integration validation: no warnings for cases where
  `tau_fluid << tau_wall` and schedules are slow.
- Pending integration validation: energy residual remains below acceptance
  criteria in packaged runs.

### WP7: Full Validation Matrix

Run a compact but meaningful validation set:

```text
configs: helical, shellntube
flow: co, counter
chemistry: finite_rate primary; frozen/equilibrium smoke
fluid_model: quasi_steady, transient_coolant
schedule: constant, smooth ramp, bang-bang helium
grid: N=40, 80, 160
dt: 5 ms, 2.5 ms, 1.25 ms for bang-bang cases
```

Acceptance:

```text
constant-BC settle-to-steady Q error < 1-2%
global energy residual < 1-2% over run
dt convergence of Q and outlet T < 2-5%
grid convergence monotonic or asymptotic
no physical clamps in nominal validation cases
outlet delay matches residence-time order of magnitude
```

### WP8: Optional Transient Hot Gas

Only implement if audits show hot-side quasi-steady is not defensible.

Additional states:

```text
h_removed[i] or T_g[i]
Yc[i] for finite-rate chemistry
possibly p_g[i] if pressure coupling is needed
```

Governing form:

```text
rho_g*A_g*cp_g*dT_g/dt + mdot_g*cp_g*dT_g/dx =
    -h_g*P_g*(T_g - T_wg) - radiation/source terms

dYc/dt + U_g*dYc/dx = omega_Yc
```

Validation:

- Hot-side pure advection delay.
- Hot-side energy conservation with wall.
- Finite-rate progress-variable transport consistency with current FPV
  quasi-steady limit.
- At `tau_hot ~ 0.6 ms`, transient-hot and quasi-steady-hot results should agree
  for the current engine-level forcing.

## 9. Validation Cases In Detail

### Case A: Zero-Flow Coolant Soak

Setup:

```text
mdot_c = 0
hot side disabled or fixed
wall and coolant start at different temperatures
```

Expected:

```text
each coolant cell and wall cell relax toward local equilibrium
no inlet boundary is injected into the domain
total wall+coolant energy conserved if hot side disabled
```

### Case B: Constant-Wall Duct Analytical Check

Setup:

```text
T_wall = constant
mdot_c = constant
no wall state integration
```

Analytical reference:

```text
T_c(x) = T_wall - (T_wall - T_in)*exp[-h*P*x/(mdot*cp)]
```

Expected:

```text
finite-volume result approaches analytical profile as N increases
```

### Case C: Step Flow Flush

Setup:

```text
period of mdot_c = 0
then step to finite mdot_c
```

Expected:

```text
stored coolant thermal field is flushed out after residence-time order
old quasi-steady solver has no comparable delayed packet
```

### Case D: Smooth High-Flow Ramp

Setup:

```text
slow helium mass-flow ramp
tau_fluid << tau_wall and tau_BC
```

Expected:

```text
transient_coolant approximates quasi_steady
```

### Case E: Counter-Flow

Setup:

```text
flow_config = "counter"
constant and bang-bang helium cases
```

Expected:

```text
coolant inlet applied at opposite end
no profile-relaxation loop
correct outlet end reported
constant-BC settle matches counter-flow steady reference
```

### Case F: Full Energy Audit

For every transient run:

```text
integral(Q_hot_total - Q_cold_total) dt =
    delta U_wall + delta U_coolant + residual
```

Acceptance:

```text
abs(residual) / max(abs(integral Q_hot dt), eps) < 1-2%
```

## 10. Start-From-Scratch Decision Criteria

Continue inside this repository unless the target expands to:

```text
compressible pressure waves
valve/tank/HX coupled mass dynamics as primary physics
two-phase LOX
hot-gas transient chemistry as mandatory first-class state
structural FEA
multi-dimensional shell-side CFD surrogate
```

If those become mandatory, a separate simulator may be cleaner. For the present
goal, a new `transient_core` package inside the existing codebase is the more
controlled path.

## 11. Definition Of Done For The First Transient Core

The first transient-core goal is complete only when:

- `fluid_model="transient_coolant"` runs for both HX configurations.
- Wall and helium coolant are both time-integrated states.
- Hot gas remains quasi-steady with a runtime audit.
- Co-flow and counter-flow both work without counter-flow profile relaxation in
  transient-coolant mode.
- Zero-flow helium is physically meaningful and stable.
- Energy accounting includes both wall and coolant internal energy.
- Constant-BC cases settle to corresponding steady solvers.
- Bang-bang helium produces delayed outlet response.
- Time-step and grid convergence checks are documented.
- Result packaging includes the new coolant state and energy diagnostics.
- Existing quasi-steady transient and steady solver paths still work.

Transient hot gas is the final optional extension after this definition of done.

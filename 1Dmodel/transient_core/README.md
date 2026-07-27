# transient_core

This package is the next-generation transient solver core for Combustor-HX.

Current implemented layer:

```text
adapters_helical.py
adapters_shelltube.py
coolant_fv.py
diagnostics.py
grid.py
integrator.py
schedules.py
state.py
wall_coolant.py
```

It provides a geometry-independent implicit-upwind finite-volume update for the
helium coolant temperature state:

```text
rho * V * cp * dT_i/dt = mdot * cp * (T_upwind - T_i) + Q_i
```

where `Q_i` is the total heat rate into coolant cell `i`.

`grid.py` provides `AxialGrid`, a geometry-neutral description of cell edges,
centers, lengths, coolant/wall volumes, perimeters, and coolant inlet/outlet
indices for either flow direction.

Important behavior:

- `flow_direction = +1`: inlet at cell 0.
- `flow_direction = -1`: inlet at the last cell.
- `abs(mdot) <= mdot_floor`: advection is disabled, so cells behave as stagnant
  thermal inventories.
- The result includes an energy residual diagnostic for constant-property
  validation.

`state.py` provides a named layout for the first wall + coolant transient state:

```text
state = [Tbar_wall[0:N], T_coolant[0:N]]
```

`wall_coolant.py` provides the first coupled wall/coolant step. It solves, per
cell, a two-equation implicit linear system:

```text
Cw_i * (Tw_i_new - Tw_i_old) / dt = Qhot_i - G_i * (Tw_i_new - Tc_i_new)

Cc_i * (Tc_i_new - Tc_i_old) / dt =
    mdot * cp_i * (Tup_i_new - Tc_i_new)
  + G_i * (Tw_i_new - Tc_i_new)
```

This layer is still geometry-neutral: helical and shell-and-tube adapters must
provide `Qhot_i`, `G_i`, cell heat capacities, coolant properties, and flow
direction.

`diagnostics.py` provides shared validation helpers:

```text
residence time = sum(rho_i * V_i) / abs(mdot)
wall tau_i     = Cw_i / (G_hot_i + G_cold_i)
energy audit   = residual / max(relevant energy terms)
timescale audit for coolant/hot-gas quasi-steady assumptions
```

`adapters_helical.py` provides the first helical bridge:

```text
legacy helical geometry -> AxialGrid
flow_config             -> flow_direction
rho, cp, grid           -> coolant/wall heat capacities
h_c, grid               -> wall-to-coolant conductance
rho, mu, k, cp, mdot    -> U, Re, Pr, f, Nu, h_c, G
Tbar, T_c, T_g, h_g,h_c -> T_wg, T_wc, Qhot_i, Qcold_i
```

The convention is explicit: `AxialGrid` stores total areas and total perimeters
across all parallel coils. The coolant-film bridge mirrors the legacy helical
solver's friction/Nusselt dispatcher calls, but receives thermophysical
properties as arrays instead of calling CoolProp internally.
The wall-flux bridge wraps the existing `fluxes_at_Tbar()` reconstruction and
converts single-tube `W/m` outputs into total `W` per transient-core cell.

`adapters_shelltube.py` provides the first shell-and-tube bridge:

```text
EchTherm shell/tube geometry -> AxialGrid
flow_config                  -> shell-side flow_direction
shell bore - tube OD volume  -> shell-side coolant hold-up
total tube annulus volume    -> wall heat capacity
h_shell, grid                -> wall-to-coolant conductance
rho, mu, k, cp, mdot         -> G_s, Re_s, Pr_s, h_shell, G
rho_g, mu_g, k_g, cp_g, mdot -> U_g, Re_g, Pr_g, f_g, Nu_g, h_g, G_hot
Tbar, T_c, T_g, h_g,h_shell  -> T_wg, T_wc, Qhot_i, Qcold_i
gas state provider           -> sequential h_removed/Yc tube march
FPV/equilibrium/oxygen       -> gas-state provider callables
assembled step inputs        -> WallCoolantStepInputs + diagnostics
adapter-level run            -> transient wall/coolant history + diagnostics
```

The current shell-side hold-up is a clear baseline: shell cylinder volume minus
displaced tube outer volume, uniformly distributed over tube length. Baffle
window/leakage path refinements belong in a later residence-volume refinement.
The shell-film bridge mirrors the maintained Bell-Delaware call convention but
receives thermophysical properties as arrays instead of calling CoolProp.
The wall-flux bridge wraps the existing `fluxes_at_Tbar()` reconstruction with
`hot_side="inner"` and converts representative-tube `W/m` outputs into total
`W` per transient-core cell using `N_tubes*dx`.
The tube-side gas-film bridge mirrors the maintained smooth/grooved tube
correlation selection using supplied gas-property arrays and divides total
hot-side mass flow by `N_tubes`.
The hot-gas march bridge advances representative-tube enthalpy removal and an
optional progress variable from an injected gas-state provider. Provider
helpers wrap FPV finite-rate manifolds, equilibrium/frozen manifolds, and
pre-ignition oxygen sensible-cooling properties.
`shelltube_step_inputs()` assembles coolant properties, Bell-Delaware shell
film, hot-gas march, wall/coolant heat capacities, and the
`WallCoolantStepInputs` consumed by the generic fixed-step integrator.
`run_shelltube_transient_core()` builds the fixed time grid with schedule
breakpoints, runs the selected wall/coolant integrator, and preserves per-step
shell-side/hot-gas diagnostics for packaging.
It accepts time-dependent hot-side provider and mass-flow callbacks, so a run
can switch from no-hot-flow or GOX sensible chilldown to FPV/equilibrium
combustion at ignition schedule breakpoints.

For production `transient_coolant` runs, the selected coolant state model is:

```text
coolant_state_model = "mass_energy"
state = [T_wall_i, m_coolant_i, U_coolant_i]
```

The shell-and-tube adapter and the helical `solve_transient_core()` reconstruct
helium `T, p, rho, h, cp, mu, k` from `(m, U, V)` with CoolProp, treat positive
scheduled helium flow as the inlet command, close the inlet face when the
schedule goes to zero, and leave internal/outlet faces pressure-driven. This
captures residual helium discharge after a bang-bang valve closure while
retaining quasi-steady momentum for now.

Because coolant advection is explicit in the conserved variables, the
mass/energy path internally refines the timestep by a cell-inventory CFL rule
even when `transientProp.max_step` is larger. It also caps face flows at a
finite valve/line-capacity scale tied to the scheduled mass-flow range; replace
that closure with real valve/supply/outlet pressure data when available.

Helical note: this path satisfies the current quasi-steady-momentum
implementation goal, but previous bang-bang momentum estimates showed helical
line inertance can be large for aggressive valve events. Revisit transient
momentum before treating high-frequency helical bang-bang predictions as final.

`integrator.py` provides a generic fixed-step driver. A geometry adapter supplies
one `WallCoolantStepInputs` object per time interval; the integrator stores wall
and coolant temperature histories, outlet temperature, heat transfer, and energy
residuals. It also provides `fixed_time_grid()` for inserting schedule
breakpoints and requested output times into a bounded-step grid.

`schedules.py` provides shared schedule helpers:

```text
interp_schedule(schedule, t, default)
schedule_times(...)
collect_transient_schedule_times(transient, names, ...)
```

The interpolation behavior matches the legacy transient solvers: linear between
points and flat-held outside the schedule range.

Shell-and-tube production dispatch is wired through:

```text
transientProp.fluid_model = "transient_coolant"
```

in `input_data.py`, then run with `python -m hps_combustor.main_transient`.
The legacy wall-only transient remains the default with
`fluid_model = "quasi_steady"`.

The full staged plan is documented in:

```text
docs/solver_design/TRANSIENT_CORE_IMPLEMENTATION_PLAN.md
```

Next layers: helical transient-coolant hot-gas wrapper and dispatch, dashboard
field polish, broader short-run validation cases, and shell-side
residence-volume refinement.

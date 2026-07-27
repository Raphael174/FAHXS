# Transient Solver Status

Last consolidated from project memory on 2026-07-08.

Implementation update on 2026-07-11:

- `1Dmodel/validation/coupled_bangbang_hx.py` now supports both
  `hx_config = "shellntube"` and `hx_config = "shellnHelicalTube"` with the
  same 0D bang-bang tank/feed/valve boundary schedule. The helical coupled
  preset applies 316L, user-level tube ID/wall/pipe-length/helix-diameter
  geometry fields, Salimpour/Mori with the existing 0.28 shell correction, and
  `transientProp.coolant_momentum_model = "low_mach"`.
- Helical `solve_transient_core()` now has a pressure-driven low-Mach coolant
  mode for coupled validation. It advances distributed wall temperature plus
  distributed helium mass/internal energy, while the helical momentum closure is
  a single implicit lumped pipe through-flow state using curved-pipe Darcy
  resistance and inertance. A fully distributed face-pressure low-Mach model was
  probed but is too stiff as an explicit line solver; it needs an implicit
  pressure-projection solve before it can be a production option.
- Coupled helical sizing skips the expensive steady-reference probe and uses
  the known helical pipe length directly for the transient grid. Normal
  standalone helical transients still keep the steady comparison unless
  `transientProp.skip_steady_reference_probe = True`.
- New sweep runner:
  `python -m hps_combustor.validation.coupled_bangbang_hx_geometry_sweep`.
  The latest helical-focused 1 s finite-rate low-Mach sweep is:
  `docs/validation/coupled_bangbang_helical_geometry_sweep_1s_pressure_metrics/SWEEP_SUMMARY.md`.
  It compared three 316L helical candidates: 13.5 mm ID / 1.0 mm wall /
  60 mm helix, 14.0 mm ID / 1.25 mm wall / 65 mm helix, and 14.5 mm ID /
  1.5 mm wall / 70 mm helix, all with 12 m pipe length. All cases met the
  system-level 30 L/s and 70 bar water-tank targets over the 1 s screen and
  gave about 0.143 kg/s solved helium flow. The flowing coolant pressure-span
  statistics were about 9.53 bar mean, 10.95 bar p95, and 11.14 bar max for all
  three candidates. This means the helical pressure loss is in the desired
  5-10 bar band on average, with the p95/max slightly above 10 bar because of
  the imposed bang-bang pressure ripple. Helical 1 s wall peaks were about
  408-443 K, and the highest early-screen helium outlet was about 467 K for the
  13.5 mm ID / 1.0 mm wall / 60 mm helix candidate. This is not yet a 100 s
  design acceptance run.

Planning update on 2026-07-09: bang-bang helium operation invalidates the
quasi-steady-coolant assumption during low-flow and zero-flow intervals. The
current maintained transient solvers remain useful as fast/reference modes, but
the next production transient architecture should be the finite-volume
`transient_core` described in
`docs/solver_design/TRANSIENT_CORE_IMPLEMENTATION_PLAN.md`: transient wall +
transient helium coolant, with hot gas kept quasi-steady initially because the
hot-side residence time is about 0.6 ms in the current regime.

Implementation update on 2026-07-10: shell-and-tube and helical
`transient_coolant` now advance helium mass and internal energy, not only
temperature. The detailed
plan/current assumptions live in
`docs/solver_design/TRANSIENT_COOLANT_MASS_ENERGY_PLAN.md`.
`1Dmodel/transient_core/compressible_coolant.py` provides conservative
mass/energy primitives, CoolProp reconstruction from `(m, U, V)`, and a
quasi-steady pressure-flow helper. `wall_compressible_coolant.py` provides the
semi-implicit wall + compressible-coolant step. `adapters_shelltube.py` now
wires these into `run_shelltube_transient_core(...,
coolant_state_model="mass_energy")`, and
`main_solve_shellntube_transient.py` selects that mode for production
`transientProp.fluid_model = "transient_coolant"` runs.
`main_solve_transient.py` now provides `transient_solver.solve_transient_core()`
for the helical configuration, reusing the existing helical `fluid_pass()` as
the quasi-steady hot-side/wall-flux evaluator against the conserved helium
temperature profile. `main_transient.py` dispatches both HX configurations to
their core path when `fluid_model="transient_coolant"`.

Implementation update on 2026-07-09: `1Dmodel/transient_core/` now contains the
first geometry-independent numerical kernels and the shell-and-tube
`transient_coolant` production dispatch:

- `coolant_fv.py`: implicit-upwind transient helium inventory and advection.
- `diagnostics.py`: energy residual scaling, residence time, wall time
  constant, and quasi-steady timescale audits.
- `grid.py`: axial cell geometry, volumes, perimeters, and inlet/outlet indices.
- `integrator.py`: generic fixed-step wall/coolant history driver using adapter
  supplied step inputs.
- `schedules.py`: shared schedule interpolation and breakpoint-time extraction.
- `state.py`: explicit state layout `[Tbar_wall, T_coolant]`.
- `wall_coolant.py`: implicit local wall/coolant exchange coupled to
  implicit-upwind coolant advection.
- `compressible_coolant.py`: conservative coolant mass/energy finite-volume
  step, CoolProp `(m,U,V)` state reconstruction, and quasi-steady face-flow
  helper. This is the new path for fixing bang-bang residual outflow.
- `wall_compressible_coolant.py`: semi-implicit wall energy plus conservative
  compressible-coolant mass/energy update. It conserves total wall+coolant
  energy and allows residual outlet flow after inlet closure.
- `adapters_helical.py`: first helical geometry/inventory/coolant-film bridge
  into `AxialGrid`, heat capacities, coolant-side conductance, and the existing
  coil friction/Nusselt dispatchers; also wraps wall face reconstruction and
  converts single-tube per-length fluxes to total per-cell heat rates.
- `adapters_shelltube.py`: first shell-and-tube geometry/inventory bridge into
  `AxialGrid`, shell-side hold-up, tube-wall heat capacity, and shell-side
  conductance; also wraps the Bell-Delaware shell-side film call using supplied
  property arrays, the tube-side smooth/grooved gas-film correlations using
  supplied gas-property arrays, plus wall face reconstruction with
  `hot_side="inner"` and representative-tube-to-total heat-rate conversion,
  plus a sequential representative-tube hot-gas march scaffold driven by an
  injected gas-state provider, with tested FPV, equilibrium/frozen, and
  pre-ignition oxygen provider wrappers; it now also assembles shell-and-tube
  `WallCoolantStepInputs` and provides adapter-level fixed-step runners for the
  temperature-only and mass/energy coolant models. The mass/energy path uses a
  scheduled inlet command, pressure-driven internal/outlet faces, a finite
  face-flow cap tied to the scheduled valve flow, and internal CFL refinement.
  The runner supports per-step hot-side property providers and hot-side
  mass-flow callbacks, so no-hot-flow, GOX chilldown, and combustion phases can
  share one transient run.
- `main_solve_shellntube_transient.py`: `solve_transient_core()` bridges the
  legacy shell-and-tube setup into `transient_core`, stores `core_result`, and
  emits the normal packaged `time_series` dictionary.
- `main_solve_transient.py`: helical `solve_transient_core()` advances
  conserved helium mass/internal energy with quasi-steady momentum while keeping
  hot gas quasi-steady through the existing helical march.
- `main_transient.py`: both shell-and-tube and helical runs use the transient
  fluid/material core when `transientProp.fluid_model = "transient_coolant"`;
  the default `"quasi_steady"` path remains the legacy wall-only transient.

Focused tests in `tests/test_transient_core_coolant_fv.py` cover zero-flow
soak, flow direction, state packing, constant-property first-law residuals,
shell-and-tube adapter assembly, adapter-level runner scheduling, CoolProp
mass/energy state round trips, coupled wall/compressible-coolant conservation,
and shell-and-tube residual outlet discharge after inlet closure. Dispatch tests
also cover helical `transient_coolant` using `solve_transient_core()` and
helical `quasi_steady` using the legacy path. Smoke/audit runs completed for
shell-and-tube and helical finite-rate production paths. The helical path still
may need transient momentum after final valve/feed data because the inertial
pressure scale can be large in the small-bore coil.

Short-run validation note:

- `docs/validation/TRANSIENT_CORE_SHORT_RUNS_2026-07-09.md` records two 2 s
  shell-and-tube `transient_coolant` runs at 20 nodes:
  finite-rate bang-bang helium, and GOX chilldown -> ignition -> hot-gas ramp.
- Both ran in about 1.5 s in memory; the packaged GOX/ignition run took about
  2.0 s and reloaded with 10 time points, 20 axial nodes, and dashboard-style
  wall/coolant/gas fields.
- Repeatable runner:
  `python -m hps_combustor.validation.transient_core_short_runs`. It writes
  `docs/validation/transient_core_short_run_results.json`.
- Latest runner result: warm-cache 2 s cases ran in about 0.8-1.2 s at 20 nodes
  and about 1.1 s at 40 nodes. Bang-bang helium timestep refinement from
  `max_step=0.25` to `0.10 s` at 20 nodes changed final `T_c_out` by about
  3.5 K and final `T_wall_max` by about 13.4 K. GOX/ignition grid refinement
  from 20 to 40 nodes changed final `T_c_out` by about 13.9 K, final
  `T_wall_max` by about 74.2 K, and final FPV progress strongly, so the
  GOX/ignition case is not yet grid converged.
- Packaged `transient_coolant` time series now include dashboard-style
  engineering fields: gas velocity/Re/Pr/Nu/friction/dp, shell mass flux/Re/Pr
  and pressure-drop estimate, coolant `rho/mu/k/cp`, gas enthalpy removed, and
  FPV progress variable. Mass/energy runs also save `coolant_mass_kg`,
  `coolant_internal_energy_J`, `p_c`, `rho_c_state`, `h_c_state`,
  `face_mdot_c`, and `coolant_mass_residual_kg`.
- Production `transient_coolant` runs now print throttled terminal progress via
  `transient_core/progress.py`: step/time, material min/max temperature,
  helium outlet temperature, pressure diagnostic, and hot-gas outlet
  temperature. In shell-and-tube quasi-steady momentum mode the pressure
  diagnostic is `dpHe [bar]`, the hydraulic shell-side pressure-drop estimate;
  in low-Mach momentum mode it is the reconstructed helium outlet pressure.
  Controls live in `transientProp.progress_print`, `progress_interval_steps`,
  and `progress_interval_time_s`.
- Transient zip folders still save one raw numeric transient dataset:
  `data/transient_timeseries.npz`. Packaging now also writes
  `HX_performance_summary.txt`, a text engineering report with multi-time-point
  duty, LMTD, coolant-capacity-referenced NTU/effectiveness, pressure drops,
  temperature extrema, and any available combustion/mechanical histories. No
  separate engineering transient `.npz` is written.
- Exact zero helium flow now uses a stagnant-shell fallback
  `h_shell = k_He / D_tube_outer`, `dp_shell = 0` instead of Bell-Delaware.
  This is a numerical/physical fallback requiring later low-flow validation.
- Bang-bang momentum audit:
  `python -m hps_combustor.validation.bangbang_momentum_audit <schedule>`.
  For the provided dummy helium file, `max |dmdot/dt| ~= 30.9 kg/s2`; the
  helical inertial pressure scale is about 193 bar, while the EchTherm-scale
  shell-and-tube shell-side scale is about 0.0097 bar. Conclusion: shell-and-
  tube can start with transient mass/energy + quasi-steady momentum; helical may
  require transient momentum/line inertance after final valve/feed data.
- First shell-and-tube low-Mach momentum mode is implemented behind
  `transientProp.coolant_momentum_model = "low_mach"`. The maintained coupled
  validation path now uses a pressure-driven lumped shell-side through-flow
  momentum state for the baffled shell path, with cell-resolved coolant
  mass/energy and wall temperatures. The 0D valve/tank model provides a
  time-local available-flow cap; the HX does not use that cap as a prescribed
  operating flow. The older cellwise face-momentum helper remains tested but is
  not the accepted coupled validation path because it produced nonphysical
  internal pressure oscillations for the bang-bang case. The dummy Excel 20 s
  validation runner is
  `python -m hps_combustor.validation.low_mach_momentum_dummy_run`. The 70 bar
  downstream test completed, but is not physically accepted: the dummy inlet
  pressure is below 70 bar for about one third of the first 20 s, so the model
  predicts weak/reversing coolant flow and hits temperature clamps. Details are
  in `docs/context/LOW_MACH_COOLANT_MOMENTUM_PLAN.md`.
- A separate dummy pressurant-system sizing surrogate exists at
  `python -m hps_combustor.validation.pressurant_bangbang_sizing`. It models a
  finite 400 bar, 100 K, 265 L helium tank with adiabatic ideal-gas blowdown,
  2/3 parallel calibrated orifice/valve branches, a nominal 80 bar pre-HX line,
  a simplified quadratic HX/feed pressure loss, and a 3000 L water tank drained
  through an exit orifice. Latest sweep output is in
  `docs/validation/pressurant_bangbang/`: best coarse candidate is 3 branches,
  3.0 mm branch orifices, 3.0 mm valve equivalent diameter,
  `Kv ~= 0.288 m3/h` per valve, 40 Hz staged bang-bang. It gives about
  69.3 bar mean water-tank pressure, 78.6 bar mean pre-HX line pressure,
  29.84 L/s water flow, and 2985 L delivered over 100 s. Final helium supply is
  about 84.5 bar and 53.7 K, so the 265 L tank is marginal rather than generous.
- A separate coupled validation sandbox now exists at
  `python -m hps_combustor.validation.coupled_bangbang_hx`. It is deliberately
  not wired into `main_transient.py`: the external system remains a cheap 0D
  pressurant/feed/water-tank surrogate, and the shell-and-tube HX remains the
  detailed 1D transient coolant/material solver. The module writes
  `docs/validation/coupled_bangbang_hx/settings_used.json`,
  `summary.json`,
  `system_timeseries.csv`, `coupled_timeseries.csv`,
  `hx_boundary_schedule.csv`, `water_outlet_orifice_sweep.json`,
  `hx_transient_timeseries.npz`, and `coupled_dashboard.html`. The default
  user-editable dummy preset is
  `inputs/coupled_bangbang_hx_dummy_validation.json`; run it with
  `python -m hps_combustor.validation.coupled_bangbang_hx --settings inputs/coupled_bangbang_hx_dummy_validation.json`.
  Previous accepted dummy baseline: 100 s, 10 HX nodes, `max_step = 0.01 s`,
  counter-flow, finite-rate chemistry, quasi-steady coolant momentum,
  EchTherm-style shell-and-tube geometry, 90 K / 265 L / 400 bar helium supply,
  70 bar water tank target, 30 L/s water discharge, three 1.75 mm helium
  valve/orifice branches at 40 Hz, and hot-gas tracking of helium demand at
  O/F = 2 capped at 72.5 g/s total hot-side flow (48.33 g/s LOX plus
  24.17 g/s diesel). The active engineering targets are wall temperature below
  1000-1200 degC, flowing helium below about 700 K, and helium-side hydraulic
  pressure loss no more than 5-10 bar.
  In quasi-steady momentum mode, use `dp_shell_total_Pa` /
  `hx_max_shell_pressure_drop_estimate_bar` for hydraulic pressure loss; the
  saved `p_c` span is a thermodynamic coolant-state reconstruction and is not a
  pressure-drop prediction in this mode.
  Accepted 100 s detailed coupled evidence is in
  `docs/validation/coupled_bangbang_hx_100s_hotgas72p5_pressure_floor/`.
  The HX portion completed in about 7.8 min and triggered no coupled sanity
  flags. System targets were met: 29.99 L/s mean water flow, 69.95 bar mean
  water-tank pressure, 0.142 kg/s mean helium flow, and 14.19 kg helium used
  over 100 s. Flowing helium outlet temperature was about 596 K mean and
  693 K max. Hydraulic pressure-drop estimates stayed low: helium shell-side
  below 0.83 bar and hot-gas tube-side below 0.14 bar. Wall-depth peaks from
  the stored fields were about 894 K hot face, 880 K mean wall, and 880 K cold
  face near `x = 0.01175 m`. Coolant temperature stayed away from the former
  60/2500 K validity limiter (`T_c_min ~= 90 K`, `T_c_max ~= 700 K`).
  Production transient coolant paths now call the CoolProp internal-energy
  bounds check in strict mode, so a future out-of-range state raises instead of
  being silently clipped.
  Current pressure-driven low-Mach coupled short validation is in
  `docs/validation/coupled_bangbang_hx_lowmach_5s_dt5ms_n5_hot70/`.
  Settings: 5 s, 5 HX nodes, `max_step = 0.005 s`, counter-flow,
  finite-rate chemistry, low-Mach pressure-driven shell-side momentum, 90 K /
  265 L / 400 bar helium supply, 70 bar water tank target, 30 L/s water
  discharge, three 1.75 mm helium valve/orifice branches at 40 Hz, and hot-gas
  tracking at O/F = 2 capped at 70 g/s. The HX portion ran in about 24 s, so
  the linear 100 s estimate is about 8 min on this machine. The 0D system hit
  69.98 bar mean water-tank pressure, 78.33 bar mean pre-HX line pressure,
  29.995 L/s water outflow, and 0.145 kg/s mean helium flow. The resolved HX
  pressure field stayed about 69.8-84.0 bar. Flowing helium outlet temperature
  peaked at about 705 K, slightly above the 700 K target; mean flowing outlet
  temperature was about 534 K. Wall peaks were about 796 K hot face, 764 K mean
  wall, and 749 K cold face. Hot-gas pressure drop stayed about 0.13 bar.
  Detailed Bell-Delaware shell-side pressure-drop estimate peaked at about
  3.13 bar, below the 5-10 bar target, so the next hydraulic task is to
  reconcile shell-side resistance/line losses before treating a long low-Mach
  run as final design evidence.
- Bang-bang coolant behavior audit:
  `python -m hps_combustor.validation.shelltube_bangbang_coolant_audit`.
  Latest counter-flow shell-and-tube result writes
  `docs/validation/shelltube_bangbang_coolant_audit.json` and passes all checks:
  inlet open before shutoff, inlet closed after shutoff, residual outlet flow
  after shutoff, and coolant inventory decrease after shutoff. The latest run
  used 15 internal steps over 0.02 s, `max_internal_dt_s ~= 1.43e-3`, and had
  final mass/energy residuals near roundoff.
- Helical bang-bang coolant behavior audit:
  `python -m hps_combustor.validation.helical_bangbang_coolant_audit`.
  Latest co-flow helical result writes
  `docs/validation/helical_bangbang_coolant_audit.json` and passes all checks:
  inlet open before shutoff, inlet closed after shutoff, residual outlet flow
  after shutoff, and coolant inventory decrease after shutoff. The latest run
  used 3 internal steps over 0.002 s, `max_internal_dt_s = 1e-3`, and had final
  mass/energy residuals near roundoff.

Performance update on 2026-07-10:

- The first strict explicit conserved-energy coolant update forced per-cell CFL
  timesteps and was far too slow for 100 s runs.
- Production now uses a bounded-cost semi-implicit thermal update with
  conserved coolant mass/face-flow diagnostics. Open-valve intervals carry the
  scheduled flow through all faces; closed-valve intervals use one-way residual
  outlet discharge with short memory from the previous through-flow.
- Shell-and-tube finite-rate counter-flow benchmark, 80 nodes, 5 s simulated,
  `max_step=0.25 s`: 5.66 s wall time, 21 steps. Linear extrapolation to 100 s
  is about 113 s, below the 10 min target.
- This is a pragmatic production model, not a fully implicit compressible-flow
  PDE. Revisit a true implicit mass/energy/momentum solve only if validation
  data show this bounded-cost model is insufficient.
- Low-Mach coolant momentum plan saved separately in
  `docs/context/LOW_MACH_COOLANT_MOMENTUM_PLAN.md`. The decision is to avoid
  full acoustic momentum for routine runs and, if needed later, add a
  low-Mach inertance/friction face-momentum model coupled to current conserved
  coolant mass/energy.

## Implemented Solvers

| Solver | File | Class | Status |
|---|---|---|---|
| Helical transient | `main_solve_transient.py` | `transient_solver` | Implemented and working. |
| Shell-and-tube transient | `main_solve_shellntube_transient.py` | `shellntube_transient_solver` | Implemented; broader validation was still in progress in the source memory. |

Both transient solvers support `flow_config in ("co", "counter")`. Counter-flow
uses a profile-relaxation pass inside each RHS evaluation; low-flow startup
falls back only in the near-stagnant regime already flagged as unreliable.

Helical counter-flow validation note:

- The legacy steady helical counter-flow march prescribes `coolantProp.T_out`
  at the gas-inlet end, while the transient solver's physical boundary is
  `coolantProp.T_in` at the gas-outlet end.
- `solve_counterflow_physical_reference()` in `main_solve.py` provides an
  opt-in steady reference for settle checks by shooting the hot-end helium
  temperature until the cold-end inlet matches `coolantProp.T_in`.
- A reduced frozen-chemistry check converged the steady reference boundary to
  0.33 K and matched transient-pass heat rate within about 0.6%; full
  finite-rate settle-to-steady validation is still open.

## Performance State

Helical transient performance improved from about 3 hours at first cut to minutes per 100 s run:

- For helical runs, BDF has been useful because RK45 takes too many small steps
  when the wall time constant remains fast at full helium flow.
- Radiation was once about 72% of runtime due to repeated WSGGM calls inside nested loops.
- Radiation was accelerated by tabulating emissivity, collapsing nested radiation/Nusselt loops into one fixed point, and using a closed-form 2x2 face-temperature solve.
- Equilibrium chemistry is now table-driven, so the transient march does not call Cantera at each node.
- Memory benchmark: 100 s equilibrium helical run was about 5.7 min and under the user's target of 10 min.

Shell-and-tube counter-flow performance decision:

- Do not use BDF as the production default for shell-and-tube counter-flow long
  runs. One RHS call includes quasi-steady tube/shell marching plus counter-flow
  profile relaxation, so BDF Jacobian probing can multiply cost severely around
  ignition or schedule discontinuities.
- `transientProp.solver_method = "fixed_step"` is now the bounded-cost
  production path. It inserts schedule breakpoints into the time grid and uses
  a linearly-implicit wall update, so it costs roughly one fluid pass per time
  step without the explicit Euler helium-film instability.
- Benchmark after this change, with cached finite-rate chemistry, shell-and-tube
  counter-flow, 80 nodes, 5 s simulated: 5.8 s wall time and 26 fluid passes.
  The prior BDF run on the same 80-node/5 s case took about 134 s and 456 RHS
  calls.
- Finite-rate FPV manifolds are cached in `cache/fpv_manifolds`. First build of
  the current reduced transient table was about 33 s; cached setup was about
  0.9 s. Do not disable this cache for normal work.

## Chemistry Status

Critical user decision from 2026-07-08: frozen chemistry is physically wrong for the main regime.

Regime: diesel-C16H34/O2 with low hot-gas mass flow, high extracted heat, and deep cooldown through recombination. Finite-rate FPV is the required default because rapid per-mass cooling can cause recombination freeze-out. Equilibrium remains a useful comparison case.

Implemented modes:

- `finite_rate`: default C1 FPV manifold using `physics/combustion_chemistry/fpv_manifold.py`; also used by steady solvers.
- `equilibrium`: C0 equilibrium manifold vs enthalpy removed. Implemented in helical transient through `_build_chem_tables()`, `_gas_at()`, `_eps_at()`, and `_h_g_rad()` style table lookups.
- `frozen`: validation-only, not the default physical model for the current user regime.

FPV details to preserve:

- Progress variable: unnormalized `Yc = Y_CO2 + Y_H2O - Y_CO`.
- Transport: `dYc/dx = omega_Yc / U_g`.
- State table: `(h, c)` built from constant-enthalpy reactor relaxations.
- The `c=1` column must be anchored to the exact equilibrium state from `gas.equilibrate("HP")` with `omega=0`; do not extrapolate from a stiff reactor trajectory at cold enthalpy.

Ignition decision:

- Use a pilot diesel/O2 low-mass-flow segment in `schedule_mass_flow_g`.
- No PLA surrogate mechanism is wanted.
- One diesel/O2 manifold covers the full 100 s mass-flow ramp because the manifold is per unit mass; only O/F or pressure changes force rebuild.
- Before ignition, scheduled LOX/GOX flow is now modeled as Oxygen with CoolProp
  properties and sensible heat exchange. Once `ignition` becomes 1, the solver
  switches to the combustion manifold.

## Numerical Fixes Already Found

- Use hot face temperature `T_wg`, not mean wall temperature, in the hot-gas Kays-Crawford correction.
- Do not verify steady transient fluxes by iterating `mean(T_wg, T_wc)`; use actual `solve_ivp` settling or a root solve on flux imbalance.
- Cache inlet Cantera `T/p/Y` once and reset from the cache on each repeated sweep. A persistent Cantera object otherwise keeps the previous cooled state.
- In shell-and-tube transient, hot-gas mass flow for tube-side velocity and enthalpy accounting must be total mass flow divided by `N_tubes`.
- Clamp shell-and-tube transient coolant temperatures during RHS probing to keep `solve_ivp` from evaluating wildly unphysical wall fields.
- `Solve1Dconduction()` must honor `hot_side="inner"` for shell-and-tube. A bug
  where the steady conduction solve always treated the outer perimeter as hot
  caused 10-25% steady/transient mismatch; after fixing the hot/cold perimeters,
  shell-and-tube transient wall reconstruction matches steady heat rate within
  about 0.4% or better for co/counter and frozen/equilibrium/finite-rate at
  16 nodes.

## Shell-And-Tube Status

Implemented and validated components from the memory:

- `physics/bell_delaware.py`: Bell-Delaware shell-side model.
- `mechanical/geometry/shelltube_geometry.py`: EchTherm-style geometry translation.
- `dispatch_nu_tube_straight()` and `dispatch_friction_tube_straight()`: blended laminar/transitional/turbulent straight-tube correlations.
- `shellTubeProp`: shell-and-tube configuration dataclass.
- `main_solve_shellntube.py`: steady sweep-iteration solver with both fluid inlets prescribed.
- `physics/combustion_chemistry/gas_manifold.py`: shared equilibrium manifold module for transient chemistry.
- `main_solve_shellntube_transient.py`: shell-and-tube transient solver.

Validation remembered:

- Steady shell-and-tube grid convergence was clean: `Q_tot` about 422.8 to 416.3 kW from `N=25` to `N=200`.
- Steady shell-and-tube energy balance error was about 0.048%.
- Co-flow vs counter-flow qualitative behavior was correct in the steady solver.
- Shell-and-tube transient energy consistency after the per-tube mass-flow fix: `Q=345.9 kW` vs coolant enthalpy estimate about `345 kW`.

Material fix already made:

- `mechanical/material_specs/material_temperature_strength.py` had 1000x fallback unit typos in the high-temperature branches of both `ElasticityModulus_316331600` and `ElasticityModulus_INCO718`. They were fixed from MPa-scale to Pa-scale values.

## Remaining Work From Memory

- Before the prior Claude session stopped, two background agents were still running:
  the low-flow freeze-out demo for the helical configuration, and the full
  shell-and-tube transient validation. Treat their results as unknown until the
  output/logs are found or the runs are repeated.
- Visually verify the dynamic HTML dashboard in `model_data_process/data_plotting_transient.py`.
- Complete broader shell-and-tube transient validation: settle-to-steady, chemistry mode comparisons, and 100 s benchmark.
- Finish full finite-rate helical counter-flow settle-to-steady validation using
  the physical steady reference helper.

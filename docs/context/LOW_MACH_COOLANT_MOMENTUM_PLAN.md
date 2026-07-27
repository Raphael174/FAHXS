# Low-Mach Coolant Momentum Plan

Context saved on 2026-07-10 before compaction.

## Decision

Do **not** implement full compressible/acoustic coolant momentum for routine
Combustor-HX transients. It would force acoustic CFL timesteps:

```text
dt_acoustic ~ dx / a_He
```

For shell-and-tube with `L ~= 0.235 m`, `N ~= 80`, and helium sound speed around
`550-600 m/s`, this gives `dt ~ 5 us`, or about 20 million steps for a 100 s
run. For the 6 m helical tube, the acoustic step is larger but still leads to
roughly 0.8-2 million steps for 100 s depending on axial resolution. That is
incompatible with the project runtime target:

```text
100 s simulated < 10 min wall time
```

The intended next upgrade is a **low-Mach inertance/friction coolant momentum
model** coupled to the existing transient coolant mass/energy and transient wall
model. It should capture valve/line inertial flow memory without resolving
pressure waves.

## Physical Model

Current production model:

```text
state = [T_wall_i, m_coolant_i, U_coolant_i]
mdot_face = quasi-steady pressure/friction closure
```

Upgrade target:

```text
state = [T_wall_i, m_coolant_i, U_coolant_i, mdot_face_j]
```

with a face or compartment momentum ODE:

```text
I_j d(mdot_j)/dt =
    p_left - p_right
  - R_j mdot_j |mdot_j|
  - K_j mdot_j |mdot_j|
  - valve_or_boundary_loss_j
```

where:

- `I_j` is hydraulic inertance, approximately `L_j / A_j` for the flow path
  segment, with units consistent with `Delta_p = I dmdot/dt`;
- `R_j mdot |mdot|` represents Darcy-Weisbach distributed friction converted to
  pressure loss;
- `K_j mdot |mdot|` represents local/minor losses, including baffle/window,
  entrance/exit, nozzle, bend, and valve losses where appropriate;
- pressures `p_i` are reconstructed from coolant cell mass/internal energy:

```text
rho_i = m_i / V_i
u_i   = U_i / m_i
p_i, T_i, h_i, cp_i, mu_i, k_i = CoolProp(rho_i, u_i)
```

This is a low-Mach model: it keeps pressure/inertance dynamics and residual
flow after valve closure, but it does not attempt to resolve acoustic waves.

## Boundary Conditions

Keep scheduled helium mass flow as the primary user input for now. For
transition compatibility:

- open-valve intervals should relax inlet-face `mdot` toward the scheduled
  command through a finite valve/feed inertance or relaxation time;
- closed-valve intervals should set inlet valve conductance to zero or near
  zero, but internal/outlet face momenta should continue evolving;
- outlet can be pressure boundary, tank/backpressure boundary, or calibrated
  restriction depending on available hardware data;
- if no upstream/downstream hardware data exist, use calibrated `I`, `R`, `K`
  values that reproduce the commanded steady-flow pressure drop at nominal
  flow.

Do not force all internal faces to the scheduled inlet flow during closed-valve
intervals. The point of this upgrade is that helium already in the exchanger can
continue to exit after a bang-bang inlet command drops to zero.

## Geometry-Specific Guidance

### Shell-And-Tube

Priority is lower. Prior momentum audit for the EchTherm-scale shell side gave:

```text
Delta_p_inertia ~ 0.0097 bar
```

for the provided aggressive dummy helium schedule. This is small compared with
typical operating pressure, so transient coolant mass/energy with quasi-steady
momentum remains defensible for first production results.

If implemented for shell-and-tube:

- lump momentum by baffle compartment or axial cell faces;
- build `R_j` from the existing Bell-Delaware pressure-drop decomposition:
  cross-flow, window, entrance/exit, leakage/bypass corrections;
- distribute the whole-shell steady pressure drop over the active faces/zones;
- keep hot gas quasi-steady.

### Helical

Priority is higher. Prior audit for a 6 m, 3.5 mm ID tube and the supplied
aggressive bang-bang derivative gave:

```text
Delta_p_inertia ~ 193 bar
```

This is large enough that quasi-steady momentum may be physically wrong if the
upstream valve/feed system really imposes such fast mass-flow changes.

If implemented for helical:

- start with a 1D face-momentum model along the coil;
- build `I_j = L_j/A_j`;
- build `R_j` from the existing curved-pipe friction model at local properties;
- add optional valve/feed-line inertance upstream if hardware data are
  available;
- validate against the prescribed bang-bang mass-flow file and pressure scale.

## Numerical Strategy

Avoid explicit acoustic timestepping. Preferred production method:

1. Evaluate schedule and boundary states at `t_n` or midpoint.
2. Reconstruct coolant thermodynamic state from `(m_i, U_i, V_i)`.
3. Advance face momenta with a linearly implicit or semi-implicit update:

```text
mdot^{n+1} = mdot^n
  + dt/I * [dp - R_eff(mdot^*) mdot^{n+1}]
```

where the nonlinear drag can be linearized with `|mdot^n|` or solved with a
small local Newton iteration per face.

4. Use `mdot_face^{n+1}` in the conservative mass/energy update.
5. Keep the wall/coolant heat exchange semi-implicit as it is now.
6. Keep hot gas quasi-steady per coolant/wall timestep.

Fallback if stiffness appears:

- use a small global implicit solve over `[p_i or m_i, mdot_j]` for the coolant
  subsystem only;
- do not let BDF probe the full combustion/wall/HX RHS repeatedly in production.

## Runtime Estimate

Expected runtime relative to current bounded-cost transient coolant model:

```text
shell-and-tube low-Mach momentum:  ~2x to 4x if face-local implicit
helical low-Mach momentum:         ~3x to 10x depending on valve/feed stiffness
```

Target:

```text
shell-and-tube 100 s finite-rate counter-flow: < 10 min
helical 100 s: aim < 10 min, accept slightly above only if physics demands it
```

Do not accept a design that requires global microsecond timesteps for routine
100 s studies.

## Code Touch Points For Later

Likely new/modified modules:

- `1Dmodel/transient_core/compressible_coolant.py`
  - add low-Mach momentum primitives;
  - keep current quasi-steady helper for validation and fallback.
- `1Dmodel/transient_core/wall_compressible_coolant.py`
  - add coupled wall + mass/energy + momentum step.
- `1Dmodel/transient_core/state.py`
  - extend state packing to optional face momentum.
- `1Dmodel/transient_core/adapters_shelltube.py`
  - build shell-side inertance/resistance by face/compartment.
- `1Dmodel/transient_core/adapters_helical.py`
  - build coil inertance/resistance by segment.
- `1Dmodel/input_data.py`
  - add a mode such as:

```text
transientProp.fluid_model = "transient_coolant_low_mach_momentum"
```

Keep `transient_coolant` as the current fast production mode until validation
shows the momentum model is needed.

## Validation Plan

Minimum validation matrix:

- constant-flow case settles to the steady solver result;
- bang-bang inlet command shows delayed outlet flow and delayed outlet
  temperature packets;
- mass and energy residuals remain small;
- momentum model tends to the quasi-steady pressure-flow solution as inertance
  tends to zero;
- pressure drops match existing steady friction/Bell-Delaware predictions at
  constant nominal flow;
- timestep convergence for bang-bang cases, e.g. `dt = 5 ms, 2.5 ms, 1.25 ms`
  where practical;
- helical inertial pressure scale matches the audit order of magnitude for the
  supplied dummy schedule.

Acceptance criteria:

- production 100 s shell-and-tube finite-rate run remains under 10 min;
- no acoustic-CFL timestep requirement in routine mode;
- no forced instantaneous zero internal flow when inlet command closes;
- saved outputs still include `face_mdot_c`, coolant mass, coolant internal
  energy, pressure, wall temperatures, and engineering summary fields.

## Implementation Update 2026-07-10

The first shell-and-tube implementation is wired as:

```python
transientProp.coolant_momentum_model = "low_mach"
transientProp.transient_coolant_outlet_pressure = 70e5
```

Internally this selects `coolant_state_model="low_mach_momentum"` in
`run_shelltube_transient_core()`. The coolant state remains conserved
`(m_i, U_i)` per cell. Face mass flows are advanced as low-Mach inertance
states:

```text
I_f (mdot_f^{n+1} - mdot_f^n) / dt =
    p_left^n - p_right^n - K_f mdot_f^{n+1} |mdot_f^{n+1}|
```

with a closed-form semi-implicit quadratic update. Boundary faces use the
scheduled helium inlet pressure and the fixed/scheduled downstream pressure.
The shell path now uses the conservative wall + compressible-coolant
mass/energy step for pressure-driven runs, rather than reconstructing coolant
energy from an intermediate temperature-only update.

Dense measured schedules can be interpolated without forcing every measurement
time into the solver grid:

```python
transientProp.insert_schedule_breakpoints = False
```

This is required for `inputs/HX_dummy_inlet.xlsx`, which is sampled every
2.5 ms.

Validation runner:

```text
python -m hps_combustor.validation.low_mach_momentum_dummy_run
```

Default case:

- Shell-and-tube, counter-flow.
- First 20 s of `inputs/HX_dummy_inlet.xlsx`.
- Helium inlet pressure/temperature from the Excel file.
- Fixed downstream pressure: 70 bar.
- LOX: 90 g/s at 100 K.
- Diesel: 30 g/s.
- Ignition always on.
- Initial wall temperature: 300 K.
- 20 axial nodes, `max_step=0.05 s`, dense schedule interpolation.

Latest 20-node run wrote:

```text
docs/validation/low_mach_momentum_dummy_run.json
```

Runtime was about 76 s. This is a completed numerical smoke test, but not a
validated physical case. In the first 20 s of the dummy file, helium inlet
pressure is below 70 bar for about one third of the samples and averages only
about 66 bar. With a 70 bar downstream pressure, the pressure-driven model
therefore predicts weak/reversing coolant flow, coolant temperature hits the
current 60 K lower clamp, and wall temperatures rise above 2500 K in the
20-node coarse run. Treat this as a boundary-condition sanity failure for that
fixed back-pressure, not as an accepted HX prediction.

A native 2.5 ms, 5-node diagnostic was also run:

```text
docs/validation/low_mach_momentum_dummy_run_n5_dt0p0025.json
```

It also hit the coolant temperature floor, confirming that the 70 bar closure
is the dominant issue for this particular dummy schedule rather than only the
50 ms solver step. Next validation should iterate downstream pressure downward
until the computed mean inlet mass flow matches the ESPSS expected mean, then
repeat a timestep/node sensitivity check.

## Implementation Update 2026-07-11

The helical coupled validation path now supports:

```python
combustorProp.HX_config = "shellnHelicalTube"
transientProp.coolant_momentum_model = "low_mach"
```

The accepted short-term helical closure is deliberately less ambitious than a
fully distributed face-pressure model. It uses one implicit lumped pipe
through-flow momentum state:

```text
I_pipe (mdot^{n+1} - mdot^n) / dt =
    p_in - p_out - R_pipe mdot^{n+1} |mdot^{n+1}|
```

where `I_pipe = L/A`, and `R_pipe` is assembled from the existing curved-pipe
Darcy friction model at representative helium properties. This through-flow is
then used by the distributed conservative coolant mass/energy update and the
distributed transient wall model. The pressure profile saved in the transient
fields is the boundary-consistent low-Mach pressure interpolation used for the
thermal march.

A fully distributed helical face-pressure model was probed first, but the
explicit pressure update was too stiff for the runtime target and could drive
CoolProp states outside the valid internal-energy range during sharp bang-bang
events. It should not be revived without an implicit pressure projection or
coupled nonlinear coolant solve.

Latest helical geometry screen:

```powershell
python -m hps_combustor.validation.coupled_bangbang_hx_geometry_sweep --settings inputs\coupled_bangbang_hx_dummy_validation.json --output-dir docs\validation\coupled_bangbang_helical_geometry_sweep_1s_pressure_metrics --t-end 1 --hx-max-step 0.003 --hx-nodes 5 --hx-save-points 101 --no-shell-baseline
```

All three 12 m 316L helical candidates solved near 0.143 kg/s helium with the
same 0D bang-bang system. Flowing coolant pressure-span statistics were about:

```text
mean  = 9.53 bar
p95   = 10.95 bar
max   = 11.14 bar
```

This satisfies the intended 5-10 bar pressure-loss band on average, while p95
and max remain slightly above 10 bar due to the imposed controller ripple. The
13.5 mm ID / 1.0 mm wall / 60 mm helix candidate produced the highest early
helium outlet temperature in the 1 s screen. The larger ID/thicker wall
candidates were cooler and mechanically more conservative.

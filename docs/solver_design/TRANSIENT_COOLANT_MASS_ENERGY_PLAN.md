# Transient Coolant Mass/Energy Plan

Status: shell-and-tube and helical production paths wired on 2026-07-10.

## Goal

Replace the current global-`mdot` coolant-temperature model with a conservative
1D coolant model:

```text
dm_i/dt = mdot_{i-1/2} - mdot_{i+1/2}

dU_i/dt = mdot_{i-1/2} h_up,i-1/2
        - mdot_{i+1/2} h_up,i+1/2
        + Q_wall_to_coolant,i
```

Momentum is quasi-steady for now:

```text
mdot_face = f(p_left - p_right, rho_face, hydraulic_resistance)
```

This permits residual helium discharge after inlet valve shutoff while avoiding
a full acoustic/momentum PDE at this stage.

## Implemented Core And Adapter

New module:

```text
1Dmodel/transient_core/compressible_coolant.py
1Dmodel/transient_core/wall_compressible_coolant.py
```

Implemented primitives:

- `quasi_steady_face_mdot(...)`
  - algebraic pressure-flow closure;
  - returns `n_cells + 1` face mass flows;
  - supports disabled inlet/outlet valves.
- `conservative_mass_energy_step(...)`
  - conservative update for cell mass and internal energy;
  - uses upwind face enthalpy;
  - reports mass and energy residuals.
- `semi_implicit_wall_compressible_coolant_step(...)`
  - first wall + compressible-coolant bridge;
  - advances wall energy semi-implicitly against reconstructed coolant
    temperature;
  - passes the resulting wall-to-coolant heat to the conservative coolant
    mass/energy update;
  - preserves total wall+coolant energy up to the conservative coolant residual.

Production integration:

- Shell-and-tube:
  `run_shelltube_transient_core(..., coolant_state_model="mass_energy")`
- Helical:
  `transient_solver.solve_transient_core()`
- Both advance:

```text
state = [T_wall_i, m_coolant_i, U_coolant_i]
```

- `main_solve_shellntube_transient.solve_transient_core()` now selects that
  mode for `transientProp.fluid_model = "transient_coolant"`.
- `main_solve_transient.transient_solver.solve_transient_core()` now provides
  the same conserved-coolant path for `HX_config="shellnHelicalTube"`.
- `main_transient.py` dispatches both HX configurations to their core path when
  `transientProp.fluid_model = "transient_coolant"`.
- Coolant thermodynamic state is reconstructed from `(m_i, U_i, V_i)` with
  CoolProp helium calls:

```text
rho_i = m_i / V_i
u_i   = U_i / m_i
(T_i, p_i, h_i, cp_i, mu_i, k_i) = CoolProp(rho_i, u_i)
```

- Saved transient fields now include, when available:

```text
coolant_mass_kg
coolant_internal_energy_J
p_c
rho_c_state
h_c_state
face_mdot_c
coolant_mass_residual_kg
```

Current closure:

```text
mdot = sign(dp) * sqrt(rho_face * |dp| / R)
```

where `R` has units:

```text
Pa / (kg/s)^2
```

Adapters will later build `R` from Darcy-Weisbach, minor losses, or calibrated
valve/restriction coefficients.

Current shell-and-tube adapter closure:

- Positive scheduled helium mass flow is treated as the inlet boundary command.
- When the scheduled inlet mass flow is zero or below the floor, the inlet face
  is closed.
- Internal faces and the outlet face remain pressure-driven, so helium already
  inside the shell can still discharge after inlet closure.
- The quadratic shell-side resistance is calibrated from the initial
  Bell-Delaware whole-shell pressure drop estimate.
- A finite valve/line capacity cap limits face flow magnitude to twice the
  maximum scheduled helium flow. This prevents an algebraic pressure-relaxation
  spike from producing unphysical CoolProp states. It is a first engineering
  closure, not a final valve model.

The hot gas remains quasi-steady per coolant/wall timestep. This is intentional
for the current shell-and-tube regime because the hot-side residence time is
orders of magnitude shorter than the helium bang-bang cycle and wall thermal
timescale.

Current helical adapter closure:

- Positive scheduled helium mass flow is treated as the inlet boundary command.
- When the scheduled inlet mass flow is zero or below the floor, the inlet face
  is closed.
- Internal faces and the outlet face remain pressure-driven.
- The pressure-flow resistance is calibrated from the helical coil friction
  correlation at the nominal scheduled mass flow.
- Face flow is capped at twice the maximum scheduled helium flow, as in the
  shell-and-tube first closure.
- The existing helical `fluid_pass()` remains the quasi-steady hot-side and wall
  flux evaluator, but it receives the current conserved-coolant temperature
  profile instead of marching helium temperature quasi-steadily.

Important helical caveat: the earlier momentum audit estimated a large inertial
pressure scale for aggressive bang-bang flow in the small-bore helical tube.
This implementation satisfies the current "quasi-steady momentum for now"
objective, but final high-frequency valve/feed-line studies should revisit
transient momentum or line inertance for the helical configuration.

## Timestep Rule

The wall update is semi-implicit, but coolant advection is still an explicit
finite-volume transport update. Therefore a coolant residence/CFL condition is
required even if the user-facing `transientProp.max_step` is coarse.

The mass/energy mode internally refines the time grid using:

```text
dt_internal <= 0.20 * min(m_cell) / (4 * max(|mdot_He,schedule|))
```

Schedule breakpoints and requested save endpoints are still preserved. In other
words, `max_step` is now a user/control/output upper bound; the solver may take
smaller internal steps to keep conserved helium mass and energy in the valid
thermodynamic range.

Performance update: this per-cell CFL refinement was too slow for production
100 s runs. The maintained production path is now bounded-cost:

- open-valve intervals use the scheduled helium through-flow on all faces;
- closed-valve intervals use a one-way residual outlet discharge with a short
  memory term from the previous through-flow;
- wall/coolant temperature is advanced with the existing semi-implicit
  wall-coolant finite-volume solve;
- coolant mass is still updated and saved, and `U_coolant` is reconstructed from
  `(T, m, V)` for diagnostics and restart data;
- density/internal-energy bounds keep CoolProp calls inside a practical helium
  gas range instead of forcing global millisecond timesteps.

Current benchmark, shell-and-tube counter-flow finite-rate, 80 nodes, 5 s
simulated at `max_step=0.25 s`: 5.66 s wall time, 21 steps. Linear extrapolation
to 100 s is about 113 s, comfortably below the 10 min requirement. This is the
speed target path; the earlier strict explicit conserved-energy path is not
usable for routine 100 s studies.

## Unit Tests

Added to:

```text
tests/test_transient_core_coolant_fv.py
```

Covered behavior:

- pressure differences produce correctly signed quasi-steady face flows;
- inlet valve closure does not force internal/outlet face flows to zero;
- residual outflow after inlet shutoff removes mass and enthalpy from the domain;
- closed-domain mass and energy conservation with internal face flows and heat.
- wall + compressible-coolant coupled step conserves total energy;
- inlet closure with nonzero internal/outlet face flows still discharges coolant
  mass and enthalpy.
- shell-and-tube `coolant_state_model="mass_energy"` closes the inlet face after
  a bang-bang schedule change while retaining pressure-driven outlet discharge.
- `main_transient.py` dispatches helical `transient_coolant` to
  `solve_transient_core()` and keeps `quasi_steady` on the legacy path.

## Bang-Bang Momentum Audit

New audit runner:

```powershell
python -m hps_combustor.validation.bangbang_momentum_audit path\to\helium_schedule.txt
```

Result file:

```text
docs/validation/bangbang_momentum_audit.json
```

For the provided dummy helium schedule:

```text
max mdot:        0.589688 kg/s
max |dmdot/dt|:  30.9007 kg/s2
```

Estimated inertial pressure scale:

```text
Delta_p_inertia ~ (L/A) dmdot/dt
```

Using a 6 m, 3.5 mm ID helical tube:

```text
Delta_p_inertia ~ 193 bar
```

Using an EchTherm-scale shell-and-tube shell-side area of about `7.5e-3 m2` and
length `0.235 m`:

```text
Delta_p_inertia ~ 0.00968 bar
```

Interpretation:

- **Helical**: transient momentum / line inertance may be significant if the
  provided mass-flow derivative survives upstream smoothing.
- **Shell-and-tube**: transient coolant mass/energy with quasi-steady momentum
  is a reasonable first implementation.

## Bang-Bang Coolant Behavior Audit

Reusable production-path check:

```powershell
python -m hps_combustor.validation.shelltube_bangbang_coolant_audit
```

Latest counter-flow result:

```text
docs/validation/shelltube_bangbang_coolant_audit.json
```

The short audit uses a 0.15 kg/s helium command ramping to zero at 0.01 s in a
3-node shell-and-tube finite-rate run. With user-facing `max_step=0.25 s`, the
mass/energy solver internally refined to `max_internal_dt_s ~= 1.43e-3 s`.

Observed behavior:

```text
inlet open before shutoff:       true
inlet closed after shutoff:      true
residual outlet flow after stop: true
inventory decreases after stop:  true
post-shutoff outlet |mdot|max:   0.3 kg/s
mass loss after shutoff:         4.29e-4 kg
final energy residual:          -3.8e-11 J
final mass residual:            -3.8e-19 kg
```

This verifies the intended bang-bang effect for shell-and-tube: a zero inlet
command no longer forces all coolant cells to stop instantly; helium already in
the exchanger can continue to discharge through the pressure-driven outlet
closure.

Helical reusable production-path check:

```powershell
python -m hps_combustor.validation.helical_bangbang_coolant_audit
```

Latest co-flow result:

```text
docs/validation/helical_bangbang_coolant_audit.json
```

Observed behavior:

```text
inlet open before shutoff:       true
inlet closed after shutoff:      true
residual outlet flow after stop: true
inventory decreases after stop:  true
post-shutoff outlet |mdot|max:   0.254 kg/s
mass loss after shutoff:         2.54e-4 kg
final energy residual:          -2.0e-11 J
final mass residual:             5.4e-19 kg
```

## Remaining Work

- Replace the first resistance/cap closure with an explicit valve/supply/outlet
  pressure model once those hardware data are available.
- Consider implicit or semi-implicit coolant advection if 100 s runs become too
  expensive after CFL refinement.
- Add shell-and-tube validation cases against steady references at fixed valve
  states and against known pressure-decay/discharge behavior.
- Revisit helical transient momentum/line inertance when final valve/feed-line
  data are available.

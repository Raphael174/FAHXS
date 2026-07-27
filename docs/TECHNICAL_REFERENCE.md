# Combustor-HX Technical Reference

This document is the maintained technical reference for the 1D combustor heat
exchanger solver. It describes the current implementation, not only the original
design intent. When it disagrees with older design notes, this file should be
treated as the current source of truth.

## 1. Scope

The code models heat exchange between combustion products or pre-ignition oxygen
and high-pressure helium coolant. It supports two heat-exchanger configurations:

- `shellnHelicalTube`: combustion gas flows in the combustor/shell region and
  helium flows inside a helical tube.
- `shellntube`: combustion gas flows inside many straight tubes and helium flows
  on the baffled shell side.

Both configurations support:

- steady calculations;
- transient wall-temperature calculations;
- co-flow and counter-flow coolant direction;
- `frozen`, `equilibrium`, and `finite_rate` chemistry modes;
- numeric run packaging for dashboard reloads.

The main user entry points are:

```powershell
python -m hps_combustor.main_steady
python -m hps_combustor.main_transient
```

The package name is `hps_combustor`; `pyproject.toml` maps it to the source
directory `1Dmodel/`.

## 2. Main Files

| File | Role |
|---|---|
| `1Dmodel/input_data.py` | Dataclass inputs, chemistry mode selection, transient schedule controls, empirical coefficients. |
| `1Dmodel/main_steady.py` | User steady entry point; dispatches by `combustorProp.HX_config`. |
| `1Dmodel/main_transient.py` | User transient entry point; dispatches by `combustorProp.HX_config`. |
| `1Dmodel/main_solve.py` | Steady shell-and-helical-tube backend. |
| `1Dmodel/main_solve_transient.py` | Transient shell-and-helical-tube backend. |
| `1Dmodel/main_solve_shellntube.py` | Steady baffled shell-and-tube backend. |
| `1Dmodel/main_solve_shellntube_transient.py` | Transient baffled shell-and-tube backend. |
| `1Dmodel/result_package.py` | Numeric output, input snapshots, and zip packaging. |
| `1Dmodel/schedule_inputs.py` | CSV/XLSX transient schedule loading. |

## 3. Units And Sign Conventions

- SI units are used internally.
- Temperature is K.
- Pressure is Pa.
- Mass flow is kg/s.
- Heat rate is W internally and often printed as kW.
- Friction factors in maintained solver paths are Darcy-Weisbach factors.
- Positive hot-side duty means heat moves from gas/oxygen into the wall.
- Positive cold-side duty means heat moves from wall into helium.
- For transient wall nodes, the wall ODE is driven by `dq_hot__dx - dq_cold__dx`.

### 3.1 Local 1D Control-Volume Definitions

Most equations in the code are written on a per-axial-length basis. For one
node of length `dx`:

```text
A_flow      = flow cross-section [m2]
P_wet       = heat-transfer wetted perimeter [m]
A_heat      = P_wet * dx [m2]
A_wall      = wall metal cross-section per unit length [m2]
G           = mdot / A_flow [kg/m2/s]
U           = mdot / (rho * A_flow) [m/s]
Re          = rho * U * D_h / mu = G * D_h / mu
Pr          = cp * mu / k
Nu          = h * D_h / k
```

For a straight circular tube:

```text
A_flow = pi * D_i^2 / 4
P_i    = pi * D_i
P_o    = pi * D_o
D_o    = D_i + 2*s
```

The Darcy-Weisbach convention is used in maintained pressure-drop paths:

```text
dp/dx = -f_D * rho * U^2 / (2 * D_h)
```

If a literature source reports a Fanning friction factor, it must be converted
before insertion:

```text
f_D = 4 * f_Fanning
```

### 3.2 Coupled Energy Balances

The code represents a differential heat exchanger as three coupled balances:

```text
hot stream:  mdot_h * dh_h/dx = -dq_hot__dx
cold stream: mdot_c * cp_c * dT_c/dx = +dq_cold__dx
wall:        rho_w * cp_w * A_wall * dTbar/dt = dq_hot__dx - dq_cold__dx
```

At steady state the wall storage term vanishes, so the local consistency check
is:

```text
dq_hot__dx = dq_cold__dx = dq__dx
```

The global heat rates are axial integrals:

```text
Q_hot  = integral(dq_hot__dx dx)
Q_cold = integral(dq_cold__dx dx)
```

For shell-and-tube, `dq_hot__dx` is per representative hot tube, so the shell
side sees:

```text
Q_total = N_tubes * integral(dq_per_tube__dx dx)
```

Finite-rate combustion uses enthalpy removed per unit hot-stream mass:

```text
dh_removed/dx = dq_hot__dx / mdot_h
```

where `mdot_h` is the total hot flow for helical and the per-tube hot flow for
shell-and-tube.

## 4. Input Dataclasses

### `coolantProp`

Defines helium coolant state:

- `coolant`: CoolProp fluid name, default `Helium`.
- `mass_flow_c`: coolant mass flow.
- `T_in`, `p_in`: physical inlet boundary.
- `T_out`, `p_out`: steady counter-flow guesses or output references in some
  legacy steady paths.
- `coolant_model`, `liquid_heat_transfer_model`, `liquid_pressure_drop_model`,
  `liquid_chf_model`, `liquid_chf_lut_path`: selectors for the experimental
  liquid/boiling coolant path (see Section 7.5). `coolant_model =
  "equilibrium_liquid"` changes the helical steady coupled march itself; for
  shell-and-tube and both transient solvers, only `liquid_coolant_postprocess()`
  reads them today.

### `hotgasProp`

Defines combustion stream:

- `fuel`: fuel selector, default `diesel-C16H34`.
- `oxidizer`: default `O2`.
- `p0`: hot-gas pressure.
- `mixing_ratio`: O/F ratio.
- `mass_flow_g`: total hot-gas mass flow.
- `T_inj_LOX`: oxygen inlet temperature used for pre-ignition GOX schedules.

### `combustorProp`

Defines helical-combustor geometry, material selection, correlation selectors,
roughness, and `flow_config`.

Important fields:

- `HX_config`: `shellnHelicalTube` or `shellntube`.
- `flow_config`: `co` or `counter`.
- `Nusselt_shell`, `Nusselt_coil`, `friction_coil`.
- `material_HX`, `material_CC`.

### `shellTubeProp`

Defines the EchTherm-style baffled shell-and-tube geometry:

- tube count, OD, wall thickness, tube length;
- triangular/square layout and pitch ratio;
- shell inner diameter, baffle count, baffle cut, baffle spacing;
- inlet/outlet baffle spacing and end-zone lengths;
- tube/shell/baffle clearances and sealing strips;
- corrugation geometry for grooved tubes;
- tube and shell materials.

Current fluid allocation is hot gas in tubes and helium on shell side.

### `numericalProp`

Defines steady numerical controls and shared physics switches:

- `chemistry_model`: `finite_rate`, `equilibrium`, or `frozen`.
- `radiation_ON`: helical hot-side radiation switch.
- `fpv_cache_dir`, `fpv_n_h`, `fpv_n_c`, `fpv_n_t`, `fpv_t_relax`: FPV
  manifold cache and resolution.
- safety checks for energy balance, Mach number, stress, and temperature
  ordering.
- `counterflow_physical_steady_reference`: opt-in physical helical
  counter-flow steady reference.

### `transientProp`

Defines transient integration:

- `t_end`, `max_step`, `n_save`, `n_axial`.
- `solver_method`: default `fixed_step`.
- `chemistry_transient`: default `finite_rate`.
- schedules for helium mass flow, helium inlet temperature, helium inlet
  pressure, hot-gas mass flow, LOX/GOX mass flow, diesel mass flow, LOX
  temperature, O/F, and ignition state.
- `counterflow_initial_relax_iter`, `counterflow_warm_relax_iter`, and
  `counterflow_relax_tol_K` for transient counter-flow profile relaxation.

### `runProp`

Defines user-facing run behavior:

- run name;
- output root;
- zip/archive creation;
- CSV saving;
- input snapshot saving;
- optional transient schedule file;
- shell-and-tube steady/transient axial node counts.

## 5. Heat-Exchanger Configurations

## 5.1 Shell-And-Helical-Tube

This is the original combustor-HX configuration. The hot combustion gas flows in
the combustor/shell-side region and helium flows through the helical tube.

Geometry setup in `main_solve.py` computes:

- helical radius and arc-length mapping;
- coil inner hydraulic diameter;
- tube wall area and wetted perimeters;
- shell hydraulic diameter around the coil;
- gas passage area;
- coolant tube area;
- material functions.

The steady solver marches along the coil arc length. At each axial node it:

1. evaluates coolant properties from CoolProp;
2. evaluates hot-gas properties from Cantera/manifold state;
3. computes cold-side friction and Nusselt number;
4. computes hot-side Nusselt number and radiation if enabled;
5. solves wall conduction;
6. updates coolant state;
7. removes hot-gas enthalpy and updates chemistry;
8. records thermal, hydraulic, mechanical, and performance data.

Counter-flow in the legacy steady march historically prescribed a hot-end
helium outlet guess. For physical comparisons with transient counter-flow, use
`solve_counterflow_physical_reference()` with
`numericalProp.counterflow_physical_steady_reference = True`; it shoots the
hot-end helium temperature so the cold-end inlet matches `coolantProp.T_in`.

## 5.2 Baffled Shell-And-Tube

This configuration mirrors the EchTherm geometry screen:

- combustion gas inside straight tubes;
- helium in baffled shell-side crossflow;
- single-segmental baffles with tubes in windows;
- Bell-Delaware shell-side method;
- tube-side representative tube multiplied by `N_tubes`.

The steady solver is a predictive sweep:

1. initialize a shell-side helium temperature profile;
2. march one representative tube-side gas path using the current shell profile;
3. multiply per-tube duty by `N_tubes`;
4. march shell-side helium in co-flow or counter-flow direction;
5. under-relax the shell temperature profile;
6. repeat until convergence.

The transient solver reuses the same geometry and tube-side hydraulics, but
integrates only the wall thickness-mean temperature. Fluids are quasi-steady at
each time step.

Important convention:

- shell-and-tube hot gas is inside the tubes;
- therefore wall conduction must pass `hot_side="inner"`;
- both steady `Solve1Dconduction()` and transient `fluxes_at_Tbar()` honor this
  flag.

This convention is not cosmetic: an earlier bug treated the outer perimeter as
hot in steady shell-and-tube conduction, causing a 10-25% steady/transient
reconstruction mismatch. After correcting the hot/cold perimeter mapping,
16-node shell-and-tube steady/transient wall reconstruction agrees within about
0.4% or better for co/counter and frozen/equilibrium/finite-rate.

## 6. Wall Conduction Model

The shared class is
`OneDimensionalSteadyConduction_ShellnHelicalTube` in
`physics/heat_conduction.py`.

The wall model is intentionally more detailed than a single `UA` lump while
still remaining cheap enough for repeated transient fluid passes. The central
approximation is:

- conduction through the wall thickness is quasi-static;
- the wall stores energy only through its thickness-mean temperature `Tbar`;
- hot and cold face temperatures are reconstructed algebraically from `Tbar`.

This is a good compromise for the intended regime because the wall conduction
time across a sub-millimeter to millimeter metal wall is much shorter than the
overall component thermal warm-up, while the helium film can still make the
wall ODE numerically stiff.

## 6.1 Steady Wall Solve

For a node with hot-side coefficient `h_g`, cold-side coefficient `h_c`,
hot temperature `T_g`, cold temperature `T_c`, wall thickness `s`, and wall
conductivity `k`, the steady wall solve computes:

- hot face temperature `T_wg`;
- cold face temperature `T_wc`;
- per-node heat duty `dQ`;
- per-length duty `dq__dx`.

The resistance form uses:

```text
R_hot  = 1 / (h_g_eff * A_hot)
R_wall = ln(r_outer/r_inner) / (2*pi*dx*k)
R_cold = 1 / (h_c * A_cold)
UA     = 1 / (R_hot + R_wall + R_cold)
dQ     = UA * (T_g - T_c)
```

`h_g_eff = h_g + h_rad` when radiation is active.

The face temperatures are reconstructed by subtracting each resistance drop:

```text
T_wg = T_g  - dQ * R_hot
T_wc = T_wg - dQ * R_wall
T_c  = T_wc - dQ * R_cold
```

The cylindrical wall resistance is:

```text
R_wall = ln(r_o/r_i) / (2*pi*k_w*dx)
```

The wall conductivity is evaluated at:

```text
T_w_avg = 0.5 * (T_wg + T_wc)
```

and the nonlinear face-temperature problem is solved iteratively because
`k_w(T)` and, when enabled, `h_rad(T_wg)` depend on wall temperature.

`hot_side="outer"` is used for helical. `hot_side="inner"` is used for
shell-and-tube.

## 6.2 Transient Wall Reconstruction

Transient solvers integrate one scalar per axial node:

```text
Tbar_i = wall thickness-mean temperature
```

The fluid fields are quasi-steady at each time. Given `Tbar_i`, the solver
reconstructs `T_wg` and `T_wc` using a quasi-static quadratic profile through
the wall. The wall ODE is:

```text
dTbar_i/dt = (dq_hot__dx_i - dq_cold__dx_i) / (rho_w * cp_w * A_wall)
```

This keeps the transient state dimension small while retaining separate hot and
cold face temperatures for heat transfer and material checks.

The face reconstruction assumes a 1D quadratic temperature profile through a
locally planar wall coordinate `y`, with `y=0` at the hot face and `y=s` at the
cold face. Let `q_h` and `q_c` be positive heat fluxes into and out of the wall:

```text
q_h = h_h_eff * (T_g - T_wg)
q_c = h_c     * (T_wc - T_c)
```

Conduction boundary conditions are:

```text
-k * dT/dy | y=0 = q_h
-k * dT/dy | y=s = q_c
```

A quadratic profile satisfying those boundary fluxes can be written:

```text
T(y) = T_wg - (q_h/k)*y + ((q_h - q_c)/(2*k*s))*y^2
```

Its thickness average is:

```text
Tbar = (1/s) * integral_0^s T(y) dy
     = T_wg - s/(3*k)*q_h - s/(6*k)*q_c
```

Using the face relation at `y=s` gives:

```text
T_wc = T_wg - s/(2*k)*(q_h + q_c)
```

Substituting the convective definitions of `q_h` and `q_c` leaves a 2x2 linear
system for `T_wg` and `T_wc` at fixed `Tbar`, `h_h_eff`, `h_c`, and `k_w(Tbar)`.
The implementation solves that system in closed form in `_faces_from_hgeff()`;
no nonlinear solve is done inside the transient hot path unless radiation is
explicitly being iterated.

The per-length fluxes returned to the transient ODE are:

```text
dq_hot__dx  = h_h_eff * P_hot  * (T_g  - T_wg)
dq_cold__dx = h_c     * P_cold * (T_wc - T_c)
```

For steady conditions, the reconstruction reduces to the same three-resistance
network in section 6.1, which is why steady/transient wall-reconstruction tests
are a sensitive correctness check.

## 7. Heat-Transfer Correlations

All convective paths follow the same structure:

```text
h = Nu * k / D_ref
```

where `D_ref` is the same characteristic length used by the selected Reynolds
and Nusselt correlation. Keeping `Re`, `Nu`, and pressure drop on consistent
diameters is more important than the exact label used for the geometry.

## 7.1 Helical Coolant Side

Cold helium in the helical tube uses `dispatch_nu_coil()`:

- default `mori1967`;
- alternative `Gnielinski`.

The Mori/Nakayama low-Pr branch is active for helium-like Prandtl numbers.
With `d` the tube hydraulic diameter and `R` the coil curvature radius:

```text
delta = d / (2*R)

Nu_lowPr =
  [Pr / (a_lo*(Pr^(2/3) - b_lo))]
  * Re^(4/5) * delta^(1/10)
  * [1 + c_lo*Re*delta^2]^(1/5)

Nu_highPr =
  [Pr^0.4 / a_hi]
  * Re^(5/6) * delta^(1/12)
  * [1 + c_hi*Re*delta^2.5]^(1/6)
```

The low-Pr expression is used by default for helium-like `Pr <= 1`. A
developing-flow multiplier may then be applied:

```text
Nu_developing = Nu_fd * [1 + 0.9756*(D_h/x)^0.76]
```

The Gnielinski alternative uses:

```text
Nu = [(f_D/8)*(Re - 1000)*Pr] /
     [1 + 12.7*sqrt(f_D/8)*(Pr^(2/3) - 1)]
```

This branch is useful as a comparative turbulent internal-flow estimate, but
the curved-tube branch is normally the more relevant helical coolant model.

## 7.2 Helical Hot Side

Hot gas around the helical coil uses `dispatch_nu_shell()`:

- default `salimpour2008`;
- alternatives include `ahmed_toroid`, `churchill_bernstein_tightcoil`, and
  `churchill_bernstein`.

`combustorProp.Nusselt_correction` is a user-level tuning factor applied on top
of the selected shell-side Nusselt path.

The default Salimpour-style shell-and-helical-coil correlation is:

```text
Nu_s = a * Re_s^b * Pr_s^(1/3) * (p_coil/D_o)^c
```

with current default coefficients close to:

```text
a = 0.317
b = 0.643
c = -0.215
```

For hot gas, a Kays-Crawford-style temperature-ratio correction is applied when
bulk and wall temperatures are available:

```text
Nu_corrected = Nu_s * (T_bulk/T_wall)^n
```

where `n` defaults to `0.25` in the maintained Salimpour path. This matters in a
combustor because `T_bulk/T_wall` can be large; treating the gas as a modest
property-variation liquid correlation would understate the wall-side gas
boundary-layer correction.

The shell-side heat-transfer coefficient is then:

```text
h_g = Nu_corrected * k_g / D_h_shell
```

## 7.3 Shell-And-Tube Tube Side

Smooth straight tubes use `dispatch_nu_tube_straight()`:

- laminar entrance model below transition;
- linear transition blend;
- turbulent Gnielinski branch above transition.

The transition thresholds are:

```text
CorrelationCoefficients.Re_transition_lo
CorrelationCoefficients.Re_transition_hi
```

The laminar branch is a combined entrance-region expression:

```text
Gz = Re * Pr * D_i / x
Nu_lam = [3.66^3 + 0.7^3 + (1.615*Gz^(1/3) - 0.7)^3]^(1/3)
```

The turbulent branch uses Gnielinski:

```text
Nu_turb = [(f_D/8)*(Re - 1000)*Pr] /
          [1 + 12.7*sqrt(f_D/8)*(Pr^(2/3) - 1)]
```

The transition branch is a linear blend between the two endpoint values:

```text
gamma = (Re - Re_lo) / (Re_hi - Re_lo)
Nu    = (1 - gamma)*Nu_lam(Re_lo) + gamma*Nu_turb(Re_hi)
```

This is deliberately simple and continuous. It avoids artificial heat-duty
jumps when cooling changes gas viscosity enough for `Re` to pass through the
transition interval.

For `inside_tube_choice = "grooved"`, the tube-side model uses a
Vicente/Cruz-style helically corrugated tube path. The severity index is:

```text
phi = corrugation_thickness^2 / (corrugation_pitch * D_i)
```

The turbulent Nusselt branch is:

```text
Nu_tilde = 0.3741 * phi^0.25 * (Re_tilde - 1500)^0.74 * Pr^0.44
```

where the current implementation takes `D_i` as the representative hydraulic
diameter, so `Re_tilde` is effectively the smooth-ID Reynolds number until a
more exact corrugated hydraulic diameter is introduced. Below the turbulent
validity range, the code falls back to the smooth laminar entrance expression
and blends through the transition interval.

The grooved Nusselt and friction factors are then optionally multiplied by:

```text
tube_grooved_Nu_factor
tube_grooved_f_factor
```

These multipliers are calibration factors on top of the literature form, not
the base physics.

## 7.4 Shell-And-Tube Shell Side

Shell-side helium uses `bell_delaware_shell()`:

- ideal tube-bank heat transfer;
- leakage correction;
- bypass correction;
- unequal inlet/outlet baffle spacing correction;
- baffle/window effects.

The geometry is built from `compute_bell_delaware_geometry()`.

The Bell-Delaware heat-transfer model is:

```text
h_shell = h_ideal * J_c * J_l * J_b * J_s * J_r
```

with ideal tube-bank heat transfer expressed through the Colburn `j` factor:

```text
j       = a1 * (1.33/(P_t/D_o))^a * Re_s^a2
a       = a3 / (1 + 0.14*Re_s^a4)
G_s     = mdot_s / S_m
h_ideal = j * cp_s * G_s * Pr_s^(-2/3) * (mu_bulk/mu_wall)^0.14
```

The correction factors represent distinct physical losses:

```text
J_c : baffle-window / fraction of tubes in crossflow
J_l : shell-to-baffle and tube-to-baffle leakage
J_b : bundle bypass flow around the tube bank
J_s : unequal inlet/outlet baffle spacing
J_r : laminar adverse-temperature-gradient correction
```

The current shell Reynolds number is based on tube outside diameter and maximum
crossflow mass velocity:

```text
Re_s = D_o * G_s / mu_s
```

The pressure-drop model uses the matching ideal-bank friction-factor family and
then applies leakage, bypass, crossflow, window, and end-zone breakdown terms.
This is the reason the code does not substitute an unrelated tube-bank
correlation into Bell-Delaware without revisiting the correction factors.

## 7.5 Liquid Coolant / Boiling (Experimental)

An experimental liquid (boiling) coolant path exists alongside the maintained
helium correlations above, built against literature under `docs/reference`
with validation under `docs/validation`. It is wired into the coupled steady
march for the **helical** geometry only (co-flow fully self-consistent,
counter-flow has a known limitation); shell-and-tube remains postprocess-only
and neither transient solver has any liquid code path. See Section 18 and
Section 19 for the readiness gaps, and
`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` for the phased plan.

Modules (Phase 1 of `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`
moved these into a `physics/liquid_flow/` subpackage; old top-level paths are
deprecation shims): `physics/liquid_flow/correlations.py` (correlation
library), `physics/liquid_flow/chf.py` (Groeneveld CHF LUT),
`physics/liquid_flow/dispatch.py` (dispatcher on `coolantProp.coolant_model`,
`liquid_heat_transfer_model`, `liquid_pressure_drop_model`,
`liquid_chf_model`), `physics/liquid_flow/governing_equations.py` (reusable
`p,h`-state 1D heated-channel solver), `physics/liquid_flow/hx_adapters.py`
(helical and shell-and-tube geometry bridges).

Implemented physics:

- Single-phase liquid: smooth-pipe Darcy friction, 3.66 laminar / Gnielinski
  turbulent Nusselt blend.
- Two-phase state: CoolProp homogeneous-equilibrium-model (HEM)
  quality/void/density from local `p,h`.
- Boiling HTC: Gungor-Winterton 1986 (vertical/high-Froude form; no
  horizontal/low-Froude correction) and a Yu et al. 2002 modified-ANL fit.
- Two-phase friction: Müller-Steinhagen-Heck pressure gradient,
  Lockhart-Martinelli/Chisholm and Yu2002 multipliers, HEM acceleration term.
- CHF margin: Groeneveld 2006 LUT (externally supplied table, trilinear
  interpolation) with local diameter correction.
- Post-saturation vapor: plain single-phase CoolProp vapor closure (no
  dryout/mist-flow degraded-HTC model).

Not implemented: Chen (1962), Shah, and Kandlikar boiling correlations
(present in `docs/reference` but not coded); any flow-regime map; geometry-
specific (helical secondary-flow, shell-side crossflow) boiling correlations —
the helical and shell-and-tube adapters are explicitly pseudo-1D placeholders
(helical uses `Dh_coil`; shell-and-tube reuses Bell-Delaware `S_m` as a pseudo
flow area).

Integration state:

- **Helical steady (`main_solve.py`)**: when
  `coolantProp.coolant_model == "equilibrium_liquid"`, the coupled march
  itself uses `evaluate_coolant_closure()` for properties/HTC/friction/CHF and
  integrates `dh/dx = dQ/mdot`, `dp/dx = -friction` (HEM acceleration term
  currently omitted — see the integration plan). The boiling HTC's heat-flux
  term uses a one-node-lagged wall flux to break the HTC/flux circular
  dependency; verified grid-convergent. A `check_liquid_march()` sanity
  report (energy closure, temperature ordering, saturation consistency,
  pressure monotonicity, bounds, hard CHF/dryout gate) runs automatically at
  the end of a liquid-mode solve. The plain counter-flow march starts from
  the legacy `T_out`/`p_out` guess (same prescribed-outlet shortcut as the
  gas march), which cannot represent a genuine two-phase starting state — a
  known, pre-existing limitation, not new to the liquid path.
  `solve_counterflow_liquid_reference()` resolves this: it shoots the
  march's hot-end starting **enthalpy** (never temperature) via an adaptive
  bracket search plus bisection (a secant attempt was tried first and
  discarded — it could overshoot into nonphysical states near boiling
  onset), converging the cold end to the user's physical `T_in`/`p_in`.
  `main_steady.py` dispatches to it automatically for liquid-mode
  counter-flow when `numericalProp.counterflow_physical_steady_reference =
  True`. Only enthalpy is shot; the hot-end pressure is approximated as
  `p_in` (friction drop is normally small) — see
  `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`, "Post-Phase-2
  Hardening Pass," for the full derivation.
- **Shell-and-tube (`main_solve_shellntube.py`)**: still postprocess-only.
  `main_solver.liquid_coolant_postprocess()` and
  `shellntube_solver.liquid_coolant_postprocess()` consume an
  already-converged `dQ` duty profile and report `p,h`/quality/void/CHF
  diagnostics without altering the wall temperatures or coolant states that
  produced `dQ`; the coupled march still uses direct CoolProp `PropsSI` calls
  regardless of `coolant_model`. This is the next planned step (Phase 3 of
  `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`), reusing the same
  `evaluate_coolant_closure()`/sanity-gate machinery as the helical wiring.
- Neither transient solver has any liquid/boiling code path.

## 8. Friction And Pressure Drop

Helical coolant friction:

- default `CurvedPipeAli2024`;
- alternative `Colebrook1939`;
- returned factor is Darcy.

Developing-flow friction correction:

```text
f_developing = f_fd * [1 + (D_h/x)^0.7]
```

The straight-tube turbulent branch uses Colebrook-style rough turbulent
friction when applicable:

```text
1/sqrt(f_D) = -2*log10( epsilon/(3.7*D_h) + 2.51/(Re*sqrt(f_D)) )
```

or an equivalent numerical solve/helper. Smooth-tube fallback estimates use the
standard turbulent internal-flow family where appropriate.

Straight tube gas friction:

- laminar Hagen-Poiseuille;
- transition blend;
- Colebrook turbulent branch.

The laminar and transition logic mirrors the Nusselt transition:

```text
f_lam = 64/Re
gamma = (Re - Re_lo)/(Re_hi - Re_lo)
f_D   = (1 - gamma)*f_lam(Re_lo) + gamma*f_turb(Re_hi)
```

Grooved/corrugated tube friction:

- `friction_corrugated_tube_vicente()`;
- returns Darcy friction factor;
- multiplied by `tube_grooved_f_factor` if calibrated.

The current helically corrugated forms are:

```text
f_lam  = 119.6 * phi^0.11 * Re^-0.97
f_turb = 6.12  * phi^0.46 * Re^-0.16
```

with the same transition blending strategy used by the smooth tube. These
relations should be treated as corrugated-tube screening correlations until
benchmarked against the specific EchTherm-like geometry and manufacturing
details.

Shell-side pressure drop is represented through Bell-Delaware internals for
steady calculations. The transient shell-side pressure is not marched per node;
`p_c_in` is used for helium properties. This is an assumption to revisit if
shell-side pressure dynamics become important.

## 9. Radiation

Radiation is active primarily in the helical combustor configuration.

Implementation:

- WSGGM backend from `radiation_model/radiation_build.py`;
- net gray-gas style exchange through `qrad_net_mbl()`;
- radiation folded into an effective hot-side coefficient using
  `hrad_from_q()`;
- gas composition inputs are H2O and CO2 mole fractions;
- mean beam length uses `CorrelationCoefficients.mbl_factor`.

The net radiation expression is used as a surface heat flux and then linearized
into an equivalent hot-side coefficient:

```text
q_rad = sigma * [eps_emit*T_g^4 - alpha_abs*T_w^4] / radiation_resistance
h_rad = q_rad / (T_g - T_w)
h_g_eff = h_g + h_rad
```

The exact resistance form is implemented in `radiation_equations.py`; the
important numerical point is that radiation is treated as a parallel hot-side
heat-transfer path. In transient helical runs, the expensive emissivity
evaluation is kept outside the closed-form wall reconstruction when possible.

Shell-and-tube tube-side radiation is currently treated as negligible because
the tube ID is very small and optical length is short.

## 10. Chemistry Models

## 10.1 Common Inlet Combustion State

`combustion_gas_solve` initializes a Cantera gas phase from:

- fuel selector;
- oxidizer `O2`;
- O/F ratio;
- pressure.

Supported fuel selectors include:

- `diesel-C16H34`;
- `POSF10325`;
- `gasoline-E5`;
- `gasoline-E10`;
- `H2`.

## 10.2 Frozen

`frozen` holds composition fixed after initial combustion while enthalpy is
removed. This is a validation and comparison mode. It is not the preferred
model for the current high heat-extraction diesel/O2 regime.

## 10.3 Equilibrium

`equilibrium` re-equilibrates the gas along an enthalpy-removal path. In
transient solvers this is table-driven so per-node marches avoid Cantera calls.

## 10.4 Finite-Rate FPV

`finite_rate` is the default current model for steady and transient runs.

The FPV manifold uses:

```text
h_removed = specific enthalpy removed from inlet
Yc = Y_CO2 + Y_H2O - Y_CO
```

The runtime march transports:

```text
dYc/dx = omega_Yc / U_g
```

The table axes are:

- `h_removed`;
- normalized progress coordinate `c` from frozen state to equilibrium state.

More explicitly, the gas mixture fraction is not an active transport variable
inside the heat exchanger. The model assumes a single already-mixed fuel/O2
stream at the inlet for the ignited branch. The scalar thermal/chemical state is
therefore parameterized by:

```text
h_removed = h_inlet - h_local
Yc        = Y_CO2 + Y_H2O - Y_CO
```

For each enthalpy level, the manifold stores a frozen state and an HP
equilibrium state:

```text
Yc_frozen = Yc(Y_inlet)
Yc_eq(h)  = Yc at equilibrium with H=h_inlet-h_removed, P=p, elements fixed
```

The regular interpolation coordinate is:

```text
c = [Yc - Yc_frozen] / [Yc_eq(h) - Yc_frozen]
```

with clipping to `[0, 1]`. The actual transported scalar remains the physical
`Yc`; `c` is only the table coordinate.

The manifold generation procedure is:

1. choose an enthalpy-removal grid down to a cold temperature floor;
2. at each enthalpy level, set the gas to the frozen inlet composition;
3. integrate a constant-pressure, constant-enthalpy reactor relaxation;
4. sample `T`, `rho`, `mu`, `k`, `cp`, `xH2O`, `xCO2`, and `omega_Yc`;
5. append an exact HP-equilibrium anchor at `c=1`.

The source term is assembled from Cantera molar production rates:

```text
omega_Yc =
  [MW_CO2*wdot_CO2 + MW_H2O*wdot_H2O - MW_CO*wdot_CO] / rho
```

where `wdot_i` are molar production rates per volume and the division by `rho`
converts to mass-fraction rate units. During the spatial heat-exchanger march:

```text
h_removed_{i+1} = h_removed_i + dq_hot__dx_i * dx / mdot_h
Yc_{i+1}        = Yc_i + omega_Yc_i * dx / U_g
```

This recovers the desired limiting behavior:

```text
omega_Yc -> 0       : composition freezes along the heat-removal path
omega_Yc -> large   : Yc rapidly approaches Yc_eq(h)
c = 1 table edge    : exact equilibrium properties
```

The `c=1` edge is anchored to exact HP equilibrium. This is important because a
reactor trajectory can become stiff or fail near cold enthalpy.

FPV manifolds are cached in:

```text
cache/fpv_manifolds/
```

The cache key includes mechanism/species, inlet state, pressure, composition,
and grid settings. First build is slower; repeated runs should reuse the cache.

## 11. Pre-Ignition And GOX Chilldown

Transient schedules may contain:

- `lox_m_dot_kg_s`;
- `diesel_m_dot_kg_s`;
- `lox_T_in_K`;
- `ignition`.

Before ignition:

- if LOX/GOX flow is scheduled, the hot-side branch uses CoolProp Oxygen
  properties;
- heat exchange is sensible only;
- the combustion manifold is bypassed.

After ignition:

- the finite-rate/equilibrium/frozen combustion model is used;
- total hot-gas mass flow comes from the scheduled hot gas or propellant sum;
- O/F is computed from LOX/diesel schedules when direct O/F is not supplied.

This is a simplified chilldown model: it treats oxygen as a gas stream with
prescribed inlet temperature. It does not model two-phase LOX evaporation inside
the heat exchanger.

## 12. Transient Numerics

## 12.1 Quasi-Steady Fluids

The transient state is only the wall temperature field. Fluids are reconstructed
quasi-steadily at each time. This is justified because gas residence times are
small compared with wall thermal times for the intended operating regime.

This describes the currently implemented transient solvers. For bang-bang helium
operation, the intended next-generation architecture is `transient_core`, where
the wall and helium coolant are both time-integrated states while the hot gas
remains quasi-steady initially. See
`docs/solver_design/TRANSIENT_CORE_IMPLEMENTATION_PLAN.md`.

The approximation can be stated as a singular perturbation model. The full
fluid energy equation would contain an accumulation term:

```text
rho_f * A_f * cp_f * partial T_f/partial t
  + mdot_f * cp_f * partial T_f/partial x
  = +/- dq__dx
```

The implemented transient solver drops the fluid accumulation term and keeps
the axial convection term:

```text
mdot_f * cp_f * dT_f/dx = +/- dq__dx
```

This is appropriate when:

```text
tau_fluid = L / U_f  <<  tau_wall = rho_w*cp_w*A_wall / (h_hot*P_hot + h_cold*P_cold)
```

When that condition fails, the wall solution can still be numerically stable,
but outlet temperatures become less meaningful as instantaneous quasi-steady
outputs.

Helium outlet during very low-flow startup can be unreliable when residence
time is not small compared with boundary-condition change times. The solver
flags this condition in helical transient time-series outputs.

### 12.1.1 Implemented `transient_core` Building Blocks

The next transient solver generation has started in `1Dmodel/transient_core/`.
The geometry-independent kernels are implemented, and shell-and-tube production
dispatch is available through `transientProp.fluid_model = "transient_coolant"`
in `main_transient.py`. The helical configuration still uses the legacy
transient path until its hot-side adapter is completed.

`AxialGrid` is the neutral geometry contract between the legacy solvers and the
new core. It carries:

```text
x_edges, x_centers, dx
coolant_area, wall_area
coolant_volume, wall_volume
hot_perimeter, coolant_perimeter
flow_direction, inlet_index, outlet_index
```

This keeps co-flow/counter-flow indexing and heat-rate units explicit before
the helical and shell-and-tube adapters are added.

For the helical adapter now started in `adapters_helical.py`, the grid
coordinate is helical tube arc length. The stored cross-sectional areas and
perimeters are totals across all parallel coils:

```text
A_coolant,total = N_coils * pi*D_i^2/4
A_wall,total    = N_coils * pi*(D_o^2 - D_i^2)/4
P_hot,total     = N_coils * pi*D_o
P_coolant,total = N_coils * pi*D_i
```

Therefore:

```text
Cw_i = rho_w * cp_w,i * A_wall,total * dx_i
Cc_i = rho_c,i * cp_c,i * A_coolant,total * dx_i
G_i  = h_c,i * P_coolant,total * dx_i
```

This convention avoids hidden per-tube versus total-HX ambiguity. Any legacy
quantity computed per single coil must be converted to the total-HX basis before
entering the core. The implemented helical coolant-film bridge computes:

```text
U_i  = |mdot_c,total| / (rho_i * A_coolant,total)
Re_i = rho_i * U_i * D_i / mu_i
Pr_i = cp_i * mu_i / k_i
f_i  = dispatch_friction_coil(...)
Nu_i = dispatch_nu_coil(...)
h_i  = Nu_i * k_i / D_i
G_i  = h_i * P_coolant,total * dx_i
```

The property arrays are supplied by the future production adapter; the helper
does not call CoolProp itself. This keeps thermodynamic backend choice outside
the finite-volume core and makes validation tests cheap.

The helical wall-flux bridge wraps the existing quasi-static wall reconstruction
`fluxes_at_Tbar()`:

```text
single-tube output: dq_hot__dx, dq_cold__dx [W/m]
core-cell output:  Qhot_i  = dq_hot__dx  * dx_i * N_coils
                   Qcold_i = dq_cold__dx * dx_i * N_coils
```

It also returns the reconstructed face temperatures `T_wg`, `T_wc`, radiation
diagnostics, and wall conductivity. This is still not a full hot-side march: it
requires `T_g`, `h_g`, `h_c`, and optional `h_g_rad` arrays from the next
adapter layer.

The generic fixed-step integrator is intentionally small. For each interval
`[t_n, t_{n+1}]`, it calls:

```text
step_inputs(t_n, T_wall^n, T_coolant^n) -> WallCoolantStepInputs
```

and then applies the implicit wall/coolant step. This keeps the core execution
path independent of helical, shell-and-tube, CoolProp, chemistry, and schedule
implementations while still providing a real transient history:

```text
T_wall[t, i]
T_coolant[t, i]
T_coolant_outlet[t]
heat_wall_to_coolant[t, i]
energy_residual[t]
```

Production dispatch still needs to call this shared grid builder and integrate
the resulting histories into the result package/dashboard format.

Schedule interpolation is centralized in `transient_core.schedules`:

```text
interp_schedule(schedule, t, default)
```

It preserves the legacy transient behavior:

```text
schedule is None       -> default
t before first point   -> first value
t after last point     -> last value
inside schedule range  -> linear interpolation
```

The same module extracts schedule breakpoint times for `fixed_time_grid()`.

For the shell-and-tube adapter in `adapters_shelltube.py`, the grid coordinate
is straight tube length. The first hold-up model uses:

```text
A_shell_free = pi*D_shell_inner^2/4 - N_tubes*pi*D_tube_outer^2/4
A_wall,total = N_tubes*pi*(D_o^2 - D_i^2)/4
P_hot,total  = N_tubes*pi*D_i
P_cold,total = N_tubes*pi*D_o
```

Therefore:

```text
Cw_i = rho_tube * cp_tube,i * A_wall,total * dx_i
Cc_i = rho_shell,i * cp_shell,i * A_shell_free * dx_i
G_i  = h_shell,i * P_cold,total * dx_i
```

This is a deliberately simple inventory baseline. Bell-Delaware still governs
the shell-side film coefficient:

```text
G_s  = |mdot_shell| / S_m
Re_s = D_o * G_s / mu_shell
Pr_s = cp_shell * mu_shell / k_shell
h_shell = bell_delaware_shell(...).h_shell
G_i = h_shell_i * P_cold,total * dx_i
```

The transient-core adapter receives thermophysical properties as arrays and
passes `rho_s` through the Bell-Delaware geometry dict for pressure-drop
diagnostics.

The tube-side gas-film bridge receives gas properties from the future production
hot-side adapter and applies the maintained representative-tube convention:

```text
mdot_tube = |mdot_hot,total| / N_tubes
U_g       = mdot_tube / (rho_g * A_tube_inner)
Re_g      = rho_g * U_g * D_i / mu_g
Pr_g      = cp_g * mu_g / k_g
```

For smooth tubes and direct intensification-factor mode, the adapter calls the
same straight-tube friction and Nusselt dispatchers used by the shell-and-tube
solver. For grooved tubes it calls the Vicente/Cruz corrugated-tube friction and
Nusselt paths with:

```text
phi = corrugation_thickness^2 / (corrugation_pitch * D_i)
```

The resulting hot-side coefficient and conductance are:

```text
h_g     = Nu_g * k_g / D_i
G_hot_i = h_g,i * P_hot,total * dx_i
```

and the pressure-drop diagnostic is:

```text
dp/dx = f_g * rho_g * U_g^2 / (2*D_i)
```

The shell-and-tube hot-gas march scaffold uses those film coefficients inside a
sequential representative-tube enthalpy march. Thermochemistry is injected
through a gas-state provider:

```text
gas_state_at(h_removed_i, Yc_i, i) -> T_g, rho_g, mu_g, k_g, cp_g, omega_Yc
```

For each cell:

```text
dq_hot__dx_i = h_g,i * P_hot,single * (T_g,i - T_wg,i)

h_removed_{i+1}
  = h_removed_i + dq_hot__dx_i * dx_i / mdot_tube

Yc_{i+1}
  = Yc_i + omega_Yc,i * dx_i / U_g,i
```

The returned wall heat rates are converted to bundle-total cell watts:

```text
Qhot_i  = dq_hot__dx_i  * dx_i * N_tubes
Qcold_i = dq_cold__dx_i * dx_i * N_tubes
```

This scaffold is testable without building chemistry tables because the
gas-state provider is injectable. The transient-core shell-and-tube adapter
provides three wrappers:

```text
fpv_gas_state_provider(fpv)
equilibrium_gas_state_provider(manifold)
oxygen_gas_state_provider(T_inlet, pressure)
```

The FPV wrapper returns `omega_Yc` as the progress source, the
equilibrium/frozen wrapper returns zero progress source, and the oxygen wrapper
uses CoolProp to map `h_removed` through `(H, P)` before evaluating gas
properties. The shell-and-tube transient-core dispatch selects these providers
per timestep, so no-hot-flow or GOX chilldown can precede FPV/equilibrium
combustion in one scheduled run.

The shell-and-tube adapter now also provides:

```text
shelltube_step_inputs(...)
```

which assembles the per-step contract consumed by the generic wall/coolant
integrator:

```text
coolant_properties_at(T_c[i], p_c) -> rho_c, mu_c, k_c, cp_c
shelltube_shell_film(...)          -> h_shell_i, G_shell_i
shelltube_hot_gas_march(...)       -> Qhot_i, T_g_i, h_g_i, h_removed_i

Cw_i = rho_wall * cp_wall_i * V_wall_i
Cc_i = rho_c_i  * cp_c_i    * V_coolant_i

WallCoolantStepInputs(
    Cw_i,
    Cc_i,
    cp_c_i,
    mdot_c,
    T_c,in,
    Qhot_i,
    G_shell_i,
    flow_direction,
)
```

This is the adapter-level production dispatch. Top-level shell-and-tube
selection is now wired through `input_data.py`/`main_transient.py`; remaining
work is dashboard semantic review, broader short-run validation, and shell-side
residence refinement.

The adapter-level runner is:

```text
run_shelltube_transient_core(...)
```

It builds a fixed-step time grid with schedule breakpoints:

```text
t_grid = fixed_time_grid(t_end, max_step, schedules, t_eval)
```

and at each interval start evaluates scheduled helium and hot-side boundary
conditions before calling `shelltube_step_inputs()`. It returns:

```text
ShellTubeTransientCoreResult(
    integration=WallCoolantIntegrationResult(...),
    step_diagnostics=(ShellTubeStepInputDiagnostics, ...),
)
```

The shell-and-tube solver bridge stores both `core_result` and a normal
`time_series` dictionary, so existing transient packaging can write the new
wall/coolant fields. A 0.25 s finite-rate combustion smoke run and a 0.25 s
cold-He no-ignition smoke run have completed through `main_transient.py`.

The shell-and-tube wall-flux bridge wraps the same quasi-static reconstruction
as the legacy transient solver, but with the shell-and-tube orientation made
explicit:

```text
hot_side = "inner"
P_hot    = pi*D_i
P_cold   = pi*D_o
```

The conduction object returns representative-tube per-length rates:

```text
dq_hot__dx, dq_cold__dx [W/m per tube]
```

The transient core stores total cell heat rates:

```text
Qhot_i  = dq_hot__dx_i  * dx_i * N_tubes
Qcold_i = dq_cold__dx_i * dx_i * N_tubes
```

This matches the maintained shell-and-tube transient convention where the gas
path is solved for one representative tube and the shell side receives the
bundle-total duty. Bell-Delaware leakage/window geometry is not yet used to
refine transient shell-side hold-up volume.

The pure coolant finite-volume kernel solves:

```text
rho_c,i * V_c,i * cp_c,i * (Tc_i^{n+1} - Tc_i^n) / dt =
    mdot_c * cp_c,i * (Tup_i^{n+1} - Tc_i^{n+1})
  + Qc_i
```

where `Tup` is the inlet boundary for the first upwind cell or the newly solved
upstream cell. At `abs(mdot_c) <= mdot_floor`, the advection term is removed and
the coolant behaves as stagnant local thermal inventory. This is the required
behavior for bang-bang helium shutoff; the inlet condition must not instantly
overwrite coolant already inside the exchanger.

The coupled wall/coolant kernel solves a local implicit 2x2 system per cell:

```text
Cw_i * (Tw_i^{n+1} - Tw_i^n) / dt =
    Qhot_i - G_i * (Tw_i^{n+1} - Tc_i^{n+1})

Cc_i * (Tc_i^{n+1} - Tc_i^n) / dt =
    mdot_c * cp_c,i * (Tup_i^{n+1} - Tc_i^{n+1})
  + G_i * (Tw_i^{n+1} - Tc_i^{n+1})
```

with:

```text
Cw_i = wall heat capacity in cell i [J/K]
Cc_i = coolant heat capacity in cell i [J/K]
G_i  = total wall-mean-to-coolant-bulk conductance [W/K]
Qhot_i = quasi-steady hot-side heat into the wall [W]
```

The global constant-property first-law residual checked by tests is:

```text
R_E =
  (U_wall^{n+1} + U_coolant^{n+1}
 - U_wall^n     - U_coolant^n)
 - (Qhot_total*dt + H_in*dt - H_out*dt)
```

For the tested linear cases, `R_E` closes near floating-point roundoff. The
implemented diagnostic scales this residual as:

```text
epsilon_E = |R_E| /
  max(|Qhot_total*dt|, |H_in*dt|, |H_out*dt|, |Delta U|, floor)
```

The shared timescale diagnostics are:

```text
tau_coolant = sum_i(rho_c,i * V_c,i) / abs(mdot_c)

tau_wall,i = Cw_i / (G_hot_i + G_cold_i)

coolant quasi-steady ratio = tau_coolant / min_i(tau_wall,i)
hot quasi-steady ratio     = tau_hot / min_i(tau_wall,i)
control-forcing ratio      = tau_residence / tau_boundary
```

These are diagnostic ratios, not solver stability limits. For the new transient
coolant model, a large coolant ratio is expected during bang-bang flow and is
handled by the coolant state. A large hot-gas ratio would challenge the
quasi-steady hot-side assumption and is the trigger for considering the optional
transient hot-gas extension.

The remaining work is not another local integrator; it is the geometry adapters
that compute `Qhot_i`, `G_i`, heat capacities, pressure/friction updates, face
temperatures, and packaged diagnostics from the existing helical and
shell-and-tube correlations.

## 12.2 Counter-Flow

Co-flow is an initial-value march.

Counter-flow is a two-boundary problem because helium enters at the opposite end
from the hot gas. The transient solvers use a warm-started coolant-temperature
profile relaxation:

1. assume a coolant profile;
2. march gas and wall fluxes forward;
3. march coolant backward from its physical inlet using the cold-side duty;
4. under-relax the profile;
5. cache the profile for the next call.

For production transient runs, this relaxation is intentionally not fully
converged on every time step. A few initial iterations and one warm-start
iteration per later step are used to keep cost bounded.

Mathematically, co-flow is an initial-value problem:

```text
T_h(0) = T_h,in
T_c(0) = T_c,in
dT_h/dx = -dq/dx / (mdot_h*cp_h)
dT_c/dx = +dq/dx / (mdot_c*cp_c)
```

Counter-flow has the coolant boundary at the opposite end:

```text
T_h(0) = T_h,in
T_c(L) = T_c,in
```

The relaxation converts the boundary-value problem into repeated initial-value
marches. At iteration `k`:

```text
given T_c^k(x), march hot side forward and compute q^k(x)
then march coolant backward from x=L using q^k(x)
T_c^{k+1}(x) = (1 - omega)*T_c^k(x) + omega*T_c,backward^k(x)
```

The warm-start cache works because, for small time steps, the converged coolant
profile at `t_n` is close to the profile at `t_{n+1}`. This is why one
relaxation pass per later fixed step can be accurate enough after an initial
settling pass.

## 12.3 Fixed-Step Wall Integrator

The production transient method is:

```text
transientProp.solver_method = "fixed_step"
```

It is linearly implicit in the local wall-film stiffness:

```text
Tbar_next = Tbar + dt * R(Tbar) / (1 + dt * lambda)
```

where:

```text
R      = dTbar/dt from the current fluid pass
lambda = (h_hot * P_hot + h_cold * P_cold) / (rho_w * cp_w * A_wall)
```

This handles the stiff helium film without BDF's hidden Jacobian-probing cost.
The time grid also inserts schedule breakpoints so ignition and ramp changes are
not stepped across.

The stability rationale is local. If fluid temperatures and coefficients are
temporarily frozen over one step, the wall equation near the current state has
the form:

```text
dTbar/dt = S - lambda*Tbar
```

where `S` contains the neighboring fluid-temperature terms. Forward Euler gives:

```text
Tbar_{n+1} = Tbar_n + dt*(S - lambda*Tbar_n)
```

which is stable only when roughly:

```text
dt * lambda < 2
```

For high-pressure helium, `h_c` can make `1/lambda` much smaller than practical
communication or output steps. Treating only the diagonal linear stiffness
implicitly gives:

```text
Tbar_{n+1}
  = Tbar_n + dt*R_n/(1 + dt*lambda_n)
```

where `R_n` is the full RHS evaluated from the current quasi-steady fluid pass.
The damping factor:

```text
1/(1 + dt*lambda_n)
```

is bounded for any positive `dt`, so the stiff local film coupling cannot cause
the explicit blow-up seen with plain Euler. This does not make time-discretized
physics magically exact; it makes the production method stable and bounded-cost
so `max_step` convergence can be assessed directly.

BDF/Radau/RK methods remain available for validation, but BDF can be very slow
for counter-flow because one RHS evaluation includes full quasi-steady fluid
marching and profile relaxation.

## 12.4 Step Size

Default:

```text
max_step = 0.25 s
```

For production, compare at least one short case with:

```text
max_step = 0.25 s
max_step = 0.10 s
```

The early helium outlet temperature may remain transient-sensitive, especially
during warm-up. Heat duty is generally a better short smoke convergence metric.

Useful numerical error metrics:

```text
energy_mismatch = abs(Q_hot - Q_cold) / max(abs(Q_hot), abs(Q_cold), eps)
step_delta_Q    = abs(Q_0p25 - Q_0p10) / max(abs(Q_0p10), eps)
settle_error_Q  = abs(Q_transient_final - Q_steady) / max(abs(Q_steady), eps)
wall_delta_Tmax = max(abs(Tbar_0p25 - Tbar_0p10))
```

For counter-flow, also compare against the physical steady counter-flow
reference, not the legacy prescribed-outlet guess path.

## 13. Steady Numerics

## 13.1 Helical Steady

The helical steady solver marches along coil arc length until one of the stop
conditions is reached:

- target/maximum HX length;
- coolant temperature safety floor;
- pressure/stress constraints in post-checks.

The wall is solved at each node with a local conduction solve. The gas state is
updated by enthalpy removal and chemistry mode.

At node `i`, the steady helical march is conceptually:

```text
state_i = {T_c, p_c, h_removed, Yc, wall guesses}
properties_i = properties(state_i)
h_c, f_c, h_g, h_rad = correlations(properties_i, geometry_i)
dQ_i = UA_i * (T_g_i - T_c_i)

T_c_{i+1} = T_c_i + dQ_i / (mdot_c*cp_c)
p_c_{i+1} = p_c_i - f_c*rho_c*U_c^2/(2*D_h)*dx
h_removed_{i+1} = h_removed_i + dQ_i/mdot_g
Yc_{i+1} = Yc_i + omega_Yc_i*dx/U_g       finite_rate only
```

This is a marching model, not a global nonlinear solve. The main exception is
physical counter-flow steady reference, where the hot-end coolant condition is
shot so the cold-end inlet equals the user-specified helium inlet.

## 13.2 Shell-And-Tube Steady

The shell-and-tube steady solver uses sweep iteration because shell-side
temperature and tube-side heat duty are coupled globally:

```text
T_shell_guess -> tube march -> duty profile -> shell march -> relaxed T_shell
```

The convergence criterion is maximum shell-profile change.

The gas-side tube march computes one representative tube:

```text
mdot_tube = mdot_hot_total / N_tubes
dQ_tube_i = UA_i * (T_g_i - T_shell_i)
```

The shell-side helium march uses the total bundle duty:

```text
dQ_shell_i = N_tubes * dQ_tube_i
T_shell_{i+1} = T_shell_i + dQ_shell_i/(mdot_shell*cp_shell)
```

For counter-flow, the shell march direction is opposite the tube-index
direction, but reported arrays keep the same physical axial index convention.
The under-relaxed iteration is:

```text
T_shell^{k+1} = (1 - omega)*T_shell^k + omega*T_shell,new^k
```

This avoids a dense coupled solve while still capturing the global coupling
introduced by shell-side flow direction.

## 14. Materials And Mechanical Checks

Material functions are dispatched in:

```text
mechanical/material_specs/material_temperature_strength.py
```

Supported maintained materials include:

- `ST316L`;
- `INCO718`.

Properties include:

- thermal conductivity;
- specific heat where needed;
- density;
- coefficient of thermal expansion;
- elastic modulus;
- yield strength.

Helical mechanical checks include:

- pressure stress in tube wall;
- thermal stress;
- stress/yield ratio;
- mass estimates.

The first-order stresses are screening estimates. Pressure stress for a thin
tube is interpreted with the usual thin-wall scaling:

```text
sigma_hoop ~ p * r / s
sigma_axial ~ p * r / (2*s)
```

Thermal stress is estimated from constrained thermal strain scaling:

```text
sigma_thermal ~ E(T) * alpha(T) * DeltaT / (1 - nu)
```

where `DeltaT` is a representative wall or through-wall temperature difference.
The checks are useful for rapid design filtering, but final tube-wall,
manifold, tube-sheet, and baffle stresses still require a dedicated structural
model.

Shell-and-tube mechanical checks include:

- external-pressure stress on tubes;
- thin-tube collapse pressure;
- thermal stress estimate;
- tube wall temperature limits.

Not all EchTherm-style structural fields are fully active yet. In particular,
`shell_thickness`, `tube_sheet_thickness`, and nozzle geometry are represented
in inputs but still need deeper mass/mechanical integration.

## 15. Meshing And Resolution

Helical:

- steady spatial step is derived from coil geometry and
  `N_arc_steps_per_turn`;
- transient grid is coarsened to `transientProp.n_axial`.

Shell-and-tube:

- steady axial nodes from `runProp.shelltube_steady_nodes`;
- transient axial nodes from `runProp.shelltube_transient_nodes`.

Typical values:

- shell-and-tube steady: 200 nodes;
- shell-and-tube transient: 80 nodes;
- shell-and-tube smoke: 16 nodes.

The transient wall field is smooth enough for moderate node counts, but
step-convergence and grid-convergence should be checked for final studies.

## 16. Outputs And Packaging

User-facing runs package:

- input dataclass snapshot;
- summary JSON;
- numeric `.npz` data;
- `HX_performance_summary.txt` engineering report;
- optional CSV tables;
- zip archive if enabled.

Transient outputs include:

- time vector;
- axial grid;
- wall mean temperature field;
- hot and cold wall-face temperatures;
- gas and helium temperatures;
- hot/cold heat fluxes;
- scalar time histories such as outlet temperatures and total duty;
- shell-and-tube `transient_coolant` engineering fields such as gas/shell
  Reynolds numbers, Nusselt number, friction, pressure-drop estimates, coolant
  properties, gas enthalpy removed, and FPV progress variable.

Wall radial/axial transient temperature evolution is stored in the same raw
transient file when produced by the solver:

```text
field_T_wg   hot-side wall face temperature [time, x]
field_Tbar   thickness-mean wall temperature [time, x]
field_T_wc   cold-side wall face temperature [time, x]
```

`HX_performance_summary.txt` is intentionally text, not another raw data file.
It summarizes heat duty, LMTD, coolant-capacity-referenced NTU/effectiveness,
temperature extrema, pressure drops, Reynolds ranges, and available
combustion/mechanical histories at selected transient points and over the run.

Packaged transient data can be reloaded by:

```powershell
hps-dashboard path\to\run_folder
hps-dashboard path\to\transient_timeseries.npz
hps-dashboard path\to\run.zip
```

## 17. Calibration

Calibration factors live in `CorrelationCoefficients`.

Important knobs:

- `ali_c_hi`: helical coolant pressure drop;
- `salimpour_a`: helical hot-side heat transfer;
- `mori_a_lo`: helical coolant heat transfer;
- `mbl_factor`, `emissivity_wall`: radiation;
- `zukauskas_C_factor`, `bell_Jl_factor`, `bell_Jb_factor`: shell-side
  Bell-Delaware calibration;
- `tube_grooved_Nu_factor`, `tube_grooved_f_factor`: grooved tube calibration.

The calibration methodology should keep published correlation shapes intact and
fit a small number of interpretable multipliers against measured data.

## 18. Validation Status

Recent targeted checks:

- Shell-and-tube steady vs transient wall reconstruction, 16 nodes:
  co/counter x frozen/equilibrium/finite-rate matched heat rate within about
  0.4% or better after the `hot_side="inner"` steady conduction fix.
- Shell-and-tube finite-rate fixed-step, 16 nodes, 5 s:
  `max_step=0.25` vs `0.10` changed `Q_hot` by about 0.5% for co and counter.
- Shell-and-tube counter-flow finite-rate fixed-step, 80 nodes, 5 s:
  about 5.8 s wall time and 26 fluid passes with cached FPV chemistry.
- Previous BDF benchmark on the same 80-node/5 s shell counter-flow case:
  about 134 s and 456 RHS calls.

Known remaining validation work:

- full helical x shell-and-tube, co x counter, 3 chemistry mode validation
  matrix on production grids;
- long 100 s transient benchmark after final schedule choice;
- GOX chilldown energy-balance check against a simple oxygen sensible-heating
  estimate;
- independent grooved/corrugated-tube reference check;
- dashboard visual verification on packaged transient outputs.

Liquid coolant / boiling (see Section 7.5) validation status, from
`docs/validation/liquid_validation_matrix.json` (`all_passed: true` for what
is validated) and `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`:

- Validated: single-phase liquid convection/friction; saturated HEM state;
  Gungor-Winterton HTC; Müller-Steinhagen-Heck friction; HEM acceleration
  gradient; Groeneveld 2006 CHF LUT + diameter correction (page-9 spot check,
  0 error); Yu et al. 2002 pressure-multiplier fit (2.7% MAE) and HTC fit
  (4.95% MAE) against representative digitized plot points; a synthetic-fixture
  postprocess audit of both solver adapters' duty conservation and counterflow
  mapping.
- **Note**: `liquid_validation_matrix.json`'s `scope.not_yet_validated` list
  predates Phase 2 of the integration plan and still literally lists "the
  fully coupled production liquid wall march" as not done; that is now true
  only for shell-and-tube. The helical steady coupled march is wired and
  covered by `tests/test_liquid_coupled_helical.py` (subcooled/boiling,
  co/counter, grid convergence, cross-check against the postprocess bridge);
  the JSON matrix itself has not been regenerated to reflect this.
- Explicitly **still not validated**: the shell-and-tube coupled production
  liquid wall march (Phase 3); the transient boiling/liquid finite-volume
  model (no transient liquid code path exists at all); geometry-specific
  shell-side or helical-coil boiling correlations (current adapters are
  pseudo-1D placeholders); dryout/post-CHF model validation. The
  counter-flow prescribed-outlet-vs-physical-inlet discrepancy (Design
  Decision 2.4) is now resolved for the plain solve path via
  `solve_counterflow_liquid_reference()` (enthalpy-shooting, adaptive
  bracket + bisection) — see the integration plan's "Post-Phase-2 Hardening
  Pass."
- Digitization risk: the Yu et al. 2002 comparison points are digitized from
  plots, not source tabular data, and are flagged as the primary remaining
  blocker before treating those fits as production-grade.

## 19. Assumptions And Limits

Core assumptions:

- one-dimensional axial or arc-length flow;
- quasi-steady fluids during transient wall integration;
- one wall temperature state per axial node;
- representative tube for shell-and-tube gas side;
- uniform distribution across tubes;
- empirical heat-transfer and pressure-drop correlations outside parts of their
  original data range;
- no shell-side transient pressure dynamics;
- no two-phase LOX evaporation model (injector/chilldown side, distinct from
  the coolant-side liquid work below);
- no detailed injector or combustion instability model;
- no full finite-volume gas momentum transient;
- the experimental liquid/boiling coolant path (Section 7.5) is wired into
  the helical steady coupled march only (counter-flow has a known
  prescribed-outlet limitation); it remains postprocess-only for
  shell-and-tube and absent from both transient solvers; and it has no
  dryout/post-CHF, flow-regime, or geometry-specific (helical/shell-side)
  boiling model yet.

Use the model for engineering screening, design iteration, calibration, and
system coupling where a fast 1D surrogate is appropriate. Do not treat it as a
replacement for CFD, structural FEA, or experimental qualification.

## 20. Practical Run Recommendations

For shell-and-tube transient production runs:

```python
combustorProp.HX_config = "shellntube"
combustorProp.flow_config = "counter"  # or "co"
transientProp.chemistry_transient = "finite_rate"
transientProp.solver_method = "fixed_step"
transientProp.max_step = 0.25
runProp.shelltube_transient_nodes = 80
```

Before a long run:

1. Confirm the FPV cache exists or accept the first-build cost.
2. Run a 5 s smoke case.
3. Compare `max_step=0.25` and `0.10` on the short case.
4. Check heat duty, outlet temperatures, max wall temperature, and wall face
   temperature difference.
5. Package numeric outputs so the dashboard can be regenerated later.

For validation runs:

- use `BDF` only when specifically comparing integrator behavior;
- use reduced grids first;
- avoid judging early helium outlet temperature during low-flow ramps without
  checking residence-time reliability.

## 21. Simulink/System Coupling Interpretation

For coupling to a larger system solver, this code can act as a component model.
At each communication time:

Inputs:

- helium inlet mass flow, temperature, and pressure;
- LOX/GOX mass flow and inlet temperature before ignition;
- diesel/LOX or total hot-gas flow after ignition;
- ignition state;
- optional O/F.

Outputs:

- helium outlet temperature;
- helium pressure estimate where available;
- gas/oxygen outlet temperature;
- heat duty;
- maximum wall/material temperature;
- stress/yield indicators where available.

Internally the component may use its own fixed-step subcycling. Simulink should
not assume the internal wall time step equals the external communication step.

## 22. Files To Read Before Modifying Core Physics

Read in this order:

1. `CLAUDE.md`
2. `docs/TECHNICAL_REFERENCE.md`
3. `docs/context/SOLVER_CONTEXT.md`
4. `docs/context/PHYSICS_CONTEXT.md`
5. `docs/context/TRANSIENT_STATUS.md`
6. `1Dmodel/input_data.py`
7. The specific solver backend being modified.

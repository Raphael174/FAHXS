# Claude Handoff - 2026-07-08

This file is the short operational handoff for the current solver state. Read
this before changing shell-and-tube geometry, transient scheduling, or chemistry
defaults.

## What Works Now

- User-facing run entry points:
  - `hps-steady`
  - `hps-transient`
  - `hps-dashboard`
- Both HX configurations are runnable:
  - `shellnHelicalTube`
  - `shellntube`
- Both steady and transient default to finite-rate chemistry.
- Packaged transient outputs can be reloaded into an HTML dashboard from:
  - a run folder
  - a `transient_timeseries.npz`
  - a zipped run archive

## Current Chemistry Truth

- Steady default:
  - `numericalProp.chemistry_model = "finite_rate"`
- Transient default:
  - `transientProp.chemistry_transient = "finite_rate"`
- The finite-rate path uses:
  - `1Dmodel/physics/combustion_chemistry/fpv_manifold.py`
- `equilibrium` and `frozen` still exist as comparison modes.

Important distinction:
- The older docs used to say steady finite-rate was not implemented.
- That is no longer true for the maintained steady helical and steady
  shell-and-tube paths.

## Schedule / Dashboard Utilities Added

- Schedule loader now accepts the human-facing single-table format:
  - `time_s`
  - `helium_m_dot_kg_s`
  - `helium_T_in_K`
  - `helium_p_in_Pa`
  - `diesel_m_dot_kg_s`
  - `lox_m_dot_kg_s`
  - `lox_T_in_K`
  - `ignition`
- Example file:
  - `inputs/example_transient_schedule.csv`
- Dashboard reload command:
  - `hps-dashboard <run_folder|npz|zip>`

Notes:
- `lox_T_in_K` and `ignition` are stored and parsed.
- `ignition` currently drives `ignition_time`.
- Cold gaseous oxygen hot-side chilldown is still not a real separate physical
  hot-side model.

## Shell-And-Tube EchTherm Alignment

The shell-and-tube input model was updated to better match the EchTherm geometry
screen.

Added / separated fields on `shellTubeProp`:
- `shell_thickness`
- `tube_sheet_thickness`
- `baffle_spacing`
- `L_inlet_spacing`
- `L_outlet_spacing`
- `L_front_end`
- `L_rear_end`
- nozzle diameter/orientation fields
- corrugation geometry fields
- `inside_tube_choice`
- `outside_tube_choice`

Important fix:
- `L_front_end` and `L_rear_end` are end zones.
- They are no longer misused as Bell-Delaware inlet/outlet spacing.
- Bell-Delaware now uses:
  - `baffle_spacing`
  - `L_inlet_spacing`
  - `L_outlet_spacing`

## What Is Physically Active vs Stored

Physically active now:
- shell-and-tube core geometry:
  - tube OD
  - tube wall thickness
  - tube length
  - tube count
  - layout
  - pitch ratio
  - shell ID
  - baffle cut
  - baffle count
  - baffle spacing
  - inlet/outlet spacing
  - clearances
  - sealing strips
- tube material
- finite-rate FPV chemistry

Active tube-side physics / calibration hooks:
- `inside_tube_choice = "grooved"`
- `tube_grooved_Nu_factor`
- `tube_grooved_f_factor`
- `tube_intensification_factor`

Stored but not yet connected to physics:
- `shell_thickness`
- `tube_sheet_thickness`
- nozzle geometry/orientation fields
- `outside_tube_choice = "low_finned"`

Do not misrepresent those last fields as already affecting h/dp/stress unless
you wire them into the model.

## Tube-Side Correlation Truth

The EchTherm screenshot had `Grooved Tube` selected.

Current implementation status:
- `inside_tube_choice = "grooved"` uses a documented helically corrugated-tube
  Nu/friction path based on the Vicente/Cruz form.
- Corrugation severity is `phi = e^2/(p*D_i)` from:
  - `corrugation_thickness`
  - `corrugation_pitch`
  - tube internal diameter
- The old multipliers remain as optional calibration factors:
  - `tube_grooved_Nu_factor`
  - `tube_grooved_f_factor`
- Defaults are `1.0`.

## Validation Runs Performed

Steady / transient smoke runs completed in this workspace during this session:

- Helical steady run packaged successfully.
- Shell-and-tube steady run packaged successfully.
- Helical transient run packaged successfully.
- Shell-and-tube transient run packaged successfully.

Finite-rate validation after the recent changes:

- Short helical steady finite-rate smoke run: passed.
- Short shell-and-tube steady finite-rate smoke run: passed.
- Short helical transient finite-rate run on normal grid: passed.
- Short shell-and-tube transient finite-rate validation run: passed.
- Short shell-and-tube counter-flow transient finite-rate smoke run: passed.

Validation artifact:
- `results/validation_runs/validation_shelltube_finite_rate_echtherm_0p1s_transient_08-07-2026_11h47.zip`

Reported final values from that short shell-and-tube transient validation:
- `T_g_out = 2589.3 K`
- `T_c_out = 306.9 K`
- `Q = 413.4 kW`
- `wall dT max = 40.8 K`

## Known Caveats

- Shell-and-tube transient now accepts `flow_config in ("co", "counter")`.
- Helical counter-flow transient support exists. A physical steady-reference
  helper is available behind
  `numericalProp.counterflow_physical_steady_reference`; it shoots the legacy
  steady solver's hot-end helium temperature so the cold-end boundary matches
  `coolantProp.T_in`. A reduced frozen-chemistry check converged that boundary
  to 0.33 K and matched transient-pass heat rate within about 0.6%, but the
  full finite-rate settle-to-steady matrix is not complete yet.
- An intentionally tiny coarse-grid helical transient smoke run produced
  unphysical helium temperatures during the startup probe. The normal-grid short
  transient run worked. Treat that as an early-ramp stiffness limitation, not a
  solved issue.
- FPV manifold build still emits Cantera warnings above 3000 K for the diesel
  mechanism.
- Finite-rate FPV manifolds are cached in `cache/fpv_manifolds`. First build of
  the current reduced transient table was about 33 s; cached setup was about
  0.9 s.
- Shell-and-tube counter-flow transient should use
  `transientProp.solver_method = "fixed_step"` for production long runs. BDF is
  retained for validation, but it was too slow because Jacobian probing calls
  the expensive counter-flow RHS many times. The fixed-step path is now
  linearly implicit in the local wall-film stiffness.
- Benchmark after switching to linearly-implicit fixed-step, cached finite-rate
  chemistry, shell-and-tube counter-flow, 80 nodes, 5 s simulated: 5.8 s wall
  time, 26 fluid passes. Earlier BDF benchmark on the same 80-node/5 s case was
  about 134 s and 456 RHS calls.
- Important correctness fix: `Solve1Dconduction()` now honors
  `hot_side="inner"` for shell-and-tube. Before this, the steady solver used the
  wrong hot/cold perimeters and disagreed with transient reconstruction by
  10-25%. After the fix, 16-node shell-and-tube steady/transient wall
  reconstruction agrees within about 0.4% or better for co/counter and
  frozen/equilibrium/finite-rate.
- The `.git` metadata in this workspace has behaved oddly before; do not assume
  `git status` is reliable until checked.

## Immediate Next Good Tasks

1. Finish the helical counter-flow settle-to-steady validation on the normal
   finite-rate setup.
2. Wire `shell_thickness` and `tube_sheet_thickness` into mass/mechanical logic.
3. Validate the new grooved/corrugated tube correlation against an independent
   reference case or experiment before calibration.
4. Validate the pre-combustion GOX chilldown branch against a simple oxygen
   sensible-heating energy balance.
5. Run a broader finite-rate validation matrix for all four solver combinations
   after the geometry/correlation changes.

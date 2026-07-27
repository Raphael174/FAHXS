# Combustor-HX User Guide

This project is a 1D heat-exchanger solver for combustor-coupled helium heating.
The intended user entry points are:

```powershell
hps-steady
hps-transient
hps-dashboard
```

The solver behavior is configured in `1Dmodel/input_data.py`. Transient boundary
conditions can also be supplied with a CSV or Excel schedule file.

## 0. Verified Quick Start (Read This First)

Every command below was run and confirmed working on 2026-07-13 against the
project's current default `input_data.py` (Helium coolant,
`combustorProp.HX_config = "shellntube"`). If `hps-steady` itself is not
found, that is almost always environmental (install/PATH) — see Section 0.4
before assuming the solver code is broken.

### 0.1 Steady, default config (shell-and-tube, Helium)

```powershell
hps-steady
```

Expected: 10-20 sweep lines converging (`max|dT_shell|` shrinking toward
`< 0.05 K`), then a `SHELL-AND-TUBE RESULTS` block with `Q_tot` a few hundred
kW, followed by `Saved steady run: folder: zip_folders\...`. Confirmed output
signature (your exact numbers may differ slightly with future correlation
changes, but the shape/order of magnitude should match):

```text
  sweep 14: max|dT_shell|=  0.037 K  T_g_out=2589.3  T_c_out=1067.6
=======================================================
SHELL-AND-TUBE RESULTS
=======================================================
  N_tubes=235  N_baffles=15  sweeps=15
  Q_tot = 595.09 kW
  ...
Saved steady run:
  folder:  zip_folders\combustor_hx_..._steady_...
  archive: ...\zip_folders\combustor_hx_..._steady_....zip
```

To run the helical-coil geometry instead (still Helium), edit
`combustorProp.HX_config = "shellnHelicalTube"` in `1Dmodel/input_data.py`, or
run it directly without touching the shared file:

```powershell
python -c "from hps_combustor.input_data import *; from hps_combustor.main_solve import main_solver; s = main_solver(coolantProp(), hotgasProp(), combustorProp(HX_config='shellnHelicalTube'), numericalProp(), system_requirements(), CorrelationCoefficients()); s.solver(); s.compute_performance(); s.print_summary()"
```

Expected: a `RESULTS` block (not `SHELL-AND-TUBE RESULTS`) with `Q_tot` a few
hundred kW and `max sigma/yield` comfortably below 1.0.

### 0.2 Transient, default config

```powershell
hps-transient
```

Expected (default `t_end=10s`):

```text
Integrating shell-and-tube wall ODE: 80 nodes, t_end=10s, method=fixed_step, max_step=0.25s ...
  final: T_g_out=2589.3 K  T_c_out=1036.0 K  Q=595.1 kW  wall dT max=94.0 K

Saved transient run:
  folder:  zip_folders\combustor_hx_..._transient_...
```

### 0.3 Dashboard

```powershell
hps-dashboard zip_folders\<the folder or .zip printed above>
```

Expected: `Transient dashboard written: ...html` (or the steady equivalent).
Open that `.html` file in a browser — no server needed.

### 0.4 If a command is not found

`hps-steady` / `hps-transient` / `hps-dashboard` are console scripts installed
by `python -m pip install -e .` into whichever Python environment you ran that
in (commonly `.venv/Scripts/` on Windows). If your shell reports "command not
found", either activate that environment first
(`.\.venv\Scripts\Activate.ps1`) or call the module form directly, which needs
no console-script PATH entry at all:

```powershell
python -m hps_combustor.main_steady
python -m hps_combustor.main_transient
```

### 0.5 Water / boiling coolant test configuration

Water (`coolant_model = "equilibrium_liquid"`) is an **experimental test
configuration**, not the project default — Helium/`single_phase_coolprop`
remains the working baseline in `input_data.py`, on purpose (do not edit the
shared `coolantProp`/`combustorProp` defaults to switch the project to water;
see Section 0.6 for why that broke things before). Use the confirmed-working,
standalone recipe instead:

```powershell
python -c "from hps_combustor.validation.water_helical_example import run_coflow; run_coflow()"
```

Expected: an ordinary `RESULTS` block, `Q_tot` around 285-300 kW, ending with
no `WARNING` lines (the sanity gates pass). This is co-flow, which needs no
guessed boundary condition. For counter-flow with water, use
`run_counterflow_physical` from the same module instead — it is much slower
(several minutes; see the module docstring) because it has to search for the
physically correct starting condition rather than just marching forward.

### 0.6 Known Gotchas (read before editing `input_data.py` coolant/combustor fields)

These are real defects this project hit and fixed; both are guarded against
now, but understanding them will save you time if you see similar symptoms
after editing shared config:

- **`combustorProp.HX_config` must match the solver class you construct.**
  `main_solver` (helical) requires `HX_config == "shellnHelicalTube"` and now
  raises a clear `ValueError` if it is anything else — it used to silently
  compute the wrong axial length instead (a ~14x error in effective coil
  length went undetected for a long time). `shellntube_solver` is chosen by
  `main_steady.py`/`main_transient.py` when `HX_config == "shellntube"`, the
  project default.
- **`coolantProp.T_out`/`p_out` are legacy single-phase-gas fields.** They
  only make sense for Helium in the current maintained solvers. If you
  experiment with `coolant_model = "equilibrium_liquid"` (water) and these
  correspond to a state past complete vaporization at that pressure (this
  combination has recurred multiple times: `coolant="Water"` with the
  Helium-tuned `T_out=650`, `p_out=13e5`), `main_solver` now raises a clear
  `ValueError` at construction time naming the exact field to check and the
  computed quality, instead of the old confusing `numpy` "zero-size array"
  crash inside `compute_performance()`. If you see that specific `ValueError`
  ("Liquid coolant march cannot start..."), the fix is exactly what it says:
  either give `T_out`/`p_out` a genuinely subcooled/saturated liquid state at
  that pressure, or (for counter-flow) use
  `solve_counterflow_liquid_reference()` instead of guessing them at all. The
  water recipe in Section 0.5 sidesteps this entirely by not depending on
  `T_out`/`p_out`.
- If you ever cannot run `hps-steady` with an unmodified `input_data.py`,
  that is a real regression — re-run Section 0.1's command verbatim and
  compare against the expected output above before debugging anything else.
  If `input_data.py`'s `coolantProp` shows `coolant="Water"` /
  `coolant_model="equilibrium_liquid"` instead of the Helium baseline, that
  by itself is not a bug — it means the shared defaults were edited for
  local experimentation — but combined with unchanged `T_out=650`/
  `p_out=13e5` it will hit the gotcha above.

## 1. Install On A Machine

Use Python 3.10 or newer. The current local environment has been tested with
Python 3.13.

From the repository root:

```powershell
python -m venv --system-site-packages .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e .
```

The `--system-site-packages` flag is useful on machines where Cantera, CoolProp,
NumPy, or SciPy are already installed globally. If you prefer an isolated venv,
install those dependencies directly into `.venv`.

Verify the install:

```powershell
python -c "import hps_combustor, numpy, scipy, cantera, CoolProp; print('ok')"
```

If PowerShell refuses to activate the venv, run:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

Then open a new terminal and activate again.

## 2. Choose The Simulation

Edit `1Dmodel/input_data.py`.

For the HX geometry:

```python
class combustorProp:
    HX_config = "shellnHelicalTube"  # or "shellntube"
    flow_config = "counter"          # or "co"
```

Shell-and-tube geometry mirrors the EchTherm geometry page. The most important
fields are:

```python
class shellTubeProp:
    D_tube_outer = 5e-3
    thickness_tube_wall = 0.75e-3
    L_tube = 235e-3
    N_tubes = 235
    layout = "triangular30"
    pitch_ratio = 1.3

    D_shell_inner = 110e-3
    shell_thickness = 8e-3
    tube_sheet_thickness = 3e-3

    N_baffles = 15
    baffle_cut = 0.20
    baffle_thickness = 3e-3
    baffle_spacing = 12e-3
    L_inlet_spacing = 8e-3
    L_outlet_spacing = 8e-3
    L_front_end = 100e-3
    L_rear_end = 10e-3

    inside_tube_choice = "grooved"
    outside_tube_choice = "smooth"
```

`L_front_end` and `L_rear_end` are end-zone/nozzle geometry. They are not used as
Bell-Delaware baffle spacing corrections. The shell-side Bell-Delaware model uses
`baffle_spacing`, `L_inlet_spacing`, and `L_outlet_spacing`.

The EchTherm grooved-tube selection is exposed through `inside_tube_choice`.
For `inside_tube_choice = "grooved"`, tube-side Nu and pressure-drop now use a
Vicente/Cruz-style helically corrugated-tube correlation with severity
`phi = e^2 / (p * D_i)`, using:

- `corrugation_thickness`
- `corrugation_pitch`
- tube internal diameter

The grooved multipliers remain as optional calibration factors on top of the
published correlation:

```python
class CorrelationCoefficients:
    tube_grooved_Nu_factor = 1.0
    tube_grooved_f_factor = 1.0
```

Leave those at `1.0` unless you have calibration data.

For the run metadata and output:

```python
class runProp:
    run_name = "my_case"
    output_root = "zip_folders"
    make_archive = True
    save_csv = True
    schedule_file = "inputs/example_transient_schedule.csv"
```

For transient time controls:

```python
class transientProp:
    t_end = 100.0
    n_save = 120
    max_step = 0.25
    solver_method = "fixed_step"
    coolant_momentum_model = "quasi_steady"  # or "low_mach"
    transient_coolant_outlet_pressure = None
    insert_schedule_breakpoints = True
    progress_print = True
    progress_interval_steps = None
    progress_interval_time_s = None
```

`fixed_step` is the default production path for shell-and-tube counter-flow
transients because its cost is bounded: roughly one quasi-steady fluid pass per
time step. It uses a linearly-implicit local wall update to damp the stiff
helium-film coupling without BDF-style Jacobian probing. `BDF` remains useful
for validation comparisons, but it can be much slower around schedule
discontinuities.

During `transient_coolant` runs, the terminal prints throttled progress lines
with step count, simulated time, material minimum/maximum temperature, helium
outlet temperature, a pressure diagnostic, and hot-gas outlet temperature. In
shell-and-tube `quasi_steady` momentum mode, that diagnostic is `dpHe [bar]`,
the shell-side hydraulic pressure-drop estimate. In pressure-resolved
`low_mach` mode, it is the reconstructed helium outlet pressure. With
`progress_interval_steps = None`, the solver prints about 20 progress lines per
run. Set `progress_interval_steps = 1` for every internal time step, set
`progress_interval_time_s` for a simulated-time interval, or set
`progress_print = False` to silence it.

The default coolant momentum closure remains `quasi_steady`. For pressure-driven
low-Mach testing, set:

```python
class transientProp:
    coolant_momentum_model = "low_mach"
    transient_coolant_outlet_pressure = 70e5
```

This uses scheduled helium inlet pressure/temperature and a fixed downstream
pressure to compute transient face mass flows. Dense measured schedules, such as
Excel files sampled every few milliseconds, should normally use:

```python
class transientProp:
    insert_schedule_breakpoints = False
```

so the solver interpolates the schedule at its own `max_step` instead of forcing
every measurement row into the integration grid.

## 3. Prepare A Transient Schedule

The simplest format is one CSV or one Excel sheet with these columns:

```text
time_s
helium_m_dot_kg_s
helium_T_in_K
helium_p_in_Pa
diesel_m_dot_kg_s
lox_m_dot_kg_s
lox_T_in_K
ignition
```

Example:

```csv
time_s,helium_m_dot_kg_s,helium_T_in_K,helium_p_in_Pa,diesel_m_dot_kg_s,lox_m_dot_kg_s,lox_T_in_K,ignition
0.0,0.001,293.15,8200000,0.0,0.0,90.0,0
2.5,0.150,303.15,8200000,0.0,0.0,90.0,0
5.0,0.150,303.15,8200000,0.001,0.003,90.0,1
6.0,0.150,303.15,8200000,0.075,0.225,90.0,1
100.0,0.150,303.15,8200000,0.075,0.225,90.0,1
```

The included template is `inputs/example_transient_schedule.csv`.

Decimal dots and decimal commas are both accepted in schedule values. If you
export a CSV with decimal commas, use semicolon or tab separators, for example
`0,0025;1,88E-07;91,25`; otherwise the comma is ambiguous with the CSV field
separator. Native Excel `.xlsx` files can use decimal-comma cells directly.

The loader computes:

```text
hot gas mass flow = diesel_m_dot_kg_s + lox_m_dot_kg_s
O/F = lox_m_dot_kg_s / diesel_m_dot_kg_s
```

The schedule is linearly interpolated between rows and held constant before the
first row and after the last row.

Two-sheet Excel files are also supported:

- Sheet `helium`: `time_s`, `m_dot_kg_s`, `T_in_K`, `p_in_Pa`
- Sheet `propellants`: `time_s`, `m_dot_lox_kg_s`, `m_dot_diesel_kg_s`

## 4. Pre-Combustion / Chilldown State

A helium-only or cold-GOX first phase is supported today:

- Set helium flow normally.
- For helium-only, set diesel and LOX mass flows to `0.0`.
- For GOX chilldown, set `lox_m_dot_kg_s > 0`, `diesel_m_dot_kg_s = 0.0`,
  `lox_T_in_K` to the oxygen inlet temperature, and `ignition = 0`.
- Set `ignition = 0`.
- At combustion start, set `ignition = 1` and add pilot/full propellant flow.

Before ignition, the hot-side combustion manifold is bypassed. If LOX/GOX flow
is scheduled, the transient solver uses CoolProp Oxygen properties and marches
that oxygen stream through the hot-side passage. Once `ignition` becomes `1`,
the finite-rate combustion manifold takes over.

## 5. Chemistry Options

Steady chemistry is controlled by:

```python
class numericalProp:
    chemistry_model = "finite_rate"
```

Accepted steady values:

- `finite_rate`: default, FPV manifold plus axial progress-variable march.
- `equilibrium`: HP re-equilibration after each heat-removal step.
- `frozen`: composition fixed after initial combustion.

Transient chemistry is controlled by:

```python
class transientProp:
    chemistry_transient = "finite_rate"
```

Accepted transient values:

- `finite_rate`: default, FPV manifold plus progress-variable transport.
- `equilibrium`: table-driven HP re-equilibration during gas cooling.
- `frozen`: no recombination during cooling; mainly a validation comparison.

The production finite-rate implementation is:

```text
1Dmodel/physics/combustion_chemistry/fpv_manifold.py
```

It is used by both steady and transient solvers when the corresponding chemistry
setting is `finite_rate`.
The standalone flamelet/FPV toolkit in:

```text
research/flamelet_kit/
```

is for methodology development and validation, not for ordinary simulation runs.

## 6. Run Simulations

Activate the environment:

```powershell
.\.venv\Scripts\Activate.ps1
```

Run steady:

```powershell
hps-steady
```

Run transient:

```powershell
hps-transient
```

Equivalent module commands:

```powershell
python -m hps_combustor.main_steady
python -m hps_combustor.main_transient
```

## 7. Output Files

Each run creates a timestamped folder and, by default, a zip archive.

Transient output structure:

```text
zip_folders/my_case_transient_08-07-2026_14h30/
  metadata.json
  summary.json
  HX_performance_summary.txt
  inputs/
    input_preset.json
    input_data_snapshot.py
  data/
    transient_timeseries.npz
    transient_scalars.csv
```

The `.npz` file is the complete low-cost numeric dataset for replotting.
It includes the transient wall temperature evolution when available:
`field_T_wg` hot wall face, `field_Tbar` mean wall temperature, and
`field_T_wc` cold wall face, each over time and axial position.
The `.csv` file contains scalar time histories only.
`HX_performance_summary.txt` is a human-readable engineering report with
multi-time-point duty, effectiveness/NTU estimates, LMTD, pressure drops,
temperature extrema, and available combustion/hydraulic metrics.

## 8. Rebuild A Dashboard From Saved Data

Generate a dashboard from a run folder:

```powershell
hps-dashboard zip_folders\my_case_transient_08-07-2026_14h30
```

Generate from the numeric data file:

```powershell
hps-dashboard zip_folders\my_case_transient_08-07-2026_14h30\data\transient_timeseries.npz
```

Generate directly from the zip archive:

```powershell
hps-dashboard zip_folders\my_case_transient_08-07-2026_14h30.zip
```

Choose an explicit output file:

```powershell
hps-dashboard zip_folders\my_case_transient_08-07-2026_14h30.zip -o dashboards\my_case.html
```

The dashboard is a self-contained HTML file. It does not need a server.

## 9. Coupled Bang-Bang Validation

The coupled bang-bang validation module is separate from the production
`hps-transient` path. It is a sandbox for checking how a 0D pressurant/feed/tank
system drives the detailed 1D HX. It can run either the EchTherm-style
shell-and-tube configuration or the helical coil configuration with the same
pressurant system.

The default dummy settings are saved in:

```text
inputs/coupled_bangbang_hx_dummy_validation.json
```

Run the saved dummy case:

```powershell
python -m hps_combustor.validation.coupled_bangbang_hx --settings inputs/coupled_bangbang_hx_dummy_validation.json
```

Select the helical HX variant explicitly:

```powershell
python -m hps_combustor.validation.coupled_bangbang_hx --settings inputs/coupled_bangbang_hx_dummy_validation.json --hx-config shellnHelicalTube
```

Fast system-only check, without running the detailed HX:

```powershell
python -m hps_combustor.validation.coupled_bangbang_hx --settings inputs/coupled_bangbang_hx_dummy_validation.json --system-only
```

Outputs are written to:

```text
docs/validation/coupled_bangbang_hx/
  settings_used.json
  summary.json
  README.md
  system_timeseries.csv
  hx_boundary_schedule.csv
  water_outlet_orifice_sweep.json
  coupled_timeseries.csv
  hx_transient_timeseries.npz
  coupled_dashboard.html
```

Open `coupled_dashboard.html` in a browser to view the dedicated bang-bang
validation dashboard. It includes system pressures/flows/valve state, HX
temperatures, heat duty, coolant face mass flows, pressure-drop diagnostics, and
validation flags.

The current dummy preset uses a 90 K, 265 L, 400 bar helium tank, a 70 bar water
tank target, 30 L/s water discharge, three 1.75 mm helium valve/orifice branches
at 40 Hz, counter-flow shell-and-tube HX, 5 axial nodes, `max_step = 0.005 s`,
finite-rate chemistry, and low-Mach pressure-driven helium momentum. The 0D
valve/tank model provides the available helium flow cap; the HX through-flow is
then solved from pre-HX line pressure, water-tank pressure, and the
Bell-Delaware shell-side resistance. Hot-gas flow tracks helium demand at
O/F = 2 and is capped at 70 g/s total flow. It uses the EchTherm-style
shell-and-tube geometry from `input_data.py`.

For `--hx-config shellnHelicalTube`, the coupled validator applies the helical
geometry fields from the validation config, maps them onto the legacy helical
input fields, uses 316L by default, and keeps the previously calibrated
Salimpour/Mori shell-side correction factor of `0.28`. The helical low-Mach
momentum model is currently a single implicit lumped pipe through-flow state
coupled to distributed coolant mass/energy and distributed wall temperature. It
is therefore appropriate for fast design screening and coupled-system sanity
checks, but it is not yet a fully distributed pressure-wave or face-pressure
solver.

Run the compact shell/helical geometry screen:

```powershell
python -m hps_combustor.validation.coupled_bangbang_hx_geometry_sweep --settings inputs/coupled_bangbang_hx_dummy_validation.json --t-end 1 --hx-max-step 0.003 --hx-nodes 5 --hx-save-points 101
```

The latest helical-focused screen is stored under
`docs/validation/coupled_bangbang_helical_geometry_sweep_1s_pressure_metrics/`.
Over the 1 s screen, the three 12 m helical candidates all solved near
0.143 kg/s helium and held the water tank near 70 bar with about 30 L/s water
outflow. Their flowing coolant pressure span was about 9.53 bar mean,
10.95 bar p95, and 11.14 bar max. Read that as mean loss inside the desired
5-10 bar band, with p95/max slightly above 10 bar due to bang-bang ripple.

In `coolant_momentum_model = "quasi_steady"`, use
`hx_max_shell_pressure_drop_estimate_bar` or `hx_dp_shell_total_Pa` for helium
pressure loss. The saved `p_c` field is a thermodynamic coolant-state
reconstruction used for transient mass/energy bookkeeping; its axial span is not
the hydraulic pressure drop in this mode. In `low_mach` mode, `p_c` is the
resolved pressure field and should be interpreted as the momentum-model
diagnostic.

The accepted 100 s detailed coupled validation is stored under
`docs/validation/coupled_bangbang_hx_100s_hotgas72p5_pressure_floor/`. It
completed in about 7.8 min and triggered no coupled sanity flags. System-level
targets were met: 29.99 L/s mean water flow, 69.95 bar mean water-tank pressure,
0.142 kg/s mean helium flow, and about 14.19 kg helium used over 100 s. The HX
hydraulic pressure drops stayed low: helium shell-side estimate below 0.83 bar,
and hot-gas tube-side estimate below 0.14 bar. Flowing helium outlet temperature
had a mean near 596 K and a peak near 693 K; the final low outlet temperature
occurs at zero outlet flow and should be interpreted as stagnant closed-valve
coolant, not delivered helium. Wall-depth peak temperatures were about 894 K at
the hot face, 880 K through the mean wall, and 880 K at the cold face, all near
`x = 0.01175 m`.

The current low-Mach pressure-driven short validation is stored under
`docs/validation/coupled_bangbang_hx_lowmach_5s_dt5ms_n5_hot70/`. It completed
5 s in about 24 s. The 0D system hit the targets: 69.98 bar mean water-tank
pressure, 78.33 bar mean pre-HX line pressure, 29.995 L/s water flow, and
0.145 kg/s mean helium flow. The HX pressure field stayed between about 69.8
and 84.0 bar. The flowing helium outlet temperature peaked near 705 K, slightly
above the 700 K target, and the wall peak stayed near 764 K. The detailed
shell-side pressure-drop estimate peaked near 3.13 bar, below the 5-10 bar
engineering target; treat a long low-Mach run as a thermal/system sanity run
until the shell-side hydraulic resistance is recalibrated or line losses are
coupled explicitly.

Production transient coolant runs use strict CoolProp internal-energy bounds:
if the coolant state leaves the configured property-validity range, the run
raises instead of silently clipping the result.

## 10. Practical Checklist

Before a serious run:

1. Set `combustorProp.HX_config`.
2. Set `combustorProp.flow_config = "counter"` or `"co"`.
3. Set `runProp.run_name`.
4. Set `runProp.schedule_file`.
5. Set `transientProp.t_end`, `n_save`, and `max_step`.
6. Keep `numericalProp.chemistry_model = "finite_rate"` and
   `transientProp.chemistry_transient = "finite_rate"` unless deliberately
   running an equilibrium/frozen comparison.
7. Confirm propellant flows are in kg/s, pressures in Pa, temperatures in K.
8. Use `hps-transient`.
9. Archive the generated zip with the design case.

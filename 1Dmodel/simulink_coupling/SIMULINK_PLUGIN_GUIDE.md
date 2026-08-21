# Simulink Plug-In Guide

Step-by-step: get the shell-and-tube transient HX running as a block inside
Simulink. Read `README.md` first if you haven't — it explains the I/O
contract referenced throughout this guide.

## 0. What you're plugging in

An FMI 2.0 **Co-Simulation** FMU. Simulink calls it once per communication
step with your commanded boundary conditions (coolant pressures/flow,
hot-gas flow, ignition state) and reads back coolant/wall/gas outputs. The
FMU internally sub-steps the wall/coolant physics as needed — Simulink's own
communication step does not have to equal the HX's internal numerical step.

**How the FMU actually runs**: this is a *source* FMU, not a statically
linked one. `pythonfmu` packages a small generic binary per platform
(`binaries/win64/*.dll`, `binaries/linux64/*.so` — both ship in every build,
Simulink picks the one matching its host OS) that embeds/calls a Python
interpreter and executes `resources/fmu_wrapper.py` inside it.
`build_fmu.py` bundles `shelltube_stepper.py` and a private vendored copy of
`hps_combustor` (`resources/simulink_coupling/_vendor/`) into the FMU
alongside it — **the target machine does NOT need this repo or a separate
`hps_combustor` install.** It does still need a Python environment with
`numpy`/`scipy`/`CoolProp`/`cantera` installed — that part can't be avoided
(there's no practical way to statically bundle Cantera into an FMU zip), but
it's the same requirement any Cantera/CoolProp-based tool has, not something
specific to this repo. Verified end-to-end in
`tests/test_simulink_coupling_fmu_wrapper.py::test_packaged_fmu_runs_standalone`:
builds the FMU, unzips it fresh, and runs it with this repo's own
`hps_combustor` install hidden.

## 1. Prerequisites

On the machine that will **build** the FMU (this repo's dev environment is
already set up for this):

```powershell
pip install pythonfmu
# or: pip install -e ".[simulink]"   (from the repo root)
```

On the machine that will **run** Simulink and import the FMU: just a Python
environment with `numpy`, `scipy`, `CoolProp`, and `cantera` installed — the
built FMU carries its own copy of the rest, see step 0. This can be the same
machine as the build machine, or a different one, as long as its Python is
compatible (this repo targets Python 3.10+; built and tested here on 3.13).

## 2. Build the FMU

```powershell
python -m hps_combustor.simulink_coupling.build_fmu path\to\output\dir
```

This produces `ShellTubeTransientFmu.fmu` in that directory. Verified in this
repo's environment — the build succeeds and produces a valid
`modelDescription.xml` with the ports listed below.

To point the FMU at your own combustor/coolant/geometry configuration
instead of the placeholder defaults, edit `_build_config()` in
`fmu_wrapper.py` before building — it returns the same config dataclasses
`main_transient.py` uses (`coolantProp`, `hotgasProp`, `shellTubeProp`,
`numericalProp`, `system_requirements`, `transientProp`,
`CorrelationCoefficients`), imported from `hps_combustor.input_data`. This is
the same file you'd edit for any other run of this solver.

## 3. Import into Simulink

In Simulink (R2018b or later, from a version that ships/has the FMU Import
block — check your MATLAB release's FMI toolbox support):

1. Open your model, add an **FMU Import** block (Simulink library:
   *Simulink > Ports & Subsystems*, or search "FMU" in the Library Browser).
2. Point it at `ShellTubeTransientFmu.fmu`.
3. Simulink reads `modelDescription.xml` and exposes the ports automatically
   — you should see 7 inputs and 10 outputs (listed in README.md's I/O
   contract table). No manual port configuration needed beyond wiring them.
4. Set the block's **communication step size**. Start conservative (a few
   milliseconds) and increase once the run is stable — see "Choosing a
   communication step" below.

## 4. Wire the ports

Minimum wiring for a useful loop:

- **Inputs you must drive every step**: `p_coolant_in`, `p_coolant_out`,
  `mdot_coolant`, `T_coolant_in`, `mdot_hot_total`, `ignited`.
- **`T_lox_in`**: only matters while `ignited=false`; wire it to your
  GOX/LOX feed line's temperature if you model pre-ignition chilldown,
  otherwise leave at its default.
- **Outputs worth scoping first**: `T_wall_max` (structural margin),
  `T_coolant_outlet` (what the downstream system sees thermally),
  `dp_coolant_hydraulic_Pa` (the model's own predicted pressure drop —
  see the warning below), `duty_W` (energy balance sanity),
  `mass_residual_kg` / `energy_residual_J` (should stay within roundoff — a
  growing residual means the communication step is too large or a boundary
  input jumped discontinuously between steps).

If your Simulink model represents the combustion chamber and feed system as
separate blocks, `mdot_hot_total` and `ignited` are naturally their outputs;
`p_coolant_in` is naturally the output of whatever upstream pressure source
feeds the HX.

**Do not wire `p_coolant_outlet` back into `p_coolant_out`.** In `low_mach`
mode `p_coolant_outlet` is mathematically forced to equal `p_coolant_out` —
it's a state reconstruction, not a prediction — so feeding it back creates a
degenerate identity loop, not a physical one. If you want the HX's own
predicted downstream pressure driving `p_coolant_out` (e.g. no separate
plenum/tank model on the Simulink side), use **`p_coolant_outlet_predicted`**
instead (`p_coolant_in - dp_coolant_hydraulic_Pa`, a genuine friction-based
prediction, safe to feed back). The more physically correct architecture,
though, is to give your downstream system its own pressure state (a
plenum/tank mass balance) driven by `mdot_coolant_achieved`, and let *that*
model's own state become the next `p_coolant_out` — see README.md's "Do not
feed `p_coolant_outlet` back into `p_coolant_out`" section for the full
explanation.

## 5. Choosing a communication step

Start small and verify before trusting a larger step:

1. Run a short scenario (a few seconds) at a conservative step (1-5 ms).
2. Check `mass_residual_kg`/`energy_residual_J` stay near zero throughout.
3. Increase the step and re-check. If residuals grow or `T_wall_max` looks
   physically implausible (a jump of hundreds of K in one step), the step is
   too large for the current grid resolution (`N_axial` in
   `_build_config()`) and boundary ramp rate — this project has documented
   the same grid/timestep stiffness interaction for the legacy
   schedule-driven solver (see `docs/context/TRANSIENT_STATUS.md`); it is
   not specific to this wrapper.
4. **Never command a step change to full hot-gas flow from a cold start.**
   Ramp `mdot_hot_total` in gradually (the project's own default schedule
   ramps over ~1 s) — a discontinuous jump to full duty at `t=0` on a fine
   grid can blow past CoolProp's valid temperature range in a couple of
   steps and raise `FloatingPointError`, exactly as it would in the legacy
   solver given an equally reckless schedule.

## 6. Deploying to another machine

Just hand over the built `.fmu` file — it's self-contained (see step 0), so
the target machine only needs `numpy`/`scipy`/`CoolProp`/`cantera` installed,
not this repo.

If instead you want to hand someone the **editable source** — e.g. so they
can also run `main_transient.py`, tweak physics, or build their own FMU with
a different `_build_config()` from source rather than using a `.fmu` you
already built — use `package_for_handoff.py`, which zips the whole
`1Dmodel/` package + `pyproject.toml` (nothing else from the repo — no git
history, docs, or accumulated validation output):

```powershell
python 1Dmodel\simulink_coupling\package_for_handoff.py path\to\output.zip
```

Or, if you want to hand over just this `simulink_coupling/` folder itself
(e.g. to someone building their own FMI integration around
`ShellTubeTransientStepper` rather than using the FMU wrapper as-is), it
already works standalone — copy the whole folder as-is, `_vendor/` included.
See README.md's "Standalone deployment" section.

## 7. Troubleshooting

| Symptom | Likely cause |
|---|---|
| `ImportError: attempted relative import with no known parent package` when building | Don't edit `fmu_wrapper.py` to use `from .shelltube_stepper import ...` — `pythonfmu` loads this file as a standalone top-level module both at build time and inside the FMU, so it must use the absolute `from hps_combustor.simulink_coupling.shelltube_stepper import ...` form (already how it's written; re-check if you've modified it). |
| `FloatingPointError: coolant internal energy left the configured CoolProp temperature range` | Boundary trajectory too aggressive for the grid/step combination — see step 5, point 4. This is an intentional strict guard (documented project-wide), not a silent-clip bug. |
| FMU imports but every output stays at its initial value | Check the master algorithm is actually calling `do_step` — some FMI masters skip `setup_experiment`; the wrapper builds the stepper lazily on first `do_step` either way, but confirm the block isn't paused/disabled. |
| Results don't match a full offline run of `main_transient.py` for the "same" scenario | Two known, documented reasons, not bugs: (1) hydraulic resistance is calibrated from `mdot_coolant_reference` (default: nominal design flow) instead of the whole run's look-ahead maximum the schedule-driven solver uses; (2) the initial pressure profile defaults to a nominal dp estimate unless you pass `p_coolant_out_initial`. See README.md's I/O contract for both. |
| Only Helium coolant behaves correctly | Expected — see README.md's "Fixed at construction" section on `coolantProp.coolant`. |

## 8. Verifying without Simulink

You don't need Simulink to sanity-check the stepper or the FMU wrapper:

```powershell
# Core stepper, pure Python:
python -c "
from hps_combustor.input_data import coolantProp, hotgasProp, shellTubeProp, numericalProp, transientProp, system_requirements, CorrelationCoefficients
from hps_combustor.simulink_coupling import ShellTubeTransientStepper, BoundaryInputs
cp, hp, stp, np_, sr, tp = coolantProp(), hotgasProp(), shellTubeProp(), numericalProp(), system_requirements(), transientProp()
tp.coolant_momentum_model = 'low_mach'
s = ShellTubeTransientStepper(cp, hp, stp, np_, sr, tp, corrCoeffs=CorrelationCoefficients(), N_axial=20, flow_config='co')
b = BoundaryInputs(mdot_coolant=cp.mass_flow_c, p_coolant_in=cp.p_in, p_coolant_out=cp.p_in-2e5, T_coolant_in=cp.T_in, mdot_hot_total=0.01, ignited=True)
out = s.step(0.01, b)
print(out.T_coolant_outlet, out.p_coolant_outlet, out.T_wall.max())
"

# FMI slave wrapper, without a real FMU import:
pytest tests/test_simulink_coupling_fmu_wrapper.py -v
```

The commands above exercise this repo's own `hps_combustor` install, not the
standalone/vendored path. To verify the standalone folder or the packaged
FMU specifically:

```powershell
pytest tests/test_simulink_coupling_standalone.py -v          # bare folder, no repo
pytest tests/test_simulink_coupling_fmu_wrapper.py::test_packaged_fmu_runs_standalone -v  # built .fmu, no repo
```

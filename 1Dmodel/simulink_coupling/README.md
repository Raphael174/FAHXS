# simulink_coupling

Decoupled, step-by-step co-simulation layer for the shell-and-tube transient
HX solver. **Nothing outside this folder is modified** to provide this
capability — see the repo root's `CLAUDE.md` / `docs/` for the maintained
solver itself; this folder only *calls into* it.

For "how do I actually wire this into Simulink", see
[`SIMULINK_PLUGIN_GUIDE.md`](SIMULINK_PLUGIN_GUIDE.md). This file is the
technical reference for the I/O contract and design limitations.

## Why this exists

Simulink's own integrator owns time-stepping in a multiphysics loop and calls
each component at fixed "communication steps," exchanging boundary values.
The existing shell-and-tube transient solver (`run_shelltube_transient_core`
in `1Dmodel/transient_core/adapters_shelltube.py`) instead owns its *entire*
run internally against a pre-known schedule — it cannot be handed to Simulink
directly. `ShellTubeTransientStepper` (`shelltube_stepper.py`) re-sequences
the same underlying physics functions into a single-call
`step(dt, boundary) -> outputs` API so an external caller can own the loop.

## Files

| File | Purpose |
|---|---|
| `shelltube_stepper.py` | The core stepper: `BoundaryInputs`, `StepOutputs`, `ShellTubeTransientStepper`. Pure Python + numpy/CoolProp/Cantera, no FMI dependency. Only absolute `hps_combustor.*` imports — never relative — so it works standalone (see below). |
| `fmu_wrapper.py` | `pythonfmu.Fmi2Slave` subclass wrapping the stepper as FMI variables. Requires `pip install pythonfmu` (or `pip install -e ".[simulink]"`) to build/run. |
| `build_fmu.py` | Packages `fmu_wrapper.py` into a `.fmu` file. |
| `_vendor/hps_combustor/` | A private, self-contained copy of exactly the `hps_combustor` files this folder needs — see "Standalone deployment" below. Regenerate with `vendor_dependencies.py`. |
| `vendor_dependencies.py` | Regenerates `_vendor/hps_combustor/` from the current `1Dmodel/` source. Run after changing `shelltube_stepper.py`'s imports. |
| `package_for_handoff.py` | Zips the *whole dev repo's* dependency closure (`1Dmodel/` + `pyproject.toml`) for someone who wants to build their own FMU with custom `_build_config()` settings from source — a different scenario from the standalone folder below. |
| `SIMULINK_PLUGIN_GUIDE.md` | Step-by-step: build the FMU, import into Simulink, wire the ports. |

## Standalone deployment — this folder really does work on its own

Earlier revisions of this module needed the rest of the repo installed as
`hps_combustor` — copying just this folder somewhere else failed on
`import hps_combustor`. That's fixed: **`1Dmodel/simulink_coupling/`, copied
by itself with nothing else from this repo, is a complete, working
deliverable.**

How: `_vendor/hps_combustor/` is a private copy of exactly the
`hps_combustor` files `shelltube_stepper.py` actually imports (traced via
`sys.modules`, not hand-guessed — see `vendor_dependencies.py`), plus the
combustion mechanism files `choose_fuel()` can select among. Both
`shelltube_stepper.py` and `fmu_wrapper.py` try a normal `import
hps_combustor` first (so the real, editable-installed package is always
preferred inside this dev repo — nothing about normal development or testing
here changes) and only fall back to `_vendor/` if that import fails, which is
exactly the situation when someone has *only* this folder.

**Verified, not just claimed**: `tests/test_simulink_coupling_standalone.py`
copies this folder to an isolated temp directory, strips this repo's
editable install from a subprocess's import machinery, and confirms the
stepper and FMI wrapper both run using nothing but the vendored copy. Run it
yourself: `pytest tests/test_simulink_coupling_standalone.py -v`.

**What you still need on the receiving machine** (this is not "the
combustor," it's generic scientific Python — the same thing any Cantera/
CoolProp-based tool needs): `numpy`, `scipy`, `CoolProp`, `cantera`, and
`pythonfmu` if building/running the FMU. None of this repo's git history,
`pyproject.toml`, docs, tests, or the rest of `1Dmodel/` are required.

**Keeping the vendor copy current**: it's a generated snapshot, not
hand-maintained — if you change what `shelltube_stepper.py`/`fmu_wrapper.py`
import, re-run `python 1Dmodel/simulink_coupling/vendor_dependencies.py`
before handing the folder to anyone, and re-run
`test_simulink_coupling_standalone.py` to confirm it still works with only
the vendor copy present.

`package_for_handoff.py` is for a different scenario: handing someone the
*real, editable* source of the whole solver (e.g. so they can also run
`main_transient.py`, edit physics, or point the FMU at a different
combustor config from source) rather than just this pre-built component.

## I/O contract

**Fixed at construction** (pass a new set of config dataclasses and build a
new `ShellTubeTransientStepper`/FMU instance to change these — see
`_build_config()` in `fmu_wrapper.py`):

- All of `shellTubeProp` (geometry, materials, correlation selection).
- `combustorProp.flow_config` (co/counter), effectively via the
  `flow_config` constructor argument.
- `coolantProp.coolant` — **this stepper's compressible coolant path is
  Helium-only**, independent of what you set here. Several private helpers
  it reuses from `adapters_shelltube.py` call CoolProp with a hardcoded
  `"Helium"` fluid string; this is an existing constraint in that code, not
  something this module adds or removes.
- Chemistry mode (`transientProp.chemistry_transient`: `finite_rate` /
  `equilibrium` / `frozen`).
- Momentum model (`transientProp.coolant_momentum_model`: `quasi_steady` /
  `low_mach`) — **this determines whether `p_coolant_out` is load-bearing**,
  see below.
- The hot-gas manifold's baked-in O/F ratio and chamber pressure
  (`hotgasProp.mixing_ratio`, `hotgasProp.p0`). Varying either mid-run would
  require rebuilding the FPV/equilibrium manifold (disk-cached under
  `cache/fpv_manifolds`, ~1 s warm / ~30 s cold) — not implemented here.
- Hydraulic resistance/inertance calibration reference flow
  (`mdot_coolant_reference` constructor arg, default: `coolantProp.mass_flow_c`).
  The existing schedule-driven solver gets to calibrate this against the
  *whole run's* maximum commanded flow (it knows the future schedule up
  front); a live Simulink-driven stepper cannot, so pass your actual design
  peak flow here if it differs materially from the nominal value, or accept
  that resistance/inertance stay calibrated at the design point rather than
  whatever peak Simulink ends up commanding.
- Initial outlet pressure (`p_coolant_out_initial` constructor arg, default:
  a nominal dp estimate). If your co-simulation starts from a known steady
  operating point, pass its actual outlet pressure here for an exact match
  to what a from-scratch steady solve would give; otherwise the stepper
  seeds a physically reasonable but approximate cold-start guess that
  self-corrects within the first few steps anyway.

**Live per-communication-step inputs** (`BoundaryInputs`, passed to every
`step()` call):

| Field | Units | Role |
|---|---|---|
| `p_coolant_in` | Pa | Coolant inlet pressure |
| `p_coolant_out` | Pa | Coolant downstream/outlet pressure. **The actual momentum-driving boundary in `low_mach` mode** — flow is solved from the two-ended Δp, not commanded. A soft reference for the hydraulic dp estimate only in `quasi_steady` mode. |
| `mdot_coolant` | kg/s | Coolant inlet mass flow. The direct flow driver in `quasi_steady` mode. In `low_mach` mode it acts only as a valve-position/cap ceiling on top of the pressure-solved flow — the true flow is `StepOutputs.face_mdot`, not this input. |
| `T_coolant_in` | K | Coolant inlet temperature |
| `mdot_hot_total` | kg/s | Total hot-gas mass flow (post-ignition: combined fuel+oxidizer at the baked-in O/F; pre-ignition: LOX/GOX-only) |
| `ignited` | bool | Switches between pre-ignition GOX/CoolProp sensible-cooling physics and the combustion manifold |
| `T_lox_in` | K | Pre-ignition LOX/GOX inlet temperature (only meaningful while `ignited=False`; cheap to vary live — pure CoolProp, no manifold) |

**Outputs** (`StepOutputs`, returned from every `step()` call):

| Field | Units | Meaning |
|---|---|---|
| `T_wall` | K (array) | Full axial wall temperature field |
| `T_coolant` | K (array) | Full axial coolant temperature field |
| `T_coolant_outlet` | K | Coolant temperature at the outlet cell |
| `p_coolant_outlet` | Pa | Reconstructed thermodynamic coolant pressure at the outlet cell. **In `low_mach` mode this is mathematically forced to equal `boundary.p_coolant_out` every step — it is not an independent prediction.** See the warning below. |
| `dp_coolant_hydraulic_Pa` | Pa | The model's own Bell-Delaware friction-based pressure-drop prediction for this step's actual mass flux/properties. Independent of `coolant_momentum_model`; this is the genuinely computed quantity, not a state echo. |
| `T_gas_outlet` | K | Hot-gas outlet temperature |
| `duty_W` | W | Heat added to the wall this step, `/dt` |
| `face_mdot` | kg/s (array, n_cells+1) | Actual solved coolant mass flow at every cell face — the ground truth flow in `low_mach` mode |
| `energy_residual_J`, `mass_residual_kg` | J, kg | Per-step conservation residuals — health diagnostics, should stay near roundoff |

## Do not feed `p_coolant_outlet` back into `p_coolant_out`

This is the single most important usage note in this file, added after exactly
this mistake caused a real Simulink coupling to fail.

`low_mach` mode builds a **linear pressure profile between `p_coolant_in` and
`p_coolant_out`** every step, solves the momentum equation for **mass flow**
from that profile, then reconstructs the coolant's thermodynamic state so it
matches that same linear profile — including at the outlet cell, where the
profile's endpoint *is*, by construction, `p_coolant_out`. So
`p_coolant_outlet` is not a prediction in this mode; it is guaranteed (to
CoolProp round-trip precision) to equal whatever you passed as
`p_coolant_out`. Wiring it straight back as the next step's `p_coolant_out`
creates a degenerate identity loop, not a physical feedback loop — and in
practice, closing that loop through an FMI master's algebraic-loop handling
can produce exactly the kind of two-value oscillation you'd expect from an
ill-posed feedback path.

**The genuinely predicted quantities** in `low_mach` mode are `face_mdot` /
`mdot_coolant_achieved` (the solved flow) and `dp_coolant_hydraulic_Pa` (the
model's own friction-based pressure-drop estimate, from the same
Bell-Delaware correlation the steady solver uses for `dp_c_bar`). Two correct
ways to close a loop with an external system model:

1. **Let your downstream model own the pressure state.** Give it
   `mdot_coolant_achieved` as an inflow to its own mass/energy balance (a
   receiving plenum, tank, etc.); its *own* evolving pressure becomes the
   next `p_coolant_out`. This is the physically correct architecture for a
   component that doesn't know what's downstream of it.
2. **Use `p_coolant_outlet_predicted`** (FMU output, `p_coolant_in -
   dp_coolant_hydraulic_Pa`) if you want a self-consistent "our own predicted
   dP" signal to feed back directly — unlike `p_coolant_outlet`, this is not
   forced to equal any input, so it's safe to close a loop around.

## Known limitation: reused private helpers

`shelltube_stepper.py` imports several underscore-prefixed helpers directly
from `adapters_shelltube.py` (hydraulic resistance/inertance calibration,
pressure-profile reconstruction, inventory/backpressure limiting for
closed-valve residual discharge). These encode calibrated physics the
project has separately validated (bang-bang coolant audits) and must be
reused rather than re-derived. They are not part of that module's public
`__all__`, so a future refactor there could rename or change them without
this folder's own tests catching it in review — re-run
`tests/test_simulink_coupling_stepper.py` after any change to
`adapters_shelltube.py`'s mass/energy path; it will fail loudly if behavior
drifts.

## Verification

- `pytest tests/test_simulink_coupling_stepper.py -v` — proves the stepper
  reproduces the existing full-schedule solver's trajectory bit-for-bit
  (within floating-point tolerance) when fed the same time-varying boundary
  values one step at a time, for both `quasi_steady` and `low_mach` momentum
  models.
- `pytest tests/test_simulink_coupling_fmu_wrapper.py -v` — confirms the FMI
  slave wrapper is wired correctly (requires `pythonfmu`).
- `pytest tests/test_transient_core_coolant_fv.py tests/test_main_transient_dispatch.py -v`
  — confirms zero regression to the existing solver (nothing here should
  ever change their result, since no existing file was edited).

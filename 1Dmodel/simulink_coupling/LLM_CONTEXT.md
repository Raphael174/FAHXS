# simulink_coupling/ LLM Context

## Scope

Parked, deprioritized-but-working feature: a decoupled, step-by-step
co-simulation layer that re-sequences the existing, validated shell-and-tube
transient solver (`run_shelltube_transient_core` in
`1Dmodel/transient_core/adapters_shelltube.py`) into a single-call
`step(dt, boundary) -> outputs` API so an external caller (Simulink, an FMU
master algorithm, a MATLAB System object) can own time-stepping instead of
this codebase owning a fixed schedule. It modifies nothing outside this
folder — it only calls into the existing solver.

## Contents

| File | Role |
|---|---|
| `shelltube_stepper.py` | Core stepper: `BoundaryInputs`, `StepOutputs`, `ShellTubeTransientStepper`. Pure Python + numpy/CoolProp/Cantera, no FMI dependency. Constructs `shellntube_transient_solver` exactly as `main_transient.py` does (geometry/materials/chemistry setup reused, not reimplemented) then re-sequences the per-step physics primitives from `transient_core/adapters_shelltube.py`. Uses only absolute `hps_combustor.*` imports (never relative) so it works standalone. |
| `fmu_wrapper.py` | `pythonfmu.Fmi2Slave` subclass wrapping the stepper as FMI 2.0 Co-Simulation variables. Requires `pip install pythonfmu` at build time. Wall temperature field is exposed as three scalars (max/mean/min over the axial grid), not the full array. Construction-time parameters (geometry, chemistry mode, momentum model, O/F, chamber pressure) are fixed at FMU build time, not live FMI variables. |
| `build_fmu.py` | Packages `fmu_wrapper.py` (+ its dependencies, including `shelltube_stepper.py` and `_vendor/`) into a `.fmu` file. |
| `package_for_handoff.py` | Zips the *whole dev repo's* dependency closure (`1Dmodel/` + `pyproject.toml`) for someone who wants to build their own FMU from source with custom `_build_config()` settings — different scenario from the standalone folder-copy path below. |
| `vendor_dependencies.py` | Regenerates `_vendor/hps_combustor/` from the current `1Dmodel/` source. File list was derived empirically (traced via `sys.modules` after constructing a stepper and calling `step()`), not hand-guessed. Re-run after changing `shelltube_stepper.py`'s or `fmu_wrapper.py`'s imports. |
| `_vendor/hps_combustor/` | Generated, vendored copy — **not source**. See "Vendored copy has diverged" below. Do not write per-file `LLM_CONTEXT.md` documentation inside its subfolders; it is a build artifact regenerated from real `1Dmodel/` source, not a place to hand-maintain. |
| `README.md` | Technical reference: I/O contract, standalone-deployment mechanism, design limitations. Read first for "how does this actually work". |
| `SIMULINK_PLUGIN_GUIDE.md` | Step-by-step: build the FMU, import into Simulink, wire the ports. |
| `ShellTubeTransientFmu.fmu` | A built FMU artifact checked into the tree (binary; treat as generated output, not source). |

## Known constraint

The compressible coolant path this stepper drives is Helium-only — several of
the reused private helpers in `transient_core/adapters_shelltube.py` call
CoolProp with a hardcoded `"Helium"` fluid string (a pre-existing gap in the
underlying solver, not introduced here — see also
`docs/solver_design/FV_CORE_REWORK_PLAN.md`'s note on the same hardcoding).
Not generalized by this folder.

## Standalone deployment

`1Dmodel/simulink_coupling/`, copied by itself with nothing else from this
repo, is meant to be a complete, working deliverable. Both
`shelltube_stepper.py` and `fmu_wrapper.py` try `import hps_combustor` first
(so the real editable-installed package is preferred in this dev repo) and
fall back to `_vendor/hps_combustor/` only if that import fails — the
situation when someone has only this folder. Verified by
`tests/test_simulink_coupling_standalone.py` (copies the folder to an
isolated temp dir, hides this repo's editable install from a subprocess, and
confirms the stepper + FMI wrapper run on the vendored copy alone) and
`tests/test_simulink_coupling_fmu_wrapper.py::test_packaged_fmu_runs_standalone`
(same technique against the actual built/unzipped `.fmu`).

## Vendored copy has diverged (as of 2026-08-19) — deliberately

`_vendor/hps_combustor/transient_core/` is a byte-identical vendored copy of
the **pre-fix** `transient_core/` modules from when it was last regenerated.
Since then, a real bug fix landed in the source `1Dmodel/transient_core/
adapters_shelltube.py`: an explicit-Euler CFL instability in
`conservative_mass_energy_step` (unconditionally unstable once a macro step
exceeds roughly one cell's residence time) was fixed by subdividing each
macro step into CFL-safe substeps for mass/energy advection — see
`docs/solver_design/FV_CORE_REWORK_PLAN.md`'s 2026-08-18 entry for the full
root-cause writeup, and `tests/test_transient_core_shelltube.py` for the
regression gate. That fix was **deliberately NOT propagated** into
`_vendor/`, because Simulink work is explicitly parked/deprioritized right
now and re-running `vendor_dependencies.py` was judged not worth doing for a
parked feature.

Practical implication: anyone building a new FMU from the current
`_vendor/` copy (or running the standalone-deployment tests, which exercise
`_vendor/` as-is) is getting the pre-fix, CFL-unstable coolant advection —
same silent under-resolution / possible `FloatingPointError` risk the fix
addressed in the real solver. This is a conscious, tracked debt, not a
surprise: flag it again before resuming Simulink work, and regenerate
`_vendor/` via `vendor_dependencies.py` at that point rather than assuming it
is current.

## How to run

```powershell
python -m pip install pythonfmu   # or: pip install -e ".[simulink]"
python 1Dmodel/simulink_coupling/vendor_dependencies.py   # regenerate _vendor/ after import changes
python 1Dmodel/simulink_coupling/build_fmu.py
pytest tests/test_simulink_coupling_stepper.py tests/test_simulink_coupling_standalone.py tests/test_simulink_coupling_fmu_wrapper.py -v
```

These three test files require `pythonfmu` and are conventionally excluded
from a normal full-suite run via `--ignore` (see `tests/LLM_CONTEXT.md`).

## TODOs

None found via grep of `TODO`/`FIXME` in this folder's own files. The one
open, evidenced item is the `_vendor/` CFL-fix divergence above — tracked in
prose in `docs/solver_design/FV_CORE_REWORK_PLAN.md`, not as an inline code
marker.

## Change history

Git history is one initial commit plus uncommitted work; not a real
changelog. Evidenced, dated events from code comments and
`docs/solver_design/FV_CORE_REWORK_PLAN.md`:

- **2026-07-27** (per project memory / repo setup context): module created,
  vendored to be genuinely standalone after an earlier revision was found not
  to be.
- **2026-08-18**: the shell-and-tube transient CFL instability was found and
  fixed in the real `transient_core/adapters_shelltube.py`; the fix was
  explicitly not propagated to `_vendor/`, creating the divergence described
  above.

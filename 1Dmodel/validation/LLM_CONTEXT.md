# validation/ LLM Context

## Scope

Standalone example/validation scripts, each runnable directly
(`python -m hps_combustor.validation.<name>`) rather than through pytest.
They exist to sanity-check specific solver paths (liquid/boiling coolant,
transient bang-bang coolant behavior, momentum closures, nozzle sizing) with
human-readable printouts and/or JSON/HTML artifacts, not as CI gates. Several
double as the canonical "how do I actually run X" recipe for a feature that
has no simpler example elsewhere in the repo.

## Contents
| File | Role | Character |
|---|---|---|
| `water_helical_example.py` | Confirmed-working Water/`equilibrium_liquid` recipe for the helical solver (`run_coflow()` fast; `run_counterflow_physical()` slower, shooting-method anchored). Explicitly the recommended way to try Water instead of editing shared defaults. | Confirmed working (own docstring says so). |
| `liquid_boiling_straight_pipe.py` | Straight-pipe liquid/boiling proof-of-concept: `p,h` state march, EOS phase/quality/void, Gungor-Winterton boiling HTC. Deliberately avoids helical/shell-and-tube geometry. | Foundational/working — feeds `liquid_validation_matrix.py`. |
| `liquid_hx_imposed_duty.py` | HX-style imposed-duty check for the liquid heated-channel adapter; mimics segment edges/geometry/`dQ` a real solver would hand it, without being one. | Working, part of validation matrix. |
| `liquid_solver_postprocess_audit.py` | Verifies the opt-in liquid postprocess hooks on `main_solver`/`shellntube_solver` consume solver-output-shaped duty fields correctly, without changing the existing steady helium march. | Working, part of validation matrix. |
| `liquid_ttse_backend_validation.py` | Accuracy/speed check of CoolProp TTSE/BICUBIC tabulated backends vs exact HEOS, dense-sampled right at the saturation boundary where CHF margin/quality are evaluated. | Working, standalone report generator. |
| `liquid_validation_matrix.py` | Aggregates the three liquid gates above (straight-pipe, imposed-duty, postprocess-hook) into one JSON readiness report. Not a replacement for the individual scripts. | Working, orchestrator. |
| `bangbang_momentum_audit.py` | Audits bang-bang helium schedule timescales against coolant momentum relevance (10-90% rise time vs line/coil transit time) for both helical and shell-and-tube. | Working, diagnostic-only — reads a tab-separated time/mdot file. |
| `helical_bangbang_coolant_audit.py` | Short audit of the production helical `transient_coolant` path under a bang-bang schedule: inlet open/close behavior, residual discharge, inventory decrease. Notes helical momentum is still quasi-steady and may need a transient-momentum upgrade for aggressive schedules. | Working, targeted audit. |
| `shelltube_bangbang_coolant_audit.py` | Same audit, shell-and-tube `transient_coolant` path. Explicitly not a physics validation against test data — checks intended numerical behavior only. | Working, targeted audit. |
| `low_mach_momentum_dummy_run.py` | Short (20 s) run of the pressure-driven low-Mach shell-side coolant momentum closure, reading `inputs/HX_dummy_inlet.xlsx` for helium inlet T/p. | Working — the input file lives at repo-root `inputs/HX_dummy_inlet.xlsx` (not `1Dmodel/inputs/`); the script's relative default path resolves correctly when run from the repo root, the project's standard invocation convention. See `inputs/LLM_CONTEXT.md`. |
| `coupled_bangbang_hx.py` | Coupled 0D tank/feed/valve pressurant model driving the detailed 1D shell-and-tube transient HX; own docstring calls it "intentionally a separate validation/coupling sandbox" ahead of a future Simulink/ESPSS integration. Emits an HTML report (embedded plotting JS in-file). | Working, larger sandbox script. |
| `coupled_bangbang_hx_geometry_sweep.py` | Small geometry sweep reusing `coupled_bangbang_hx.py`'s tank/feed/valve model and HX runner, to screen candidates before committing to 100 s runs. | Working, sweep wrapper. |
| `pressurant_bangbang_sizing.py` | System-level (not detailed-ESPSS) sizing of a parallel valve/orifice helium bang-bang pressurization feed from a 400 bar/100 K tank into a water-tank ullage. | Working, standalone sizing script. |
| `transient_core_short_runs.py` | Short validation matrix for the shell-and-tube `transient_core`/`transient_coolant` path — this is the "documented validation matrix" referenced by `transient_core`'s CFL-stability fix (see `transient_core/LLM_CONTEXT.md`). Was crashing on every case pre-fix (2026-08-18); now passes with energy residuals ~1e-8 J instead of the pre-fix ~30-93 J that was actually unresolved CFL error, not real non-closure. | Working (post-2026-08-18 fix); pre-fix archived numbers at `docs/validation/transient_core_short_run_results_PRE_CFL_FIX_2026-07-10.json` (gitignored, not always visible via search). |
| `nozzle_c2h4_o2_bartz_example.py` | First-pass C2H4/O2 conical-nozzle Bartz sizing example exercising `core/geometry/nozzle_contour.py`/`core/hotgas/nozzle_gas.py`. Own docstring flags it standalone, not wired into the coupled FV core, uses an ASSUMED uniform hot-wall temperature (no solved wall/coolant coupling), and documents a found mass-flow/throat-diameter inconsistency in the stated design point (resolved by treating mdot+p0 as authoritative and deriving throat diameter). Chemistry uses Cantera's `gri30.yaml`, tuned for natural gas, not C2H4/O2. | Exploratory — explicitly a first pass with stated caveats. |

## Pitfalls

- These scripts are NOT covered by `pytest`; do not assume `pytest` failures
  or passes say anything about them. Run each with
  `python -m hps_combustor.validation.<module>` from the repo root after an
  editable install.
- Several write JSON/HTML artifacts under `docs/validation/` (gitignored —
  see `/CLAUDE.md`'s gitignore note) — a clean git status does not mean these
  scripts haven't been run, and missing output files does not mean they
  don't work.
- `low_mach_momentum_dummy_run.py` depends on an Excel input file at
  repo-root `inputs/HX_dummy_inlet.xlsx` (confirmed present) — the path is
  relative to the repo root, not to this `validation/` folder, so run it
  from the repo root.

## TODOs

None found as explicit TODO/FIXME markers in this folder.

## Change history

No usable git history for most of this folder (single initial commit; the
folder shows largely uncommitted/untracked in `git status`, except
`nozzle_c2h4_o2_bartz_example.py` which is explicitly untracked per the
snapshot at session start). Dated notes found in-file: `nozzle_c2h4_o2_bartz_example.py`
records a design-point discrepancy finding dated 2026-07-31 in its own
docstring. `transient_core_short_runs.py`'s role in the 2026-08-18 CFL fix is
documented in `docs/solver_design/FV_CORE_REWORK_PLAN.md`, not in this file
itself.

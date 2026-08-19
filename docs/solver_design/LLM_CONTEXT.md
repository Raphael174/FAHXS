# docs/solver_design/ LLM Context

## Scope

Solver architecture design documents, ordered roughly by when they were
superseded. This folder holds both the **current plan of record** for the
active FV-core rework and several historical/completed design docs that fed
into it. Reading the wrong one as "current" is the main risk here — use the
status table below, not filename recency, to decide what's live.

## Contents

| File | Role |
|---|---|
| `FV_CORE_REWORK_PLAN.md` | **CURRENT PLAN OF RECORD.** The active, still-being-executed design for the fluid-agnostic quasi-1D FV core rework (`1Dmodel/core/`) plus two brand-new rocket-nozzle regen configs. Large, updated in place with dated notes. Read this first for any solver-architecture question. |
| `FLUID_AGNOSTIC_CLOSURES_AND_SUPERCRITICAL_PLAN.md` | **COMPLETE for Phases 0–1** (closure registry, regime detection, supercritical closure family — done 2026-07-19). Its own Phase 2 (generic FV core) is explicitly deferred and is now **superseded by `FV_CORE_REWORK_PLAN.md`**, which says so in its own status banner ("Supersedes 'Phase 2' of `FLUID_AGNOSTIC_CLOSURES_AND_SUPERCRITICAL_PLAN.md`"). Still the authoritative source for how the closure registry (`physics/liquid_flow/registry.py`), regime detection (`regime.py`), and supercritical closures (`supercritical.py`) work and why (Cheng2020 vs. McCarthy-Wolf/Taylor ranking decision, PDF minus-sign-drop verification discipline, LN2 validation case). |
| `LIQUID_COOLANT_INTEGRATION_PLAN.md` | **Phases 0–2 (helical) DONE; Phase 3 (shell-and-tube shell-side coupling) is the next open phase; Phases 4–6 (module renaming, TOML presets, transient liquid) not started.** Authoritative source for exactly what's wired vs. postprocess-only for liquid/boiling coolant — matches root `CLAUDE.md`'s liquid-coolant status line. Contains the detailed 2026-07-13 incident writeup (the `HX_config` axial-bookkeeping bug) and the "Post-Phase-2 Hardening Pass" for `solve_counterflow_liquid_reference()`. |
| `water_coolant_conversion_plan.md` | Original from-scratch Helium→Water conversion plan (roadmap with week-level effort estimates, phase 1–5). **Superseded as a status source** by `STEADY_LIQUID_BOILING_REQUIREMENTS.md`, which explicitly says it "updates `water_coolant_conversion_plan.md` against what is actually implemented." Still useful for the original physical/numerical reasoning (HEM vs. drift-flux vs. two-fluid, why `(p,h)` state, correlation stack table) but do not read it for current status. |
| `STEADY_LIQUID_BOILING_REQUIREMENTS.md` | **Status-tracked, mostly DONE as of 2026-07-16.** Item-by-item ("required, done" / "required, partially done") checklist against `water_coolant_conversion_plan.md`, covering state variables, momentum (incl. shell-and-tube shell-side pressure march), correlation stack, sound speed, numerical robustness. Remaining open items: literature CHF/wall-temperature benchmark validation (Bennett et al. 1967 — not attempted, data unavailable) and shell-side cross-flow boiling enhancement (deliberately unmodeled, flagged as a real gap). This is the file to read for "is liquid/boiling physics actually correct" status. |
| `STEADY_LIQUID_BOILING_REQUIREMENTS_V2.md` | A generic, **status-free** "build list" covering the same six requirement categories (state variables, momentum, correlation chain, sound speed, numerical robustness, validation) but framed as a from-scratch checklist for *any* HX design in this architecture, not tied to what's implemented in this repo. Despite the "V2" name it does not supersede the non-V2 file — the non-V2 file is the one with actual dated implementation status; V2 reads as an abstracted/generalized companion with no per-item status markers. Treat V2 as a reusable requirements checklist, and the non-V2 file as the status report against it. |
| `TRANSIENT_CORE_IMPLEMENTATION_PLAN.md` | Design + running implementation-status log for `1Dmodel/transient_core/` (the finite-volume transient coolant kernel: `coolant_fv.py`, `wall_coolant.py`, `compressible_coolant.py`, adapters, integrator). Per its own per-work-package status tables: WP1/WP1b (coolant FV kernel, axial grid) and the geometry/inventory halves of WP3/WP4 (helical and shell-and-tube adapters) are implemented; WP5 (unified dispatch) has shell-and-tube done but helical transient-coolant dispatch not yet wired at the time this doc was last touched; WP2 nonlinear Picard coupling and WP8 (transient hot gas) are not started. **`docs/context/TRANSIENT_STATUS.md` is more current** than this file for exact wiring status — read that first, use this file for the original governing-equation derivations and numbered work-package structure. |
| `TRANSIENT_COOLANT_MASS_ENERGY_PLAN.md` | Design + status log for the conservative `(m, U)` coolant mass/energy model in `transient_core/compressible_coolant.py` / `wall_compressible_coolant.py`. Status line at top: "shell-and-tube and helical production paths wired on 2026-07-10." **Important**: `FV_CORE_REWORK_PLAN.md`'s 2026-08-18 note found and fixed a real CFL-instability bug in exactly this code path (`adapters_shelltube.py`'s use of `conservative_mass_energy_step`) — this plan document predates that fix and does not reflect it; the "Timestep Rule" and benchmark numbers here are from before the fix and are now known to have been running inside the CFL-violated regime (see `FV_CORE_REWORK_PLAN.md`'s "Stages B/C found already built" note and `docs/validation/LLM_CONTEXT.md`). |
| `DESIGN_PLAN_shellntube_transient.md` | Original architecture design for the shell-and-tube config and its transient solver (WP1–WP4, plus chemistry WP-C1–C3). **Explicitly flagged superseded at its own top**: "the implementation has advanced beyond the original plan... For current truth, read `docs/TECHNICAL_REFERENCE.md`, `docs/context/TRANSIENT_STATUS.md`, and `CLAUDE.md`." Still the source for the lumped-wall-ODE-vs-resolved-PDE validation reasoning (the three-way A/B/C comparison) that justifies the current thin-wall transient wall model — that reasoning is not restated anywhere newer. |
| `check_wall_quasi_static_validity.py` | The guardrail script backing `DESIGN_PLAN_shellntube_transient.md` §4.1 and referenced as WP3 acceptance test #2. Compares three wall-model candidates (resolved radial PDE truth, lumped-ODE + quadratic reconstruction, zero-thermal-mass strawman) during a fast He ramp; the lumped model matches truth to <2 K and is what the maintained transient solvers use. Still runnable as a standalone check. |

## Status of each design doc (quick reference)

| Doc | Status |
|---|---|
| `FV_CORE_REWORK_PLAN.md` | **ACTIVE — current plan of record.** Stage A complete; Stages B/C found already built; Stage D in progress (split D1–D4, D1 done); Stages E–G not started. |
| `FLUID_AGNOSTIC_CLOSURES_AND_SUPERCRITICAL_PLAN.md` | Phases 0–1 complete (2026-07-19); Phase 2 superseded by `FV_CORE_REWORK_PLAN.md`. |
| `LIQUID_COOLANT_INTEGRATION_PLAN.md` | Phases 0–2 done; Phase 3 next; 4–6 not started. |
| `water_coolant_conversion_plan.md` | Historical/superseded as a status source; keep for physical reasoning. |
| `STEADY_LIQUID_BOILING_REQUIREMENTS.md` | Status-tracked, mostly done (2026-07-16); two open gaps flagged. |
| `STEADY_LIQUID_BOILING_REQUIREMENTS_V2.md` | Reusable requirements checklist, not a status doc. |
| `TRANSIENT_CORE_IMPLEMENTATION_PLAN.md` | Mostly implemented per its own WP tables; `docs/context/TRANSIENT_STATUS.md` is more current. |
| `TRANSIENT_COOLANT_MASS_ENERGY_PLAN.md` | Wired 2026-07-10; **stale** re: the CFL-instability fix documented later in `FV_CORE_REWORK_PLAN.md`. |
| `DESIGN_PLAN_shellntube_transient.md` | Superseded by implementation; keep for wall-model validation reasoning. |

## The FV-core rework's current state (from `FV_CORE_REWORK_PLAN.md`, as of its last update)

- **Stage A (core/thermo.py)**: complete. All four legacy solvers route
  through `IdealGasBackend`/thermo backends; baseline regressions bit-identical.
- **Stages B (core/mesh.py) and C (core/wall.py)**: discovered already built
  (46 passing tests) while starting Stage D — the plan document's own text
  describing them as "doesn't exist yet" was stale and has been corrected.
- **Stage D (core/residual.py + drivers/transient.py)**: in progress, split
  into slices D1–D4. D1 (closure-registry unification for shell-and-tube's
  tube-side gas closures) is done. D2–D4 not started.
- **A real, fixed bug found along the way**: `transient_core/
  compressible_coolant.py::conservative_mass_energy_step` was unconditionally
  CFL-unstable (explicit forward Euler in coolant conserved variables) at the
  documented validation matrix's step sizes — fixed via CFL-safe substepping
  in `adapters_shelltube.py`. The pre-fix numbers are archived at
  `docs/validation/transient_core_short_run_results_PRE_CFL_FIX_2026-07-10.json`
  for comparison against the current `docs/validation/
  transient_core_short_run_results.json`. See `docs/validation/LLM_CONTEXT.md`.
  **Also fixed, 2026-08-19**: the helical transient core had the same bug
  (confirmed same root cause), fixed with the same mechanism. See
  `tests/test_transient_core_helical.py`.
- Stages E (steady-solver migration), F (nozzle configs), G (cleanup) not
  started.

## TODOs (only evidenced items)

- `FV_CORE_REWORK_PLAN.md`: D2 must decide implicit-vs-subcycled coolant
  advection; D3/D4 must decide whether to reproduce two legacy shell-and-tube
  quirks (`dq_cold`/`T_wc` discarded in favor of `G·(Tbar−Tc)`; `mdot_effective`
  as a mean-of-faces scalar) or fix them; the helical transient core's CFL
  bug is unaddressed; `coolprop_state_from_mass_energy` and related call
  sites hardcode `"Helium"` rather than reading `coolantProp.coolant`.
- `LIQUID_COOLANT_INTEGRATION_PLAN.md`: Phase 3 (shell-and-tube shell-side
  coupled liquid march) is the next open phase; Phases 4–6 (module rename,
  TOML presets, transient liquid) not started.
- `STEADY_LIQUID_BOILING_REQUIREMENTS.md`: literature CHF/wall-temperature
  benchmark validation (Bennett et al. 1967) not attempted; shell-side
  cross-flow boiling enhancement deliberately unmodeled.
- `TRANSIENT_CORE_IMPLEMENTATION_PLAN.md`: WP2 nonlinear Picard coupling of
  `G_i` from `h_c(T, mdot, p)` not done; WP8 transient hot gas not started
  (optional extension, only needed if audits show hot-side quasi-steady is
  not defensible).

## Change history

No meaningful git history (effectively one initial commit plus large
uncommitted work). Dated notes evidenced within the files: `FV_CORE_REWORK_PLAN.md`
carries entries through 2026-08-18 (today's context date is 2026-08-19);
`FLUID_AGNOSTIC_CLOSURES_AND_SUPERCRITICAL_PLAN.md` through 2026-07-19;
`LIQUID_COOLANT_INTEGRATION_PLAN.md` through 2026-07-13;
`STEADY_LIQUID_BOILING_REQUIREMENTS.md` through 2026-07-16;
`TRANSIENT_COOLANT_MASS_ENERGY_PLAN.md` and `TRANSIENT_CORE_IMPLEMENTATION_PLAN.md`
through 2026-07-10/11; `DESIGN_PLAN_shellntube_transient.md` written 2026-07-07.

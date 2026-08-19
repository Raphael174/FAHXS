# docs/context/ LLM Context

## Scope

Canonical per-topic context files, explicitly named in the repo root
`CLAUDE.md`'s "Canonical Context Files" list, meant to be read by future LLM
sessions in place of re-deriving solver state from source. Most are living
documents updated in place with dated notes rather than rewritten.

## Contents

| File | Role |
|---|---|
| `SOLVER_CONTEXT.md` | Maintained solver architecture: entry points, `input_data.py` dataclasses, helical/shell-and-tube/transient solver flow, liquid-coolant wiring status, common pitfalls. Listed as canonical in root `CLAUDE.md`. |
| `PHYSICS_CONTEXT.md` | Correlations, chemistry (FPV/equilibrium/frozen), heat transfer/friction dispatchers, liquid/boiling coolant physics and wiring status, wall conduction, radiation, materials/stress, calibration knobs. Listed as canonical in root `CLAUDE.md`. |
| `TRANSIENT_STATUS.md` | Transient solver implementation status and numerical lessons — largest and most frequently updated file here; a running log of `transient_core` implementation progress, benchmark numbers, and validation-run results dated through 2026-07-11. Listed as canonical in root `CLAUDE.md`. |
| `1Dmodel_CLAUDE.md` | Older, helical-focused project guide (module-by-module key-file table, supported fuels/materials, correlation selectors). Explicitly flagged at its own top as superseded for current truth by `TECHNICAL_REFERENCE.md`, root `CLAUDE.md`, and `TRANSIENT_STATUS.md` — still useful as a quick module map since it lists physics/mechanical/data-module file responsibilities that the newer docs don't restate. |
| `CLAUDE_HANDOFF_2026-07-08.md` | Point-in-time operational handoff (2026-07-08): finite-rate chemistry becoming default, EchTherm shell-and-tube field alignment, what's physically active vs. merely stored on `shellTubeProp`, validation runs performed that session. Historical snapshot — superseded piecemeal by `SOLVER_CONTEXT.md`/`PHYSICS_CONTEXT.md`/`TRANSIENT_STATUS.md` for current truth, but still the best record of *why* certain `shellTubeProp` fields exist and which are wired vs. cosmetic as of that date. |
| `LOW_MACH_COOLANT_MOMENTUM_PLAN.md` | Design + implementation log for pluggable low-Mach coolant momentum (vs. full acoustic momentum, rejected as too expensive). Documents the decision, the physical model (inertance/friction face-momentum ODE), shell-and-tube (`coolant_momentum_model="low_mach"`) and helical (lumped pipe through-flow) implementations, and validation runs through 2026-07-11. Superseded in spirit by the FV-core rework's `momentum.py` design (`solver_design/FV_CORE_REWORK_PLAN.md` §3.6), which generalizes this into a pluggable `MomentumModel` protocol — but this file is still the primary source for *why* low-Mach was chosen over full acoustic and for the legacy solvers' actual current momentum behavior. |
| `COOLANT_COMPARISON_HELIUM_VS_WATER.md` | Focused analysis: why swapping Helium↔Water coolant barely moves absorbed duty `Q` in the helical-coil-in-shell geometry specifically (gas-side resistance dominates, ~97.5–99.6% of total). Explicit scope-of-validity section: does NOT generalize to shell-and-tube (untested, and liquid coupling isn't wired into shell-and-tube's coupled march yet) or to nozzle regen cooling (opposite regime — coolant-side often dominant there). Explains the resistance-network sensitivity law `d(ln UA)/d(ln R_coolant) = -R_coolant/R_total` used to justify the conclusion. |
| `shell_and_tube_architecture_target.png` | Reference geometry screenshot (EchTherm-style shell-and-tube target architecture). Not read as text; see `EchTherm_correlation_docs/` for the extracted correlation content this geometry motivated. |

## Reading order for solver work

1. `SOLVER_CONTEXT.md` and `PHYSICS_CONTEXT.md` first (both are current,
   actively maintained, and cross-reference each other).
2. `TRANSIENT_STATUS.md` if the task touches any transient solver — it is the
   single source of truth for `transient_core` implementation status, and is
   more current than `TRANSIENT_CORE_IMPLEMENTATION_PLAN.md` in
   `docs/solver_design/` (that file is the original plan; this one is the log
   of what actually landed against it).
3. `LOW_MACH_COOLANT_MOMENTUM_PLAN.md` only if the task touches coolant
   momentum/pressure-driven transient behavior specifically.
4. `COOLANT_COMPARISON_HELIUM_VS_WATER.md` only if asked to reason about why
   a fluid swap did or didn't change a result — check its scope-of-validity
   section before applying its conclusion outside the helical-coil-in-shell
   geometry.
5. `1Dmodel_CLAUDE.md` and `CLAUDE_HANDOFF_2026-07-08.md` are historical —
   consult only for module-map detail or point-in-time rationale not restated
   elsewhere; do not treat either as current status.

## TODOs (only evidenced items)

- `TRANSIENT_STATUS.md` itself records, as of its last update, two
  incomplete items: full finite-rate helical counter-flow settle-to-steady
  validation, and broader shell-and-tube transient validation (settle-to-
  steady, chemistry-mode comparison, 100 s benchmark). Root `CLAUDE.md`
  separately notes two background validations (helical low-flow freeze-out
  demo, full shell-and-tube transient validation) were still in flight at
  the prior session's stop — do not assume either completed.

## Change history

No meaningful git history (effectively one initial commit). Dated notes
evidenced within the files themselves: `TRANSIENT_STATUS.md` carries update
entries through 2026-07-11; `LOW_MACH_COOLANT_MOMENTUM_PLAN.md` through
2026-07-11; `CLAUDE_HANDOFF_2026-07-08.md` is a fixed 2026-07-08 snapshot;
`COOLANT_COMPARISON_HELIUM_VS_WATER.md` is undated internally but references
"the user's own solver runs" and the 2026-07-13 liquid-coolant integration
plan, placing it after that date.

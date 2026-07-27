# Claude Project Context

Load this file first when working in Claude Code. It is intentionally short; use the linked files only when the task needs that area.

## Canonical Context Files

- `AGENTS.md`: general agent instructions and read order.
- `CODEBASE_MAP.md`: repository structure, package import model, and entry points.
- `TESTING_CONTEXT.md`: setup, smoke tests, and verification gaps.
- `docs/context/SOLVER_CONTEXT.md`: maintained solver architecture and common pitfalls.
- `docs/context/PHYSICS_CONTEXT.md`: correlations, chemistry, radiation, materials, and calibration knobs.
- `docs/context/TRANSIENT_STATUS.md`: transient solver implementation status and numerical lessons.
- `docs/TECHNICAL_REFERENCE.md`: master technical reference for physics, numerics, materials, assumptions, and validation status.
- `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`: liquid/boiling coolant implementation status. Wired into the helical steady coupled march (co-flow validated); shell-and-tube and both transient solvers still postprocess-only/unwired. See `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` for phase status.
- `optimization/LLM_CONTEXT.md`: calibration and optimization context.
- `research/flamelet_kit/LLM_CONTEXT.md`: standalone flamelet/PFR context.

## Current Project Truths

- Package imports go through `hps_combustor`; `pyproject.toml` maps it to `1Dmodel/`.
- User steady entry point: `1Dmodel/main_steady.py`.
- User transient entry point: `1Dmodel/main_transient.py`.
- Maintained backend steady helical solver: `1Dmodel/main_solve.py`.
- Maintained backend steady shell-and-tube solver: `1Dmodel/main_solve_shellntube.py`.
- Maintained backend helical transient solver: `1Dmodel/main_solve_transient.py`.
- Maintained backend shell-and-tube transient solver: `1Dmodel/main_solve_shellntube_transient.py`.
- Every user-facing run should archive numeric data and input presets through `result_package.py`.
- `CorrelationCoefficients` in `1Dmodel/input_data.py` is the source of calibration knobs; do not rename fields casually.
- `data_master` must be fresh per solver instance; use `make_solver_data()`.
- Use SI units unless a file explicitly says otherwise.
- Liquid/boiling coolant physics (`1Dmodel/physics/liquid_flow/`) is real and literature-validated for its own scope. It IS wired into the helical steady coupled march (`main_solve.py`, co-flow validated; plain counter-flow has a known prescribed-outlet limitation, resolved via `solve_counterflow_liquid_reference()` — see the integration plan). It is still **postprocess-only** for the shell-and-tube steady solver and has zero presence in either transient solver. Do not assume `coolantProp.coolant_model="equilibrium_liquid"` changes solved results for shell-and-tube or transient runs.
- `coolantProp`/`combustorProp` defaults are the Helium/`single_phase_coolprop`/`shellntube` working baseline (confirmed with the user 2026-07-13: fluid-agnostic, Helium is the baseline, Water is a test configuration — do not hardcode a fluid's properties, like a molar mass, in shared defaults). For a working Water recipe, use `1Dmodel/validation/water_helical_example.py` rather than editing the shared defaults.

## High-Priority Technical Memory

- For steady and transient chemistry, finite-rate FPV is the required default for the current diesel/O2 high heat-extraction regime. Frozen chemistry is validation-only and is physically wrong for the user's main regime.
- Finite-rate chemistry is implemented through the FPV manifold in `physics/combustion_chemistry/fpv_manifold.py`; equilibrium remains a comparison mode.
- Equilibrium transient chemistry is implemented through precomputed enthalpy-removed manifolds/table lookup, so transient marches should not make per-node Cantera calls.
- Pilot ignition is represented as a low-diesel/O2 mass-flow segment in `schedule_mass_flow_g`; no PLA ignition surrogate is wanted.
- The helical and shell-and-tube transient solvers support `flow_config in ("co", "counter")`.
- Long shell-and-tube counter-flow transients should use
  `transientProp.solver_method = "fixed_step"`: it is linearly implicit in the
  local wall-film stiffness and avoids BDF Jacobian-probing cost.
- For helical counter-flow settle checks, enable
  `numericalProp.counterflow_physical_steady_reference` so the steady reference
  uses the same physical cold-end inlet boundary as the transient solver.
- Shell-and-tube steady/transient wall conduction must respect
  `hot_side="inner"` because hot gas is inside the tubes. Do not regress the
  corrected hot/cold perimeter mapping in `physics/heat_conduction.py`.
- Pre-ignition GOX chilldown is modeled with CoolProp Oxygen properties when
  `ignition=0` and LOX/GOX mass flow is scheduled.
- Finite-rate FPV manifolds are cached under `cache/fpv_manifolds`; keep this
  cache enabled for normal solver work.
- Hot gas transient modeling intentionally uses continuity plus Cantera/manifold enthalpy balance, not a quasi-1D momentum equation, because Mach_g is designed below about 0.3.

## Known Sharp Edges

- Persistent Cantera objects must be reset from cached inlet `T/p/Y` state before repeated sweeps or table builds. Do not read a prior sweep's cooled state as a new inlet.
- In shell-and-tube calculations, tube-side gas quantities are per tube. Divide total hot-gas mass flow by `N_tubes` for per-tube velocity and enthalpy accounting.
- Darcy friction convention is used in maintained solver paths. Do not apply a Fanning-to-Darcy factor.
- `inside_tube_choice="grooved"` uses the corrugated-tube Nu/friction path plus
  optional calibration multipliers; it is no longer a smooth-tube placeholder.
- `main_solver` (helical) requires `combustorProp.HX_config == "shellnHelicalTube"`
  and raises otherwise (added 2026-07-13). Its axial-length bookkeeping
  (`_advance_state()`) silently used a wrong linear-`dx` approximation for any
  other `HX_config` value — this bit every Phase 0-2 liquid-coolant test
  before the guard existed. The real helical coil is ~1378 arc-length nodes
  (not ~100) for this combustor's geometry; duty is ~150-300 kW, not ~20 kW.

## Verification Bias

For solver work, prefer focused scalar checks over dashboard inspection. For transient work, verify energy closure, settle-to-steady behavior, and chemistry mode behavior before tuning performance.

At the prior Claude session stop, two background validations were still in flight:
the helical low-flow freeze-out demo and the full shell-and-tube transient
validation. See `docs/context/TRANSIENT_STATUS.md`; do not assume either completed.

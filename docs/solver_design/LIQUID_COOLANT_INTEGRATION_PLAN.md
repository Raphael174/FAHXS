# Liquid Coolant Integration And Architecture Plan

Status: Phase 0, Phase 1, and Phase 2 (helical) DONE (2026-07-13, implemented
in Sonnet from the Fable-written plan below). Phase 3 (shell-and-tube shell
side) is next.

## Progress Log

- **Phase 0 (done)**: `tests/test_steady_baseline_regression.py` freezes
  helical + shell-and-tube helium scalars (co/counter, `chemistry_model =
  "frozen"` for speed/determinism) as a regression gate. 4/4 passing,
  ~6 s total, bit-identical across independent runs.
- **Phase 1 (done)**: created `1Dmodel/physics/liquid_flow/` with
  `correlations.py` (properties/HTC/friction), `chf.py` (Groeneveld LUT —
  split out, no fluid-property calls of its own), `dispatch.py` (the former
  `coolant_models.py`), `governing_equations.py` (the former
  `heated_liquid_channel.py`), and `hx_adapters.py` (the former
  `liquid_hx_adapters.py`). Also created `1Dmodel/physics/gas_flow/` and moved
  the ideal-gas `governing_equations.py` there unchanged, since it is a
  completely different equation set from the liquid one and the shared
  filename was the exact ambiguity this phase was meant to remove. Old paths
  (`physics/liquid_coolant.py`, `coolant_models.py`, `heated_liquid_channel.py`,
  `liquid_hx_adapters.py`, `governing_equations.py`) are now thin re-export
  shims that emit `DeprecationWarning` (to be deleted in Phase 4). Every real
  importer in the repo (`main_solve.py`, `main_solve_shellntube.py`, all
  `tests/`, all `1Dmodel/validation/*.py`) was updated to the new canonical
  paths; only `tests/test_liquid_flow_shims.py` deliberately still imports the
  old paths, to keep the shims themselves under regression coverage.
  `Helium_thermodynamics.py` was audited (dead, unimported, `#%%`-cell
  exploratory script) and left in place with a clarifying docstring rather
  than moved — moving dead code would have been unnecessary churn.
  Verification: full `pytest tests/ -q` — 122/122 passing, including the
  Phase 0 baseline unchanged and 5 new shim tests
  (`tests/test_liquid_flow_shims.py`).
- **Phase 2 (helical, done)**: `main_solver` (`main_solve.py`) now has a real
  `(p,h)`-state liquid branch gated on
  `coolantProp.coolant_model == "equilibrium_liquid"`, for co- and
  counter-flow. Per-node closure goes through
  `physics/liquid_flow/dispatch.evaluate_coolant_closure()` (state, HTC,
  friction, CHF margin); governing equations are `dh/dx = dQ/mdot` and
  `dp/dx = -friction` (see "Known simplification" below); the loop guard,
  gas-only-diagnostic nan-ing, and `p_c`/`enthalpy_c` state advancement all
  follow Design Decisions 1/6/7. New
  `physics/liquid_flow/sanity_checks.check_liquid_march()` runs automatically
  in `_check_global()` for liquid mode (energy closure, temperature ordering,
  saturation consistency, pressure monotonicity, bounds, and a **hard** CHF/
  dryout gate — no post-CHF closure exists). `coolantProp.liquid_chf_lut_path`
  added (CHF margin is `None` until a LUT path is supplied). New
  `data_master` keys: `enthalpy_c`, `dh_c__dx`, `quality_c`, `void_c`,
  `chf_margin_c`.
  Verified with `tests/test_liquid_coupled_helical.py` (6 tests): subcooled
  co/counter, boiling reached via heating (co-flow), a near-saturated
  counter-flow case, grid convergence (50 vs 100 arc-steps/turn, outlet
  quality changes <0.05), and a cross-check against the pre-existing
  postprocess bridge fed the identical converged duty profile (co-flow only,
  see limitation below) — matches to <1%, confirming the new coupled wiring
  and the already-validated postprocess bridge use the same closure
  self-consistently. Full suite: 128/128 passing, Phase 0 baseline still
  bit-identical.
  - **Known simplification**: the HEM acceleration pressure-gradient term is
    currently omitted from `dp/dx` (friction only). It requires a local
    quality *gradient*, which needs either a lookahead or a lagged finite
    difference; deferred rather than adding another lag on top of the
    heat-flux lag below. Small relative to friction pre-CHF for the validated
    water cases exercised so far, but should be added before treating results
    near CHF as quantitative.
  - **Known simplification**: the boiling HTC's heat-flux (boiling-number)
    term uses the *previous* node's wall flux (`self.q_w`, already available
    from `_unpack_heat_transfer_node()` before this node's own conduction
    solve) — a one-node-lagged closure. This breaks the circular dependency
    (HTC needs flux; flux needs HTC) at first order; accuracy improves as
    `dx -> 0`, confirmed by the grid-convergence test above.
  - **Discovered limitation, since resolved**: counter-flow's coupled march
    used the same legacy "prescribed-outlet" shortcut the gas coolant march
    already uses — it started from `coolantProp.T_out`/`p_out` as a
    single-phase `(T,P)` state, which cannot represent a genuine two-phase
    state and therefore acted as a hard ceiling on enthalpy. Design Decision
    2.4 is now implemented (see "Post-Phase-2 hardening pass" below):
    `solve_counterflow_liquid_reference()` converges the march's hot-end
    starting enthalpy against the user's physical `T_in`/`p_in`, removing the
    need to guess `T_out`/`p_out` at all for the liquid counter-flow case.

## Configuration-Drift Incident And HX_config Fix (2026-07-13)

Triggered by the user reporting they could not run the steady helical solver
with water themselves. Root-caused through several layers, all fixed:

1. **Immediate crash**: `coolantProp.T_out=650`/`p_out=13e5` are leftover
   helium-era numbers. At 13 bar, 650 K is superheated steam past complete
   vaporization (`quality ~= 1.21`), so the liquid-mode loop guard
   (`_coolant_flow_continues()`) correctly refused to run even one node —
   `main_solver` never touches these fields for gas mode's sake, but nothing
   validated them for liquid mode's very different validity range.
2. **Deeper bug, found while fixing #1**: `main_solver`'s axial-length
   bookkeeping (`_advance_state()`) branches on
   `combustorProp.HX_config == "shellnHelicalTube"` — that branch (true
   helical arc-length-to-axial-position mapping) vs. the `else` branch (naive
   linear `dx` accumulation, wrong for a wound coil) was silently selected by
   an ambient dataclass field with **no validation**. Every Phase 0-2 test
   and probe in this document constructed `combustorProp(flow_config=...)`
   without pinning `HX_config`, so all of it ran against whatever the file's
   *current* default happened to be — "shellntube" at the time (main_solver
   doesn't handle "shellntube"; that's `shellntube_solver`'s own label), so
   the WRONG (linear) branch was silently used throughout. `main_solver`
   now raises `ValueError` in `__init__` if `HX_config != "shellnHelicalTube"`.
3. **Consequence of #2**: the real helical coil is **~1378 arc-length nodes**,
   not the ~98 every earlier Phase 0-2 case assumed — about 14x longer,
   because a tightly-wound coil needs much more arc length than axial length
   to traverse the same combustor length. Duty is correspondingly ~150-300 kW
   for this geometry, not ~20 kW. **Every Phase 0-2 baseline value and test
   case parameter was recaptured/retuned** against the corrected geometry:
   `tests/test_steady_baseline_regression.py` (new hardcoded values, plus
   `HX_config` now pinned explicitly rather than relying on the file
   default — see that file's module docstring for why), and both
   `tests/test_liquid_coupled_helical.py` and
   `tests/test_liquid_counterflow_reference.py` (T_in/p_in/mass_flow_c
   rescaled to the real duty; `HX_config` pinned in every `combustorProp(...)`
   construction).
4. **Project decision, asked and answered**: `coolantProp`'s default fluid
   had also been changed to Water/`equilibrium_liquid` (and `combustorProp`
   to `shellnHelicalTube`/counter) directly in `input_data.py`, which is
   exactly what let #1-#3 surface undetected in shared defaults. Asked the
   user: keep Water as the permanent default, or revert to the historical
   Helium baseline? Answer: **fluid-agnostic; Helium is the baseline, Water
   is a test configuration — and no fluid's properties should be hardcoded
   in the shared defaults.** Reverted `coolantProp`
   (`coolant="Helium"`, `coolant_model="single_phase_coolprop"`,
   `mass_flow_c=150e-3`) and `combustorProp.HX_config="shellntube"` to the
   historical baseline. Removed the hardcoded `molar_mass=4.002602` field
   (helium's molar mass) from `coolantProp` entirely; `main_solver` now looks
   up molar mass from CoolProp for whatever `coolant` is actually configured
   (`self._coolant_molar_mass_g_mol`, computed once in `__init__`), so no
   fluid-specific constant needs to track whatever fluid the field says.
5. **Working water recipe**: `1Dmodel/validation/water_helical_example.py`
   is the confirmed-working, standalone recipe for testing water (does not
   touch shared defaults). `run_coflow()` is fast (starts directly at
   physical `T_in`/`p_in`, no guessing) — default `mass_flow_c=0.2` gives a
   safe boiling case (CHF margin ~3.7) *under production `finite_rate`
   chemistry specifically* (a `mass_flow_c=0.1` case validated only under
   `frozen` chemistry earlier in this document sits at CHF margin ~0.01 —
   effectively dryout — under `finite_rate`; chemistry model measurably
   changes duty in this regime, so a margin validated under one is not
   automatically safe under the other). `run_counterflow_physical()` uses
   `solve_counterflow_liquid_reference()` for a physically-anchored counter-
   flow result instead of guessing `T_out`/`p_out`.
6. **Full test suite** (134+ tests, several minutes for the slow shooting-
   reference cases) re-verified passing end to end after all of the above.

## Post-Phase-2 Hardening Pass (2026-07-13)

A code review after Phase 2 surfaced three findings and one explicit user
request; all four are resolved:

1. **Sanity-gate energy-balance off-by-one (bug, fixed)**:
   `data_master` records each node's state BEFORE that node's step advances
   it, so the last recorded node's `dQ` advances to an *unrecorded* final
   state. The energy-closure check was summing `dQ` over all nodes but
   comparing against the recorded (one-node-short) enthalpy span, producing a
   spurious ~1/N (~1%) "imbalance" that silently consumed the whole 2%
   tolerance. Fixed in `physics/liquid_flow/sanity_checks.py`: sum
   `dQ[:-1]`, matching what the recorded span actually reflects. Now closes
   to float precision (confirmed: 0.00% in the Phase 2 boiling case).
2. **Groeneveld CHF LUT re-read from disk every node (perf, fixed)**:
   `load_groeneveld_2006_lut()` in `physics/liquid_flow/chf.py` is now
   `functools.lru_cache`-wrapped (per resolved path string). Was ~98 file
   parses per helical run; would have been ~thousands per shell-and-tube
   sweep in Phase 3.
3. **`counterflow_physical_steady_reference` was gas-only and silently wrong
   for liquid (bug, fixed)**: `main_steady.py`'s `run_steady()` now dispatches
   liquid-mode counter-flow to the new `solve_counterflow_liquid_reference()`
   instead of the gas-only, temperature-shooting
   `solve_counterflow_physical_reference()` (invalid inside the two-phase
   dome).
4. **User request: converge counter-flow on the physical T_in/p_in inlet
   instead of guessing T_out/p_out.** Implemented as
   `solve_counterflow_liquid_reference()` in `main_solve.py`:
   - Shoots on the hot-end starting **enthalpy** (never temperature — Design
     Decision 1), targeting the enthalpy of the user's physical
     `coolantProp.T_in`/`p_in`. `T_out`/`p_out` are not used by this path.
   - New `main_solver.__init__(..., _liquid_enthalpy_hot_end_override=...)`
     hook (private, not a public `coolantProp` field) lets the shooting
     helper inject a starting enthalpy directly, bypassing the `(T,P)`
     conversion that cannot represent a two-phase state.
   - **HEM check (user-requested)**: empirically confirmed
     `residual(h_hot_end)` is monotonic (`tests/test_liquid_counterflow_reference.py
     ::test_hem_closure_is_monotonic_for_shooting`), which is what makes a
     single-bracket root-find valid.
   - **Root-finding method — two hardening iterations were needed**:
     - First attempt used secant iteration; a non-smooth region (near
       boiling onset) produced a wild overshoot that crashed CoolProp
       mid-march. Switched to **bracketed bisection** (a guess can never
       leave a bracket with confirmed opposite-sign residuals), per the
       plan's own "Known Risks" section written before Phase 2 started.
     - Bisection with a *fixed-size* initial bracket (a flat "150 K worth of
       cp" guess) still failed: the true root's distance from the cold-inlet
       enthalpy varies over orders of magnitude with mass flow (a high-flow
       case needed ~15-20 K worth of margin; a low-flow boiling case needed
       >150 K worth), so a fixed bracket missed the root entirely for the
       high-flow case (both ends landed on the same side). Replaced with an
       **adaptive bracket search**: start from a small physical probe
       (~2 K worth of margin above the cold inlet) and geometrically expand
       or shrink until the residual's sign flips, then bisect.
   - **Enthalpy floor safety guard (bug, fixed)**: the adaptive search's
     small initial probe could legitimately drive the march's `(p,h)` state
     past CoolProp's valid range in a single arc-length step, crashing
     `PropsSI` deep inside `_advance_state()` before the loop guard's
     quality/pressure checks ever ran. Added `self._liquid_min_enthalpy_J_kg`
     (computed and validated once at `__init__`, generously below the
     physical cold-inlet enthalpy — never legitimately reached), clamped in
     `_advance_state()` exactly like the existing pressure floor, and added
     to `_coolant_flow_continues()` as a third stop condition.
   - **Also fixed while wiring this**: `_advance_state()` previously only
     updated `enthalpy_c`/`p_c` in the liquid branch, leaving `T_c`/
     `quality_c`/`void_c` one node stale after the loop exited (each
     per-node closure eval happens at the TOP of the next iteration, which
     never runs after the last step). Now refreshed via `equilibrium_state_ph()`
     immediately after every advance, so `solver.T_c` etc. are correct
     post-`solve()` — this is what the shooting residual actually depends on.
   - Verified: `tests/test_liquid_counterflow_reference.py` (6 tests) —
     convergence for both a high-flow subcooled case and a low-flow boiling
     case, sanity-gate pass, HEM monotonicity, co-flow/helium fallback, and
     the `main_steady.py` dispatch routing.
   - **Documented residual limitation**: only enthalpy is shot; the hot-end
     starting pressure is approximated as `p_in` (friction drop is normally
     small, ~0.01-1 bar in validated water cases). Near saturation, where
     `T = Tsat(p)`, this can shift the converged cold-end temperature by more
     than the raw enthalpy tolerance suggests — a documented one-variable
     approximation, not a coupled `(h, p)` 2D root-find.

## Goals

Priority 1 — physical capability:

- Run water (liquid, with boiling) as the coolant in the **helical tube**
  (shell-and-helical-tube config), co-flow and counter-flow, in the coupled
  steady march — not just the current postprocess bridge.
- Run water as the coolant on the **shell side** of the shell-and-tube config,
  co-flow and counter-flow, in the coupled steady sweep.
- Enforce thermodynamic sanity gates on every liquid run so results are
  trusted for engineering decisions: energy closure, temperature ordering,
  saturation consistency, pressure monotonicity, quality/void bounds, CHF
  margin.

Priority 2 — architecture and rigour:

- Liquid physics lives in a clearly named subpackage
  (`physics/liquid_flow/` with `correlations.py` and
  `governing_equations.py`), not in four loosely named top-level modules.
- At most one steady and one transient user-runnable script; backend solver
  classes move to a `solvers/` package and lose the `main_*` naming.
- Simulation presets move to `.toml` files; `input_data.py` remains the
  schema/defaults, a loader applies a chosen preset per run.

Scope: **steady solvers only** (confirmed by the user 2026-07-13). Transient
liquid coolant is explicitly deferred; the transient solvers stay on the
gas/compressible path. A later plan should reuse the same `p,h` closure
inside `transient_core` (see
`docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`, Remaining Limits).

Fluid-generality (user note, 2026-07-13): the gas coolant side should be
gas-agnostic (any non-reacting gas, properties from CoolProp via
`coolantProp`), and the liquid side fluid-agnostic to the extent the
correlations allow — correlation validity ranges (Re, Pr, reduced pressure)
are the honest limit and must be surfaced, not hidden. See Design Decision 6.

## Current State (updated post-Phase-2, 2026-07-13)

- Liquid physics exists and is literature-validated standalone:
  `physics/liquid_flow/correlations.py` (properties/HTC/friction),
  `physics/liquid_flow/chf.py` (Groeneveld CHF LUT),
  `physics/liquid_flow/dispatch.py` (coolant state/closure dispatcher),
  `physics/liquid_flow/governing_equations.py` (`p,h` channel march),
  `physics/liquid_flow/hx_adapters.py` (geometry bridges),
  `physics/liquid_flow/sanity_checks.py` (thermodynamic gates). See
  `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`.
- **Helical steady (`main_solve.py`) now has a real coupled liquid march**,
  not just the postprocess bridge, for `coolantProp.coolant_model ==
  "equilibrium_liquid"`. Co-flow is fully self-consistent (verified against
  the postprocess bridge). Counter-flow works numerically and passes all
  sanity gates but inherits the same prescribed-outlet-vs-physical-inlet
  discrepancy the gas march already has (see Phase 2 progress log above).
- **Shell-and-tube (`main_solve_shellntube.py`) is still postprocess-only** —
  Phase 3, not yet done.
- It is postprocess-only. The helical coupled march
  (`main_solve.py::solver()`, ~lines 328-530) computes coolant properties with
  direct `PropsSI(T,p)` calls, `dispatch_nu_coil()` /
  `dispatch_friction_coil()`, and ideal-gas governing equations
  (`dT__dx_IdealGas`, `dp__dx_IdealGas_logical`, line ~478-480). The
  shell-and-tube shell side (`main_solve_shellntube.py::_shell_side_march`,
  ~line 291) is a temperature-only `cp` march, and `_shell_h_at` (~line 273)
  evaluates Bell-Delaware with single-phase `PropsSI(T,p)` properties.
- `coolantProp` already carries the liquid selectors (`coolant_model`,
  `liquid_heat_transfer_model`, `liquid_pressure_drop_model`,
  `liquid_chf_model`), unused by the march.
- Tests: `tests/test_liquid_boiling_poc.py`, `tests/test_liquid_hx_adapters.py`.
  No coupled-solve liquid test exists yet.

## Design Decisions (fixed up front)

1. **State variables for liquid mode are `(p, h)`, never `(T, p)`.**
   Temperature is not a valid state inside the two-phase dome (it plateaus at
   `Tsat(p)` while enthalpy keeps rising). All marching, convergence checks,
   and under-relaxation in liquid mode must operate on enthalpy (and
   pressure), with `T`, quality, void, density derived via the HEM closure in
   `coolant_models` / CoolProp. This is the single most important numerical
   rule of the whole integration; converging or relaxing on `T` will silently
   stall at saturation.
2. **Helium behavior must be bit-for-bit unaffected.** Every liquid branch is
   gated on `coolantProp.coolant_model == "equilibrium_liquid"`. The default
   `"single_phase_coolprop"` path keeps the existing code untouched. A frozen
   helium regression harness (Phase 0) guards this.
3. **The coupled liquid march reuses the already-validated closure**, not a
   re-implementation: property/HTC/friction evaluation goes through
   `evaluate_coolant_closure()` (dispatcher) and the correlation library, the
   same code path the standalone channel solver and the validation matrix
   exercise. The coupled solver supplies geometry, local `dQ` coupling, and
   the wall solve; the physics stays in one place.
4. **Boiling HTC on non-straight geometries is an explicit interim model.**
   Helical: single-phase liquid uses the existing curvature-aware coil path
   (Gnielinski-based; Mori 1967 is a low-Pr gas correlation and must not be
   the water default); saturated boiling uses Gungor-Winterton with an
   optional `CorrelationCoefficients` multiplier (`liquid_boil_htc_factor`,
   default 1.0) as the calibration hook. Shell side: single-phase keeps
   Bell-Delaware; boiling uses `max(h_bell_delaware_single_phase,
   h_gungor_winterton_pseudo1D)` and the result dict flags
   `shell_boiling_model = "interim_pseudo1d"`. Geometry-specific boiling
   closures (helical secondary-flow enhancement, shell-side bundle boiling —
   Jeong 2025, Kumar 2022, Moharana 2022 in `docs/reference`) are a later,
   separate correlation task; the plan keeps the dispatch seam ready for them.
5. **CHF is a hard gate, not just a diagnostic.** If local heat flux exceeds
   the Groeneveld 2006 CHF (with diameter correction), the run is marked
   `dryout_risk = True` and the sanity gate FAILS the run summary (the model
   has no post-CHF closure, so results beyond that point are not physical).
   Quality reaching 1.0 (complete vaporization) is likewise a hard flag.
6. **Fluid-agnostic on both sides.** The existing gas coolant path must work
   for any non-reacting CoolProp gas, not just helium: no hardcoded
   `"Helium"` strings or helium constants in solver/physics code — fluid name
   and molar mass come from `coolantProp`. The compressible governing
   equations (`dT__dx_IdealGas` etc.) assume a single-component ideal gas;
   that assumption (and the known Z ≈ 1.04-1.06 supercritical-He deviation
   noted in `main_solve.py`) must be stated where the equations live, not
   silently tied to helium. The liquid path is likewise any pure CoolProp
   fluid — the `(p,h)` HEM closure and governing equations are
   fluid-agnostic; only the correlations carry validity ranges (Re, Pr,
   reduced pressure, geometry). The dispatcher should emit a soft warning
   (not a failure) when a correlation is evaluated outside its published
   range, and the validation evidence remains water-anchored until other
   fluids are validated.
7. **Gas-only quantities are skipped in liquid mode.** `Mach_c`, `gamma_c`,
   compressibility `Z`, and the CoolProp minimum-temperature while-loop guard
   in the helical march are gas-path logic. Liquid mode replaces the loop
   guards with: pressure above floor, quality below 1, length below max.
   Stress checks are geometry/wall based and stay active.

## Phase 0 — Helium Baseline Regression Harness

Before touching solver code, freeze current behavior:

- New test `tests/test_steady_baseline_regression.py`: run both steady
  solvers (helical + shell-and-tube) x (co, counter) on the default
  `input_data.py` presets with frozen chemistry or the cheapest stable
  chemistry mode that avoids long FPV builds in CI-less local runs, and
  assert scalar outputs (`Q_tot`, `dp_c`, `T_c_out`, max wall temperature)
  against recorded values with tight relative tolerance (1e-6 where
  deterministic; document any nondeterminism found).
- These four cases are re-run at the end of every phase. Any drift on the
  helium path is a defect, full stop.

Acceptance: test passes twice in a row from a clean interpreter.

## Phase 1 — `physics/liquid_flow/` Subpackage (moves only, no behavior change)

Create `1Dmodel/physics/liquid_flow/` and move content:

| New module | Content from | Role |
|---|---|---|
| `liquid_flow/correlations.py` | `physics/liquid_coolant.py` | Saturation/HEM state, single-phase Nu/friction, Gungor-Winterton, Yu2002, Müller-Steinhagen-Heck, multipliers |
| `liquid_flow/chf.py` | CHF portion of `liquid_coolant.py` | Groeneveld 2006 LUT load/interp/diameter correction |
| `liquid_flow/governing_equations.py` | `physics/heated_liquid_channel.py` | `p,h` channel march, profile/HX-grid solvers, cell/node field mapping, diagnostics summary |
| `liquid_flow/dispatch.py` | `physics/coolant_models.py` | `CoolantState`, `evaluate_coolant_closure`, model-name validation |
| `liquid_flow/hx_adapters.py` | `physics/liquid_hx_adapters.py` | Helical + shell-and-tube duty bridges |
| `liquid_flow/sanity_checks.py` | new (Phase 2) | Thermodynamic gate functions |

Rules:

- **Audit every file's content before moving it — do not move by filename.**
  In particular, the existing top-level `physics/governing_equations.py` is
  the **compressible ideal-gas** equation set used by the helium/gas march
  (`dT__dx_IdealGas`, `dp__dx_IdealGas_logical`, ...). It must NOT be moved
  into or confused with `liquid_flow/governing_equations.py`; the two are
  unrelated equation sets. To remove the ambiguity, this phase also renames
  the gas-side file to `physics/gas_flow/governing_equations.py` (new
  `gas_flow/` subpackage, with a shim at the old path) and adds a module
  docstring stating its assumptions: single-component non-reacting ideal gas,
  quasi-1D. `physics/Helium_thermodynamics.py` is audited in the same pass:
  if its content is generic ideal/real-gas helpers, it moves to
  `gas_flow/gas_thermodynamics.py` with helium-specific constants replaced by
  `coolantProp`-driven values (Design Decision 6); if it is truly
  helium-specific, it keeps its name but gains a docstring saying so.
- Old module paths become 3-line deprecation shims (`from
  .liquid_flow.correlations import *` plus a `DeprecationWarning`) so
  existing tests/validation runners keep working; shims are deleted in
  Phase 4 when all importers are updated.
- Grep-and-update all importers (`tests/`, `1Dmodel/validation/`,
  `main_solve*.py`) to the new paths in the same change.
- No function bodies change in this phase.

Acceptance: full liquid test suite (`test_liquid_boiling_poc.py`,
`test_liquid_hx_adapters.py`) and the Phase 0 baseline pass; the
`liquid_validation_matrix` runner reproduces `all_passed: true`.

## Phase 2 — Coupled Steady Helical Liquid March (priority 1a)

Target: `main_solve.py::main_solver` (or its Phase-4 successor
`solvers/helical_steady.py` if phases are reordered — do NOT do both moves at
once; finish this phase on the current file layout).

2.1 State initialization: when `coolant_model == "equilibrium_liquid"`, the
coolant state is `(p_c, h_c)`; initialize from `T_in, p_in` (co-flow) via
`coolant_inlet_state()`. Counter-flow: see 2.4.

2.2 Per-node property/closure evaluation: replace the `PropsSI(T,p)` block and
`dispatch_friction_coil`/`dispatch_nu_coil` calls with one call to
`evaluate_coolant_closure()` given local `(p_c, h_c)`, mass flux, `Dh_ch`, and
heated perimeter. It returns `rho, T, quality, void, htc, dp/dx components`.
Single-phase liquid Nu keeps the curvature-aware coil enhancement
(Gnielinski + Dean-number factor already available in the coil dispatcher);
two-phase uses Gungor-Winterton per Design Decision 4. The wall conduction
call (`OneDimensionalSteadyConduction_ShellnHelicalTube`) is unchanged — it
receives `h_c` (film coefficient) and `T_c` exactly as today.

2.3 Governing equations: replace `dT__dx_IdealGas` / `dp__dx_IdealGas_logical`
with, per channel:

```text
dh/dx = (dq/dx) / (mdot / N_ch)                    # energy, exact for 1D
dp/dx = -(dp/dx)_friction - (dp/dx)_acceleration   # MSH + HEM acceleration
```

These already exist in `governing_equations.py` (liquid_flow); the solver just
integrates them with the existing node stepping and `sign` convention.

2.3b Gas-agnostic check (Design Decision 6): while touching the helical march,
verify the gas-mode branch has no hardcoded helium constants — sound speed,
molar mass, and property calls must come from `coolantProp.coolant` /
`coolantProp.molar_mass`. Fix any found; the Phase 0 helium baseline guards
against behavior change.

2.4 Counter-flow: physical boundary is `T_in` at the cold end while the march
starts at the hot end. Implement `solve_counterflow_liquid_reference()`
mirroring the existing `solve_counterflow_physical_reference()` shooting
pattern, but shooting on hot-end **enthalpy** until the cold-end enthalpy
matches `h(T_in, p_in)` (bisection/secant on a bracketed h-guess). Reuse of
the existing helper's loop structure is encouraged; the shot variable must be
h, not T (Design Decision 1).

2.5 Loop guards and data recording: liquid-mode while-loop guard per Design
Decision 7. Record liquid fields into `data_master` under new keys
(`h_c`, `quality_c`, `void_c`, `chf_margin`, `T_sat_c`) — do not overload
existing helium keys; extend `make_solver_data()` and the plotting/sanity
paths deliberately (this was an explicit prior warning in the status doc).
`Mach_c`-dependent sanity checks are skipped in liquid mode.

2.6 Sanity gates (`liquid_flow/sanity_checks.py`), run automatically at end of
`solver()` in liquid mode and reported in `HX_sizing_brief()`:

- Global energy closure: `|Σ dQ − mdot·(h_out − h_in)| / Σ dQ < 1e-3`.
- Temperature ordering at every node: `T_g > T_wg > T_wc > T_c` (heating).
- Saturation consistency: where `0 < x < 1`, `|T_c − Tsat(p_c)| < 0.1 K`.
- Pressure strictly non-increasing along flow; total Δp positive and finite.
- Bounds: `0 ≤ void ≤ 1`; equilibrium quality within `[-0.5, 1.0]`.
- CHF: `q_wall_local < CHF_local` everywhere, else FAIL with location.
- Cross-check: run the existing postprocess bridge on the converged `dQ` and
  assert the coupled march's `p,h` profile matches within 1% (same closure,
  same duty ⇒ must agree; this ties the new wiring to the validated path).

2.7 Validation cases (new runner
`1Dmodel/validation/liquid_coupled_helical.py`, artifacts under
`docs/validation/liquid_coupled_helical/`):

- Water, subcooled throughout (low duty): coupled result ≈ single-phase
  analytic `mdot·cp·ΔT` estimate; co and counter.
- Water reaching saturated boiling mid-length: all sanity gates pass; outlet
  quality grid-converges (dx halving changes `x_out` < 2%); co and counter.
- Helium baseline unchanged (Phase 0 harness).

Acceptance: all of 2.7 plus new pytest coverage
(`tests/test_liquid_coupled_helical.py`) exercising a real `main_solver` run
with `coolant_model="equilibrium_liquid"` for co and counter.

## Phase 3 — Coupled Steady Shell-And-Tube Shell-Side Liquid (priority 1b)

Target: `main_solve_shellntube.py::shellntube_solver`.

3.1 Replace `_shell_side_march(dQ_profile)` in liquid mode with a `p,h` march
on the HX grid using `solve_steady_heated_channel_on_hx_grid(...)` with real
Bell-Delaware geometry (`S_m` pseudo flow area, total outer tube perimeter as
heated perimeter — same mapping the validated adapter already uses) and
`coolant_enters_at` derived from `flow_config`. Return both the `T_shell`
array (derived from `h`) for the tube-side march and the full liquid fields.

3.2 Sweep convergence: keep the existing tube/shell sweep structure but
under-relax and measure convergence on the **enthalpy profile**
(`max|Δh|/h_fg_scale`), not `max|ΔT_shell|` (Design Decision 1). The
tube-side march continues to consume the derived `T_shell` profile unchanged.

3.3 `_shell_h_at` equivalent in liquid mode: per Design Decision 4 —
single-phase Bell-Delaware from `(p,h)`-derived properties; boiling region
uses the interim `max(...)` model with an explicit flag in results.

3.4 Sanity gates: identical list to 2.6 (shared `sanity_checks.py`), plus the
per-tube duty scaling check (`dQ_total = N_tubes · dQ_per_tube` — a known
sharp edge in this solver).

3.5 Validation runner `1Dmodel/validation/liquid_coupled_shelltube.py` with
the same case structure as 2.7 (subcooled, boiling, co/counter, grid
convergence, helium baseline unchanged) plus a sweep-convergence robustness
check: the sweep must converge for a case where boiling onset moves between
sweeps (the historically fragile scenario).

Acceptance: mirrors Phase 2, with `tests/test_liquid_coupled_shelltube.py`.

## Phase 4 — Solver Layout And Entry-Point Consolidation (priority 2)

Only after priority 1 is green.

4.1 New layout:

| New path | From |
|---|---|
| `1Dmodel/solvers/helical_steady.py` | `main_solve.py` |
| `1Dmodel/solvers/shelltube_steady.py` | `main_solve_shellntube.py` |
| `1Dmodel/solvers/helical_transient.py` | `main_solve_transient.py` |
| `1Dmodel/solvers/shelltube_transient.py` | `main_solve_shellntube_transient.py` |
| `1Dmodel/steady.py` | `main_steady.py` (single user steady script) |
| `1Dmodel/transient.py` | `main_transient.py` (single user transient script) |

4.2 Before moving, grep the whole repo for importers (`optimization/`,
`studies/`, `1Dmodel/validation/`, `tests/`, `research/`) and update them.
Leave `main_steady.py`/`main_transient.py`/`main_solve*.py` as one-cycle
deprecation shims re-exporting from the new locations (several external
memories/docs reference them); delete the Phase 1 physics shims here too.

4.3 Class names: keep `main_solver`/`shellntube_solver` importable under old
names via the shims, but rename canonically to `HelicalSteadySolver`,
`ShellTubeSteadySolver`, `HelicalTransientSolver`, `ShellTubeTransientSolver`.

4.4 Update all context markdowns (`CLAUDE.md`, `AGENTS.md`,
`CODEBASE_MAP.md`, `docs/context/*.md`, `docs/TECHNICAL_REFERENCE.md`) in the
same change — they are load-bearing for future agent sessions.

Acceptance: Phase 0 harness, full pytest suite, and
`python -m hps_combustor.steady` / `python -m hps_combustor.transient` smoke
runs all pass; no repo-wide references to the old module paths outside the
shims.

## Phase 5 — TOML Presets

5.1 `input_data.py` remains the single schema and source of defaults (the
dataclasses, including `CorrelationCoefficients` — calibration depends on
those field names). Add `1Dmodel/input_loader.py`:

- `load_preset(path) -> SimulationInputs` (a small container of all dataclass
  instances).
- TOML sections map 1:1 to dataclasses: `[coolant]`, `[hotgas]`,
  `[combustor]`, `[shelltube]`, `[numerical]`, `[transient]`, `[run]`,
  `[correlation_coefficients]`, `[system_requirements]`.
- Unknown section or key ⇒ hard error listing valid names (rigour: presets
  must never silently ignore a typo).
- Values are validated against the dataclass field types; units stay SI, and
  the preset header comment must say so.
- Use stdlib `tomllib` (Python ≥ 3.11); add `tomli` fallback only if the
  environment requires it.

5.2 CLI on both user scripts:

```powershell
python -m hps_combustor.steady --preset inputs/presets/water_helical_counter.toml
python -m hps_combustor.transient --preset inputs/presets/helium_shelltube_bangbang.toml
```

No `--preset` ⇒ current `input_data.py` defaults (unchanged behavior).

5.3 Starter presets under `inputs/presets/`:

- `helium_helical_co.toml`, `helium_helical_counter.toml`
- `helium_shelltube_co.toml`, `helium_shelltube_counter.toml`
- `water_helical_co.toml`, `water_helical_counter.toml`
- `water_shelltube_co.toml`, `water_shelltube_counter.toml`

The helium presets must reproduce the Phase 0 baseline scalars exactly (this
doubles as the loader's correctness test). `result_package.py` archives the
resolved preset (both the source `.toml` and the fully-resolved field dump)
with every run.

Acceptance: loader unit tests (round-trip, unknown-key rejection, type
validation), baseline-equivalence test, and one water preset run end-to-end
through `steady.py` with all sanity gates green.

## Phase 6 (deferred) — Transient Liquid

Not in this plan. Prerequisites established here: the `p,h` closure and
sanity gates are geometry- and time-scheme-agnostic, and `transient_core`
already has the finite-volume scaffolding. A future plan should add a
mass/energy `(p,h)` coolant model to `transient_core` mirroring
`compressible_coolant.py`.

## Known Risks And Mitigations

- **Saturation-dome property calls**: CoolProp `PropsSI('...','T',T,'P',p)`
  fails or is ill-defined inside the dome. All liquid-mode property access
  must go through the `(p,h)` closure (which uses quality-based calls). Any
  remaining `(T,p)` call in a liquid branch is a bug.
- **Boiling-onset sweep instability (shell-and-tube)**: onset location can
  oscillate between sweeps. Mitigation: enthalpy-based under-relaxation
  (3.2), and if needed a smaller `omega` in liquid mode.
- **Counter-flow shooting with boiling**: `h_cold_end(h_hot_end_guess)` can be
  non-smooth across boiling onset. Use a bracketing method (bisection) rather
  than pure secant/Newton.
- **Helium regression**: guarded by Phase 0; run it at every phase boundary.
- **Interim boiling closures on real geometry**: results in the boiling
  region carry an explicit interim-model flag until geometry-specific
  correlations are adopted; decision-making should weight CHF margin and
  energy balance (robust) over local HTC detail (interim).

## Documentation Deliverables Per Phase

Each phase updates: `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`
(status + new evidence), `docs/context/PHYSICS_CONTEXT.md` and
`docs/context/SOLVER_CONTEXT.md` (wiring status), `TESTING_CONTEXT.md` (new
tests/runners), and `docs/TECHNICAL_REFERENCE.md` Section 7.5/18/19. The
"postprocess-only" warnings written on 2026-07-13 must be revised as each
phase lands — stale safety warnings are as harmful as missing ones.

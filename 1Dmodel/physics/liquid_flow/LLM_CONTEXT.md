# physics/liquid_flow/ LLM Context

## Scope

The `(p,h)`-state real-fluid coolant path: single-phase, boiling/two-phase, and
supercritical closures, a fluid/geometry/regime-tagged correlation registry, CHF/dryout
lookup, sanity gates, and HX-integration adapters. Despite the package name, it is not
"liquid-only" — it also owns the supercritical family (McCarthy-Wolf, Taylor, Cheng2020,
Wang2023, Krasnoshchekov-Protopopov) and, as of this session, the gas-side
forced-convection closures used by shell-and-tube's tube-side hot gas. `docs/solver_design/FV_CORE_REWORK_PLAN.md`
flags the name itself as a documented misnomer, planned for rename to `physics/fluid/`
at Stage G with a deprecation shim (same pattern as the existing top-level
`physics/liquid_coolant.py` etc. shims). This is the largest, most actively developed
package in `physics/` — treat it as the reference implementation the new `core/` FV
rework is built on top of, not something being replaced.

## Contents

| File | Role |
|---|---|
| `registry.py` | `ClosureContext`, `ClosureRecord`, `ExtrapolationReport`, `REGISTRY`/`register()`/`get_record()`, `select_supercritical()` — the fluid-agnostic closure-selection mechanism: rank candidate closures by in-validity-range first, then tier, fluid-specificity, geometry/orientation match, and a final `priority` tiebreak. |
| `gas_closures.py` | **New this session.** Registers shell-and-tube's four tube-side gas Nu/friction correlations (`dispatch_nu_tube_straight`, `dispatch_friction_tube_straight`, `nu_corrugated_tube_vicente`, `friction_corrugated_tube_vicente`) into the same registry as delegating wrappers — no physics reimplemented, only registry-driven selection. Introduces two new regime tags, `"gas_forced_convection"` (HTC-returning) and `"gas_forced_convection_friction"` (Darcy-factor-returning). |
| `regime.py` | Fluid-agnostic regime classification for ANY CoolProp fluid at any `(p,h)`: subcritical (`subcooled_liquid`/`two_phase`/`superheated_vapor`) vs. supercritical (`supercritical_liquid_like`/`pseudo_critical`/`supercritical_gas_like`, split on `T` vs. `T_pc(p)`). Also `buoyancy_parameter` (McEligot-Jackson Bo*) and `acceleration_parameter` (McEligot Kv) — detection-only HTD (heat-transfer-deterioration) precursor flags, `real_fluid_state_ph()` (dome-based below `p_crit`, real-EOS flash above it). |
| `correlations.py` | Single-phase and two-phase property/HTC/friction correlations: `saturation_state`, `equilibrium_state_ph` (HEM two-phase closure), `liquid_single_phase_nusselt` (3.66 laminar / Gnielinski turbulent blend), `gungor_winterton_boiling_htc`, `muller_steinhagen_heck_friction_gradient`, `bergles_rohsenow_onb_wall_superheat`, `post_chf_dispersed_flow_htc`, `two_phase_sound_speed` (Wood's equation), plus the registered cross-check-only closure `rpe_dittus_boelter_8_24`. |
| `chf.py` | Critical heat flux / dryout: `groeneveld_2006_chf()` (trilinear LUT interpolation + diameter correction), `chf_regime()` (DNB vs. dryout label). The full 2006 LUT text table is supplied externally (`docs/reference/external/2006LUTdata.txt`), not encoded in the repo. |
| `coolprop_state_cache.py` | `get_cached_state()` — persistent low-level `CP.AbstractState` objects per (backend, fluid) pair, replacing repeated high-level `PropsSI` calls (~10x faster for HEOS; TTSE/BICUBIC opt-in, ~300-1000x faster but interpolation-approximate, off by default). |
| `dispatch.py` | `evaluate_coolant_closure()` — the integration entry point: closes state from `(p,h)`, branches subcritical vs. supercritical, returns a `CoolantClosureResult` (state, HTC, dp/dz, CHF margin, ONB margin, HTD risk, extrapolation report, cross-check HTC). Also re-exports `ThermoState`/`coolant_state_from_Tp`/`coolant_state_from_ph`/`coolant_inlet_state` from `hps_combustor.core.thermo` (moved there in Stage A of the FV core rework; this is a pure re-export, proven identical-object by `tests/test_core_thermo.py`). |
| `supercritical.py` | Property-ratio-corrected supercritical Nusselt closures (`mccarthy_wolf_nu`, `taylor_nu`, `wang2023_eckert_split_nu`, `cheng2020_supercritical_nu`, `krasnoshchekov_protopopov_nu`, `_conservative_bound_htc`), all registered at import. |
| `governing_equations.py` | `(p,h)`-state 1D heated-channel solver (`solve_steady_heated_channel`, `solve_steady_heated_channel_on_hx_grid`, `HeatedChannelCase`/`HeatedChannelResult`) — geometry-agnostic; state is always `(p,h)`, never `T`, because `T` plateaus at `Tsat(p)` inside the two-phase dome and is not a valid marching/convergence variable there. |
| `hx_adapters.py` | `solve_helical_coil_liquid_from_duty/_from_data_master()`, `solve_shelltube_shellside_liquid_from_duty/_from_tube_result()` — convert existing solver geometry + a converged `dQ` duty profile into the generic heated-channel interface. Explicitly self-documented as pseudo-1D integration bridges, not final geometry-specific boiling models (no centrifugal/secondary-flow enhancement for the coil; Bell-Delaware `S_m` as a placeholder shell-side flow area). |
| `sanity_checks.py` | `check_liquid_march()` — end-of-solve engineering gates (energy closure, temperature ordering, saturation consistency, pressure monotonicity, quality/void bounds, CHF margin, two-phase choking Mach). CHF margin and dryout are **hard failures**, not diagnostics-only, because there is no post-CHF/mist-flow closure. |

## Wiring status — READ BEFORE ASSUMING THIS CHANGES A RESULT

This is real, literature-validated physics for its own scope, but its integration into
the maintained solvers is uneven and must not be assumed uniform:

- **Helical steady (`main_solve.py`): wired into the coupled march.** When
  `coolantProp.coolant_model == "equilibrium_liquid"`, the march uses
  `evaluate_coolant_closure()` for properties/HTC/friction/CHF and `dh/dx = dQ/mdot`,
  `dp/dx = -friction` as the governing equations. Co-flow is fully self-consistent.
  Plain counter-flow has a known prescribed-outlet limitation (the legacy
  `T_out`/`p_out` guess is a single-phase `(T,P)` state and cannot represent a genuine
  two-phase starting point); `solve_counterflow_liquid_reference()` resolves this by
  shooting the hot-end starting enthalpy against the physical `T_in`/`p_in` instead.
  `check_liquid_march()` runs automatically at the end of a liquid-mode solve.
- **Shell-and-tube steady (`main_solve_shellntube.py`): corrected 2026-08-19, was stale.**
  `self._liquid_mode` (set from `coolant_model == "equilibrium_liquid"`) couples
  `evaluate_coolant_closure` directly into `_shell_h_at`/`_shell_side_march` — the same
  coupled-march pattern as the helical solver, not a postprocess-only bridge. Verified by
  running Water and LN2/supercritical-N2 cases through it (see
  `1Dmodel/validation/friday_shelltube_water.py`/`friday_shelltube_n2.py`), not just
  reading the code. `liquid_coolant_postprocess()` is a *separate*, additional opt-in
  diagnostic bridge on top of an already-converged gas-mode `dQ` profile — not the only
  path, and not what `coolant_model="equilibrium_liquid"` actually uses.
- **Both transient solvers (`main_solve_transient.py`, `main_solve_shellntube_transient.py`):
  zero liquid/boiling code path** — not even a postprocess bridge.
- Do not assume `coolant_model="equilibrium_liquid"` changes solved results for
  shell-and-tube or transient runs (CLAUDE.md is explicit on this point).
- `coolantProp`/`combustorProp` defaults remain the Helium/`single_phase_coolprop`/
  `shellntube` baseline; use `1Dmodel/validation/water_helical_example.py` for a working
  Water/boiling recipe rather than editing shared defaults.

See `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md` and
`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` for the authoritative status
writeup; `docs/context/PHYSICS_CONTEXT.md` for the correlation-level detail.

## The new closure-registry mechanism (this session, 2026-08-18)

`gas_closures.py` is the newest file in this package and is Slice D1 of the active
`FV_CORE_REWORK_PLAN.md` "Stage D" work. Context that matters for anyone touching it
next:

- It registers exactly the four correlation functions shell-and-tube's tube-side hot gas
  actually uses (confirmed by reading `transient_core/adapters_shelltube.py`'s
  `shelltube_tube_gas_film`, the only function that computes that film today), selected
  by `shellTubeProp.inside_tube_choice ∈ {"smooth", "grooved"}`. It deliberately does
  **not** register `dispatch_nu_coil`/`dispatch_nu_shell`/`dispatch_friction_coil`
  (helical-only) or Bell-Delaware (returns `(h, dp)` from one call against a whole
  geometry dict, not a single scalar from a `ClosureContext` — needs its own two-output
  closure protocol, deferred).
- Every registered closure in `gas_closures.py` is a thin wrapper that **delegates** to
  the existing, validated function in `heat_transfer_correlations.py`/
  `friction_correlations.py` — no physics is reimplemented. Bit-identical equivalence to
  the direct legacy calls is asserted in `tests/test_core_closures.py` (14 tests, `==`
  not `approx`).
- `ClosureContext` gained two new optional fields to support this:
  `corrCoeffs` (a `CorrelationCoefficients` instance — needed because the 21
  `CorrelationCoefficients` fields are load-bearing for `optimization/calibrate.py` and
  the gas closures need calibration knobs like `Re_transition_lo`/`_hi`) and `extra`
  (a dict escape hatch for closure-specific scalars with no natural home in the common
  bulk-property fields — raw axial position `x_m`, `roughness_m`,
  `corrugation_thickness_m`/`corrugation_pitch_m`). Both are optional and backward
  compatible with every existing `ClosureContext` construction site.
- Two new regime tags exist specifically so an HTC-returning and a friction-factor-
  returning `ClosureRecord` can never be ranked against each other by `select_supercritical`
  or `get_record`: `"gas_forced_convection"` and `"gas_forced_convection_friction"`. This
  is a deliberate, documented broadening of `ClosureRecord.callable`'s "always returns an
  HTC" contract — see `ClosureRecord`'s own docstring in `registry.py`.
- **Not yet wired into any production solver path** — additive only, zero existing call
  sites touched as of this session (per the plan doc's Stage D1 closure note). The next
  slice (D2) is `core/state.py` + `core/momentum.py` + a standalone coolant mass/energy
  kernel; D1 is a prerequisite building block, not itself a behavior change.

## Sharp edges

- **Groeneveld CHF LUT and Bergles-Rohsenow ONB are water-only fits** — the registry
  surfaces this via `validity_report_for`/`check_validity`, but the subcritical dispatch
  path itself does not hard-block a non-water fluid; check the extrapolation report.
- **Gungor-Winterton here is vertical/high-Froude only** — the horizontal/low-Froude
  correction is not implemented, so horizontal shell-and-tube liquid-side orientation
  effects are unvalidated.
- **No flow-regime map and no post-CHF/dryout degraded-HTC model** — post-saturation
  vapor is treated as plain single-phase vapor (`post_chf_dispersed_flow_htc` is
  explicitly a conservative simplification, not a validated dispersed-flow correlation).
- **Chen (1962), Shah, and Kandlikar boiling correlations are referenced in
  `docs/reference` but not implemented** — only Gungor-Winterton and the Yu2002 fit are
  coded (`yu2002_modified_anl_boiling_htc` exists in `correlations.py` but is not the
  active default; `liquid_heat_transfer_model` only accepts `"gungor_winterton"` today in
  `dispatch.py`).
- **`select_supercritical`'s in-range-first ranking** (changed 2026-07-17): an in-range
  generic closure now beats an out-of-range fluid-specific one — a fluid-specific
  correlation used far outside its own validated envelope can be less trustworthy than an
  in-range generic property-ratio form (Locke & Landrum 2008's bulk-reference correlations
  overpredict near the pseudo-critical line when extrapolated). Tier still outranks
  geometry/orientation tag matching (2026-07-19 fix — a broadly-tagged always-in-range
  fallback was previously out-scoring a genuinely validated-in-range correlation).
- **`x_over_D` vs. the new `extra["x_m"]`** are two distinct fields on `ClosureContext`
  consumed differently — `x_over_D` is Taylor's supercritical entrance-effect ratio;
  `extra["x_m"]` (added for `gas_closures.py`) is the raw axial position for the
  developing-length correction in `dispatch_nu_tube_straight`. Do not conflate them.
- **PDF text extraction in this repo drops minus signs** — every hardcoded supercritical
  correlation coefficient in `supercritical.py` is sourced from a rendered page image,
  not text extraction, after this bit Cheng2020, Krasnoshchekov-Protopopov, and
  McEligot-Jackson during development. Follow the same discipline for any new closure.

## Test coverage

`tests/test_liquid_boiling_poc.py`, `tests/test_liquid_hx_adapters.py`,
`tests/test_liquid_flow_shims.py`, `tests/test_liquid_coupled_helical.py`,
`tests/test_liquid_counterflow_reference.py`, `tests/test_core_closures.py` (the new
`gas_closures.py` equivalence suite), `tests/test_core_thermo.py` (the `dispatch.py`
re-export identity), `tests/test_coolant_models.py` (top-level shim). No test file
directly targets `chf.py`'s LUT interpolation against a known reference point beyond
what `test_liquid_boiling_poc.py` exercises.

## TODOs

No explicit `TODO`/`FIXME` code comments found in this folder; open items are tracked as
prose gaps in module docstrings (surfaced above) and in
`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` (Phase 3: shell-and-tube coupled
wiring, not started) and `docs/solver_design/FV_CORE_REWORK_PLAN.md` (Stage D2 onward).

## Change history

- Phase 1 (2026-07-13 per CLAUDE.md): moved from top-level `physics/liquid_coolant.py`
  etc. into this package; old locations left as deprecation shims.
  2026-07-17/2026-07-19: `select_supercritical` ranking fixes (in-range-first; tier
  before geometry/orientation) documented in `registry.py`/`supercritical.py` docstrings.
- 2026-07-31: `core/thermo.py` extraction (Stage A of the FV core rework) — `CoolantState`
  renamed `ThermoState` and moved to `hps_combustor.core.thermo`; `dispatch.py` now
  re-exports it unchanged (proven identical-object in `tests/test_core_thermo.py`).
  Same session: `rpe_dittus_boelter_8_24` cross-check closure added to `correlations.py`.
- 2026-08-18 (this session, Stage D Slice 1): `gas_closures.py` created; `registry.py`
  gained `ClosureContext.corrCoeffs`/`.extra` and the two new gas-forced-convection regime
  tags (all detailed above).

## Shell-side model corrections (2026-08-20)

Four changes to `main_solve_shellntube.py` / `bell_delaware.py` /
`shelltube_geometry.py` / `liquid_flow/correlations.py`, in the order they
matter:

1. **Bell-Delaware Sieder-Tate term enabled.** `(mu_b/mu_w)^0.14` had been left
   at its neutral default of 1.0 — the correlation's own property-variation
   correction was switched off. It is now evaluated from the lagged
   coolant-side wall temperature (`_shell_Tw_lagged`), clamped to [0.25, 4].
   Lowered the water design point's peak wall temperature by 67 K and resolved
   a convergence stall; moves the Helium shell-and-tube baseline ~0.1%.

2. **Regime-based closure dispatch.** Which shell-side closure a node used was
   previously decided by the `coolantProp.coolant_model` STRING, not by the
   coolant's state. At supercritical pressure the property-ratio closures
   (McCarthy-Wolf/Taylor) are now used only where the bulk->wall interval
   actually reaches the pseudo-critical band around `T_pc(p)`; elsewhere
   Bell-Delaware supplies `h`. Helium at 80 bar has `T_pc = 11.4 K` against a
   300-1400 K march, so it never needed a property-ratio closure; N2 at 88 bar
   has `T_pc = 147.8 K` with bulk 100-124 K and wall ~164 K, so it does.
   Latched one-way per node so the choice cannot chatter.

3. **Pressure march moved onto Bell-Delaware, plus the momentum term.** The
   closure gradients (MSH, supercritical registry) are straight-TUBE
   correlations: axial flow along one `L_tube`-long channel with skin friction.
   The real path crosses the bundle `N_baffles+1` times through `N_tcc` rows
   each — ~7.5x the path length here, with form drag — so they under-predict by
   ~25x. Two-phase nodes now take the ALL-LIQUID Bell-Delaware drop scaled by
   the **Chisholm two-phase multiplier** (`chisholm_two_phase_multiplier`);
   supercritical/single-phase nodes take Bell-Delaware at local properties.
   The **momentum (acceleration) term** `dp_acc = G^2 * delta(1/rho)` was
   omitted entirely and is worth **26% of the water design point's total drop**
   (2.99 of 11.05 bar) and **23% of N2's** (6.96 of 29.67 bar).

4. **`S_sb` factor-2 bug fixed** in `shelltube_geometry.py`. Shell-to-baffle
   leakage area is `Ds*Lsb*(pi - theta_ds/2)` (gap x arc length); the code
   carried a spurious extra `/2`. Confirmed against first principles and
   Hellborg (2017) eq. 47. `S_sb` 1.218 -> 2.436 cm2, `r_lm` 6.51 -> 6.93.

### Open items / known extrapolations

- **`r_lm = 6.93` is far outside Bell-Delaware's fitted range (`r_lm <~ 1`)** on
  this geometry: the tube-to-baffle clearance area (17.7 cm2) is ~6.5x the
  cross-flow area (2.9 cm2), driven by `clearance_tube_baffle = 1.0 mm` on a
  5 mm tube (2.5x TEMA practice) and a 12 mm baffle spacing. Both leakage
  corrections are therefore extrapolated — `J_l = 0.39` on heat transfer,
  `R_l ~ 0.01` on pressure drop. **This is the dominant uncertainty on shell-side
  dp, larger than any boiling-correlation choice.** The VDI/"modified Delaware"
  method normalizes leakage differently (`R_Q = (S_sb+S_tb)/(B*L_e)`) and lands
  at `f_L ~ 0.37` — independent corroboration that ~60% of ideal heat transfer
  is lost to leakage, so this is not an artifact of one correlation's algebra.
  Reported at runtime by `print_summary()`.
- **`chisholm_B` is reconstructed, not transcribed.** Chisholm (1973) /
  Grant & Chisholm (1979) were unavailable. Corroborated at exactly one point
  (Hellborg eq. 137 hardcodes `B = 21/Gamma`, which falls in the
  `9.5 < Gamma < 28` branch). Verify against the primary source before relying
  on a dp that depends on it. Note Hellborg's hardcoded branch is WRONG for this
  exchanger: water sits at `Gamma = 3.60`, `G ~ 3000 kg/m2s`, giving `B ~ 1.01`
  rather than 5.83.
- **Boiling HTC is still Gungor-Winterton, a tube-flow correlation.** A bundle
  correlation (Shah 2017, or Doo 2005 for shellside evaporation specifically)
  would be more defensible, but the measured leverage is small: the coolant film
  carries only **10.0%** of total thermal resistance in the dome (water) and
  **0.58%** (N2), so a +/-30% change in the boiling correlation moves UA by ~3%.
  Also note our shell-side mass flux (~3000 kg/m2s water, ~6600 N2) is far above
  the kettle-reboiler conditions those bundle correlations are fitted on.
- **Bell-Delaware's `_pick_row` switches j/f coefficients discontinuously** at
  Re = 10/100/1000/10^4. Hellborg (2017) flags this as unsuitable for iterative
  solving and smooths it. Verified NOT to bite at present: both design points sit
  wholly in the Re > 10^4 band (water 1.7e4 -> 7.6e5, N2 3.5e5 -> 1.9e6).

# Fluid-Agnostic Closure Architecture + Supercritical Flow Plan

> **STATUS 2026-07-19: Phase 0 and Phase 1 COMPLETE, closure-selection decision
> RESOLVED.** New modules: `physics/liquid_flow/regime.py` (regime detection,
> T_pc, McEligot-Jackson buoyancy/acceleration), `registry.py` (validity/fluid/
> geometry-tagged closures + in-range-first selection + extrapolation
> reporting), `supercritical.py` (McCarthy-Wolf and Taylor -- rocket-
> regenerative-cooling property-ratio correlations, Locke & Landrum 2008 -- now
> the default in-range pick; Cheng2020 demoted to its own narrow low-Re niche;
> Krasnoshchekov-Protopopov and Gnielinski as generic fallbacks). All
> coefficients visually verified from rendered source pages (minus-sign drops
> caught in Cheng2020, K-P, AND the McEligot-Jackson buoyancy formula). Bit-
> identical acceptance PASSED exactly, twice (before and after the closure-
> selection fix): helical helium 182.14810060417753, helical water
> 184.0092196238135, shell-and-tube water 323.60702108508826. LN2 case (100 K
> in, 80 bar, helical, co-flow) with the corrected closures: N2 100->1104 K
> (was ->1231 K with Cheng), Q=127.6 kW (was 142.9 kW), energy closure 0.056%,
> all sanity gates pass, 91 nodes HTD-flagged by the corrected (McEligot-
> Jackson) criterion (was 27, with a since-replaced, less literature-faithful
> formula -- not directly comparable).
> **UPDATE 2026-07-19 (later same day): shell-and-tube wiring + Wang2023 done.**
> `_shell_h_at`/`_tube_side_march`/`_shell_side_march`/`solve()` in
> `main_solve_shellntube.py` now route supercritical states through
> `evaluate_coolant_closure` exactly like the two-phase dome (merged branch:
> `quality_local` is NaN supercritically, which is now a trigger condition
> alongside `0<=quality<=1`); added a lagged shell-side wall-temperature array
> (`_shell_Tw_lagged`, same one-sweep-lag pattern as `_shell_qw_lagged`/
> `_shell_dp_lagged`) and a per-node shell-side `x_over_D` (measured from the
> coolant's own entrance, co- or counter-flow aware) for the property-ratio
> closures. All three `equilibrium_state_ph` crash sites fixed to
> `real_fluid_state_ph`. `print_summary()` reports closure name, regimes
> traversed, HTD-risk node count, and an explicit warning that shell-side
> cross-flow has no validated supercritical closure (McCarthy-Wolf/Taylor are
> straight-tube fits, applied here as an extrapolation) -- same honesty
> pattern as the pre-existing shell-side boiling gap.
>
> **Found and fixed a real registry ranking bug while validating this**: the
> shell-and-tube LN2 run initially selected `gnielinski_bulk_bound` instead of
> the in-range `mccarthy_wolf`, because the score tuple compared geometry-tag
> match BEFORE tier, and Gnielinski's broad `shell_crossflow` tag let a crude
> always-in-range fallback outrank a genuinely validated-in-range correlation
> that simply wasn't geometry-tagged for cross-flow. Fixed by moving tier
> ahead of fluid/geometry/orientation matching in `select_supercritical`'s
> score tuple (right after the `in_range` primary key). Effect on the LN2
> shell-and-tube case: T_wg max rose from 1878 K (buggy Gnielinski pick, no
> property-ratio correction) to 2328 K (corrected McCarthy-Wolf pick) -- the
> fix moved the result in the conservative direction, consistent with
> McCarthy-Wolf's `(T_s/T_b)^-0.55` term correctly reducing HTC for a heated
> wall. Verified against four hand-built scenarios (shell-side in-range,
> Cheng's own niche, helical high-Re, and a genuine out-of-range fallback) --
> all four now pick the intended closure.
>
> `wang2023_eckert_split_nu` is now registered too (Eq. 8-10, verified from
> the rendered page -- signs matched the original text extraction exactly this
> time, no drops found). Validity `p_reduced in (1.0, 1.1)` — confirmed it
> wins only its own near-critical niche (P/Pc~1.03) and never displaces
> McCarthy-Wolf at our actual 80 bar target (P/Pc~2.36).
>
> Bit-identical acceptance re-confirmed exactly after ALL of the above (three
> separate re-runs across this update): helical helium 182.14810060417753,
> helical water 184.0092196238135, shell-and-tube water 323.60702108508826.
>
> Phase 2 (FV core) not started -- deferred per user instruction.

## Closure selection policy -- RESOLVED 2026-07-19

Measured at a representative LN2 helical node (100 g/s, d=15.2 mm) with the
ORIGINAL registry ranking (fluid-specific-out-of-range beat generic-in-range):

| T_bulk [K] | Re | Cheng2020 | K-P | Gnielinski |
|---|---|---|---|---|
| 110 | 1.2e5 | 5064 | 1882 | 2116 |
| 145.7 (pseudo-crit) | 3.4e5 | 15639 | 1945 | 2605 |
| 250 | 4.8e5 | 6082 | 1146 | 1227 |
| 600 | 2.8e5 | 5145 | 1167 | 1241 |

Cheng's high fitted Pr exponent (1.3873) amplified HTC near the pseudo-critical
Pr spike, and the original ranking picked Cheng everywhere despite a ~300-770%
Re extrapolation. The user supplied three more references
(`2008_Locke_Landrum_...`, `Pizzarelli2018`, `JTHT2013`/`Urbano2013`) that
settled this:

- **Locke & Landrum (2008)** compared six H2 rocket-regen-cooling correlations
  against 2992 measurements at Re up to >1e6 -- our 1e5-6e5 sits inside that
  envelope, not 4-20x outside it like Cheng2020's 7000-27000. Their verified
  conclusion: bulk-reference property-ratio correlations (McCarthy-Wolf,
  Taylor) OVERPREDICT HTC around the pseudo-critical line when extrapolated --
  exactly the mechanism behind Cheng's spike above. This directly motivated
  switching the registry's ranking to **in-range-first** (was fluid-specific-
  first) in `select_supercritical()`: an in-range generic closure now beats an
  out-of-range fluid-specific one, while an in-range fluid-specific closure
  (Cheng, inside its own 7000-27000 window) still wins outright.
- **McCarthy-Wolf** (`Nu_b=0.025 Re_b^0.8 Pr_b^0.4 (T_s/T_b)^-0.55`) and
  **Taylor** (`Nu_b=0.023 Re_b^0.8 Pr_b^0.4 (T_s/T_b)^-(0.57-1.59/(x/D))`) are
  now registered, `FLUID_ANY`, validity `1e4<Re<1e7`, priority 2/1 respectively
  (McCarthy-Wolf wins ties as the simpler baseline; Taylor is the entrance-
  corrected alternative, better for P_r>4 which is not our case). Cheng2020's
  priority was raised to 3 so it still wins inside its own in-range niche.
- **Urbano & Nasuti (2013)** supplied the verified McEligot-Jackson buoyancy
  parameter (`Bo* = Gr*/(Re_b^3.425 Pr_b^0.8)`, `Gr* = gβ_b q_w D⁴/(k_b ν_b²)`,
  significant when `Bo*>6e-7`) and flow-acceleration parameter (`Kv =
  4q+/(Re_b^0.625 Pr_b^0.4)`, significant when `Kv>2.9e-5`) -- both used
  `regime.py`'s isobaric thermal expansion coefficient β, NOT the wall-bulk
  density-DEFECT Grashof the module previously used (a different, less
  literature-faithful formula). `htd_risk()` now flags either threshold. This
  is still detection-only: no deteriorated-HTC magnitude model exists for N2
  (Urbano-Nasuti's own `q_w/G_tr` threshold fit, Eq. 19-21, is methane-specific
  via CH4's `cp/β` behavior and was NOT adapted for nitrogen).
- **Pizzarelli (2018)** review's conclusion -- "no single heat transfer
  correlation has been found to properly describe" supercritical HTD --
  validates keeping the tagged multi-closure registry rather than committing to
  one fit.

Alternatives that were on the table before this: (a) keep Cheng, accept the
extrapolation (fluid-specificity captures N2 property behavior even out of Re
range); (b) prefer the in-range generic K-P above some extrapolation threshold;
(c) size the coil/flow
so Re lands in Cheng's 7000-27000 window (would need much lower mass flux or a
smaller channel). This is a real 3-8x wall-temperature-relevant choice.

---

Written 2026-07-17 (planning session). Three workstreams, in execution order:

- **Phase 0** — fluid-agnostic closure architecture (regime detection + correlation
  registry). Pure refactor, bit-identical acceptance.
- **Phase 1** — supercritical closure family, first client of the registry.
  Delivers the LN2 (90-100 K in, 80 bar, ~700 K out) coolant test case on the
  EXISTING steady solvers.
- **Phase 2** — solver unification: generic transient quasi-1D FV core with
  settle-to-steady as the canonical strongly-coupled path and a fast steady
  screening mode. Design-level here; own detailed design doc when reached.

Ordering rationale: supercritical N2 would otherwise become the *third* ad-hoc
fluid path (legacy helium (T,p) gas, water `equilibrium_liquid`, N2
supercritical). Building the registry first makes N2 the first client that
*proves* the fluid-agnostic design instead of another retrofit target. But the
**full solver unification (Phase 2) must NOT gate the N2 case** — Phases 0-1
ride on the existing solvers and deliver LN2 results in days, not weeks.

---

## Steady vs. transient-to-steady: decision

Question raised: should the codebase default to a fully transient model advanced
to steady state, dropping the dedicated steady solvers?

**Decision: unify the physics, not the drivers. Transient-settle becomes the
canonical path for strongly-coupled/BVP problems; the steady space-march is kept
as a fast screening mode and as the transient initializer.**

| | Steady space-march | Transient settle-to-steady |
|---|---|---|
| Co-flow (both states known at inlet end) | IVP — single pass, **~1.4 s** (helical, frozen chem) | needs O(100 s) physical settle, **~minutes** (memory: ~2.4 min per 100 s) |
| Counter-flow / prescribed-outlet | BVP — needs shooting (`solve_counterflow_liquid_reference`), fragile with two-phase states | **natural** — both physical inlets prescribed, integrate to steady |
| Shell-and-tube coupled sweep | under-relaxed sweep iteration; can oscillate without diagnosis (observed: water 50 g/s case stuck at ±7 K after 25 sweeps) | oscillation becomes *information*: either settles, or reveals a genuine limit cycle (density-wave-type instability) |
| Two-phase/supercritical instability physics | invisible / masked as "non-convergence" | captured (if time-accurate) |
| Design sweeps & calibration loops (dozens-hundreds of calls) | **indispensable** — 50-100x cheaper today | too slow as the only path |
| Physics duplication risk |今 duplicated across 4 solver files — recurring drift pain in this repo | single residual = single physics |

Key insight: the expensive, drift-prone thing to deduplicate is the **closure
and residual layer** (property dispatch, HTC/friction/CHF/regime logic, wall
conduction, chemistry coupling) — not the ~hundred lines of driver loop around
it. Phases 0-1 deduplicate the closures; Phase 2 deduplicates the residual;
drivers stay plural on purpose:

1. **Space-march driver** (steady, co-flow IVP): quick design checks, parametric
   sweeps, calibration. Also the default **initializer** for transient settles
   (start the settle from the marched steady field instead of a cold guess —
   large speedup, already conceptually proven by
   `counterflow_physical_steady_reference`).
2. **Transient driver, time-accurate**: transient studies proper + instability
   detection near suspicious operating points.
3. **Settle-to-steady driver** (transient integrator + steadiness criterion,
   later optionally pseudo-transient acceleration: implicit large-Δt once
   transients decay, wall-heat-capacity scaling): canonical steady answer for
   counter-flow, strong coupling, and final design verification.

Practical rule after Phase 2: *screen with the march, certify with the settle.*
A disagreement between the two is itself a finding (usually multiple steady
states or an instability), not a bug.

---

## Current-state inventory (what Phase 0 must preserve)

- `coolantProp.coolant` = any CoolProp name; `coolant_model` ∈
  {`single_phase_coolprop`, `equilibrium_liquid`}; `liquid_property_backend` ∈
  {HEOS, TTSE, BICUBIC} (validated 2026-07-16).
- `physics/liquid_flow/`: (p,h) state closure (`equilibrium_state_ph`),
  correlations (Gungor-Winterton, Müller-Steinhagen-Heck, Groeneveld 2006 LUT,
  Bergles-Rohsenow ONB, post-CHF vapor bound, Wood sound speed), dispatcher
  (`evaluate_coolant_closure`: quality-branched, CHF-gated, ONB-diagnosed,
  blend-windowed), `coolprop_state_cache`, `sanity_checks.check_liquid_march`.
- **Hard subcritical assumption**: `saturation_state` raises for p ≥ Pc;
  quality/void/dome machinery presumes a dome. Confirmed: `equilibrium_state_ph`
  crashes immediately for supercritical N2 at 80 bar.
- **Silently water-only**: Groeneveld LUT (no fluid parameter — raw water fit),
  Bergles-Rohsenow constants (1964 water fit). The registry must make these
  restrictions *explicit and machine-checked* instead of implicit.
- Materials: 316L property tables clamp flat below 27 °C — no cryogenic data.
- Four maintained solver files each carrying partially duplicated physics.

---

## Phase 0 — closure registry & regime detection (pure refactor)

New modules under `physics/liquid_flow/` (package name kept to avoid import
churn; it is now a documented misnomer — rename/alias deferred to Phase 2
cleanup):

### 0.1 `regime.py` — thermodynamic + heat-transfer regime detection
- `thermo_regime(fluid, p, h) -> ThermoRegime`:
  - p < Pc: `subcooled_liquid` / `two_phase` / `superheated_vapor` (delegates to
    existing `equilibrium_state_ph`).
  - p ≥ Pc: `supercritical_liquid_like` (T < T_pc), `pseudo_critical`
    (|T − T_pc| inside a band), `supercritical_gas_like` (T > T_pc).
- `pseudo_critical_state(fluid, p)`: T_pc(p) by bounded cp-maximization
  (bracket [T_crit, ~2·T_crit], cached per (fluid, rounded p) — measured for N2
  at 80 bar this session: T_pc ≈ 145.7 K, cp_peak ≈ 3479 J/kg·K).
- Buoyancy/HTD indicators (used by Phase 1, defined here):
  `Bu = Gr̄_b/Re_b^2.7` with threshold 1e-5 (vertical) and `Gr/Re² < 1e-3`
  (horizontal) — Hall-Jackson criterion as reproduced with coefficients in
  Ting Part II (docs/reference/Ting_part2_2024.md, Eq. 31-32).
  Eckert number `E = (T_pc − T_b)/(T_w − T_b)` (Wang2023/Yamagata regime split).

### 0.2 `registry.py` — correlation registry + selection policy
Each closure is a record: `name`, `callable`, `regime_tags`, `geometry_tags`
(straight_tube / helical_coil / shell_crossflow), `orientation_tags`,
`validity` (ranges over P/Pc, Re, G, q, D, quality, and an explicit fluid list
or "any"), `provenance` (docs/reference stem), `tier`
(`validated_in_range` / `structural_extrapolation` / `conservative_bound`).

Selection policy given (state, regime, geometry, operating point):
1. filter by regime + geometry (+ orientation when tagged);
2. rank: fluid-specific & in-validity > fluid-specific out-of-validity >
   generic in-validity > generic out;
3. **always** return the pick plus an `extrapolation_report` (which limits are
   violated, by how much) — this feeds the same honest-warning channel as the
   existing CHF/ONB messages. `coolantProp` gets one optional override knob
   (`closure_override: name`) instead of per-model string fields; existing
   fields (`liquid_heat_transfer_model` etc.) keep working as forced overrides.

Registered at Phase 0 (existing physics, no behavior change): Gungor-Winterton
(any fluid, two-phase, tube), MSH friction (any), Groeneveld LUT (**Water
only** — now explicit), Bergles-Rohsenow ONB (**Water only** — now explicit),
single-phase Gnielinski/Darcy (any), post-CHF vapor bound (any,
`conservative_bound`).

### 0.3 `dispatch.py` becomes an orchestrator
`evaluate_coolant_closure` keeps its exact signature; internally: state →
regime → registry → assemble `CoolantClosureResult` (gains `regime`,
`closure_name`, `extrapolation_report` fields, all defaulted so existing
consumers are untouched). `CoolantState` gains `T_pc`, `p_reduced`,
`is_supercritical`; `quality`/`void_fraction` become NaN/0 in supercritical
regimes (never raise).

### 0.4 `sanity_checks.py` regime-aware
Dome-specific gates (saturation consistency, quality bounds, CHF, dryout) run
only when any node is subcritical; supercritical nodes get placeholder HTD/
Mach/energy gates (filled in Phase 1).

**Acceptance gate (hard):** helium gas mode + water `equilibrium_liquid`, both
solvers, bit-identical outputs vs. pre-refactor (same discipline as the
2026-07-16 cache refactor); full liquid validation matrix unchanged; new unit
tests for regime detection and registry ranking/extrapolation-reporting.

---

## Phase 1 — supercritical closure family (LN2 test case)

### 1.1 Verify-before-hardcode checklist (blocking, first task)
PDF text extraction **drops minus signs** (confirmed on cheng2020.md). Before
implementation, render the equation pages to images and visually confirm:
- Cheng2020 Eq. 18-19: signs of the `Gr*` exponent (extracted "0.0013",
  almost certainly −0.0013) and the `(μ_b/μ_w)` exponent (extracted "0.8709",
  sign uncertain), and the exact Pr definition used (plain bulk Pr vs.
  enthalpy-averaged).
- Jackson property-ratio correlation coefficients as reproduced in Ting Part II
  (use it only with the exact form verified there; cite Ting as the source
  since the 1979 primary is not in the repo).

### 1.2 New closures (`physics/liquid_flow/supercritical.py`)
- `cheng2020_supercritical_nu` — **primary for the N2 case**. Fluid: Nitrogen.
  Validity: 7000 < Re < 27,000; P/Pc ≈ 2.06-2.65 (7-9 MPa) — brackets our
  80 bar (P/Pc ≈ 2.36); vertical 20 mm tube; wall-viscosity ratio term needs
  lagged wall temperature (same one-node/one-sweep lag pattern as the boiling
  Bo term; shell-and-tube adds `_shell_Tw_lagged` beside `_shell_qw_lagged`).
- `jackson_property_ratio_nu` — generic fallback, any fluid, tier
  `structural_extrapolation` outside its water/CO2 fit envelope.
- `wang2023_eckert_split_nu` — N2, validity 1.0 < P/Pc < 1.1 (will be
  auto-flagged out-of-range at 80 bar; registered as cross-check only).
- Plain Gnielinski at bulk properties — tier `conservative_bound` envelope.
- Friction: single-phase Darcy with local (p,h) properties (existing path);
  property-ratio friction correction only if coefficients are verifiable in
  Ting Part II, else skip.
- Sound speed: real-EOS `speed_sound` (already valid supercritically); no dome
  → Wood's equation never selected. Mach/choking gate carries over unchanged.

### 1.3 HTD (heat-transfer deterioration): detect, don't fabricate
Report `htd_risk` (buoyancy parameter vs. threshold, orientation-aware) in the
closure result and solver summaries — mirror of the ONB/CHF warning pattern.
We have **no validated closure for deteriorated-HTC magnitude at our
conditions**; the plan is to flag onset honestly and treat flagged runs as
outside design envelope, not to invent a degraded-HTC model. (Wire-insert
mitigation exists in the literature — Wang2020 — if HTD turns out to bite.)

### 1.4 Solver integration
- `coolant_model="equilibrium_liquid"` semantics generalize: the (p,h) march
  engages for supercritical inlets too (march mechanics — enthalpy march,
  lagged q_w, pressure march — are already regime-agnostic; the crashes were
  all in closure calls, fixed by Phase 0). Alias name `real_fluid_ph` accepted;
  old name kept working.
- Both maintained steady solvers (helical + shell-and-tube) via the shared
  dispatcher; expected solver-file diff is small (lag array + summary lines).
- Materials: print a warning when computed wall temperature falls below the
  316L table floor (27 °C) instead of silently clamping — cryogenic material
  data acquisition stays an open item.

### 1.5 Validation gates
1. Closure-level: reproduce Cheng2020 HTC magnitudes at their stated conditions
   (their Fig. 6/8 operating points); Dimitrov1989 qualitative two-peak HTC
   shape at 4 MPa.
2. March-level: energy-balance closure to machine precision, both solvers,
   supercritical N2.
3. **The LN2 case**: 90-100 K in, 80 bar, hot-gas side per current baseline.
   Cross-check duty against this session's hand calc (~0.90-0.93 MW for 30 L/s
   at 700 K out). **Early checkpoint: compute the actual channel Re** — if it
   falls outside Cheng's 7000-27,000 window the extrapolation report must say
   so and Jackson-form becomes the ranked pick; find this out first, not last.
4. Non-regression: water + helium suites still bit-identical.

### Known gaps carried openly (not blockers, must stay visible)
- Horizontal-tube and shell-side cross-flow supercritical HT: **no literature
  in hand for any supercritical cryogen**. Same treatment as shell-side boiling
  enhancement: documented gap, straight-tube closure applied, flagged.
- Plain-circular helical-coil supercritical enhancement: Fan2024 is cruciform
  geometry (per explicit user caveat) — evidence secondary-flow enhancement is
  real (+4-13%), **not** a transplantable coefficient. Apply Cheng unmodified,
  flag the geometry extrapolation.
- HTD magnitude closure (see 1.3). Cryo material properties (see 1.4).

---

## Phase 2 — generic quasi-1D FV core (design outline only)

Prereq: Phases 0-1 merged; all four solvers consuming the shared closure layer.
Detailed design doc to be written when this phase starts; constraints already
learned in this repo (see docs/context/TRANSIENT_STATUS.md): linearly-implicit
fixed-step beats BDF for wall-film stiffness; per-node Cantera calls forbidden
in marches (manifold/tabulated chemistry only); tabulated radiation.

Staged migration, each stage shippable:
1. `core/`: state vector + mesh + **single spatial residual** (mass/energy now,
   momentum per current low-Mach design; quasi-1D momentum only if/when the
   generic-NS backlog item demands it) with pluggable closure layer from
   Phase 0.
2. Drivers over that residual: transient (time-accurate), settle-to-steady
   (with steadiness criterion + march-initialized start), space-march steady
   (co-flow IVP fast path).
3. Migrate solvers one at a time (shell-and-tube transient first — newest, best
   documented), acceptance = reproduce that solver's validated results within
   stated tolerance before deleting its private physics.
4. Cleanup: rename/alias `liquid_flow` → fluid-agnostic package name; retire
   duplicated physics; single closure provenance table in docs.

---

## Effort & session plan

| Phase | Estimate | Model |
|---|---|---|
| 0 — registry/regime refactor | ~1 session | Opus (implement) |
| 1 — supercritical closures + LN2 validation | 1-2 sessions (solver runs dominate) | Opus (implement) |
| 2 — FV core | multi-week, staged; own design doc first | Fable (design) → Opus |

Risks: T_pc search robustness very close to Pc (bracket + cache, fail loud);
Cheng Re-window miss on our geometry (checked early, 1.5.3); sweep-iteration
convergence with a second lagged field (wall T) on the shell side (same
under-relaxation lever available); PDF sign verification is blocking for any
hardcoded coefficient (1.1).

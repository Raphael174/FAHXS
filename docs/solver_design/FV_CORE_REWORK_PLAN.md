# FV Core Rework — Fluid-Agnostic Quasi-1D Finite Volume

**Status: Stage A COMPLETE (closed 2026-08-18). Stages B and C turn out to be
ALREADY LARGELY BUILT (undocumented until now — see the 2026-08-18 note
below). Stage D IN PROGRESS, split into slices D1-D4 (see the second
2026-08-18 note below); D1 (closure registry unification for shell-and-tube's
tube-side gas closures) is done. Supersedes "Phase 2" of
`FLUID_AGNOSTIC_CLOSURES_AND_SUPERCRITICAL_PLAN.md`, which is the prerequisite
(Phases 0–1 complete: regime detection, closure registry, supercritical
family).**

> **2026-08-18 — Stage D split into slices; D1 done (tube-side gas closures
> unified into the registry).**
>
> Stage D as specced (`core/residual.py` + `drivers/transient.py` in one pass)
> was too large to execute and review as one unit, especially once "unify the
> closure mechanism into the registry" (decided over the facade alternative —
> real risk to `optimization/calibrate.py`'s calibration paths otherwise) was
> in scope. Split into a roadmap:
>
> | Slice | Deliverable | Status |
> |---|---|---|
> | D1 | `core/closures.py` + registry extension; shell-and-tube tube-side gas closures | **DONE** |
> | D2 | `core/state.py` + `core/momentum.py` + standalone coolant mass/energy kernel | not started |
> | D3 | `core/residual.py` — wall + coolant + closures + hot-gas march, one time step | not started |
> | D4 | `core/drivers/transient.py` + repoint `main_solve_shellntube_transient.py` | not started |
>
> **D1, done**: confirmed by reading `shelltube_tube_gas_film`
> (`transient_core/adapters_shelltube.py:405-517` — the only function that
> computes shell-and-tube's hot-gas film today) that shell-and-tube uses
> exactly four correlation functions, selected by `shellTubeProp.
> inside_tube_choice ∈ {"smooth", "grooved"}`: `dispatch_nu_tube_straight` /
> `dispatch_friction_tube_straight` (smooth) and `nu_corrugated_tube_vicente` /
> `friction_corrugated_tube_vicente` (grooved). NOT `dispatch_nu_coil`/
> `dispatch_nu_shell`/`dispatch_friction_coil` — those are helical-only and
> shell-and-tube never calls them, narrowing "unify the registry" from the
> originally-sketched 6 correlations down to these 4. Bell-Delaware
> (shell-side) stays a direct call — it returns `(h, dp)` from one call
> against a whole geometry dict, not one scalar from a `ClosureContext`;
> forcing it into "one closure = one scalar" needs its own two-output closure
> protocol, deferred rather than shoehorned into D1.
>
> Added to `physics/liquid_flow/registry.py`: `ClosureContext.corrCoeffs`
> (calibration knobs) and `.extra` (dict escape hatch for closure-specific
> scalars — raw axial position, roughness, corrugation geometry — with no
> natural home in the common bulk-property fields), both optional and
> backward compatible. Two new regime tags: `"gas_forced_convection"`
> (HTC-returning) and `"gas_forced_convection_friction"` (Darcy-factor-
> returning — a deliberate, documented broadening of `ClosureRecord.callable`'s
> "returns an HTC" contract, kept as a separate regime tag specifically so an
> h-returning and f-returning record can never be ranked against each other).
> Four new `ClosureRecord`s in the new `physics/liquid_flow/gas_closures.py`,
> each **delegating** to the existing validated function (no physics
> reimplemented), selected by name (not inferred ranking) via the new
> `core/closures.py::tube_htc_closure`/`tube_friction_closure`.
>
> Caught while writing the equivalence tests
> (`tests/test_core_closures.py`, 14 tests, all bit-identical `==` not
> `approx`): the grooved (corrugated) functions' own default `Re_lo`/`Re_hi`
> (2000/4000) differ from `CorrelationCoefficients`' (2300/4000), which the
> legacy adapter always passes explicitly
> (`adapters_shelltube.py:460-461,467-480`) — a naive wrapper using the
> functions' own defaults would have silently shifted the laminar/turbulent
> blend band. Locked in with a dedicated regression test
> (`test_grooved_re_thresholds_come_from_corrcoeffs_not_function_defaults`).
>
> Full suite after D1: 233 passed (was 219), same 4 pre-existing frozen-
> chemistry failures, unchanged. Zero existing call sites touched — the new
> closures are additive and not yet wired into any production solver path.
>
> D2 must decide whether to carry the legacy's CFL-subcycling pattern forward
> or make coolant advection genuinely implicit. D3/D4 must decide, explicitly,
> whether to reproduce two legacy shell-and-tube quirks found investigating
> the CFL fix above (`dq_cold`/`T_wc` discarded in favor of `G·(Tbar−Tc)`;
> `mdot_effective` a mean-of-faces scalar fed to per-cell closures) or fix
> them — default is reproduce-then-flag, since the stage gate is "reproduce
> ... within stated tolerance", not redesign.

> **2026-08-18 — Stages B/C found already built; Stage D blocked on a
> pre-existing regression, now fixed.**
>
> Starting Stage D work surfaced that `core/mesh.py` (`FlowPath`,
> `StreamCoupling`, `HXAssembly`), `core/wall.py` (`CylindricalWall` + analytic
> rib/fin), and both `core/geometry/shell_and_*.py` builders already exist,
> with 46 passing tests including the Stage B conservative-overlap gate and
> Stage C's `rel=1e-12` cross-validation against `physics/heat_conduction.py`.
> This section previously said `core/mesh.py` "doesn't exist yet" — it was
> stale. The one Stage B item still outstanding is "legacy solvers consume
> geometry from the builders", explicitly deferred to Stage E by
> `shell_and_helical_tube.py`'s own docstring.
>
> **Stage D's own migration target was broken.** Its gate is "reproduce the
> validated shell-and-tube transient" (`solve_transient_core` via
> `transient_core/adapters_shelltube.py`) — that path crashed with
> `FloatingPointError: coolant internal energy left the configured CoolProp
> temperature range (60-2500 K)` on every case in its own documented
> validation matrix (`1Dmodel/validation/transient_core_short_runs.py`),
> confirmed unrelated to the Stage A rewire (calls none of the edited methods).
>
> **Root cause**: `transient_core/compressible_coolant.py::
> conservative_mass_energy_step` is explicit forward Euler in the coolant
> conserved variables — unconditionally unstable once a macro step exceeds
> roughly one cell's residence time (`mass/mdot`). The documented validation
> matrix runs at `max_step`/`tau` ≈ 50–200 (confirmed by direct probing: an
> amplifying single-cell temperature spike originating at the inlet cell,
> traveling downstream step by step — the textbook explicit-upwind CFL
> signature). Reproduces identically at `max_step` from 0.25 down to 0.005 s
> (50× reduction), ruling out "just use the documented step size" as a fix.
> The 2026-07-10 stored validation numbers this file's matrix once passed with
> were themselves already inside the CFL-violated regime (archived at
> `docs/validation/transient_core_short_run_results_PRE_CFL_FIX_2026-07-10.json`
> for reference) — the guard's reported `energy_residual_abs_max` of 30–93 J
> per case was the unresolved CFL error, not real energy non-closure; it just
> hadn't yet pushed a cell property past the (also independently wrong —
> 2500 K exceeds CoolProp's real 2000 K ceiling for He/N₂/Water, so it was
> silently extrapolating rather than protecting) guard bound.
>
> **Fix**: `adapters_shelltube.py::_cfl_stable_substep_count` subdivides a
> macro step into CFL-safe substeps for the mass/energy advection only —
> momentum, the hot-gas march, and wall/coolant conductance stay frozen across
> substeps (the same quasi-steady-per-macro-step assumption already documented
> for the hot side). The outer time grid callers see is unchanged; diagnostics
> sum over substeps. Fixed for both `quasi_steady` and `low_mach` momentum
> models (same shared kernel). Regenerated the validation matrix: energy
> residual drops from O(10–90) J to O(1e-8) J on every case; headline
> `T_c_out_peak` moves +0.5% (bangbang cases) to −7% (GOX cases — larger
> because lower flow there meant the pre-fix CFL violation was *more* severe,
> not less). **Performance cost is real**: 45–140 s/case now vs <1 s/case
> before, because CFL now honestly resolves a step that used to silently (and
> wrongly) skip most of the physics. Test: `tests/test_transient_core_shelltube.py`.
>
> **Update, 2026-08-19: the helical transient core is now fixed too.**
> `main_solve_transient.py`'s own inline mass/energy loop (separate code
> from `adapters_shelltube.py`, not the file the original fix touched) had
> the identical `FloatingPointError` signature — confirmed to be the exact
> same root cause (not just "almost certainly"): running it at a CFL-safe
> `max_step` produced sane, converged results, and the guard trips
> immediately at the config that previously crashed. Fixed with the same
> mechanism: `_cfl_stable_substep_count` (still defined once, in
> `transient_core/adapters_shelltube.py`) is now also imported into
> `main_solve_transient.py` — that file already imported other private
> helpers from `adapters_shelltube.py`, so this follows existing precedent
> rather than introducing a new cross-module dependency. Verified for both
> `quasi_steady` and `low_mach` momentum models; energy residual at true
> machine precision (~1e-10–1e-9 J) on both. Test:
> `tests/test_transient_core_helical.py` (6 tests, mirrors
> `test_transient_core_shelltube.py`'s structure).
>
> **Implication for Stage D's own design**: an explicit CFL limit is real and
> currently invisible anywhere in `1Dmodel/` (confirmed by search — zero
> `cfl`/`courant` hits). Subcycling was the right minimal fix for legacy code
> being retired, but `core/residual.py`/`drivers/transient.py` should not
> inherit "no CFL awareness" as a default — either make the coolant advection
> genuinely implicit (consistent with §3.5's existing "linearly implicit per
> cell" direction for the wall/coolant coupling) or give the transient driver
> an explicit, config-independent CFL cap. Silent 50–200× under-resolution
> should not be reachable from any valid config.
>
> **Also found, deferred**: `coolprop_state_from_mass_energy` and the
> mass/energy update in `adapters_shelltube.py` hardcode `"Helium"` as the
> fluid string in several call sites (not read from `coolantProp.coolant`) —
> a fluid-agnosticity gap Stage D must not inherit, separate from the CFL fix
> and not touched here. `1Dmodel/simulink_coupling/_vendor/hps_combustor/transient_core/`
> holds byte-identical vendored copies of the pre-fix `transient_core/`
> modules — now diverged from the fixed originals. Simulink work is parked, so
> deliberately not touched; flagging so the divergence is a conscious debt,
> not a surprise later.
>
> Ready to start Stage D proper. Findings gathered while investigating this
> regression, not yet acted on:
>
> 1. **Factor-of-N trap, §3.1.** `transient_core.AxialGrid` stores TOTALS
>    (×`n_parallel`); `FlowPath` deliberately stores PER-CHANNEL. Critically,
>    `AxialGrid.coolant_volume` ≡ `FlowPath.volume_total`, **not**
>    `volume_per_channel`. Wiring `volume_per_channel` into
>    `initial_mass_energy_from_TP` would make every coolant inventory,
>    residence time, and thermal mass wrong by `N_ch` — and the symptom would
>    look exactly like the blow-up just fixed above. Shell-side is the reverse
>    trap: its area is already a single shared passage (`n_parallel=1`) while
>    the tube stream carries `n_parallel=N_tubes`.
> 2. **`FlowPath` is not a strict superset of `AxialGrid`, §3.1/§3.4.**
>    `wall_area`/`wall_volume` have no home in `core/mesh.py`, and
>    `core/wall.py`'s `CylindricalWall` stores no `A_wall`/`rho`/`cp` despite
>    its own docstring saying the driver integrates
>    `(rho*cp*A_wall) dT_bar/dt`. Add it as a **total** (×`n_parallel`),
>    matching how `WallCoolantStepInputs.wall_heat_capacity` is consumed today.
> 3. **§3.3's "thin closure adapter" undersells a real gap.** The
>    liquid/supercritical registry (`physics/liquid_flow/registry.py`) and the
>    gas-side `dispatch_nu_*`/`dispatch_friction_*` chains are entirely
>    separate mechanisms — different selection keys (inferred vs.
>    caller-forced string), different return types (h vs. Nu), different
>    failure semantics (`LookupError` + `ExtrapolationReport` vs.
>    `warnings.warn` + silent fallback). `ClosureContext` has **no channel for
>    `corrCoeffs`** — the concrete blocker, since the 21 `CorrelationCoefficients`
>    fields are load-bearing for `optimization/calibrate.py`. The adapter must
>    bridge **three** populations (hot gas, ideal-gas coolant, real-fluid
>    coolant), not two. Only the supercritical branch is actually
>    registry-driven today; subcritical branches are inline, so
>    `extrapolation_report` is `None` there — decide the adapter's contract for
>    that case rather than assuming every closure reports one.
> 4. **Two places a "clean" generic residual would silently NOT reproduce the
>    validated shell-and-tube transient, §3.5**: (a) `dq_cold`/`T_wc` are
>    computed by the wall solve and **discarded** — actual coolant heating is
>    `G·(Tbar_new − Tc)`, driven by the mean wall temperature, not the
>    reconstructed cold face; (b) `mdot_effective = mean(|faces|)`, a scalar,
>    is fed to per-cell closures rather than the local per-cell face flow.
>    Reproducing the legacy numbers requires copying these choices exactly;
>    "fixing" them changes the answer and is a separate decision, not a
>    migration detail.
> 5. **The 2×2-not-global-implicit lesson (§2, §3.5) is already being violated
>    in the code being migrated**: `semi_implicit_wall_compressible_coolant_step`
>    downgrades the proven fully-implicit 2×2 (`wall_coolant.py`, temperature-only)
>    to wall-implicit/coolant-explicit, because `T(m,U)` needs a CoolProp
>    inversion that can't be inlined into a 2×2 the way `T` alone can. That
>    explicit half is exactly what needed the CFL fix above.

> **Stage A closure note (2026-08-18):** all four legacy solvers
> (`main_solve.py`, `main_solve_shellntube.py`, `main_solve_transient.py`,
> `main_solve_shellntube_transient.py`) now route every gas-mode/GOX raw
> `CP.PropsSI` call through `IdealGasBackend` (`self._thermo`, instantiated
> once in each `__init__`) instead of calling `CoolProp.CoolProp.PropsSI`
> inline — 47 call sites converted 1:1, no expression rewriting. Added
> `IdealGasBackend.molar_mass(fluid)` (the one call shape, `PropsSI('MOLAR_MASS',
> fluid)`, with no prior getter) for `main_solve.py`'s coolant-molar-mass
> lookup. Liquid/two-phase/supercritical branches were already on the
> dispatch-routed path from the original Stage A extraction and were untouched.
>
> Gate result: the three steady solvers' baseline-regression tests
> (`tests/test_steady_baseline_regression.py`,
> `tests/test_steady_baseline_regression_finite_rate.py`) reproduce
> bit-identical obtained values before/after, including the same 4
> pre-existing unrelated frozen-chemistry failures (see memory
> `baseline-regression-failures-characterized` — not re-litigated here, and
> confirmed not to have drifted further). No pytest test previously exercised
> `main_solve_transient.py`/`main_solve_shellntube_transient.py` numerically
> (`test_main_transient_dispatch.py` mocks the solver classes out entirely) —
> added `tests/test_transient_baseline_regression.py` to close that gap and
> serve as this stage's gate for the two transient files; it will also be
> useful at Stage D's "reproduce the validated shell-and-tube transient" gate.
>
> Two findings worth carrying forward, both pre-existing and unrelated to this
> change:
> 1. `solve_transient_core` (the `transient_coolant` fluid_model path) raises a
>    pre-existing `FloatingPointError` ("coolant internal energy left the
>    configured CoolProp temperature range") on short/coarse smoke
>    configurations, for both helical and shell-and-tube, independent of
>    schedule choice. Not investigated further here (out of scope for a Stage A
>    rewire); the new baseline test routes around it by testing the
>    `quasi_steady` path end-to-end and the handful of `solve_transient_core`-
>    exclusive helper methods directly. Worth a dedicated look before Stage D
>    leans on `solve_transient_core` for real work.
> 2. `main_solve_shellntube_transient.py`'s `solve_transient_core` delegates
>    entirely to `transient_core.adapters_shelltube.run_shelltube_transient_core`
>    (a separate, already dispatch-routed module) — none of this file's own 9
>    raw call sites were ever reachable from that path, only from the legacy
>    `solve_transient()`. Asymmetric with the helical file, where
>    `solve_transient_core` still calls back into several of the class's own
>    methods. Relevant context for Stage D's residual-migration ordering.

> **Stage A progress note (2026-07-31):** `core/thermo.py` created —
> `ThermoState` (formerly `CoolantState`), `coolant_state_from_Tp`,
> `coolant_state_from_ph`, `coolant_inlet_state`, `ThermoBackend` protocol,
> `RealFluidBackend`, `IdealGasBackend`, `ReactingGasBackend` (thin adapter
> around the existing `transient_core/adapters_shelltube.py` gas-state
> provider pattern). `physics/liquid_flow/dispatch.py` now re-exports these
> as the SAME objects (pure relocation, not a reimplementation) — proven by
> `tests/test_core_thermo.py::test_dispatch_reexports_are_identical_objects`.
> Full test suite re-run after the move: 141 passed, only the same 4
> pre-existing baseline-regression failures (see
> [[baseline-regression-failures-characterized]] in memory), with **identical
> obtained values** to before the change — zero regression.
>
> **Not yet done** (the harder, riskier part of Stage A's gate): the four
> `main_solve*.py` files still call `CP.PropsSI` inline rather than through
> `IdealGasBackend`/`RealFluidBackend`. `IdealGasBackend`'s per-property
> getters were built to match each inline call site's exact granularity
> (proven equivalent in `test_core_thermo.py`) specifically so that rewiring
> is low-risk when it happens, but the actual solver-file edits are deferred
> to their own pass — touching a 1300-line coupled march deserves its own
> careful, fully-tested turn rather than being bundled in with the backend
> extraction.
>
> **Stage F groundwork started early (2026-07-31), standalone, ahead of
> Stages B-E** — justified because these pieces have zero dependency on
> `core/mesh.py`/`core/residual.py` (which don't exist yet):
> `core/geometry/nozzle_contour.py` (conical contour builder, throat+
> expansion-ratio parameterized) and `core/hotgas/nozzle_gas.py`
> (constant-gamma frozen quasi-1D area-Mach expansion, adiabatic wall
> temperature, and Bartz-Cornelisse film-property HTC — the form the user
> chose over the RPE stagnation-referenced variant). 20 unit tests, all
> passing, including exact closed-form isentropic checks per the Stage F
> acceptance item "Area-Mach and T_aw against closed-form isentropic
> relations". Exercised end-to-end by
> `1Dmodel/validation/nozzle_c2h4_o2_bartz_example.py` for the user's first
> real design point:
>
> | | |
> |---|---|
> | Propellants | C2H4 / O2, O/F = 2.3 (mass), 50 bar chamber |
> | Contour | conical, ER = 10, convergent/divergent half-angles and contraction ratio are ASSUMED defaults, not this engine's real geometry — see `nozzle_contour.py` |
> | Chamber (Cantera `gri30.yaml` equilibrium — caveat: tuned for natural gas, not purpose-built for C2H4/O2; fine for a first pass, not for final numbers) | T0 = 3760 K, gamma = 1.227 |
> | Throat (145.8 mm, see finding below) | M=1, T=3377 K, p=27.97 bar, h_g=11651 W/m²K, T_aw=3691 K, **q_w ≈ 33.7 MW/m²** at T_wall=800 K |
> | Exit | M=3.36, T=1651 K, p=0.58 bar, h_g=1126 W/m²K, q_w ≈ 2.9 MW/m² |
>
> Peak heat flux lands exactly at the throat, consistent with the RPE
> excerpt the user supplied ("the largest convective heat flux can be
> expected at the throat").
>
> **Geometry-consistency finding, same session:** the user's three numbers
> (D_throat=120mm, p0=50 bar, mdot~45 kg/s) are NOT mutually consistent —
> choking ties throat area, chamber pressure, and mass flow together
> (`mdot = G* · A_t`, `G*` fixed by chamber `T0`/`gamma`/`p0`). A 120mm
> throat at 50 bar with this chamber chemistry only chokes **~30.5 kg/s**,
> 68% of the stated 45 kg/s. Resolved by treating mass flow as authoritative
> (it is normally the thrust/Isp-driven quantity) and deriving throat
> diameter from it: **D_throat = 145.8 mm**, not 120 mm. New reusable
> functions `choked_mass_flux`/`throat_diameter_for_mass_flow` in
> `core/hotgas/nozzle_gas.py` do this check/conversion; the example script
> reports both the stated-vs-implied mismatch and the corrected geometry so
> nothing is silently overridden. Effect on peak `q_w`: small (~4% lower at
> the larger throat, `h_g ∝ D_t^-0.2` at fixed chamber conditions along the
> choked-flow family) — the resize is NOT what was making the number "seem
> high."
>
> **Wall-temperature sensitivity** (the user flagged q_w as high; this is
> the actual dominant lever right now, more than the geometry fix): swept
> the placeholder `T_wall_guess_K` 600–1400 K at the corrected throat —
> `q_w` ranges from 36.6 MW/m² (600 K) down to 25.6 MW/m² (1400 K), a ~40%
> swing, while `h_g` itself only moves ~5.5% (the film-property ratios
> partially self-compensate). **No wall/coolant coupling exists yet** — the
> reported `q_w` uses an ASSUMED uniform wall temperature, not a solved one;
> a real answer needs `core/wall.py` generalization (Stage C) and the
> coupled residual (Stage D/E) extended to this config, still gated behind
> Stages B-E per the original ordering. Even at the hot end of the sweep
> (1200–1400 K, closer to a real regen wall under high flux), q_w stays in
> the 25–28 MW/m² range — a rough cross-check against a Merlin-class engine
> (~97 bar, throat flux order 40–50 MW/m²) scaled by Bartz's `p0^0.8`
> dependence to 50 bar lands in a similar 23–28 MW/m² band, so the
> magnitude itself is not obviously wrong for this chamber pressure class;
> the honest uncertainty is the wall-temperature assumption and the gri30
> chamber chemistry, not (so far) a code defect.
>
> Side item, same session, user-requested: registered RPE (Sutton & Biblarz,
> *Rocket Propulsion Elements*) Eq. 8-24 — the classic Dittus-Boelter/Colburn
> liquid-film correlation — as a **cross-check-only** closure
> (`rpe_dittus_boelter_8_24` in the registry) for the subcooled-liquid regime.
> It is reported via new `CoolantClosureResult.cross_check_closure_name` /
> `cross_check_htc_W_m2_K` fields alongside the active Gnielinski HTC, never
> substituted for it. Measured divergence: the two agree within ~5% near
> Pr≈0.9, but Gnielinski runs up to **+56% higher** than Eq. 8-24 at Pr=6,
> Re=1e6 — worth remembering when the nozzle regen channel liquid-side film
> coefficient (Stage F) is implemented, since that is exactly this same
> physical regime.

This document is the executable design for the full solver rework. It is written
to be implemented by a model that has not seen the exploration session — every
stage names concrete files, concrete acceptance gates, and the existing code it
must preserve or absorb.

---

## 1. Scope

### 1.1 Configurations to support

| Config key | Hot side | Cold side | Status today |
|---|---|---|---|
| `shellntube` | combustion gas **inside** straight tubes | coolant in baffled shell cross-flow | steady + transient solvers exist |
| `shellnHelicalTube` | combustion gas in shell | coolant **inside** helical coil | steady + transient solvers exist |
| `nozzle_axial_channels` | nozzle gas (sub→supersonic) | coolant in axial milled channels | **NEW — nothing exists** |
| `nozzle_helical_channels` | nozzle gas (sub→supersonic) | coolant in N helically-wrapped channels | **NEW — nothing exists** |

All four × `flow_config ∈ {co, counter}`.

### 1.2 Fluids

Primary: **GHe**, **N2 / LN2** (incl. supercritical, 80 bar), **water / steam**
(incl. subcooled, boiling, superheated). Secondary already in the tree: O2/GOX
(pre-ignition chilldown), combustion products (Cantera/FPV).

Design rule, non-negotiable: **no fluid name appears in core code.** Fluid
identity enters only through `coolantProp.coolant` (a CoolProp name) and is
consumed by the thermo backend and the closure registry's `fluid_scope` tags.

### 1.3 Decisions taken (user, 2026-07-30)

1. **Staged migration onto a new `core/` package.** Legacy solvers stay live and
   are retired one at a time behind acceptance gates.
2. **Pluggable momentum.** One conservative mass+energy FV core; momentum is a
   swappable term, `quasi_steady` default, `low_mach` opt-in, interface shaped so
   full compressible drops in later.
3. **Nozzle hot side: internal quasi-1D solve + prescribed-table override.**
4. **Wall: 1D radial + analytic fin/rib efficiency.**

### 1.4 Explicit non-goals

- Film / transpiration cooling.
- 2D or 3D conduction (rib model is analytic).
- Combustion-chamber CFD; chamber state stays Cantera/FPV as today.
- Two-phase momentum beyond the existing homogeneous/drift closures.

---

## 2. What already exists and must be reused, not rewritten

Do **not** re-derive these. They are validated and are the reason the rework is
tractable.

| Asset | Location | Role in new core |
|---|---|---|
| Closure registry, validity/extrapolation reporting | `physics/liquid_flow/registry.py` | becomes `core`'s only closure-selection mechanism; widen tags only |
| Regime detection, `T_pc`, McEligot-Jackson HTD | `physics/liquid_flow/regime.py` | unchanged |
| Supercritical closures (McCarthy-Wolf, Taylor, Cheng2020, K-P, Wang2023) | `physics/liquid_flow/supercritical.py` | **directly reusable for nozzle regen** — McCarthy-Wolf/Taylor *are* rocket-regen correlations (Locke & Landrum 2008) |
| Boiling / CHF / ONB (Gungor-Winterton, MSH, Groeneveld LUT, Bergles-Rohsenow) | `physics/liquid_flow/correlations.py`, `chf.py` | unchanged |
| `(p,h)` real-fluid state + CoolProp cache | `physics/liquid_flow/coolprop_state_cache.py` | becomes the core thermo backend |
| Quasi-static quadratic wall reconstruction `fluxes_at_Tbar()` | `physics/heat_conduction.py` | generalized into `core/wall.py`; validated <2 K vs resolved PDE |
| Fin-efficiency helpers | `physics/heat_conduction.py:43-84` | already present, currently unused — the rib model's basis |
| FPV manifold + equilibrium manifold (tabulated, no per-node Cantera) | `physics/combustion_chemistry/` | unchanged; hot-gas provider interface |
| WSGGM tabulated radiation | `physics/radiation_model/` | unchanged |
| Bell-Delaware shell-side | `physics/bell_delaware.py` | unchanged |
| `AxialGrid`, fixed-step linearly-implicit integrator, schedules, diagnostics | `transient_core/` | **absorbed and generalized** into `core/` (see §4.1) |
| Result packaging | `result_package.py` | unchanged; core writes the same `data_master` contract |

**Hard-won numerical lessons that constrain the design** (from
`docs/context/TRANSIENT_STATUS.md` and the CLAUDE.md memory):

- Linearly-implicit fixed-step beats BDF for wall-film stiffness. BDF's Jacobian
  probing is prohibitively expensive when the RHS contains a profile relaxation.
- **No per-node Cantera calls in a march.** Manifold/tabulated only.
- Radiation must stay tabulated.
- Persistent Cantera objects must be reset from cached inlet `T/p/Y` before
  repeated sweeps.
- Shell-and-tube tube-side gas quantities are **per tube** (divide by `N_tubes`).
- Darcy friction convention throughout. No Fanning conversion.

---

## 3. Architecture

New package `1Dmodel/core/` → `hps_combustor.core`.

```
core/
  mesh.py              FlowPath, WallStack geometry containers, shared axial map
  thermo.py            fluid-agnostic property backends (the fluid-agnostic crux)
  closures.py          registry adapter: (stream, geometry, regime) -> h, f, CHF
  wall.py              radial resistance stack + analytic fin/rib + transient C
  state.py             conservative state vector packing, primitive recovery
  momentum.py          pluggable momentum terms
  residual.py          THE single spatial residual
  assembly.py          HXAssembly: streams + wall + coupling map + BCs
  geometry/
    shell_and_tube.py
    shell_and_helical_tube.py
    nozzle_contour.py          contour generation / import
    nozzle_axial_channels.py
    nozzle_helical_channels.py
  hotgas/
    combustor.py       chamber Cantera/FPV provider (wraps existing)
    nozzle_gas.py      NEW quasi-1D area-Mach expansion + Bartz-family HTC
    prescribed.py      table-driven h_g(x), T_aw(x), p(x) override
  drivers/
    march.py           steady space-march (co-flow IVP; shooting for counter)
    transient.py       time-accurate
    settle.py          settle-to-steady (canonical steady answer)
  diagnostics.py       energy closure, sanity gates, extrapolation roll-up
```

### 3.1 `mesh.py` — the geometry generalization

The single idea that makes all four configs one solver: **a stream is a 1D path
with its own arc-length coordinate `s`, plus a monotonic map to a shared axial
coordinate `z`.**

```python
@dataclass(frozen=True)
class FlowPath:
    name: str                    # "hot" | "cold"
    s_edges: np.ndarray          # [m] arc length along THIS stream's own path
    z_of_s: np.ndarray           # [m] shared axial coordinate at each s edge
    A_flow: np.ndarray           # [m^2] per-cell flow area (per channel)
    Dh: np.ndarray               # [m] hydraulic diameter
    P_wetted: np.ndarray         # [m] friction perimeter (per channel)
    P_heated: np.ndarray         # [m] heat-transfer perimeter (per channel)
    R_curv: np.ndarray | None    # [m] curvature radius (helix/bend); None = straight
    inclination: np.ndarray      # [rad] vs gravity — buoyancy/HTD criteria need it
    roughness: float
    n_parallel: int              # channels/tubes/coil starts in parallel
    geometry_tag: str            # registry geometry tag
    aspect_ratio: np.ndarray|None  # h_ch/w_ch for rectangular channels
    flow_direction: int          # +1 / -1 relative to increasing z
```

This subsumes every existing case:

| Config | hot path | cold path |
|---|---|---|
| shell-and-tube | `s = z`, `n_parallel = N_tubes`, straight | shell-side, `s = z`, cross-flow tag, Bell-Delaware `Dh` |
| shell-and-helical | shell-side, `s = z` | `s` = coil arc length, `z_of_s` = axial advance per turn, `R_curv = D_coil/2` |
| nozzle axial | `s = z` along contour, `A_flow = A_nozzle(z)` | `s ≈ z`, `n_parallel = N_channels`, rectangular |
| nozzle helical | `s = z` along contour | `s` = wrap arc length, `n_parallel = N_starts`, `R_curv` from local nozzle radius + helix angle |

`z_of_s` is what replaces the fragile `_advance_state()` axial bookkeeping in
`main_solve.py` — the bug the `HX_config == "shellnHelicalTube"` guard was added
to catch (CLAUDE.md, 2026-07-13). Here it is data, computed once by the geometry
builder, not re-derived in the march.

Coupling: `assembly.py` builds a **conservative interpolation operator** between
the hot and cold cell partitions of the shared `z` axis. For the helical case one
hot cell overlaps many coil cells; the operator distributes flux by overlap
length so energy is conserved exactly. Validate with an energy-closure unit test
before any physics rides on it.

### 3.2 `thermo.py` — the fluid-agnostic crux

One protocol, three backends:

```python
class ThermoBackend(Protocol):
    def state_ph(self, p: float, h: float) -> ThermoState: ...
    def state_pT(self, p: float, T: float) -> ThermoState: ...
    @property
    def p_crit(self) -> float: ...
```

`ThermoState`: `T, p, h, rho, cp, mu, k, a (sound speed), Pr, quality,
void_fraction, phase, is_supercritical, p_reduced, T_pc`. This is deliberately
`CoolantState` from `physics/liquid_flow/dispatch.py` — reuse that dataclass,
move it to `core/thermo.py`, leave a shim.

Backends:

- `RealFluidBackend` — CoolProp, `(p,h)` primary, backend selectable
  `HEOS|TTSE|BICUBIC` via `coolantProp.liquid_property_backend`. Wraps the
  existing `coolprop_state_cache`. **Covers GHe, N2/LN2, water/steam, O2 — every
  target fluid, single phase, two phase, and supercritical, with no branching in
  core code.**
- `ReactingGasBackend` — FPV / equilibrium manifold lookup for combustion
  products, keyed on `(h_removed, progress_variable, p)`. Tabulated only.
- `IdealGasBackend` — legacy fast path. Exists **solely** so the helium
  bit-identical regression can be reproduced during migration. Delete at Stage G
  unless it proves to be a needed speed lever.

`(p,h)` is the primary state pair everywhere. This is what makes boiling,
supercritical, and single-phase one code path — it is the recommendation
`water_coolant_conversion_plan.md` already reached, now made global.

### 3.3 `closures.py` — registry adapter

Thin. Builds a `ClosureContext` from a `FlowPath` cell + `ThermoState` +
lagged wall temperature, calls the registry, returns `h`, `f`, `CHF`, plus the
`ExtrapolationReport` and `htd_risk`. No physics here.

**Registry work required (small, additive):**

New geometry tags: `nozzle_channel_rect`, `nozzle_channel_helical`,
`annular_gap`. New closures to register:

| Closure | Fluid scope | Purpose | Source |
|---|---|---|---|
| Gnielinski + rectangular aspect-ratio correction | any | regen channel single phase | standard; cite in registry `provenance` |
| Curved-channel Dean enhancement (existing `mori1967`, `dispatch_nu_coil`) | any | helical wrap | already in `heat_transfer_correlations.py` |
| **Bartz** (+ σ correction) | combustion gas | nozzle hot side | see §5.2 |
| **Bartz local-property variant / Cinjarew** | combustion gas | cross-check | see §5.2 |

McCarthy-Wolf and Taylor are **already registered** and already tagged
`FLUID_ANY`, `1e4 < Re < 1e7`, with the `x_over_D` entrance term. They were fitted
on H2 rocket regenerative cooling. They are the correct default for the nozzle
coolant side with no new work beyond adding the geometry tags.

### 3.4 `wall.py` — radial stack + analytic rib

Generalizes `OneDimensionalSteadyConduction_ShellnHelicalTube` and
`fluxes_at_Tbar()`. Keeps:

- `hot_side ∈ {"inner", "outer"}` perimeter mapping — **do not regress this**,
  shell-and-tube needs `"inner"` (hot gas inside tubes).
- The quadratic quasi-static radial profile (`a2 = δ/2k`, `a6 = δ/6k`) and the
  closed-form 2×2 face solve. Validated <2 K vs a resolved PDE.
- Temperature-dependent `k_w` via `f_kw_at_T`.
- Tabulated-radiation fast path (`h_g_rad` supplied) vs. slow iterating path.

Adds, for `nozzle_*` and any ribbed/multi-channel geometry:

```
Effective cold-side heated perimeter per channel:
    P_c,eff = w_ch + 2 * eta_fin * h_ch          (+ closeout term if attached)
    m       = sqrt(2 * h_c / (k_w * t_rib))
    eta_fin = tanh(m * h_ch) / (m * h_ch)        (adiabatic tip)
```

Use `compute_fin_efficiency_ch` / `compute_eta_fin_rectangular` in
`heat_conduction.py:43,72` — already written, currently unused; check the factor-
of-2 convention against the rib being wetted on **both** faces before wiring
(`compute_fin_efficiency_ch` uses ×1, `compute_eta_fin_rectangular` uses ×2 —
for a rib between two channels the correct one is ×2).

Reported per node: `T_wg` (hot gas wall), `T_wc_root` (channel root), `T_rib_tip`,
`T_bar`. `T_wg` is the sizing quantity for a regen nozzle and is what the
1D-radial-only model would under-predict.

`WallStack` also carries per-cell `rho*cp*V` for the transient heat capacity, and
supports N material layers (base + optional coating/liner) as a series resistance
stack — needed if a thermal barrier or copper liner is ever specified.

### 3.5 `residual.py` — the single spatial residual

```
State per stream, per cell:
    quasi_steady momentum:  y = [ m_i = rho*A*dx ,  U_i = rho*e*A*dx ]
    low_mach momentum:      y = [ m_i , U_i , (rho*u*A)_face ]
Wall:                       y = [ Tbar_i ]  (per wall layer if multilayer)
```

Residual assembly, per stream:

1. Recover primitives: `v = m/V` → `rho`; `e = U/m` → `h = e + p/rho` (iterate
   once on `p` or carry `p` from momentum). `thermo.state_ph` → full state.
2. Upwind FV advective fluxes at faces from the signed mass flow.
3. Closures at each cell → `h_conv`, `f`, regime, extrapolation report.
4. Wall coupling → `Q_hot,i`, `Q_cold,i` via `wall.fluxes_at_Tbar`, mapped
   between stream partitions by the `assembly` overlap operator.
5. Momentum term (pluggable, §3.6) → `dp/ds`.
6. Sources: radiation (tabulated), chemistry (`ReactingGasBackend`).

The **same** residual is consumed by all three drivers. That is the whole point:
one place where physics lives.

Wall/coolant coupling stays **linearly implicit** per cell — reuse the 2×2 solve
pattern in `transient_core/wall_coolant.py`. This is the lesson that made the
shell-and-tube transient tractable; do not replace it with a global implicit
solve without measuring.

### 3.6 `momentum.py` — pluggable

```python
class MomentumModel(Protocol):
    def dpds(self, cell, state, closure) -> float: ...
    def extra_state(self) -> int:   # 0 for quasi-steady, n_faces for low-mach
```

- `QuasiSteadyMomentum` (default): `dp/ds = -(f/Dh)(rho u²/2) - rho u du/ds -
  rho g sinθ`. The acceleration term is **not** optional here — a regen channel
  with N2 going 100 K → 700 K has a large `du/ds`, and the existing helical solver
  already carries it (`dU__dx_IdealGas`). Generalize to real-fluid via
  `du/ds = -(u/rho)(drho/ds) + ...` from continuity, with `drho/ds` from the
  `(p,h)` EOS derivatives.
- `LowMachMomentum`: transient face momentum with acoustics filtered. Captures
  line inertance — the known open item for helical bang-bang valve events
  (`transient_core/README.md`).
- `CompressibleMomentum`: **not implemented now.** The protocol exists so it can
  be added without touching `residual.py`. Needed only if coolant-side choking
  becomes real.

Choking gate: keep the existing Mach check as a hard diagnostic on both streams,
using real-EOS `speed_sound` (already valid supercritically; Wood's two-phase
sound speed for the dome, already implemented).

### 3.7 `drivers/` — three drivers, one residual

1. **`march.py`** — steady space-march. Co-flow is an IVP: one pass, ~1.4 s
   today. Counter-flow needs shooting; reuse the working
   `solve_counterflow_liquid_reference` (bisection on hot-end enthalpy) and
   `solve_counterflow_physical_reference` (secant) logic from `main_solve.py`.
   Purpose: design sweeps, calibration loops, and **initializing the settle**.
2. **`transient.py`** — time-accurate. `fixed_step` linearly-implicit default;
   BDF/Radau retained as validation options with the sparse-Jacobian flag.
   Absorbs `transient_core/integrator.py` and `schedules.py` unchanged in spirit.
3. **`settle.py`** — transient integrator + steadiness criterion, march-
   initialized. **This is the canonical steady answer** for counter-flow, strong
   coupling, and final design verification.

Operating rule, unchanged from the Phase-2 outline: **screen with the march,
certify with the settle.** A disagreement between them is a finding (multiple
steady states, or a density-wave instability), not a bug.

---

## 4. What happens to the existing code

### 4.1 `transient_core/` is absorbed, not kept

`AxialGrid` → `FlowPath` (superset). `integrator.py`, `schedules.py`,
`diagnostics.py`, `wall_coolant.py` move into `core/` largely as-is.
`adapters_helical.py` and `adapters_shelltube.py` (2195 lines) are the parts that
**dissolve**: their job was bridging legacy geometry into a generic core, and in
the new design the geometry builders emit `FlowPath` directly. Expect the
adapters' physics content to end up in `core/geometry/*.py` and
`core/closures.py`, and the file count to drop sharply.

### 4.2 `physics/liquid_flow/` → `physics/fluid/`

The package name has been a documented misnomer since Phase 0 — it holds the
supercritical and single-phase closures too. Rename at Stage G with a
deprecation shim, same pattern as the existing
`physics/liquid_coolant.py` / `coolant_models.py` shims.

### 4.3 The four `main_solve*.py` files

Retired one at a time, each behind its own acceptance gate (§6). Each retirement
is a separate commit that deletes the file's private physics only after the
corresponding `core` path reproduces its validated numbers. `main_steady.py` and
`main_transient.py` keep their user-facing signatures and dispatch to `core`.

---

## 5. New physics: the nozzle configurations

Everything in this section is greenfield. There is no nozzle contour, no
supersonic gas model, and no Bartz correlation anywhere in the repo today
(verified: the existing "nozzle" hits are shell-and-tube inlet/outlet pipe
connections and `combustorProp.exhaust_diameter`).

### 5.1 `geometry/nozzle_contour.py`

Contour `r(z)` from either:
- **Import**: user CSV `(z, r)` table, monotonic in `z`, spline-smoothed. The
  path most users will take (from RPA/CEA/CAD).
- **Parametric**: conical (`theta_conv`, `theta_div`), or Rao/80 % bell
  (parabolic approximation with `theta_n`, `theta_e`).

Derived per cell: `A(z)`, `A/A_t`, wall curvature radius (needed by Bartz's
`(D_t/R_c)^0.1` throat-curvature term), and the throat index.

Sanity gates: monotonic converging section, single throat, `A/A_t ≥ 1`
everywhere, `dr/dz` continuous.

### 5.2 `hotgas/nozzle_gas.py`

**Expansion.** Given the chamber state from the existing Cantera/FPV path
(`p_c`, `T_c`, composition), solve the quasi-1D area-Mach relation for `M(z)`:
subsonic branch upstream of the throat, supersonic downstream, sonic at `A_t`.
Two chemistry modes:

- `frozen` — `gamma` fixed at chamber composition. Cheap, the validation default.
- `equilibrium` — shifting equilibrium, **tabulated** as a 1D table vs `p/p_c`
  built once per chamber state and cached alongside the FPV manifolds under
  `cache/`. Per-node Cantera calls are forbidden (§2).

**Adiabatic wall temperature.** `T_aw = T_0 * (1 + r*(γ-1)/2*M²) / (1 +
(γ-1)/2*M²)`, recovery factor `r = Pr^(1/3)` turbulent, `Pr^(1/2)` laminar.
`T_aw`, not `T_0` and not `T_static`, is the driving temperature for the hot-side
flux. **This is the most common modelling error in regen cooling — get it right
and assert it in a unit test.**

**Hot-side HTC.** Bartz. **Two forms of the correlation exist in the
literature and are NOT interchangeable line for line** — confirmed 2026-07-31
from two independent sources supplied by the user (rendered page images,
visually verified per the discipline below), plus a third source for the
liquid-side check:

*Form 1 — local free-stream properties + film-property ratio* (Cornelisse,
*Rocket Propulsion and Spaceflight Dynamics*, 1979, Eq. 8.3-3):

```
h_c = 0.026 * (mu^0.2 * cp / Pr^0.6) * (rho*V)^0.8 / D^0.2 * (rho_f/rho) * (mu_f/mu)
```

evaluated at LOCAL free-stream conditions (`rho`, `V`, `D`, `mu`, `Pr`, `cp` at
that axial station); `f` subscript = film temperature (arithmetic mean of
wall and free-stream static temperature).

*Form 2 — same core group, different correction term* (Sutton & Biblarz,
*Rocket Propulsion Elements*, latest ed., Eq. 8-22):

```
h_g = (0.026 / D^0.2) * (cp * mu^0.2 / Pr^0.6) * (rho*v)^0.8
      * (rho_am / rho') * (mu_am / mu_0)^0.2
```

`am` = arithmetic mean of local free-stream static T and wall T (same film
concept as Cornelisse's `f`); `rho'` = local free-stream density; `mu_0` =
**stagnation/chamber** viscosity (NOT local free-stream, unlike Cornelisse's
`mu`).

**The core group is now confirmed identical across both sources**: coefficient
`0.026`, `mu^0.2`, `cp` (linear), `Pr^-0.6`, mass-flux exponent `0.8` on
`(rho*V)`, `D^-0.2`. This is strong independent cross-validation of those six
numbers specifically.

**The correction term is a real, unresolved discrepancy, not a relabeling**:
Cornelisse's `(mu_f/mu)` is linear and referenced to local free-stream
viscosity; RPE's `(mu_am/mu_0)^0.2` is raised to the 0.2 power and referenced
to stagnation viscosity. Both density-ratio terms are linear and structurally
similar (`rho_f/rho` vs. `rho_am/rho'`) but use slightly different reference
points (film-vs-local vs. arithmetic-mean-vs-local). **Do not silently merge
these into one hybrid form.** When Stage F starts:

1. Register BOTH as separate closures (`bartz_cornelisse_local` and
   `bartz_rpe_8_22` or similar), same pattern as McCarthy-Wolf vs. Taylor for
   supercritical N2 — let the registry's extrapolation/comparison machinery
   make the difference visible rather than picking one by assumption.
2. The more commonly-cited "chamber-stagnation + `c*` + area-ratio `(A_t/A)^0.9`
   + throat-curvature `(D_t/R_curv)^0.1` + sigma correction" parameterization
   (the form originally sketched in this doc's first draft, before the two
   images above were verified) is a further reformulation built from this same
   core group via continuity/isentropic substitution — it needs its OWN
   independent source verification before being hardcoded; do not assume it
   from memory even though it is the most commonly quoted "the Bartz
   equation." Ask the user for that page specifically if/when it's needed.
3. `sigma`-based or film-property-based, the correction depends on `T_wg`
   either way, so this is a fixed point with the wall solve regardless of
   which form is used — use the same **lagged wall temperature** pattern
   already used for the boiling `Bo` term and the supercritical property-ratio
   closures (`_shell_Tw_lagged`). Do not add a new iteration scheme.

Tier both Bartz forms `structural_extrapolation` outside their fitted
envelope; neither source page gave an explicit validity range in the excerpt
seen so far — check for one when the full page/chapter is available.

> **VERIFY-BEFORE-HARDCODE, blocking — status: PARTIALLY DONE.** Phase 1
> established that PDF text extraction in this repo **drops minus signs** —
> caught in Cheng2020, Krasnoshchekov-Protopopov, *and* the McEligot-Jackson
> formula. The core Bartz group above is now verified from two rendered-page
> images (not text extraction) and cross-checked against each other. Still
> open before any Bartz closure is registered: (a) exact numeric value/units
> for the RPE viscosity-ratio exponent `0.2` reconfirmed directly from the
> image (recorded here from the same conversation, but re-verify against the
> saved image file before hardcoding, not from this doc's transcription of
> it); (b) a source for the chamber-stagnation/area-ratio/sigma
> reformulation if that variant is wanted instead of/alongside the two above.

**Prescribed override** (`hotgas/prescribed.py`): user table of `h_g(z)`,
`T_aw(z)`, `p_g(z)` from CEA/RPA, consumed through the identical provider
interface. Gives a cheap cross-check against the internal solve and unblocks
users who already have RPA output.

### 5.3 `geometry/nozzle_axial_channels.py`

`N_channels` rectangular channels, wall thickness `t_w` (hot wall), rib thickness
`t_rib`, channel height `h_ch(z)` and width `w_ch(z)` — both allowed to vary
axially (channels are typically narrowed at the throat). Per-channel
`mdot = mdot_total / N_channels`.

`Dh = 4*A/P`. Aspect ratio `AR = h_ch/w_ch` drives both the Nusselt correction
and the fin model. High-AR channels (AR > 4) are standard practice and are
exactly where the 1D-radial-only wall model fails.

### 5.4 `geometry/nozzle_helical_channels.py`

`N_starts` channels wrapped at helix angle `alpha(z)`. Arc length relation
`ds = dz / sin(alpha)`; local curvature radius from the nozzle radius `r(z)` and
`alpha`. Reuse `dispatch_nu_coil` / `dispatch_friction_coil` Dean-number
enhancement (`heat_transfer_correlations.py:431`,
`friction_correlations.py:158`) — already fluid-agnostic.

Single-channel (`N_starts = 1`) is just the `n_parallel = 1` case; no separate
code path.

### 5.5 Known gaps to carry openly

Same honesty discipline as the existing CHF/ONB/HTD warnings — flag, do not
fabricate:

- **Curvature-induced HTC asymmetry near the throat.** Real regen channels see
  concave/convex wall differences and secondary flows in the throat region. Not
  modelled. Flag when `|dr/dz|` is large.
- **Supercritical + strong curvature + strong heating simultaneously.** No
  validated closure. McCarthy-Wolf/Taylor applied with the curvature enhancement
  multiplied on is a `structural_extrapolation`, and must report as such.
- **HTD magnitude.** Detection only, unchanged from Phase 1. No degraded-HTC
  magnitude model exists for N2 at these conditions.
- **Cryogenic material data.** `ST316L` conductivity/cp/strength tables clamp
  flat below **27 °C** (`material_temperature_strength.py:34,58`). LN2 at 100 K
  is 173 °C below the floor. `INCO718` has data to −240 °C and is usable.
  Regen nozzles usually want **CuCrZr / NARloy-Z**, which is not in the tree at
  all. Action: warn loudly on clamp (do not silently extrapolate), and treat
  acquiring CuCrZr + cryogenic 316L data as a prerequisite for trusting any
  LN2-cooled nozzle wall temperature.
- **Shell-side cross-flow supercritical.** Existing documented gap, unchanged.

---

## 6. Staged implementation

Every stage is independently shippable and has a hard acceptance gate. **Do not
start stage N+1 until stage N's gate passes.**

The three bit-identical reference values from Phase 0/1, to be re-asserted at
every stage that touches the existing configs:

```
helical helium         182.14810060417753
helical water          184.0092196238135
shell-and-tube water   323.60702108508826
```

Plus `tests/test_steady_baseline_regression*.py`. Note: 4 baseline-regression
failures pre-exist and are unrelated to this work (see the simulink-coupling
session) — characterize them **before** Stage A so they are not mistaken for
rework regressions.

### Stage A — `core/thermo.py`
Extract `CoolantState` and the CoolProp cache into the backend protocol. Point
all four legacy solvers at it.
**Gate:** all three reference values bit-identical. No behavior change.

### Stage B — `core/mesh.py` + geometry builders for the two existing configs
`FlowPath`, `z_of_s`, the conservative overlap operator, and builders for
`shellntube` and `shellnHelicalTube`. Legacy solvers consume geometry from the
builders instead of computing it inline.
**Gate:** bit-identical, **plus** a unit test proving the overlap operator
conserves energy exactly on a non-uniform helical↔shell mapping.

### Stage C — `core/wall.py`
Generalize the radial stack; add the analytic rib/fin path (unused by the two
existing configs, so it is dead code until Stage F — test it standalone against
a hand calculation and against a 2D FEM/analytic reference for one high-AR
channel).
**Gate:** bit-identical for existing configs; rib model verified independently.

### Stage D — `core/residual.py` + `drivers/transient.py`
Migrate the **shell-and-tube transient** first — newest, best documented, has
`DESIGN_PLAN_shellntube_transient.md`.
**Gate:** reproduce the validated shell-and-tube transient within stated
tolerance; energy closure to machine precision; settle-to-steady agrees with the
steady solver.

### Stage E — `drivers/settle.py` + `drivers/march.py`; migrate the steady solvers
Shell-and-tube steady, then helical steady (helical last — it is the largest and
carries the counter-flow shooting logic).
**Gate:** all three reference values reproduced within tolerance by the march
driver; the settle driver agrees with the march on the co-flow cases; counter-flow
settle agrees with `solve_counterflow_liquid_reference`.

### Stage F — nozzle configurations
Contour, `nozzle_gas.py`, both channel geometries. **No legacy to match** — this
is where validation must come from outside the repo:
1. Bartz reproduced against a published worked example (verify coefficients from
   a rendered page first, §5.2).
2. Area-Mach and `T_aw` against closed-form isentropic relations — exact, so
   assert tightly.
3. Whole-nozzle heat flux and coolant temperature rise cross-checked against an
   RPA/CEA case via the prescribed-table path.
4. Energy closure on the coupled solve to machine precision.

### Stage G — cleanup
Delete legacy physics from the four `main_solve*.py`. Rename
`physics/liquid_flow` → `physics/fluid` with a shim. Delete `IdealGasBackend` if
unused. Single closure-provenance table in `docs/`. Update `CODEBASE_MAP.md`,
`CLAUDE.md`, `AGENTS.md`, `TESTING_CONTEXT.md`.

---

## 7. Risks

| Risk | Mitigation |
|---|---|
| Bit-identical acceptance breaks on a floating-point reassociation during extraction | Extract *without* rewriting expressions; if a value moves, bisect to the expression and restore associativity before proceeding |
| Overlap operator leaks energy on helical↔shell mapping | Dedicated conservation unit test at Stage B, before any physics rides on it |
| Nozzle gas + wall + coolant fixed point converges poorly at the throat | Lagged-`T_wg` pattern (already proven for boiling `Bo` and supercritical property ratios); under-relax; if it oscillates, that is the settle driver's job to reveal |
| Bartz coefficient transcription error | Blocking visual verification from a rendered page (§5.2). This has bitten three times already in this repo |
| Cryogenic/copper material data absent | Warn on table clamp; treat as an explicit prerequisite, not a silent extrapolation |
| Scope creep into full compressible NS | `CompressibleMomentum` is a protocol slot only. Do not implement in this rework |
| Losing calibrated correlation knobs | `CorrelationCoefficients` field names are load-bearing for the optimizer — do not rename during migration |

---

## 8. Immediate next actions

1. Characterize the 4 pre-existing `test_steady_baseline_regression*` failures so
   Stage A has a clean baseline.
2. Stage A: `core/thermo.py`.
3. In parallel (independent of A–E): acquire and verify the Bartz source page,
   and decide on CuCrZr / cryogenic-316L material data.

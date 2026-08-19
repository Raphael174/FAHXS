# transient_core/ LLM Context

## Scope

Next-generation, geometry-independent transient solver core. Provides a
finite-volume wall+coolant time-stepping kernel that helical and
shell-and-tube "adapter" modules bridge into legacy solver geometry/physics.
This is the maintained path for `transientProp.fluid_model =
"transient_coolant"` production dispatch (`main_solve_transient.py`,
`main_solve_shellntube_transient.py`), and is also the intended landing site
for most of the upcoming fluid-agnostic FV core rework — see
`docs/solver_design/FV_CORE_REWORK_PLAN.md` and `docs/context/TRANSIENT_STATUS.md`.

## Contents
| File | Role |
|---|---|
| `grid.py` | `AxialGrid` — geometry-neutral cell edges/centers/lengths/coolant+wall volumes/perimeters, inlet/outlet indices for either flow direction. Stores TOTALS across all parallel coils/tubes (`n_parallel`-scaled), unlike `core/mesh.py`'s `FlowPath` which stores PER-CHANNEL — a documented "factor-of-N trap" between the two (see FV_CORE_REWORK_PLAN.md's 2026-08-18 note). |
| `state.py` | `TransientState`/`TransientStateLayout` — named layout `[Tbar_wall[0:N], T_coolant[0:N]]` for the first wall+coolant state vector. |
| `coolant_fv.py` | Geometry-independent implicit-upwind FV update for coolant TEMPERATURE only (`rho*V*cp*dT/dt = mdot*cp*(T_up-T) + Q`). The earliest/simplest coolant layer; superseded for production by the mass/energy path below. |
| `wall_coolant.py` | First coupled wall/coolant step: per-cell 2x2 fully-implicit linear system (temperature-only coolant state). Proven/stable — referenced as the model for what "genuinely implicit" coolant advection should look like. |
| `compressible_coolant.py` | Conservative mass/energy coolant primitives: `dm/dt = mdot_in - mdot_out`, `dU/dt = mdot*h_upwind + Q`, reconstructed to `p,T,rho,h,cp,mu,k` via CoolProp (`coolprop_state_from_mass_energy`). **`conservative_mass_energy_step` is explicit forward-Euler and unconditionally unstable past roughly one cell residence time — see Sharp Edges below.** |
| `wall_compressible_coolant.py` | `semi_implicit_wall_compressible_coolant_step` — wall-implicit/coolant-explicit coupling (downgraded from `wall_coolant.py`'s fully-implicit 2x2 because inverting `T(m,U)` needs a CoolProp call that can't be inlined into a small linear solve). The explicit half here is exactly what needed the CFL fix. |
| `integrator.py` | Generic fixed-step driver: `integrate_wall_coolant_fixed_step`, `fixed_time_grid` (inserts schedule breakpoints + requested output times into a bounded-step grid). Geometry-neutral; adapters supply one `WallCoolantStepInputs` per interval. |
| `schedules.py` | Shared schedule helpers (`interp_schedule`, `schedule_times`, `collect_transient_schedule_times`) — linear interpolation, flat-held outside range, matching legacy transient solver behavior. |
| `diagnostics.py` | `energy_audit`, `residence_time_s`, `wall_time_constant_s`, `timescale_audit` — shared validation helpers (residence time = `sum(rho*V)/|mdot|`, energy residual normalized by relevant energy terms). |
| `progress.py` | `TransientProgressPrinter` — step/time/material min-max temperature progress reporting during long runs. |
| `adapters_helical.py` | Bridges legacy helical geometry/correlations into `AxialGrid` + the wall/coolant step. Mirrors the legacy friction/Nusselt dispatcher calls but takes thermophysical properties as arrays instead of calling CoolProp internally. |
| `adapters_shelltube.py` | Bridges EchTherm shell/tube geometry, Bell-Delaware shell film, hot-gas march (FPV/equilibrium/oxygen gas-state providers), and wall flux (`hot_side="inner"`) into the generic integrator. ~2200 lines; owns the CFL-subcycling fix (see below) and `run_shelltube_transient_core()`, the production entry point for shell-and-tube `transient_coolant` dispatch. |

## `__init__.py` public surface (two-tier)

`transient_core/__init__.py` re-exports essentially everything from every
module in the list above as the package's flat public API (`AxialGrid`,
`TransientState`, `CompressibleCoolantStepResult`, `run_shelltube_transient_core`,
`build_helical_core_geometry`, `energy_audit`, etc. — see the file for the
full `__all__`). Per `docs/solver_design/FV_CORE_REWORK_PLAN.md`, this surface
splits into two forward-looking categories:

- **Pure infrastructure — planned to move into `1Dmodel/core/` largely as-is**:
  `integrator.py`, `wall_coolant.py`, `wall_compressible_coolant.py`,
  `compressible_coolant.py`, `schedules.py`, `diagnostics.py`, `grid.py`
  (`AxialGrid`), `state.py`, `progress.py`.
- **Adapters — planned to dissolve** once FV-core-rework Stages D/E land:
  `adapters_shelltube.py` and `adapters_helical.py`. Their job (translate
  legacy geometry/correlation calls into the generic step-input shape) is
  meant to be replaced by `core/geometry/`'s `FlowPath`/`HXAssembly` builders
  feeding the generic core directly, without a bridging-adapter layer.

Treat this as forward-looking, not yet-done: as of 2026-08-18 the adapters
are still the only way production code reaches this package.

## Sharp edges

- **`conservative_mass_energy_step` is explicit forward-Euler in the
  conserved variables and unconditionally unstable once a step advances a
  cell by more than roughly its residence time (`mass/mdot`).** This bit
  every case in the package's own documented validation matrix
  (`1Dmodel/validation/transient_core_short_runs.py`), which had been running
  at `max_step`/`tau` ~ 50-200, crashing with `FloatingPointError` from
  `enforce_internal_energy_bounds`. Confirmed empirically 2026-08-18: stable
  at `dt/tau <= ~0.2`, a fast-growing single-cell spike at `dt/tau ~0.4`,
  crash within a handful of steps above `dt/tau ~1` — the textbook explicit-
  upwind CFL signature (traced to an amplifying spike originating at the
  inlet cell and traveling downstream step by step). Reproduced identically
  across a 50x range of `max_step` (0.25 s down to 0.005 s), ruling out "just
  use a smaller documented step" as a fix.
  - The archived pre-fix validation numbers
    (`docs/validation/transient_core_short_run_results_PRE_CFL_FIX_2026-07-10.json`,
    gitignored, on-disk but may not surface via normal search) had an
    `energy_residual_abs_max` of 30-93 J per case — this was the unresolved
    CFL error itself, not real energy non-closure, and hadn't yet pushed a
    cell property past the also-independently-wrong 2500 K guard bound
    (CoolProp's real ceiling for He/N2/Water is 2000 K, so the guard was
    silently extrapolating rather than protecting).
  - **Fix (2026-08-18, `adapters_shelltube.py::_cfl_stable_substep_count`)**:
    subdivides a macro step into CFL-safe substeps (`safety=0.25` of
    `min(cell mass)/max(|face mdot|)`) for the mass/energy advection only —
    `faces`/`hot_heat_W`/conductance stay frozen across substeps, matching
    the quasi-steady-per-macro-step assumption already used for the hot-gas
    march. Applied **only in `adapters_shelltube.py`**, and only wraps the
    explicit half of `semi_implicit_wall_compressible_coolant_step`. Test:
    `tests/test_transient_core_shelltube.py`
    (`test_grooved_re_thresholds_come_from_corrcoeffs_not_function_defaults`
    is a related but distinct regression test from the same investigation).
    Regenerated validation matrix: energy residual O(10-90) J -> O(1e-8) J;
    headline `T_c_out_peak` moves +0.5% (bang-bang cases) to -7% (GOX cases,
    because lower flow there meant the pre-fix violation was worse). Real
    performance cost: 45-140 s/case now vs <1 s/case before, because CFL
    resolution now honestly runs a step that used to silently skip most of
    the physics.
  - **Also fixed, 2026-08-19**: the separate helical transient-core path in
    `main_solve_transient.py` had the same explicit forward-Euler crash
    signature (confirmed identical root cause, not just similar) — fixed
    with the same `_cfl_stable_substep_count` helper, imported from this
    package into `main_solve_transient.py` alongside the other private
    helpers it already borrows from `adapters_shelltube.py`. Verified for
    both `quasi_steady` and `low_mach` momentum models. Test:
    `tests/test_transient_core_helical.py`.
  - Zero `cfl`/`courant` string hits exist anywhere else in `1Dmodel/`
    (confirmed by search during the investigation) — CFL awareness is
    currently invisible outside this one adapter function. The design doc
    explicitly warns the future `core/residual.py`/`drivers/transient.py`
    should not inherit "no CFL awareness" as a default: either make coolant
    advection genuinely implicit (consistent with `wall_coolant.py`'s
    existing fully-implicit 2x2 pattern) or give the transient driver an
    explicit, config-independent CFL cap.
- **Vendored copy has diverged.** `1Dmodel/simulink_coupling/_vendor/hps_combustor/transient_core/`
  holds a vendored copy of parts of this package for a separate,
  currently-parked FMU-export feature. As of 2026-08-19 that vendored copy
  is byte-identical to the PRE-fix originals — the CFL fix was deliberately
  not applied there. Do not assume the vendor copy and this package are in
  sync; if the Simulink work resumes, the CFL fix needs to be ported over.
- `coolprop_state_from_mass_energy` and the mass/energy update in
  `adapters_shelltube.py` hardcode `"Helium"` as the CoolProp fluid string in
  several call sites rather than reading `coolantProp.coolant` — a known
  fluid-agnosticity gap, found but deliberately not touched during the CFL
  fix. Any future Stage D work must not inherit this.
- `AxialGrid` (`grid.py`) stores per-`n_parallel` TOTALS, while `core/mesh.py`'s
  `FlowPath` stores PER-CHANNEL quantities — mixing the two conventions
  without an explicit conversion is the "factor-of-N trap" flagged in
  FV_CORE_REWORK_PLAN.md; relevant to anyone wiring `core/geometry/` builders
  into this package's integrator.
- Shell-and-tube wall-flux bridging in `adapters_shelltube.py` uses
  `hot_side="inner"` (hot gas inside tubes) — do not regress this convention;
  it mirrors the repo-wide steady-solver convention in
  `physics/heat_conduction.py` (see `/CLAUDE.md`).
- Two legacy shell-and-tube quirks were found (not yet decided whether to fix
  or keep) during the CFL investigation: `dq_cold`/`T_wc` are discarded in
  favor of `G*(Tbar-Tc)`, and `mdot_effective` is a mean-of-faces scalar fed
  into per-cell closures. Default stance per the design doc is
  reproduce-then-flag, not silently fix, since the current stage gate is
  "reproduce within stated tolerance."

## TODOs

- Explicitly stated in `README.md`: "Next layers: helical transient-coolant
  hot-gas wrapper and dispatch, dashboard field polish, broader short-run
  validation cases, and shell-side residence-volume refinement."
- `README.md` also flags the shell-side coolant hold-up model (shell
  cylinder volume minus displaced tube volume, uniformly distributed) as "a
  clear baseline" pending "baffle window/leakage path refinements... in a
  later residence-volume refinement."
- Helical momentum: README notes bang-bang momentum audits showed "real
  helical line inertance can be large for aggressive valve events" and that
  transient momentum should be revisited before treating high-frequency
  helical bang-bang predictions as final.
- The `"Helium"`-hardcoding fluid-agnosticity gap and the unfixed helical
  transient-core CFL crash (both above) are explicit, evidenced open items
  from the 2026-08-18 investigation, not yet scheduled to a stage.

## Change history

Sparse dated evidence (single initial commit; this folder is currently
uncommitted/untracked work per `git status`). Dated notes recovered from
`docs/solver_design/FV_CORE_REWORK_PLAN.md` and `docs/context/TRANSIENT_STATUS.md`:

- 2026-07-09: `transient_core/` package first appears with the initial
  coolant/wall FV layer (per `TRANSIENT_STATUS.md`).
- 2026-08-18: `_cfl_stable_substep_count` added to `adapters_shelltube.py`,
  fixing the forward-Euler mass/energy-advection instability described above
  under Sharp Edges. Same session found (but did not fix) the helical
  transient-core equivalent crash and the `"Helium"`-hardcoding gap, and
  confirmed the Simulink vendor copy has diverged.

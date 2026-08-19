# physics/ LLM Context

## Scope

Top-level physics package for `hps_combustor`. Holds the wall-conduction solver, the
helical-coil and shell-and-tube Nusselt/friction correlation dispatchers, the
Bell-Delaware shell-side method, and a standalone prescribed-temperature pressure-drop
study. The real fluid/chemistry/radiation physics live in the sibling subpackages
(`combustion_chemistry/`, `gas_flow/`, `liquid_flow/`, `radiation_model/` — each has its
own `LLM_CONTEXT.md`); several files directly in this folder are deprecation shims left
over from the Phase 1 restructure described in
`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` and now just re-export from the
subpackages with a `DeprecationWarning`.

## Contents

| File | Role |
|---|---|
| `heat_conduction.py` | `OneDimensionalSteadyConduction_ShellnHelicalTube` — the maintained 1D radial wall-conduction solver (steady `Solve1Dconduction()` + transient `fluxes_at_Tbar()` companion), shared by the helical and shell-and-tube solvers. Also holds unused-but-real fin-efficiency helpers (`compute_fin_efficiency_ch`, `compute_eta_fin_rectangular`, `compute_eta_fin_triangular`) flagged by `docs/solver_design/FV_CORE_REWORK_PLAN.md` as the basis for the future nozzle rib model. |
| `heat_transfer_correlations.py` | Nusselt correlations and the `dispatch_nu_coil` / `dispatch_nu_shell` / `dispatch_nu_tube_straight` selectors; also `nu_corrugated_tube_vicente` for grooved tubes and several standalone correlations (Salimpour2008, Ahmed1997 toroid, Mori-Nakayama curved tube, Churchill-Bernstein) not all of which are wired into a dispatcher. |
| `friction_correlations.py` | Darcy friction correlations and the `dispatch_friction_coil` / `dispatch_friction_tube_straight` selectors, plus `friction_corrugated_tube_vicente` for grooved tubes. |
| `bell_delaware.py` | Bell-Delaware shell-side h/dp method for the baffled shell-and-tube config (`bell_delaware_shell()`), with the ideal-tube-bank j/f coefficient tables and J_c/J_l/J_b/J_s/J_r correction factors from Serth (2007). |
| `pressure_drop_prescribedT.py` | Standalone helium-side compressible pressure-drop march under a prescribed linear coil temperature profile — decouples pressure-drop confidence from the unvalidated hot-gas Nusselt number. Reads `dispatch_friction_coil` only; does not touch the coupled solver. |
| `Helium_thermodynamics.py` | Standalone exploratory Jupyter-style script plotting He compressibility factor Z vs (T,p). Confirmed no importers (audited 2026-07-13); not part of any solver path. |
| `coolant_models.py`, `liquid_coolant.py`, `heated_liquid_channel.py`, `liquid_hx_adapters.py`, `governing_equations.py` | Deprecation shims (Phase 1 restructure). Each re-exports the real implementation from `physics/liquid_flow/*` (or `physics/gas_flow/governing_equations.py` for the ideal-gas one) and raises `DeprecationWarning` on import. Do not add new code here — import the target module directly. |

## Key correctness points

- **Darcy convention.** All friction dispatchers in this folder return Darcy-Weisbach
  factors; the momentum closures in `physics/gas_flow/governing_equations.py` and the
  hot-gas `dp` form assume Darcy. Do not apply a Fanning→Darcy ×4 correction anywhere in
  this package (see `dispatch_friction_coil`'s own docstring, and CLAUDE.md's "Known
  Sharp Edges").
- **`hot_side` orientation in `heat_conduction.py`.** `OneDimensionalSteadyConduction_ShellnHelicalTube`
  takes `hot_side="outer"` (helical coil: hot gas in the shell, coolant inside the tube —
  the default) or `hot_side="inner"` (shell-and-tube: hot combustion gas inside the
  tubes). This flag swaps which perimeter (`P_inner` vs `P_outer`) is used for the hot vs
  cold flux in both `Solve1Dconduction()` and `fluxes_at_Tbar()`. Shell-and-tube callers
  MUST pass `hot_side="inner"`; regressing this silently swaps which side sees which
  perimeter.
- **Per-tube gas quantities.** Anything in this folder called from the shell-and-tube
  solver (tube-side Nusselt/friction, Bell-Delaware shell-side h) operates on per-tube
  hot-gas quantities; the caller is responsible for dividing total hot-gas mass flow by
  `N_tubes` before calling in.
- **`inside_tube_choice="grooved"` is real physics**, not a smooth-tube placeholder: it
  routes through `nu_corrugated_tube_vicente()` / `friction_corrugated_tube_vicente()`
  (Vicente/Cruz-style corrugated-tube correlations) with `tube_grooved_Nu_factor` /
  `tube_grooved_f_factor` as calibration multipliers on top of the literature form.
- **Grooved Re thresholds must come from `CorrelationCoefficients`, not the grooved
  functions' own defaults.** `nu_corrugated_tube_vicente`/`friction_corrugated_tube_vicente`
  default to `Re_lo=2000/Re_hi=4000`, but the legacy shell-and-tube adapter (and the new
  `physics/liquid_flow/gas_closures.py` wrappers) always pass `Re_transition_lo`/`_hi`
  from `corrCoeffs` explicitly (usually 2300/4000) — using the functions' own defaults
  silently shifts the laminar/turbulent blend band (caught while writing
  `tests/test_core_closures.py`, see `docs/solver_design/FV_CORE_REWORK_PLAN.md`'s
  2026-08-18 note).
- **`dispatch_nu_shell` uses a different characteristic length per selector** (`Dh_cc`
  for `salimpour2008`, `D_tube_outer` for both Churchill-Bernstein variants, an `Asqrt`
  toroid length for the Ahmed1997 fallback) — the dispatcher handles this internally so
  callers always get the matching `h_g`; do not reuse a cached `Dh` across selectors.
- Radiation coupling in `heat_conduction.py` is a parallel-path addition to `h_g`
  (`h_g_eff = h_g + h_g_rad`), computed via `radiation_model/radiation_equations.py`'s
  `qrad_net_mbl`/`hrad_from_q` — see `physics/radiation_model/LLM_CONTEXT.md`.

## Test coverage

No dedicated unit test file targets `heat_conduction.py`, `bell_delaware.py`,
`heat_transfer_correlations.py`, or `friction_correlations.py` directly by name; they are
exercised indirectly through `tests/test_steady_baseline_regression*.py` and
`tests/test_core_closures.py` (which specifically asserts bit-identical equivalence
between the new `gas_closures.py` registry wrappers and these dispatch functions — see
`physics/liquid_flow/LLM_CONTEXT.md`). `pressure_drop_prescribedT.py` and
`Helium_thermodynamics.py` have no test coverage and are standalone/exploratory.

## TODOs

- `heat_conduction.py`'s own comment flags a known gap one level up (not in this file):
  `main_solve.py`'s He supercritical Z-correction TODO (Z ~ 1.04-1.06), referenced from
  `physics/gas_flow/governing_equations.py`'s docstring.

## Change history

No dated changelog is evidenced for most files in this folder beyond the Phase 1
restructure that produced the deprecation shims (`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`,
2026-07-13 per CLAUDE.md). `docs/solver_design/FV_CORE_REWORK_PLAN.md` (2026-08-18 note)
records `heat_conduction.py`'s `fluxes_at_Tbar()` as cross-validated at `rel=1e-12`
against the new `core/wall.py`, and its fin-efficiency helpers as the intended basis for
Stage C/F's rib model.

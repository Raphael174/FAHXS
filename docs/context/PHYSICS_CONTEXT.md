# Physics Context

## Units And Conventions

- SI units throughout unless comments state otherwise.
- Pressure is Pa, temperature is K, mass flow is kg/s.
- Friction factors in maintained solver paths are Darcy-Weisbach factors.
- Heat transfer coefficients are W/m2/K, heat rates are W internally, and some summaries print kW.

## Chemistry

Combustion utilities live in `physics/combustion_chemistry/`.

Main path:

- `choose_fuel(fuel)` maps fuel names to mechanism files and fuel composition.
- `combustion_gas_solve` initializes a Cantera phase and can remove enthalpy during the HX march.

Supported fuel selectors in current docs include:

| Fuel selector | Mechanism |
|---|---|
| `diesel-C16H34` | `RenKokjohn_surrogate.yaml` |
| `POSF10325` | `A2highT.yaml` |
| `gasoline-E5`, `gasoline-E10` | `llnl_gasoline_323.yaml` |
| `H2` | `H2-O2_Burke2012.yaml` |

Steady `numericalProp.chemistry_model` and transient
`transientProp.chemistry_transient` both default to `finite_rate`.

The `finite_rate` path uses the FPV manifold in
`physics/combustion_chemistry/fpv_manifold.py` and marches progress variable
through the HX cooling leg. The `equilibrium` and `frozen` modes remain available
for comparison and validation.

## Heat Transfer

Core module: `physics/heat_transfer_correlations.py`.

Helical-coil coolant side:

- `dispatch_nu_coil()`
- Selectors: `mori1967`, `Gnielinski`
- Mori/Nakayama low-Pr branch is the default helium curvature correction.

Helical-coil hot-gas shell side:

- `dispatch_nu_shell()`
- Selectors: `salimpour2008`, `ahmed_toroid`, `churchill_bernstein_tightcoil`, `churchill_bernstein`
- Default in `input_data.py`: `salimpour2008` with `combustorProp.Nusselt_correction`.

Straight tube side for shell-and-tube:

- `dispatch_nu_tube_straight()`
- Current selector: `gnielinski_blended`
- Blends laminar entrance, transitional, and turbulent Gnielinski regimes.
- For `shellTubeProp.inside_tube_choice = "grooved"`, tube-side Nu and friction
  use a Vicente/Cruz-style corrugated-tube path:
  `phi = corrugation_thickness^2 / (corrugation_pitch * D_i)`.
  `tube_grooved_Nu_factor` and `tube_grooved_f_factor` are calibration factors
  on top of that literature form, not the base physics.

## Friction

Core module: `physics/friction_correlations.py`.

Helical-coil coolant side:

- `dispatch_friction_coil()`
- Selectors: `CurvedPipeAli2024`, `Colebrook1939`
- Ali et al. high-curvature branch is usually active at design conditions.

Straight tube side:

- `dispatch_friction_tube_straight()`
- Blends Hagen-Poiseuille laminar, transitional, and Colebrook turbulent friction.
- `friction_corrugated_tube_vicente()` is used for grooved/corrugated tubes and
  returns a Darcy friction factor.

## Liquid Coolant / Boiling (Experimental)

Core modules (moved into `physics/liquid_flow/` in Phase 1 of the integration
plan; old top-level paths are deprecation shims, see
`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`):
`physics/liquid_flow/correlations.py`, `physics/liquid_flow/chf.py`,
`physics/liquid_flow/dispatch.py`, `physics/liquid_flow/governing_equations.py`,
`physics/liquid_flow/hx_adapters.py`, `physics/liquid_flow/sanity_checks.py`.

This is a `p,h`-state liquid (boiling) coolant path, developed alongside
literature in `docs/reference`. Wiring status (Phase 2 of the integration
plan):

- **Helical steady (`main_solve.py`) is wired in**: when
  `coolantProp.coolant_model == "equilibrium_liquid"`, the coupled march uses
  `evaluate_coolant_closure()` for properties/HTC/friction/CHF and
  `dh/dx = dQ/mdot`, `dp/dx = -friction` as the governing equations, gated
  against the legacy ideal-gas/PropsSI(T,p) path. Co-flow is fully
  self-consistent (cross-checked against the postprocess bridge). The plain
  counter-flow march starts from the legacy `T_out`/`p_out` guess (a
  single-phase `(T,P)` state, so it cannot represent a genuine two-phase
  starting point) — the same prescribed-outlet discrepancy already
  documented for the gas march below. `solve_counterflow_liquid_reference()`
  resolves this by shooting the hot-end starting enthalpy to match the
  user's physical `T_in`/`p_in` instead (see
  `docs/context/SOLVER_CONTEXT.md`). A `check_liquid_march()` sanity report
  (energy closure, temperature ordering, saturation consistency, pressure
  monotonicity, bounds, and a hard CHF/dryout gate) runs automatically at the
  end of a liquid-mode solve.
- **Shell-and-tube (`main_solve_shellntube.py`) is still postprocess-only**:
  it computes the coolant march with direct CoolProp `PropsSI` calls
  regardless of `coolant_model`; liquid physics only runs as an opt-in
  post-process step (`liquid_coolant_postprocess()`) that consumes an
  already-converged `dQ` duty profile to report `p,h`/quality/void/CHF
  diagnostics without feeding back into wall temperature. This is Phase 3 of
  the integration plan, not yet done.
- Both transient solvers have no liquid coolant path at all.

See `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md` and
`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` for the authoritative
status writeup.

Implemented correlations in `liquid_flow/correlations.py` (CHF in
`liquid_flow/chf.py`):

- Single-phase: smooth-pipe Darcy friction blend, 3.66 laminar / Gnielinski
  turbulent Nusselt blend.
- Two-phase state: CoolProp-backed homogeneous-equilibrium-model (HEM)
  saturation/quality/void fraction from `p,h`.
- Boiling HTC: Gungor-Winterton 1986 (vertical/high-Froude form only; the
  horizontal/low-Froude correction is not implemented, so horizontal
  shell-and-tube use is unvalidated for orientation effects), plus a Yu et al.
  2002 modified-ANL fit.
- Two-phase friction: Müller-Steinhagen-Heck pressure-gradient correlation,
  Lockhart-Martinelli/Chisholm and Yu2002 multipliers, HEM acceleration
  pressure-gradient term.
- CHF: Groeneveld 2006 LUT lookup (`docs/reference/external/2006LUTdata.txt`,
  trilinear interpolation) with local diameter correction. The full 2006 table
  is supplied externally, not encoded in the repo; only spot-checked at the
  LUT's page-9 reference values.
- **Chen (1962), Shah, and Kandlikar boiling correlations are referenced in
  `docs/reference` but not implemented** — only Gungor-Winterton and the
  Yu2002 fit are coded.
- No flow-regime map (bubbly/slug/annular/mist) and no post-CHF/dryout
  degraded-HTC model; post-saturation vapor is treated as plain single-phase
  vapor.

Integration bridges (`liquid_flow/hx_adapters.py`) exist for both geometries but are
explicitly pseudo-1D placeholders:

- `solve_helical_coil_liquid_from_duty/_from_data_master()` maps
  `data_master["dQ"]`/`L_ch` onto the heated-channel solver using `Dh_coil`
  and per-coil flow area; no centrifugal/secondary-flow boiling enhancement.
- `solve_shelltube_shellside_liquid_from_duty/_from_tube_result()` uses
  Bell-Delaware `S_m` as a pseudo shell-side flow area; self-documented as
  "not yet a final shell-side boiling model."

`liquid_flow/dispatch.py` is the dispatcher entry point for future wiring
(`coolantProp.coolant_model in ("single_phase_coolprop", "equilibrium_liquid")`,
`liquid_heat_transfer_model`, `liquid_pressure_drop_model`, `liquid_chf_model`);
today only the postprocess path calls it.

Transient solvers (`main_solve_transient.py`,
`main_solve_shellntube_transient.py`) have no liquid/boiling code path at all —
not even a postprocess bridge. Transient liquid coolant is unimplemented.

## Wall Conduction

Core module: `physics/heat_conduction.py`.

`OneDimensionalSteadyConduction_ShellnHelicalTube` solves a radial wall conduction balance with convection on both sides and optional radiation. It is reused by both helical-coil and shell-and-tube solvers; shell-and-tube passes `hot_side="inner"` because hot gas is inside the tubes. Both steady conduction and transient `fluxes_at_Tbar()` must honor this flag so hot/cold perimeters match.

Transient wall integration uses a single thickness-mean wall temperature per
axial node. `fluxes_at_Tbar()` reconstructs hot and cold face temperatures from
a quasi-static quadratic profile and returns hot/cold heat fluxes.

Pre-ignition GOX chilldown:

- If `ignition=0` and LOX/GOX flow is scheduled, the transient hot-side branch
  uses CoolProp Oxygen properties and sensible heat exchange.
- When `ignition=1`, the solver switches to the combustion gas manifold.

## Radiation

Core modules:

- `physics/radiation_model/radiation_build.py`
- `physics/radiation_model/radiation_equations.py`

The helical solver builds a WSGGM backend from `ehlme2025_mixture.json` when `numericalProp.radiation_ON` is true. It uses H2O and CO2 mole fractions from the Cantera gas state and a mean beam length based on geometry and `CorrelationCoefficients.mbl_factor`.

## Materials And Stress

Material property dispatch:

- `mechanical/material_specs/material_temperature_strength.py`
- Supported material names include `ST316L` and `INCO718`.

Loads:

- `mechanical/loads.py`
- Helical coil uses pressure stress plus thermal stress through the tube wall.
- Shell-and-tube checks external-pressure stress and thin-tube collapse.

The maintained helical solver evaluates yield at the hot wall by default when `numericalProp.yield_at_hot_wall` is true.

## Calibration Knobs

`CorrelationCoefficients` in `input_data.py` is the single source of tunable empirical factors. The most important current knobs are:

| Field | Primary observable |
|---|---|
| `ali_c_hi` | Helium pressure drop `dp_c`. |
| `salimpour_a` | Hot-gas outlet temperature if measured. |
| `mori_a_lo` | Secondary coil-side heat-transfer sensitivity. |
| `mbl_factor`, `emissivity_wall` | Radiation contribution. |
| `n_tube_gas`, `zukauskas_C_factor`, `bell_Jl_factor`, `bell_Jb_factor` | Shell-and-tube calibration. |
| `tube_grooved_Nu_factor`, `tube_grooved_f_factor` | Calibration on top of corrugated-tube Nu/friction. |

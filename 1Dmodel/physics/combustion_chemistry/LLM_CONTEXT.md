# physics/combustion_chemistry/ LLM Context

## Scope

Cantera-backed combustion chemistry for the hot-gas side: initial-equilibrium
combustion state, the fuel/mechanism registry, and the two "no per-node Cantera call"
manifold tables (finite-rate FPV and equilibrium/frozen) that the maintained steady and
transient solvers actually march against. This is where `numericalProp.chemistry_model` /
`transientProp.chemistry_transient` resolve to real code.

## Contents

| File | Role |
|---|---|
| `combustion_gas.py` | `choose_fuel(fuel)` — maps a fuel selector string to its Cantera mechanism file, fuel mass-fraction composition, and heat of vaporization. `combustion_gas_solve` — builds the initial Cantera `Solution`, solves LOX→GOX injection enthalpy loss plus fuel HP-equilibrium to get the initial combustion state, and exposes `remove_energy()` to march enthalpy removal (with or without re-equilibration) for validation/dQ studies. Imported by all four `main_solve*.py` files — this is the real, wired entry point. |
| `fpv_manifold.py` | `build_fpv_manifold()` / `FPVManifold` — the finite-rate FPV (h, progress-variable Yc) manifold: offline constant-enthalpy reactor relaxations sampled onto a regular (h, c) grid, cached to disk under `cache/fpv_manifolds/` (hashed on gas/inlet/grid params). Runtime interpolation is table lookup only, so the transient march makes zero Cantera calls. This is the required default chemistry path for the project's main diesel/O2 regime (per CLAUDE.md). |
| `gas_manifold.py` | `build_equilibrium_manifold()` / `EquilibriumManifold` — the simpler C0 slice: a 1D table of gas state vs. enthalpy removed at fixed inlet composition, either re-equilibrated at each level (`mode="equilibrium"`) or frozen (`mode="frozen"`). Shared by both transient solvers so neither duplicates the "no Cantera in the march" tabulation logic. |
| `stoichiometry.py` | `compute_OF_st()` — small standalone script computing mass-based stoichiometric O/F for a pure-O2 mixture from a mechanism's elemental fuel composition. Not imported elsewhere. |
| `A2highT.yaml`/`.cti`, `H2-O2_Burke2012.yaml`, `RenKokjohn_surrogate.yaml`, `llnl_gasoline_323.yaml`, `llnl_gasoline_Detailed.yaml` | Cantera mechanism files referenced by `choose_fuel()`. |

## Data subfolders (not code)

`Pei2015_Diesel_surrogate/`, `RenKokjohn_Diesel-JetA-Gasoline/`, and `gasoline_surrogate/`
hold raw Chemkin-format mechanism data (`mech.inp`, `therm.dat`, `trans.dat`) and, in
`gasoline_surrogate/`, a handful of one-off build/fix scripts (`build_SansPlomb95.py`,
`fix_ck_numbers.py`, `latent_heats_gasoline.py`) used to produce the `llnl_gasoline_*.yaml`
mechanisms consumed by `choose_fuel()`. Treat these as mechanism source data, not a
maintained code path.

## Wiring / fuel selectors

Per `docs/context/PHYSICS_CONTEXT.md`, currently supported `choose_fuel()` selectors:

| Fuel selector | Mechanism |
|---|---|
| `diesel-C16H34` | `RenKokjohn_surrogate.yaml` |
| `POSF10325` | `A2highT.yaml` |
| `gasoline-E5`, `gasoline-E10` | `llnl_gasoline_323.yaml` |
| `H2` | `H2-O2_Burke2012.yaml` |

`finite_rate` (via `fpv_manifold.py`) is the required default for steady and transient
runs in the project's main diesel/O2 high-heat-extraction regime; `equilibrium`
(`gas_manifold.py`, `mode="equilibrium"`) and `frozen` (`mode="frozen"`) remain available
for comparison/validation only — frozen chemistry is explicitly flagged in CLAUDE.md as
physically wrong for the main regime.

## Sharp edges

- **Persistent Cantera objects must be reset from cached inlet T/p/Y state before
  repeated sweeps or manifold builds.** Reading a prior sweep's already-cooled state as a
  new inlet is a known failure mode (CLAUDE.md). `combustion_gas_solve` and the manifold
  builders both re-establish `gas.TPY` from the caller-supplied inlet state at the start
  of their build/solve routines for this reason.
- **No per-node Cantera calls in a march** — both `fpv_manifold.FPVManifold.state()` and
  `gas_manifold.EquilibriumManifold.at()` are pure table interpolation
  (`np.interp`/bilinear), by design, so the transient integrators never call into Cantera
  per timestep/node.
- `fpv_manifold.py`'s CVODE relaxation can hit stiffness failures at cold enthalpy levels
  (chemistry frozen); the code explicitly tolerates a truncated trajectory there because
  the `c=1` end is always separately anchored to the exact equilibrium state, not
  extrapolated from the failed trajectory.

## Test coverage

No test imports `combustion_gas.py`/`fpv_manifold.py`/`gas_manifold.py` by name directly;
they are exercised through `tests/test_steady_baseline_regression_finite_rate.py`,
`tests/test_transient_baseline_regression.py`, and `tests/test_transient_core_coolant_fv.py`.
The three data subfolders have zero test coverage (expected — mechanism data, not code).

## Change history

No dated changelog evidenced in this folder's own files beyond the FPV manifold's own
docstring reference to "DESIGN_PLAN_shellntube_transient.md section 4b/C1" as its design
origin. `docs/context/PHYSICS_CONTEXT.md` documents finite-rate FPV as the required
default; no other dated history found for this folder specifically.

- **2026-08-19**: deleted five dead, zero-importer exploratory/legacy scripts
  (`finite_rate_solver.py`, `combustion_sizing_proto.py` — hardcoded a stale
  absolute Windows path, would not have run as-is —, `script1_sizingChemistry.py`
  and `study_energy_removal_combustion.py` — undocumented near-duplicates of
  `combustion_gas.py`'s `choose_fuel`/dQ-study content —, and `other.py` — a
  single print helper). Confirmed zero importers via grep before removal.

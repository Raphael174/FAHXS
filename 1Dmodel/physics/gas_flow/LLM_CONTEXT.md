# physics/gas_flow/ LLM Context

## Scope

The compressible, single-component ideal-gas quasi-1D governing equations used by the
gas/helium coolant march and the hot-gas side of the maintained solvers. This is a
small, single-file package — it is the real home of what `physics/governing_equations.py`
(top-level) now just re-exports as a deprecation shim. Distinct in both name and physics
from `physics/liquid_flow/governing_equations.py`, which is a `p,h`-state real-fluid
(boiling/liquid) equation set; the two happened to share a filename before the Phase 1
restructure and must not be confused.

## Contents

| File | Role |
|---|---|
| `governing_equations.py` | `dT__dx_IdealGas`, `dU__dx_IdealGas`, `dp__dx_IdealGas` / `dp__dx_IdealGas_logical`, `drho__dx_IdealGas` / `drho__dx_IdealGas_logical`, `dT_g__dx` — the coolant temperature/velocity/pressure/density gradients along a channel for compressible ideal-gas flow, plus the hot-gas temperature evolution equation (`dT_g__dx`, no pressure loss, pure UA/(mdot·cp) relaxation). |
| `__init__.py` | Package docstring only (no re-exports); callers import `governing_equations` directly. |

## Key correctness points

- Valid for a single-component ideal gas (fluid name/molar mass supplied by the caller
  via `coolantProp`; nothing here hardcodes a specific gas), quasi-1D, no real-gas
  compressibility correction. The module docstring flags this explicitly and points to
  `main_solve.py`'s known TODO for the helium supercritical Z-correction (Z ~ 1.04-1.06).
- `dU__dx_IdealGas`'s denominator (`m_dot - A*p/U - A*p*U/(T*cp)`) is the compressible
  choking term — it can approach zero near Mach 1; the maintained solvers keep coolant
  Mach numbers low by design rather than handling the singularity here.
- `dp__dx_IdealGas` and `dp__dx_IdealGas_logical` are two algebraically different forms
  of the same pressure gradient (one via T/U/A, one via T/rho) — check which one a given
  solver actually calls before assuming they are interchangeable in a partially-updated
  state.
- Friction enters through the caller-supplied Darcy factor `f` inside `dU__dx_IdealGas`
  — this module does not compute friction itself (`physics/friction_correlations.py`
  does); the convention throughout maintained solver paths is Darcy, not Fanning.

## Test coverage

No dedicated test file imports `physics.gas_flow.governing_equations` by name; it is
exercised indirectly through `tests/test_steady_baseline_regression*.py` and
`tests/test_transient_baseline_regression.py` via the four `main_solve*.py` solvers that
call into it.

## TODOs

None stated directly in this file; the Z-correction TODO it references lives in
`main_solve.py`, outside this folder's scope.

## Change history

Per `CODEBASE_MAP.md`: this module was moved here from the top-level
`physics/governing_equations.py` during the Phase 1 restructure
(`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`); the old location is now a
deprecation shim. No other dated history evidenced for this folder.

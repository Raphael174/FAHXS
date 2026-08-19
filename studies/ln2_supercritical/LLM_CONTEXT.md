# studies/ln2_supercritical/ LLM Context

## Scope

Single-script study: first supercritical-LN2 (liquid nitrogen) coolant run on the maintained helical steady solver, exercising `coolant_model="equilibrium_liquid"` with Nitrogen above its critical pressure (Pc = 33.96 bar) so there is no boiling dome, only a pseudo-critical transition near 145.7 K.

## Contents

| File | Role |
|---|---|
| `run_ln2_helical.py` | Runs `main_solver` with `coolantProp(coolant="Nitrogen", coolant_model="equilibrium_liquid")`, T_in=100 K, p_in=80 bar, mdot_c=0.10 kg/s, co-flow, `chemistry_model="frozen"` (frozen is used deliberately: coolant-side supercritical physics is independent of hot-gas chemistry, and frozen is fast/sufficient for this first check). Imports `L_HX_MAX_6M` from `studies/water_vs_helium/run_case.py` via a manual `sys.path` insert, and calls `solver._check_global()` for an energy-balance sanity check against an independent CoolProp enthalpy-delta calculation. |

## Pitfalls

- Depends on `studies/water_vs_helium/run_case.py` being importable (path-inserted, not package-relative) for the `L_HX_MAX_6M` constant (coil arc length calibrated to exactly 6.000 m). If that file moves or its constant changes, this script breaks silently until run.
- Run with `python -m studies.ln2_supercritical.run_ln2_helical` from repo root (or directly; both are documented in the module docstring).

## TODOs (only evidenced items -- do not invent)

None evidenced.

## Change history

No dated change history evidenced. Docstring calls this the "First supercritical-N2 (LN2 coolant) run," indicating this is exploratory/first-pass rather than a settled validation case.

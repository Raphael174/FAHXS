# mechanical/material_specs/ LLM Context

## Scope

Temperature-dependent material property tables (yield strength, elastic
modulus, CTE, thermal conductivity, specific heat) for the wall materials the
solvers support, plus one dispatch function that returns all of them as
callables for a given material name. Feeds `mechanical/loads.py`'s stress
checks and the transient solvers' wall-energy ODE.

## Contents
| File | Role |
|---|---|
| `material_temperature_strength.py` | Piecewise-linear (`scipy.interpolate.interp1d`) property tables for `"ST316L"` (316L stainless) and `"INCO718"` (Inconel 718): `YieldStrengthR02_*`, `ElasticityModulus_*`, `CTE_*`, `compute_*_conductivity`, `compute_*_cp`. `init_material_temperature_properties(material)` is the single entry point, returning `(CTE, E, Yield, Lambda, density, poisson, Cp)`. |

## Key correctness points

- All property functions take temperature in **degrees Celsius**, not Kelvin
  — callers must convert. Each function clamps outside its tabulated range
  (returns the boundary value) rather than extrapolating, except
  `compute_PER718_conductivity` which prints nothing (the analogous
  out-of-range print statement is commented out) but still clamps.
- `density` is returned as a constant scalar per material (not
  temperature-dependent) — documented in the function's own docstring as a
  deliberate approximation ("metals' rho varies <2% up to 800 degC, negligible
  for wall thermal mass"), not an oversight.
- `Cp` (specific heat) is appended as the 7th/last return value specifically
  so any older 6-tuple-unpacking call site keeps working; only the transient
  wall-energy ODE (`rho*cp*delta*dTbar/dt`) actually needs it — irrelevant to
  the steady solver. Both `Cp` tables cite `DESIGN_PLAN_shellntube_transient.md`
  section 4.7 as the consumer.
- **`init_material_temperature_properties` only supports `"ST316L"` and
  `"INCO718"`** — passing anything else raises `ValueError`. Note the
  `_POISSON` dict separately has a `"CuCrZr"` entry (0.30, sourced from ASM
  Handbook) that is dead/unused: there is no `CuCrZr` branch in
  `init_material_temperature_properties`, so requesting that material still
  raises `ValueError` despite its Poisson ratio being tabulated. This looks
  like groundwork for a copper-alloy material that was never finished — do
  not assume CuCrZr is a supported wall material.
- 316L and Inconel 718 tables come from different sources: 316L from AISI
  "High Temperature Characteristics of Stainless Steels" (Designer's Handbook
  N°9004); Inconel 718 yield/E/CTE credited in-file to "Abdulah from Mech
  team", conductivity table unattributed (`PER718`), Cp from the Special
  Metals INCONEL 718 datasheet (SMC-045).

## TODOs

None found as explicit TODO/FIXME markers. The unused `CuCrZr` Poisson entry
(see above) is a strong implicit gap but is not a stated TODO — flagged above
as a sharp edge rather than invented as a TODO.

## Change history

No usable git history (single initial commit; folder currently
uncommitted/untracked). No dated comments in-file.

# mechanical/ LLM Context

## Scope

Structural-check formulas for the pressurized coolant tube (helical coil or
shell-and-tube), consumed by the maintained steady solvers to fill the
`stress_*`/`CTE`/`Modulus`/`Yield` columns of `data_master`. This package
does not size anything or run its own march — it is called per-node from the
solver's coupled loop and its outputs are purely diagnostic (margin
reporting), not fed back into the thermal/flow solve.

## Contents
| File | Role |
|---|---|
| `loads.py` | Roark's-formula stress functions: internal-pressure hoop stress, external-pressure hoop stress, thermal (radial gradient) stress at inner/outer wall, and elastic collapse-pressure buckling check. |
| `geometry/` | Coil and shell-and-tube geometry builders — see `geometry/LLM_CONTEXT.md`. |
| `material_specs/` | Temperature-dependent material property tables — see `material_specs/LLM_CONTEXT.md`. |

## Key correctness points

- `stress_pressure_tube`/`stress_external_pressure_tube` are the SAME thin-wall
  hoop-stress formula (`P*Dh/(2*t)`), differing only in sign convention:
  external pressure returns a negative (compressive) stress. Which one a
  solver calls depends on loading direction — the helical coil sees internal
  pressure (coolant inside), the shell-and-tube tubes see the opposite
  (hot gas at ~1-5 bar inside, ~90 bar He shell-side outside) — this mirrors
  the repo-wide `hot_side="inner"` convention documented in `/CLAUDE.md`.
- `collapse_pressure_thin_tube` is the Euler/Bresse-Bryan elastic-buckling
  critical pressure for a long unsupported thin tube (no stiffening rings).
  Compare `|P_ext| / P_cr` as a margin; only meaningful for the
  external-pressure (shell-and-tube tube-wall) loading case.
- `stress_thermal_tube` returns `[inner, outer]` stresses for a radial
  temperature gradient through the wall; sign convention documented inline
  (cold-inside/hot-outside gives +tension inside).
- All formulas are from Roark's Formulas for Stress and Strain, 7th ed.
  (page numbers cited in each docstring) — treat page references as the
  primary source of truth if a formula looks off.

## TODOs

None found (no TODO/FIXME markers, no stated gaps in this folder's files).

## Change history

No usable git history (single initial commit; this folder shows as
uncommitted/untracked in the working tree per current `git status`). No
dated comments in-file to build a changelog from.

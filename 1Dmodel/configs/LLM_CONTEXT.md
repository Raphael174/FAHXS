# configs/ LLM Context

## Scope

Reference/example `input_data.py` dataclass presets — data, not code paths
the maintained solvers import automatically. Each file is a standalone
Python module that builds `coolantProp`/`hotgasProp`/`combustorProp`/
`numericalProp`/`system_requirements` instances for one named design point,
meant to be imported explicitly by a user script or VS Code interactive cell
(see each file's own "Usage" docstring, e.g.
`from hps_combustor.configs import preset_500kW_INCO718 as cfg`) — never
imported by `main_solve*.py` or any other production solver on its own.

## Contents
| File | Role |
|---|---|
| `input_data - 500kW Inconel 718.py` | 500 kW shell-and-helical-coil design point, JET-A (POSF10325)/O2, O/F=4.12, Inconel 718 coil. |
| `input_data - 500kW ShellnCoil - ToroidalConvection.py` | 500 kW shell-and-coil variant exercising a toroidal-convection correlation option. |
| `input_data - CFD validation study - corrected correlations.py` | Preset used for comparison against a CFD validation study, with corrected correlation coefficients. |
| `__init__.py` | Package marker only (empty). |

## Notes

- Filenames contain spaces and a leading `input_data - ` prefix, matching how
  the user names these presets outside the codebase — not a Python-importable
  module name pattern by itself; consumers import via the package
  (`hps_combustor.configs`) using whatever `__init__.py`/explicit-path
  mechanism is set up, per each file's own docstring.
- Treat these as calibration/reference recipes to copy from, not shared
  defaults — `/CLAUDE.md` is explicit that `coolantProp`/`combustorProp`
  defaults elsewhere stay fluid-agnostic (Helium baseline); these presets are
  where fluid- and design-point-specific values belong instead.

## TODOs

None found.

## Change history

No usable git history (single initial commit; folder currently
uncommitted/untracked per `git status`). No dated comments in-file.

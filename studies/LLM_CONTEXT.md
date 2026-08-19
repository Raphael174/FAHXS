# studies/ LLM Context

## Scope

Parametric-sweep and calibration scripts that exercise the maintained solvers (`main_solver`, `shellntube_solver`) from outside the `1Dmodel/` package surface. Per `CODEBASE_MAP.md` these live here instead of in the package root so the package stays free of one-off experiments. Two subfolders (`ln2_supercritical/`, `water_vs_helium/`) hold self-contained studies with their own outputs; see their own `LLM_CONTEXT.md` files.

## Contents

| File | Role |
|---|---|
| `README.md` | States the production entry points (`python -m hps_combustor.main_steady` / `main_transient`) and warns that these scripts may need import-path fixes since they were moved out of `1Dmodel/`. |
| `run_correlation_study.py` | Sweeps shell x coil Nusselt-correlation combinations against a 200 kW experimental diesel/He case (chamber 136 mm, coil 7/2.4 mm, co-flow shellnHelicalTube); prints an error-sorted comparison table, resistance breakdown, and proposed `CorrelationCoefficients` tuning. |
| `run_grid_study.py` | 2D grid sweep of `N_arc_steps_per_turn` (50-120) x `Nusselt_correction` (0.01-0.26) against the same experimental case; targets T_c_out ~ 420degC / Q_He ~ 200 kW. Has a `FAST_MODE` flag to use frozen chemistry (~10x faster, <=5% error on absolute Q). |
| `run_diameter_pressure_study.py` | Helium hydraulic-diameter -> pressure-drop sweep under a *prescribed* linear coolant temperature ramp (90 K -> 650 K), deliberately sidestepping the unvalidated hot-gas Nusselt number to answer a purely mechanical question. CLI flags `--mdot`, `--friction-error`. Standalone; does not modify the core solver. |

## Import mechanism

Each script begins with an identical "package bootstrap" block (`if __name__ == "__main__" and __package__ is None: ...`) that synthesizes a fake package alias (`_hps`) around the script's own directory and re-invokes itself via `runpy.run_module`, purely because the folder name `1Dmodel` is not a valid Python identifier and these scripts use package-relative imports (`from .main_solve import main_solver`). This is a workaround, not an example to imitate elsewhere — new scripts outside `1Dmodel/` should just `import hps_combustor` after an editable install, as the two subfolder studies do.

## TODOs (only evidenced items -- do not invent)

- `README.md` itself flags that these scripts "may need small import-path updates before rerunning."

## Change history

No dated change history is evidenced in these files beyond the runtime-estimate note in `run_grid_study.py`'s header comment (not a changelog entry).

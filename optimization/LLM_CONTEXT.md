# Optimization And Calibration Context

## Scope

The `optimization/` directory wraps the `hps_combustor` steady helical-coil solver for parameter calibration and design optimization.

Run from repository root after editable install:

```powershell
python -m pip install -e .
python -m optimization.quick_optimize
```

## Main Files

| File | Role |
|---|---|
| `calibrate.py` | Least-squares, MAP, MCMC, sensitivity, identifiability, and posterior predictive utilities. |
| `quick_optimize.py` | Mass minimization with wall-temperature, pressure-drop, and Mach constraints. |
| `optimizer_desVar.py` | Design variable definitions, ranges, and activation flags. |
| `optimizer_dataMap.py` | Converts optimizer vectors into input dataclass instances. |

## Calibration Model

`calibrate.py` calibrates fields in `hps_combustor.input_data.CorrelationCoefficients`. It uses:

- `CalibrationRecord` for one test condition and measured observables.
- `compute_sensitivities()` for normalized finite-difference sensitivities.
- `identifiability_check()` for Fisher information conditioning.
- `calibrate_ls()` for fast weighted least squares.
- `calibrate_map()` for prior-regularized point estimates.
- `calibrate_mcmc()` for full posteriors; requires `emcee`.

Important identifiability guidance from `docs/project_docs/CALIBRATION_METHODOLOGY.md`:

- With helium pressure drop only, calibrate `ali_c_hi`.
- With hot-gas outlet temperature, also calibrate `salimpour_a`.
- `Q_He` can be weakly sensitive to Nusselt number in saturated counter-flow regimes and should often be treated as an energy-balance check.

## Optimization Model

`quick_optimize.py` currently minimizes HX plus combustor mass using `scipy.optimize.minimize` and constraints:

- `T_wg_max < 950 K`
- `dP_c < 5 bar`
- `Mach_g_max < 0.5`

The objective runs a full `main_solver` evaluation, so optimization is expensive. Cache or surrogate only if the task explicitly calls for performance work.

## Pitfalls

- Calibration imports `hps_combustor`, so editable install is required unless `PYTHONPATH` is adjusted.
- MCMC can be slow because each likelihood call runs the full solver.
- Do not rename `CorrelationCoefficients` fields casually; calibration prior specs and design maps depend on exact names.

## Change history

- **2026-08-19**: deleted `optimizer_solve.py` (an older optimization execution
  path referencing a nonexistent `ToasterProp`/`toaster` object — stale/broken,
  zero importers, confirmed via grep before removal).

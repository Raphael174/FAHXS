# studies/water_vs_helium/ LLM Context

## Scope

Full comparative study: Helical (6 m coil) vs. Shell-and-tube geometry, x Helium vs. Water coolant, run at 50/100/150 g/s coolant mass flow (12 base runs), plus a secondary `steam_design_study/` feasibility sub-study. This study's Part 3 work (shell-side `(p,h)` liquid coupling, Phase 3 of `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`) actually changed `1Dmodel/main_solve_shellntube.py` — see `STUDY_REPORT.md` section 0.

## Contents

| File | Role |
|---|---|
| `STUDY_REPORT.md` | The full write-up: methodology table, results grid, and section 0's account of the shell-side `(p,h)` enthalpy-based coolant march added to `main_solve_shellntube.py::_shell_side_march`/`_shell_h_at` (previously temperature-only with zero latent-heat accounting). Also documents a pre-existing Helium baseline regression failure (both helical and shell-and-tube) found and root-caused as `input_data.py` default drift unrelated to this study's changes -- confirmed by byte-for-byte reversion testing. |
| `run_case.py` | Single-case runner (CLI: `python run_case.py --geometry {helical,shellntube} --fluid {water,helium} --mdot <kg/s>`). Fixed conditions: T_in=303.15 K, p_in=80 bar, co-flow, `finite_rate` chemistry, `hotgasProp()` defaults. Defines `run_helical()`/`run_shellntube()`, `L_HX_MAX_6M = 0.47922164350748064` (bisection-calibrated so the helical coil's arc length, not axial `L_HX_max`, is exactly 6.000 m -- imported by `studies/ln2_supercritical/run_ln2_helical.py`), `make_plot()` (4-panel: wall T, coolant T, pressures, thermal-resistance breakdown), and `save_case()` which zips `.npz` + JSON summary + PNG into `raw_data/`. |
| `make_comparison_plots.py` | Loads all `raw_data/*_summary.json` files and builds `plots/comparison_summary.png` across the full (geometry, fluid) grid. |
| `steam_design_study/calibrate_hotgas_scaling.py` | Calibrates how `Q_tot` scales with hot-gas mass flow (shell-and-tube + helium, sweeping `mdot_hotgas` in [0.10, 0.30, 0.60, 1.0] kg/s) to estimate the hot-gas flow needed for a 30 L/s steam design-feasibility target; writes `hotgas_scaling_calibration.json`. Imports `run_shellntube` from `run_case.py` via a `sys.path` insert (not package-relative). |

## Generated-output subfolders (not code)

`logs/`, `plots/`, `raw_data/`, and `steam_design_study/` (partly) hold generated artifacts from running the above scripts -- per-run stdout logs, PNG plots, zipped `.npz`+JSON result bundles, and a calibration JSON. Treat these as data, not something to hand-edit; regenerate via the scripts above if stale.

## Pitfalls

- Coolant CHF checks for water use `liquid_chf_lut_path = "docs/reference/external/2006LUTdata.txt"` (Groeneveld 2006 LUT), a path relative to repo root.
- `run_case.py` imports are package-style (`from hps_combustor...`), unlike the top-level `studies/*.py` scripts' bootstrap-shim pattern -- this folder assumes an editable install (`pip install -e .`), not the `_hps` alias trick.

## TODOs (only evidenced items -- do not invent)

- `STUDY_REPORT.md` section 0.1 flags the pre-existing Helium baseline regression-test failure (both geometries) as confirmed unrelated to this study but explicitly "not fixed here ... flagged for separate follow-up."
- Documented known limitation: shell-side cross-flow enhancement of boiling HTC is not captured by the Gungor-Winterton correlation used (a tube-flow correlation), called out in `STUDY_REPORT.md` as a deliberate simplification, not an oversight.

## Change history

`STUDY_REPORT.md` documents, as part of this study, the addition of shell-side `(p,h)` liquid coupling to `main_solve_shellntube.py` (Phase 3 of the liquid-coolant integration plan) plus an EOS-ceiling safety net for superheated steam exceeding CoolProp's Water-backend validity range (~3000 K), freezing state at the last EOS-valid point instead of crashing. No specific date is given in-file.

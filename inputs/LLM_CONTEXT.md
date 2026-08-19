# inputs/ LLM Context

## Scope

Example/dummy input files consumed by specific validation scripts and the transient schedule-file mechanism. Small, standalone folder — not part of the `hps_combustor` package import surface.

## Contents

| File | Role |
|---|---|
| `example_transient_schedule.csv` | Example schedule for `runProp.schedule_file` / `1Dmodel/schedule_inputs.py::apply_schedule_file()`, loaded by `main_transient.py`. Columns: `time_s, helium_m_dot_kg_s, helium_T_in_K, helium_p_in_Pa, diesel_m_dot_kg_s, lox_m_dot_kg_s, lox_T_in_K, ignition` — one of the "descriptive column name" layouts the loader accepts (it also accepts direct `m_dot_g_kg_s`/`OF` hot-gas columns, or a two-sheet Excel layout with `helium`/`propellants` sheets). Encodes a pilot-ignition ramp: near-zero He flow to t=2.5s, then He flow up with propellants still off, then a low-flow diesel/LOX pilot segment at t=5s (`ignition=1`), then full flow from t=6s to t=100s. |
| `HX_dummy_inlet.xlsx` | Helium inlet pressure/schedule input consumed by `1Dmodel/validation/low_mach_momentum_dummy_run.py` (default `schedule_file` and CLI `--schedule-file` value) and referenced by `bangbang_momentum_audit.py`. |
| `coupled_bangbang_hx_dummy_validation.json` | Default settings payload (`DEFAULT_SETTINGS_FILE`) for `1Dmodel/validation/coupled_bangbang_hx.py`; loaded via `json.loads`, drives the coupled bang-bang HX validation run and gets echoed back out as `settings_used.json` in that script's output directory. |

## Pitfalls

- `apply_schedule_file()` (`1Dmodel/schedule_inputs.py`) accepts decimal-comma numeric cells too, but CSV schedules using decimal commas must use semicolon or tab delimiters — `example_transient_schedule.csv` itself uses plain commas/dots.

## TODOs (only evidenced items -- do not invent)

None evidenced.

## Change history

No dated change history evidenced.

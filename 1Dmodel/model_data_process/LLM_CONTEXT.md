# model_data_process/ LLM Context

## Scope

Defines the `data_master` contract (the flat per-node dictionary every
maintained steady solver appends into during its march) and the two
dashboard classes that turn `data_master`/transient time-series data into
plots. This is the shared data layer between "the solver computed something"
and "the user can see it" — see `CODEBASE_MAP.md`'s "Primary Data Flow"
(steps 5-8): node quantities are appended to `data_master`, scalar summaries
are derived separately, `result_package.py` archives it, and this folder's
two plotting modules visualize it.

## Contents
| File | Role |
|---|---|
| `data_processing.py` | `_SOLVER_DATA_KEYS`: the full flat list of field names every steady solver populates (compressibility, thermal resistances, heat transfer, wall temps, coolant flow incl. liquid/boiling fields, hot gas flow, geometry, Biot numbers, mechanical stress). `make_solver_data()` returns a fresh `{key: []}` dict — the actual `data_master` factory. |
| `data_plotting.py` | `HXDashboard` — static multi-panel matplotlib figures for VS Code interactive use against a finished steady `data_master`: `.thermal()`, `.helium()`, `.combustion()`, `.mechanical()`, `.radiation()`, `.boiling()`, `.phase_change()`, `.mega()`, `.all()`. |
| `data_plotting_transient.py` | `TransientDashboard` — self-contained single-file HTML dashboard (embedded JSON + inline Canvas, no server/external assets) for time-resolved transient runs, consuming the `time_series` dict from `transient_solver._build_time_series`. Time scrubber + play control, wall-temperature x-t heatmap as the headline view. |

## Key correctness points

- `make_solver_data()` must be called fresh per solver instance —
  `_SOLVER_DATA_KEYS` entries are all mutable lists, so reusing a dict across
  solver runs silently accumulates data across runs. This matches the
  `/CLAUDE.md` "Current Project Truths" instruction to always call
  `make_solver_data()` per instance.
- The liquid/boiling fields (`enthalpy_c`, `dh_c__dx`, `quality_c`, `void_c`,
  `chf_margin_c`, `dp_c__dx_accel`) are only populated when
  `coolantProp.coolant_model == "equilibrium_liquid"`; for the default
  single-phase gas march they stay empty lists. `HXDashboard.boiling()`/
  `.phase_change()` are no-ops (with a message) when that data isn't present,
  and `.all()` only calls them if data is present.
- `HXDashboard.helium()`'s docstring explains a Z-vs-sound-speed panel swap:
  compressibility factor `Z` is gas-only (no standard liquid/two-phase
  definition), so the coolant panel plots sound speed instead, which the
  liquid march computes via Wood's equation inside the two-phase dome and
  the real EOS outside it — meaningful in both gas and liquid coolant modes.
- `TransientDashboard.to_html()` writes a fully self-contained `.html` file
  (default path: `1Dmodel/transient_dashboard.html`) — no server, no external
  JS/CSS dependencies; recursively converts numpy arrays/scalars to plain
  Python via `_jsonable()` before embedding as inline JSON. It follows the
  repo's `dataviz` skill conventions: categorical hues assigned by role in a
  fixed order, one sequential blue ramp for the heatmap, theme-aware
  light/dark CSS variables, always-present legends, no dual axes.

## TODOs

None found (no TODO/FIXME markers in this folder).

## Change history

No usable git history (single initial commit; folder currently
uncommitted/untracked per `git status`). No dated comments in-file beyond
the design-note style docstrings above.

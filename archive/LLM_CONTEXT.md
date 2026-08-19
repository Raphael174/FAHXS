# archive/ LLM Context

## Scope

Legacy solver variants, old input presets, injector/exploratory side work, and build artifacts, preserved for reference only. Per `CODEBASE_MAP.md` this is "usually safe to ignore unless recovering old behavior or presets." Do not build new workflows on archive files without first promoting the relevant logic back into the maintained package (`archive/README.md`'s own instruction).

## Contents

| Path | Role |
|---|---|
| `legacy_solvers/main_solve_CHT.py`, `main_solve_CHT_real_gas_localQ.py` | Older conjugate-heat-transfer solver variants, superseded by the maintained `1Dmodel/main_solve*.py` entry points. Not production entry points. |
| `input_presets/` | Old copied `input_data.py` snapshots (e.g. `input_data_shock_example.py`, `input_data_repoducingtest2_22.06_0.8Nuc.py`, `input_data_reproducingtest3_22.06_0.45Nuc_70arcs.py`) — historical parameter sets, useful only for recovering a specific old configuration. |
| `other_work/DeepFryerBench.py/` | `DeepFriedCoil_model_0Dtransient.py` — unrelated 0D transient exploratory model, not part of the combustor-HX solver chain. |
| `other_work/Injector_design/` (+ `Injector_design.zip`) | Standalone injector sizing scripts: straight/swirl/impinging orifice correlations, gasoline vapour pressure, combined gas-liquid jet sizing. Independent side work, no dependency on `hps_combustor`. |
| `build_artifacts/` | Generated `__pycache__` trees and `hps_combustor.egg-info*` metadata snapshots moved out of the working tree surface (several dated variants, e.g. `_after_venv_repair`, `_after_fpv_updates`, `_after_smoke_runs`) — pure build byproducts, not source. |
| `combustion_chemistry.zip` | Archived chemistry bundle. |

## TODOs (only evidenced items -- do not invent)

None evidenced.

## Change history

No dated change history evidenced beyond the descriptive `build_artifacts` subfolder names themselves (e.g. "after_venv_repair", "after_fpv_updates"), which mark snapshot points but carry no dates.

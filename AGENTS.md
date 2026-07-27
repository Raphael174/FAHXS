# LLM Working Context

This repository is a private Python engineering model for combustor/heat-exchanger co-design. Load this file first, then only open the focused context files needed for the task.

## Fast Orientation

- Main package name: `hps_combustor`.
- Source directory: `1Dmodel/`, mapped through `pyproject.toml` because `1Dmodel` is not a valid Python package identifier.
- User steady entry point: `1Dmodel/main_steady.py`.
- User transient entry point: `1Dmodel/main_transient.py`.
- Backend helical steady solver: `1Dmodel/main_solve.py`.
- Backend shell-and-tube steady solver: `1Dmodel/main_solve_shellntube.py`.
- Input dataclasses and calibration coefficients: `1Dmodel/input_data.py`.
- Standalone finite-rate flamelet/PFR research code: `research/flamelet_kit/`.
- In-progress liquid/boiling coolant physics:
  `1Dmodel/physics/liquid_flow/correlations.py`, `chf.py`, `dispatch.py`,
  `governing_equations.py`, `hx_adapters.py`, `sanity_checks.py`. Wired into
  the coupled steady march for the helical solver (co-flow validated;
  counter-flow has a known prescribed-outlet limitation); shell-and-tube and
  both transient solvers remain postprocess-only/unwired. Status:
  `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md` and
  `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`.
- Optimization and calibration scripts: `optimization/`.
- Existing design/calibration docs: `docs/`.

## Read These Before Editing

- `CODEBASE_MAP.md`: repository-wide module map and entry points.
- `docs/context/SOLVER_CONTEXT.md`: solver architecture, state flow, physics couplings, and common pitfalls.
- `docs/context/PHYSICS_CONTEXT.md`: correlations, chemistry, radiation, materials, and units.
- `docs/context/TRANSIENT_STATUS.md`: transient implementation status, chemistry decisions, and known fixes.
- `TESTING_CONTEXT.md`: setup and verification commands.
- `research/flamelet_kit/LLM_CONTEXT.md`: only for flamelet/PFR work.
- `optimization/LLM_CONTEXT.md`: only for calibration/optimization work.
- `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`: only for liquid/
  boiling coolant work; this feature is not yet wired into the coupled march.

## Development Rules For This Codebase

- Prefer existing dataclasses in `1Dmodel/input_data.py`; do not add global constants in solver files unless they are truly local implementation details.
- Keep unit conventions explicit. The code is SI unless comments say otherwise; pressure is Pa, temperature is K, mass flow is kg/s.
- Preserve the package bootstrap blocks in runnable scripts under `1Dmodel/`; they allow direct execution from a directory whose name starts with a digit.
- Treat `main_solve.py` as the maintained helical-coil backend. Older solver variants live in `archive/legacy_solvers/`.
- Avoid editing generated/cache files, images, zip archives, and office documents unless the task explicitly targets them.
- If changing correlations or solver behavior, add or update a narrow verification path rather than relying only on dashboard plots.

## Usual Commands

From the repository root after installing dependencies:

```powershell
python -m pip install -e .
python -m hps_combustor.main_steady
python -m hps_combustor.main_transient
python -m pytest research/flamelet_kit/tests/ -q
```

`1Dmodel/requirements.txt` lists the scientific runtime dependencies. Cantera can be installed by pip or conda depending on platform.

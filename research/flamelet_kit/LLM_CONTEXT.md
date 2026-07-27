# flamelet_kit LLM Context

## Scope

`flamelet_kit` is a standalone extraction of a Representative Interactive Flamelet style workflow plus a downstream cooling plug-flow reactor. It deliberately has no dependency on `hps_combustor`.

Use this directory when the task is about finite-rate chemistry reduction, flamelet-in-mixture-fraction space, steady-cache reuse, or cooling-PFR behavior. For ordinary combustor-HX sizing, start in `1Dmodel/` instead.

## Main Files

| File | Role |
|---|---|
| `flamelet.py` | `Flamelet` class, Z grid, Bilger stoichiometric mixture fraction, scalar dissipation profile, Strang splitting. |
| `cooling_pfr.py` | `CoolingPFR` class for a single already-mixed gas stream with wall heat loss. |
| `steady_cache.py` | Tolerance-gated cache for avoiding unnecessary re-advancement. |
| `flamelet_bank.py` | Optional multi-condition manager; most work does not need it. |
| `example_run.py` | Runnable flamelet demo. |
| `example_cooling.py` | Runnable cooling-PFR demo. |
| `tests/test_flamelet_kit.py` | Current pytest coverage. |

## Docs To Read

- `README.md`: quickstart and file map.
- `METHODOLOGY.md`: equations and modeling rationale.
- `ADAPTATION_GUIDE.md`: retargeting to different mechanisms/fuels/geometries.
- `REPRODUCE_SPEC.md`: algorithmic spec for reimplementation.

## Verification

From repository root:

```powershell
python -m pytest flamelet_kit/tests/ -q
python flamelet_kit/example_run.py
python flamelet_kit/example_cooling.py
```

The tests use Cantera `gri30.yaml`, so they do not require repository-specific mechanism files.

## Important Modeling Boundary

The flamelet model is for non-premixed two-feed-stream systems where conserved mixture fraction is meaningful. If the inlet is premixed or single-stream, use a PSR/PFR style model instead; `CoolingPFR` is already the local pattern for that.

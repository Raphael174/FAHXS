# Combustor-HX 1D Model — Project Guide

Status note: this is an older helical-focused project guide. For current solver
truth, including shell-and-tube, transient counter-flow, finite-rate FPV caching,
GOX chilldown, and the fixed-step transient integrator, read
`docs/TECHNICAL_REFERENCE.md`, `CLAUDE.md`, and
`docs/context/TRANSIENT_STATUS.md` first.

## What this code does

1D coupled heat-exchanger simulation for a **shell-and-helical-tube combustor**.  
Hot side: combustion gas from liquid-propellant combustion (diesel / gasoline / H₂ in O₂).  
Cold side: Helium coolant at high pressure (~90 bar) flowing through helical coils wound around the combustor.

The solver marches along the combustor axis, resolving at each step:
- Convective and radiative heat transfer from gas to coil wall
- 1D cylindrical conduction through the coil tube wall
- Coolant thermodynamics updated via CoolProp (real-gas Helium)
- Hot-gas chemistry updated via Cantera (equilibrium or frozen composition)
- Hoop and thermal stresses in the coil tube

## Running

```bash
cd 1Dmodel          # or whatever the folder was renamed to
pip install -r requirements.txt
python main_solve.py
```

> **Note on Cantera**: `pip install cantera` works on most platforms.  
> If it fails, use conda: `conda install -c cantera cantera`.

The entry point prints a results summary and opens matplotlib dashboard figures.

## Package bootstrap (important for contributors)

The folder `1Dmodel` is **not a valid Python identifier** (starts with a digit), so the standard `python -m package` mechanism cannot be used directly. `main_solve.py` contains a self-contained bootstrap (lines 5–22) that:

1. Registers the folder as package `_hps` in `sys.modules` using `importlib.util`.
2. Re-runs `main_solve.py` via `runpy.run_module("_hps.main_solve")` so all relative imports resolve.

This is transparent to the user: `python main_solve.py` just works.

## Key files

| File | Purpose |
|---|---|
| `main_solve.py` | Main solver class `main_solver` + `__main__` entry point |
| `input_data.py` | All configuration dataclasses (`coolantProp`, `hotgasProp`, `combustorProp`, `numericalProp`, `system_requirements`, `CorrelationCoefficients`) |
| `configs/` | Alternative input configurations (e.g. 500 kW, CFD validation study) |

## Physics modules

| Module | Key contents |
|---|---|
| `physics/combustion_chemistry/combustion_gas.py` | `combustion_gas_solve` — Cantera HP-equilibrium combustion; `choose_fuel` — YAML mechanism selector |
| `physics/heat_transfer_correlations.py` | Nu correlations: Mori & Nakayama 1967 (coil), Salimpour 2008 / Ahmed 1997 / Churchill-Bernstein (shell); dispatchers `dispatch_nu_coil`, `dispatch_nu_shell` |
| `physics/friction_correlations.py` | Friction: Ali et al. 2024 (curved pipe), Colebrook 1939; dispatcher `dispatch_friction_coil` |
| `physics/heat_conduction.py` | `OneDimensionalSteadyConduction_ShellnHelicalTube` — cylindrical wall conduction + radiation via fsolve |
| `physics/governing_equations.py` | Quasi-1D ideal-gas equations: dT/dx, dU/dx, dp/dx, dρ/dx |
| `physics/radiation_model/radiation_build.py` | `make_ehlme_backend` — builds WSGGM radiation backend from Ehlmé 2025 JSON coefficients |
| `physics/radiation_model/radiation_equations.py` | `qrad_net_mbl`, `hrad_from_q` — Stefan-Boltzmann net flux with mean beam length |

## Mechanical modules

| Module | Key contents |
|---|---|
| `mechanical/geometry/helix_geometry.py` | `HelixGeometryRadiusCST` — arc-length ↔ axial-position mapping for helix; `compute_Dh_shell` |
| `mechanical/loads.py` | `stress_pressure_tube`, `stress_thermal_tube` — Roark's formulas |
| `mechanical/material_specs/material_temperature_strength.py` | Temperature-dependent CTE, E, yield, λ for ST316L and INCO718; `init_material_temperature_properties` dispatcher |

## Data / output modules

| Module | Key contents |
|---|---|
| `model_data_process/data_processing.py` | `make_solver_data` — creates the per-run dict of all logged quantities |
| `model_data_process/data_plotting.py` | `HXDashboard` — five thematic matplotlib figures (thermal, helium, combustion, mechanical, radiation) |

## Supported fuels

| `hotgasProp.fuel` | Mechanism YAML | Notes |
|---|---|---|
| `"diesel-C16H34"` | `RenKokjohn_surrogate.yaml` | Default |
| `"POSF10325"` | `A2highT.yaml` | JET-A surrogate |
| `"gasoline-E5"` / `"gasoline-E10"` | `llnl_gasoline_323.yaml` | LLNL 323-species |
| `"H2"` | `H2-O2_Burke2012.yaml` | Hydrogen |

## Supported materials

`material_HX` / `material_CC` in `combustorProp`:
- `"ST316L"` — stainless steel 316L
- `"INCO718"` — Inconel 718

## Key design parameters (default in `input_data.py`)

- Combustor inner diameter: 115 mm
- Coil inner diameter: 14.07 mm, wall 0.5 mm
- Helium: 90 bar inlet, 0.15 kg/s, 120 K → ~750 K
- Hot gas: diesel O/F = 2.9, 70 g/s, 5 bar
- HX length: up to 691 mm axial
- Radiation: WSGGM Ehlmé 2025, enabled by default

## Chemistry model

`numericalProp.chemistry_model`:
- `"equilibrium"` — Cantera HP re-equilibrate at each step (default, most accurate)
- `"frozen"` — composition fixed after initial combustion (faster)
- `"finite_rate"` - default FPV manifold + progress-variable transport

Current chemistry override:
- `numericalProp.chemistry_model = "finite_rate"` is the default steady chemistry model.
- `transientProp.chemistry_transient = "finite_rate"` is the default transient chemistry model.
- Both finite-rate paths use the FPV manifold in
  `physics/combustion_chemistry/fpv_manifold.py`.

## Correlation selection

`combustorProp.Nusselt_shell`: `"ahmed_toroid"` | `"salimpour2008"` | `"churchill_bernstein_tightcoil"` | `"churchill_bernstein"`  
`combustorProp.Nusselt_coil`: `"mori1967"` | `"Gnielinski"`  
`combustorProp.friction_coil`: `"CurvedPipeAli2024"` | `"Colebrook1939"`

## Calibration knobs

`CorrelationCoefficients` dataclass in `input_data.py` exposes empirical prefactors for each active correlation (Salimpour a/b/c, Mori a_lo/b_lo/c_lo, Ali c_lo/c_hi, radiation mbl_factor, emissivity_wall). Pass a customised instance to `main_solver(..., corrCoeffs=CorrelationCoefficients(...))`.

## Flow configuration

`combustorProp.flow_config`:
- `"co"` — co-flow (He enters at gas inlet end)
- `"counter"` — counter-flow (He enters at gas outlet end, default is `"co"` in `input_data.py`)

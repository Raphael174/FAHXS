"""Confirmed-working water/boiling coolant example for the helical solver.

Helium/`single_phase_coolprop` remains the project's working baseline
(`input_data.py` defaults). Water/`equilibrium_liquid` is an experimental
test configuration — this module is the recommended recipe for trying it,
built as a standalone example rather than by editing the shared
`coolantProp`/`combustorProp` defaults (see
docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md).

Two entry points:

- `run_coflow()`: fast, no guessing required. Co-flow starts the coolant
  march at the real physical inlet (`T_in`, `p_in`) directly.
- `run_counterflow_physical()`: slower (up to ~40 full coil solves), but
  physically anchored on the same `T_in`/`p_in` via
  `solve_counterflow_liquid_reference()` instead of a guessed
  `T_out`/`p_out`. Prefer this over the plain (guess-based) counter-flow path
  for water — see that function's docstring for why.

Scale note: this combustor's real helical coil is ~1378 arc-length nodes
(main_solver enforces `combustorProp.HX_config == "shellnHelicalTube"` for
correct axial-length accounting — see its `__init__`), so duty is on the
order of 100-200 kW, not the ~20 kW a shorter/flat-duct assumption would
suggest. `T_in`/`p_in`/`mass_flow_c` below are picked for this real length.

Run directly:
    python -m hps_combustor.validation.water_helical_example
"""

from __future__ import annotations

from hps_combustor.input_data import (
    CorrelationCoefficients,
    combustorProp,
    coolantProp,
    hotgasProp,
    numericalProp,
    system_requirements,
)
from hps_combustor.main_solve import main_solver, solve_counterflow_liquid_reference

LUT_PATH = "docs/reference/external/2006LUTdata.txt"


def _water_coolant(mass_flow_c: float = 0.2) -> coolantProp:
    return coolantProp(
        coolant="Water",
        coolant_model="equilibrium_liquid",
        mass_flow_c=mass_flow_c,
        T_in=303.15,   # 30 degC
        p_in=82e5,     # 82 bar
        liquid_chf_lut_path=LUT_PATH,
    )


def run_coflow(mass_flow_c: float = 0.2):
    """Fast: starts directly at the physical T_in/p_in, no guessing needed.

    At the default mass_flow_c=0.2 kg/s, under the production
    (finite_rate) chemistry model, this reaches saturated boiling (outlet
    quality ~0.18) with a safe CHF margin (~3.7); it does not require
    T_out/p_out at all in co-flow. Chemistry model matters here: finite_rate
    and frozen give meaningfully different duty for this regime (per
    docs/context/TRANSIENT_STATUS.md), so a margin validated under one is not
    automatically safe under the other — mass_flow_c=0.1 (used in earlier
    frozen-chemistry validation) sits at CHF margin ~0.01 (dryout risk) under
    finite_rate.
    """
    solver = main_solver(
        coolantProp=_water_coolant(mass_flow_c),
        hotgasProp=hotgasProp(),
        combustorProp=combustorProp(HX_config="shellnHelicalTube", flow_config="co"),
        numericalProp=numericalProp(),  # production default: finite_rate chemistry
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
    )
    solver.solver()
    solver.compute_performance()
    solver._check_global()
    solver.print_summary()
    return solver


def run_counterflow_physical(mass_flow_c: float = 1.0):
    """Physically-anchored counter-flow: converges the hot-end starting
    enthalpy so the cold end matches T_in/p_in exactly, instead of requiring
    a guessed T_out/p_out (which cannot even represent a two-phase state as a
    single (T,P) pair — see solve_counterflow_liquid_reference's docstring).

    mass_flow_c=1.0 keeps this fast (smaller relative duty -> smaller
    required search range); lower it toward 0.1 for a boiling case, at the
    cost of a much longer shooting search (up to several minutes with
    finite_rate chemistry).
    """
    solver = solve_counterflow_liquid_reference(
        _water_coolant(mass_flow_c),
        hotgasProp(),
        combustorProp(HX_config="shellnHelicalTube", flow_config="counter"),
        numericalProp(),
        system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
    )
    solver.compute_performance()
    solver._check_global()
    solver.print_summary()
    print(f"Counter-flow shooting residual: {solver.counterflow_reference_residual_J_kg:.1f} J/kg")
    return solver


if __name__ == "__main__":
    run_coflow()

"""Tests for solve_counterflow_liquid_reference() (main_solve.py).

This is a follow-up to Phase 2 of
docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md: the plain liquid
counter-flow march starts from the legacy ``T_out``/``p_out`` guess, which
cannot express a genuine two-phase starting state. This shooting helper
converges the march's hot-end starting enthalpy so the cold end matches the
user's physically supplied ``T_in``/``p_in`` instead.

Scale note: this combustor's real helical coil is ~1378 arc-length nodes (see
docs/context/PHYSICS_CONTEXT.md / main_solver's HX_config guard), so full
adaptive-bracket-search-then-bisection shooting calls are slow even with
"frozen" chemistry (each full evaluation marches the whole coil, and the
search needs up to ~40 evaluations in the worst case). A single full
shooting case (``mass_flow_c=1.0``, small relative duty, fast convergence) is
used for the end-to-end test; the HEM-monotonicity check instead calls
``main_solver`` directly with ``_liquid_enthalpy_hot_end_override`` (no
search), which is fast regardless of scale.
"""
from __future__ import annotations

import numpy as np
import pytest
from CoolProp.CoolProp import PropsSI

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


def _water_coolant(T_in, p_in, mass_flow_c):
    return coolantProp(
        coolant="Water",
        coolant_model="equilibrium_liquid",
        mass_flow_c=mass_flow_c,
        T_in=T_in,
        p_in=p_in,
        liquid_chf_lut_path=LUT_PATH,
    )


def _solve(coolant, mass_flow_c=None):
    return solve_counterflow_liquid_reference(
        coolant,
        hotgasProp(),
        combustorProp(HX_config="shellnHelicalTube", flow_config="counter"),
        numericalProp(chemistry_model="frozen"),
        system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
    )


def test_shooting_converges_to_physical_cold_inlet():
    T_in, p_in, mass_flow_c = 303.15, 82e5, 1.0  # subcooled, small relative duty -> fast
    coolant = _water_coolant(T_in, p_in, mass_flow_c)
    solver = _solve(coolant)

    tol_h = 2.0e3  # default counterflow_reference_tol_J_kg in solve_counterflow_liquid_reference
    assert abs(solver.counterflow_reference_residual_J_kg) <= tol_h

    # The live post-march state (solver.enthalpy_c) is the true converged
    # endpoint (see main_solve.py's note on data_master being one node stale)
    # and should be close to the user's physical cold inlet in enthalpy terms.
    target_h = PropsSI("H", "T", T_in, "P", p_in, "Water")
    assert abs(solver.enthalpy_c - target_h) <= tol_h


def test_shooting_result_passes_liquid_sanity_gates():
    coolant = _water_coolant(303.15, 82e5, 1.0)
    solver = _solve(coolant)
    solver.compute_performance()
    solver._check_global()
    report = solver.liquid_sanity_report
    assert report.passed, report.messages


def test_hem_closure_is_monotonic_for_shooting():
    """Necessary condition for secant-method convergence: as the hot-end
    starting enthalpy increases, the resulting cold-end enthalpy must also
    increase (no non-monotonic HEM behavior that would break shooting).

    Uses ``mass_flow_c=1.0`` (small relative duty, ~213 kJ/kg of coolant
    enthalpy rise per the co-flow reference run) so the chosen guess span
    produces genuine variation instead of every guess saturating against the
    enthalpy floor (a large mass_flow_c needs only a modest starting margin).
    """
    coolant = _water_coolant(303.15, 82e5, 1.0)
    target_h = PropsSI("H", "T", 303.15, "P", 82e5, "Water")
    cp_ref = PropsSI("C", "T", 303.15, "P", 82e5, "Water")
    span = cp_ref * 150.0

    residuals = []
    for factor in (0.3, 0.6, 1.0, 1.4, 1.8, 2.2):
        solver = main_solver(
            coolant, hotgasProp(), combustorProp(HX_config="shellnHelicalTube", flow_config="counter"),
            numericalProp(chemistry_model="frozen"), system_requirements(),
            corrCoeffs=CorrelationCoefficients(),
            _liquid_enthalpy_hot_end_override=float(target_h + factor * span),
        )
        solver.solver()
        residuals.append(float(solver.enthalpy_c) - target_h)

    diffs = np.diff(residuals)
    assert np.all(diffs > 0) or np.all(diffs < 0), (
        f"HEM closure is not monotonic for shooting: residuals={residuals}"
    )


def test_non_counterflow_or_non_liquid_falls_back_to_plain_solve():
    """Co-flow or helium mode should just run the plain solver, no shooting."""
    coolant = coolantProp()  # default helium, single_phase_coolprop
    solver = solve_counterflow_liquid_reference(
        coolant, hotgasProp(), combustorProp(HX_config="shellnHelicalTube", flow_config="co"),
        numericalProp(chemistry_model="frozen"), system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
    )
    assert not hasattr(solver, "counterflow_reference_residual_J_kg")


def test_main_steady_dispatches_liquid_counterflow_to_liquid_reference():
    """main_steady.run_steady() must route liquid-mode counter-flow with
    counterflow_physical_steady_reference=True to
    solve_counterflow_liquid_reference(), not the gas-only (temperature-shot)
    solve_counterflow_physical_reference() — the latter is ill-posed inside
    the two-phase dome.
    """
    from hps_combustor.main_steady import run_steady
    from hps_combustor.input_data import shellTubeProp, runProp

    coolant = _water_coolant(303.15, 82e5, 1.0)
    inputs = {
        "coolant": coolant,
        "hotgas": hotgasProp(),
        "combustor": combustorProp(flow_config="counter"),
        "shelltube": shellTubeProp(),
        "numerical": numericalProp(
            chemistry_model="frozen", counterflow_physical_steady_reference=True
        ),
        "system": system_requirements(),
        "correlations": CorrelationCoefficients(),
        "run": runProp(),
    }
    inputs["combustor"].HX_config = "shellnHelicalTube"
    solver, _summary = run_steady(inputs)
    assert hasattr(solver, "counterflow_reference_residual_J_kg")

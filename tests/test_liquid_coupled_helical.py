"""Phase 2 of docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md.

Exercises the coupled (p,h) liquid-coolant march wired directly into
``main_solver`` (helical steady solver), as opposed to the pre-existing
postprocess-only bridge. Covers subcooled and boiling regimes, co- and
counter-flow, a grid-convergence check, and a cross-check against the
already-validated postprocess bridge fed the same converged duty profile.

Test case scale note: ``combustorProp.HX_config`` must be pinned to
"shellnHelicalTube" (main_solver now raises otherwise — see its __init__).
The real helical coil at this combustor's geometry is about 1378 arc-length
nodes (not ~100): main_solver's axial-length bookkeeping maps arc length to
axial position through the coil's actual helix geometry, and a tightly-wound
coil needs much more arc length than axial length to traverse the same
combustor length. All case parameters (T_in, p_in, mass_flow_c) below are
picked for THIS real coil length; do not reuse them assuming a short/flat
duct's worth of duty.
"""
from __future__ import annotations

import numpy as np
import pytest

from hps_combustor.input_data import (
    CorrelationCoefficients,
    combustorProp,
    coolantProp,
    hotgasProp,
    numericalProp,
    system_requirements,
)
from hps_combustor.main_solve import main_solver

LUT_PATH = "docs/reference/external/2006LUTdata.txt"


def _run(flow_config, *, T_in, p_in, T_out, p_out, mass_flow_c, n_arc_steps_per_turn=50):
    coolant = coolantProp(
        coolant="Water",
        coolant_model="equilibrium_liquid",
        mass_flow_c=mass_flow_c,
        T_in=T_in,
        p_in=p_in,
        T_out=T_out,
        p_out=p_out,
        liquid_chf_lut_path=LUT_PATH,
    )
    solver = main_solver(
        coolantProp=coolant,
        hotgasProp=hotgasProp(),
        combustorProp=combustorProp(HX_config="shellnHelicalTube", flow_config=flow_config),
        numericalProp=numericalProp(chemistry_model="frozen", N_arc_steps_per_turn=n_arc_steps_per_turn),
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
    )
    solver.solver()
    solver.compute_performance()
    solver._check_global()
    return solver


@pytest.mark.parametrize("flow_config", ["co", "counter"])
def test_subcooled_liquid_march_passes_sanity_gates(flow_config):
    solver = _run(
        flow_config,
        T_in=303.15, p_in=82e5, T_out=310.0, p_out=82e5, mass_flow_c=0.3,
    )
    report = solver.liquid_sanity_report
    assert report.passed, report.messages

    d = solver.data_master
    quality = np.asarray(d["quality_c"], dtype=float)
    assert np.all(quality < 0.0), "case is tuned to stay subcooled throughout"

    p_c = np.asarray(d["p_c"], dtype=float)
    dp = np.diff(p_c)
    assert np.all(dp <= 1.0e-6) or np.all(dp >= -1.0e-6), "pressure must be monotonic along the march"


def test_boiling_liquid_march_reaches_saturation_and_passes_gates_coflow():
    """Co-flow: the march starts subcooled at T_in and heats up, so it should
    genuinely cross into saturated boiling partway along the length.
    """
    solver = _run(
        "co", T_in=303.15, p_in=82e5, T_out=330.0, p_out=82e5, mass_flow_c=0.1,
    )
    report = solver.liquid_sanity_report
    assert report.passed, report.messages
    assert not report.dryout_risk

    d = solver.data_master
    quality = np.asarray(d["quality_c"], dtype=float)
    assert np.max(quality) > 0.0, "case is tuned to reach saturated boiling"
    assert np.max(quality) < 1.0, "case is tuned to stay short of complete vaporization"

    void = np.asarray(d["void_c"], dtype=float)
    assert np.all((void >= 0.0) & (void <= 1.0))


def test_tractable_counterflow_march_passes_gates():
    """Counter-flow structural note: the march starts from the ``T_out``/
    ``p_out`` guess as a single-phase (T,P) state (the same legacy
    "prescribed-outlet" shortcut the gas coolant march already uses — see
    docs/context/SOLVER_CONTEXT.md). A (T,P) pair cannot represent a genuine
    two-phase state, so that starting guess is a hard ceiling on enthalpy:
    the march can only move away from it (colder), never past it. Reaching
    saturated boiling *via heating* along a counter-flow march with a
    specific, physically-anchored cold inlet now uses
    ``solve_counterflow_liquid_reference()``
    (tests/test_liquid_counterflow_reference.py) instead of this MVP guess
    path. This test instead checks the achievable, meaningful behavior for
    the plain path: a modest, tractable single-phase starting guess stays
    numerically well behaved for the full coil length and passes every
    sanity gate.
    """
    solver = _run(
        "counter", T_in=303.15, p_in=82e5, T_out=330.0, p_out=82e5, mass_flow_c=1.0,
    )
    report = solver.liquid_sanity_report
    assert report.passed, report.messages
    assert not report.dryout_risk

    d = solver.data_master
    quality = np.asarray(d["quality_c"], dtype=float)
    assert np.all(quality < 0.0), "case is tuned to stay subcooled throughout"

    void = np.asarray(d["void_c"], dtype=float)
    assert np.all((void >= 0.0) & (void <= 1.0))


def test_boiling_liquid_march_grid_convergence():
    """Outlet quality should not swing wildly under arc-step refinement."""
    coarse = _run(
        "co", T_in=303.15, p_in=82e5, T_out=330.0, p_out=82e5, mass_flow_c=0.1,
        n_arc_steps_per_turn=50,
    )
    fine = _run(
        "co", T_in=303.15, p_in=82e5, T_out=330.0, p_out=82e5, mass_flow_c=0.1,
        n_arc_steps_per_turn=100,
    )
    x_out_coarse = coarse.data_master["quality_c"][-1]
    x_out_fine = fine.data_master["quality_c"][-1]
    assert abs(x_out_fine - x_out_coarse) < 0.05, (
        f"outlet quality not grid-converged: coarse={x_out_coarse:.4f}, fine={x_out_fine:.4f}"
    )


def test_coupled_march_matches_postprocess_bridge_on_same_duty_coflow():
    """The coupled march and the pre-existing postprocess bridge share the same
    (p,h) closure; fed the SAME converged duty profile, they must agree.

    Co-flow only: both the coupled march and the postprocess bridge start
    from ``coolantProp.T_in``/``p_in`` in co-flow, so they share the same
    physical anchor. Plain counter-flow is excluded here — it still uses the
    legacy "prescribed-outlet" shortcut, starting from the ``T_out``/``p_out``
    *guess*, while the postprocess bridge always anchors on the physically
    correct ``T_in``. That is the same pre-existing discrepancy already
    documented for the gas march in docs/context/SOLVER_CONTEXT.md
    ("Helical counter-flow validation note"), not a new liquid-path defect;
    ``solve_counterflow_liquid_reference()``
    (tests/test_liquid_counterflow_reference.py) resolves it for counter-flow
    by shooting on the physical T_in/p_in instead.
    """
    solver = _run(
        "co",
        T_in=303.15, p_in=82e5, T_out=330.0, p_out=82e5, mass_flow_c=0.1,
    )
    bridge_result = solver.liquid_coolant_postprocess(lut_path=LUT_PATH)

    coupled_p = np.asarray(solver.data_master["p_c"], dtype=float)
    coupled_h = np.asarray(solver.data_master["enthalpy_c"], dtype=float)
    # The coupled march records the state at the START of each step (N points
    # for N steps); the postprocess bridge is edge-based (N+1 points, adding
    # the final outlet edge the march advances to but never re-appends).
    # Compare the coupled march against the bridge's matching leading edges.
    bridge_p = bridge_result.node_fields_hx_order["p_Pa"][:-1]
    bridge_h = bridge_result.node_fields_hx_order["h_J_kg"][:-1]

    assert coupled_p.shape == bridge_p.shape
    np.testing.assert_allclose(coupled_p, bridge_p, rtol=1.0e-2)
    np.testing.assert_allclose(coupled_h, bridge_h, rtol=1.0e-2)

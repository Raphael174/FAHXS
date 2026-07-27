"""Phase 0 of docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md.

Freezes current helium/gas steady-solver behavior (helical + shell-and-tube,
co + counter flow) BEFORE any liquid-coolant coupling work touches
`main_solve.py` / `main_solve_shellntube.py`. Every later phase of that plan
must re-run this file unchanged; any drift on the gas path is a defect.

`numericalProp.chemistry_model="frozen"` is used deliberately, not as the
production default (which is `finite_rate`): frozen mode needs no FPV
manifold build/cache and no per-node Cantera equilibrate calls, so this
fixture is fast (~1 s/case) and its results are exactly reproducible run to
run. It still exercises the full coupled march (coolant properties,
correlations, wall conduction, governing-equation integration) that the
liquid-coolant work will modify.

Shell-and-tube uses a reduced `N_axial=40` (default production grid is 200)
purely so this regression fixture stays fast; this is not a grid-convergence
claim, only a frozen baseline for catching unintended behavior changes.

IMPORTANT: every fluid/geometry field this fixture depends on is pinned
explicitly (``coolant``, ``coolant_model``, ``mass_flow_c``, ``T_in``/``p_in``,
``T_out``/``p_out``, ``combustorProp.HX_config``) rather than left to
dataclass defaults in ``input_data.py``. Those defaults are shared, mutable
project state — they get changed for day-to-day experimentation (e.g.
switching the default coolant to water for a specific run), and a baseline
fixture that silently inherits "whatever the current default happens to be"
is not a baseline at all. This file exists specifically to catch drift in the
gas/helium march; it must not itself drift when someone edits
``input_data.py`` defaults for an unrelated reason.
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
    shellTubeProp,
    system_requirements,
)
from hps_combustor.main_solve import main_solver
from hps_combustor.main_solve_shellntube import shellntube_solver

# Recorded 2026-07-13 on the pre-liquid-coupling gas/helium march, WITH
# combustorProp.HX_config explicitly pinned to "shellnHelicalTube". An earlier
# version of this fixture never pinned HX_config and silently inherited
# whatever combustorProp's dataclass default happened to be; at the time that
# default was "shellntube" (main_solve_shellntube.py's own label, not a valid
# main_solver configuration), so main_solver's axial-length bookkeeping
# (_advance_state) used the wrong linear-dx branch instead of the true
# helical arc-length-to-axial-position mapping, giving a ~14x-too-short
# effective coil length (98 nodes instead of the correct 1378). main_solver
# now raises if HX_config != "shellnHelicalTube" (see its __init__), so this
# mistake cannot recur silently. These are the values under the CORRECT
# mapping; confirmed reproducible across two independent runs before being
# hardcoded.
HELICAL_BASELINE = {
    "co": {
        "Q_tot_kW": 199.29188128198106,
        "dp_c_bar": 8.539605328052556,
        "T_wg_max": 573.7031995127154,
        "T_wc_max": 561.5004140617833,
        "T_c_end0": 303.15,
        "T_c_end1": 556.8407305781116,
        "n_nodes": 1378,
    },
    "counter": {
        "Q_tot_kW": 184.49205653434237,
        "dp_c_bar": 33.99384935631357,
        "T_wg_max": 725.357704728729,
        "T_wc_max": 699.5563916496604,
        "T_c_end0": 650.0,
        "T_c_end1": 522.4148448999435,
        "n_nodes": 1378,
    },
}

SHELLTUBE_BASELINE = {
    "co": {
        "Q_tot_kW": 383.3146228927451,
        "T_g_out": 1066.0845376319512,
        "T_c_out": 795.5142039230026,
        "T_wg_max": 845.3530512957477,
        "T_wc_max": 811.3511750979103,
        "n_sweeps": 13,
    },
    "counter": {
        "Q_tot_kW": 425.77199172931074,
        "T_g_out": 800.6837238189005,
        "T_c_out": 850.097768422657,
        "T_wg_max": 1172.6332668374746,
        "T_wc_max": 1088.9553057263352,
        "n_sweeps": 14,
    },
}

REL_TOL = 1e-6


def _helium_coolant():
    """Pinned helium coolant matching the values this fixture's baseline
    numbers were recorded against — see the module docstring on why this must
    not be left to ``coolantProp()`` defaults.
    """
    return coolantProp(
        coolant="Helium",
        coolant_model="single_phase_coolprop",
        mass_flow_c=150e-3,
        T_in=30 + 273.15,
        T_out=650,
        p_in=82e5,
        p_out=13e5,
    )


def _run_helical(flow_config):
    solver = main_solver(
        coolantProp=_helium_coolant(),
        hotgasProp=hotgasProp(),
        combustorProp=combustorProp(HX_config="shellnHelicalTube", flow_config=flow_config),
        numericalProp=numericalProp(chemistry_model="frozen"),
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
    )
    solver.solver()
    summary = solver.compute_performance()
    return solver, summary


def _run_shelltube(flow_config):
    solver = shellntube_solver(
        coolantProp=_helium_coolant(),
        hotgasProp=hotgasProp(),
        shellTubeProp=shellTubeProp(),
        numericalProp=numericalProp(chemistry_model="frozen"),
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
        N_axial=40,
        flow_config=flow_config,
    )
    solver.solve(verbose=False)
    return solver


@pytest.mark.parametrize("flow_config", ["co", "counter"])
def test_helical_helium_baseline_unchanged(flow_config):
    solver, summary = _run_helical(flow_config)
    expected = HELICAL_BASELINE[flow_config]

    assert summary["Q_tot_kW"] == pytest.approx(expected["Q_tot_kW"], rel=REL_TOL)
    assert summary["dp_c_bar"] == pytest.approx(expected["dp_c_bar"], rel=REL_TOL)
    assert summary["T_wg_max"] == pytest.approx(expected["T_wg_max"], rel=REL_TOL)
    assert summary["T_wc_max"] == pytest.approx(expected["T_wc_max"], rel=REL_TOL)
    assert float(solver.data_master["T_c"][0]) == pytest.approx(expected["T_c_end0"], rel=REL_TOL)
    assert float(solver.data_master["T_c"][-1]) == pytest.approx(expected["T_c_end1"], rel=REL_TOL)
    assert len(solver.data_master["T_c"]) == expected["n_nodes"]


@pytest.mark.parametrize("flow_config", ["co", "counter"])
def test_shelltube_helium_baseline_unchanged(flow_config):
    solver = _run_shelltube(flow_config)
    expected = SHELLTUBE_BASELINE[flow_config]

    assert float(solver.Q_tot / 1e3) == pytest.approx(expected["Q_tot_kW"], rel=REL_TOL)
    assert float(solver.tube["T_g_out"]) == pytest.approx(expected["T_g_out"], rel=REL_TOL)
    assert float(solver.T_c_out) == pytest.approx(expected["T_c_out"], rel=REL_TOL)
    assert float(np.max(solver.tube["T_wg"])) == pytest.approx(expected["T_wg_max"], rel=REL_TOL)
    assert float(np.max(solver.tube["T_wc"])) == pytest.approx(expected["T_wc_max"], rel=REL_TOL)
    assert int(solver.n_sweeps) == expected["n_sweeps"]

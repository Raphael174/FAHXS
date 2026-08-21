"""Finite-rate companion to `test_steady_baseline_regression.py`.

That file deliberately uses `numericalProp.chemistry_model="frozen"` for
speed/determinism (no FPV manifold build, no per-node Cantera calls) and
freezes the pre-liquid-coupling gas/helium march for
`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` Phase 0. Frozen
chemistry is validation-only, though: per `CLAUDE.md`, finite-rate FPV is the
required default for the project's actual diesel/O2 high-heat-extraction
regime (recombination freeze-out during rapid cooldown makes frozen chemistry
physically wrong there). This file freezes the same four cases (helical +
shell-and-tube, co + counter) under `chemistry_model="finite_rate"` instead,
so a change that only breaks the production chemistry path doesn't slip
through unnoticed just because the frozen fixture still passes.

Recorded 2026-07-27, same pinned helium coolant/geometry as the frozen
fixture (see `_helium_coolant()` there for why every field is pinned
explicitly rather than left to `input_data.py` defaults). Confirmed
bit-identical across two independent runs before being hardcoded. Uses the
project's disk-cached FPV manifold (`cache/fpv_manifolds/`) — first run after
clearing that cache will be slower (~30 s) than subsequent runs (~1 s), per
`docs/context/TRANSIENT_STATUS.md`; both give the same results since the
manifold build is deterministic for fixed inlet composition/pressure.

If this file starts failing after intentional physics changes, re-derive the
expected values the same way (`git diff` to confirm no accidental default
drift, run twice, confirm bit-identical, update the numbers) rather than
loosening the tolerance.
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

HELICAL_BASELINE = {
    "co": {
        "Q_tot_kW": 254.93630294392477,
        "dp_c_bar": 9.410063950597877,
        "T_wg_max": 678.3332649967667,
        "T_wc_max": 660.4150140972256,
        "T_c_end0": 303.15,
        "T_c_end1": 627.4265850497876,
        "n_nodes": 1378,
    },
    "counter": {
        "Q_tot_kW": 242.41325666218154,
        "dp_c_bar": 33.38192981570719,
        "T_wg_max": 765.0102351902142,
        "T_wc_max": 738.3580037115677,
        "T_c_end0": 650.0,
        "T_c_end1": 449.46325585437785,
        "n_nodes": 1378,
    },
}

# Re-derived 2026-08-20 after three shell-side corrections, all of which move
# the shell-and-tube result and none of which touch the helical solver:
#   1. Bell-Delaware's Sieder-Tate (mu_b/mu_w)^0.14 term enabled -- it had been
#      left at its neutral default of 1.0, i.e. the correlation's own
#      property-variation correction was switched off.
#   2. Shell-to-baffle leakage area S_sb: a spurious factor of 1/2 removed
#      (S_sb = Ds*Lsb*(pi - theta_ds/2), the gap times the arc the baffle edge
#      follows). J_l moves 0.412 -> 0.387.
#   3. The pressure march moved onto Bell-Delaware plus a Chisholm two-phase
#      multiplier and a momentum (acceleration) term. Inert for this
#      single-phase helium fixture, but recorded here for provenance.
# Net effect on this fixture: Q_tot -0.51 % (co), -0.31 % (counter).
#
# Separately, the shared `input_data.py` default `mass_flow_g` was edited
# 0.100 -> 0.075 kg/s on 2026-08-20 for an unrelated LN2 run, which had shifted
# every value here by ~24 %. `_hotgas()` now pins it so these baselines cannot
# be moved again by a default edit. Confirmed bit-identical across two
# independent runs before being hardcoded, per this file's docstring.
SHELLTUBE_BASELINE = {
    "co": {
        "Q_tot_kW": 499.47644015210034,
        "T_g_out": 1479.9909108593774,
        "T_c_out": 944.7583890730842,
        "T_wg_max": 983.2378822091334,
        "T_wc_max": 973.4573800514054,
        "n_sweeps": 13,
    },
    "counter": {
        "Q_tot_kW": 548.6658765837284,
        "T_g_out": 1336.962142972538,
        "T_c_out": 1007.9889338275057,
        "T_wg_max": 1332.2504482812515,
        "T_wc_max": 1251.2492620816995,
        "n_sweeps": 14,
    },
}

REL_TOL = 1e-6


def _helium_coolant():
    """Same pinned helium coolant as test_steady_baseline_regression.py -
    keep these in sync if that fixture's pin ever changes, so the two
    baselines stay directly comparable case-for-case."""
    return coolantProp(
        coolant="Helium",
        coolant_model="single_phase_coolprop",
        mass_flow_c=150e-3,
        T_in=30 + 273.15,
        T_out=650,
        p_in=82e5,
        p_out=13e5,
    )


def _hotgas():
    """Pin the hot-gas inlet explicitly, for the same reason `_helium_coolant()`
    pins the coolant: these baselines must be immune to edits of the shared
    `input_data.py` defaults. Leaving this to `hotgasProp()` meant a change to
    the default `mass_flow_g` silently moved every recorded value (this happened
    2026-08-20, when the default was edited 0.100 -> 0.075 kg/s for an unrelated
    LN2 run and shifted the shelltube duty by -24%).
    """
    return hotgasProp(mass_flow_g=100e-3)


def _run_helical(flow_config):
    solver = main_solver(
        coolantProp=_helium_coolant(),
        hotgasProp=_hotgas(),
        combustorProp=combustorProp(HX_config="shellnHelicalTube", flow_config=flow_config),
        numericalProp=numericalProp(chemistry_model="finite_rate"),
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
    )
    solver.solver()
    summary = solver.compute_performance()
    return solver, summary


def _run_shelltube(flow_config):
    solver = shellntube_solver(
        coolantProp=_helium_coolant(),
        hotgasProp=_hotgas(),
        shellTubeProp=shellTubeProp(),
        numericalProp=numericalProp(chemistry_model="finite_rate"),
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
        N_axial=40,
        flow_config=flow_config,
    )
    solver.solve(verbose=False)
    return solver


@pytest.mark.parametrize("flow_config", ["co", "counter"])
def test_helical_helium_finite_rate_baseline_unchanged(flow_config):
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
def test_shelltube_helium_finite_rate_baseline_unchanged(flow_config):
    solver = _run_shelltube(flow_config)
    expected = SHELLTUBE_BASELINE[flow_config]

    assert float(solver.Q_tot / 1e3) == pytest.approx(expected["Q_tot_kW"], rel=REL_TOL)
    assert float(solver.tube["T_g_out"]) == pytest.approx(expected["T_g_out"], rel=REL_TOL)
    assert float(solver.T_c_out) == pytest.approx(expected["T_c_out"], rel=REL_TOL)
    assert float(np.max(solver.tube["T_wg"])) == pytest.approx(expected["T_wg_max"], rel=REL_TOL)
    assert float(np.max(solver.tube["T_wc"])) == pytest.approx(expected["T_wc_max"], rel=REL_TOL)
    assert int(solver.n_sweeps) == expected["n_sweeps"]

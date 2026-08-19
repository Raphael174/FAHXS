"""Stage A of docs/solver_design/FV_CORE_REWORK_PLAN.md.

Freezes current transient-solver behavior BEFORE `main_solve_transient.py` /
`main_solve_shellntube_transient.py` are rewired off inline `CP.PropsSI` calls
onto the `core/thermo.py` `IdealGasBackend`. Unlike the steady solvers, no
pre-existing pytest test numerically exercises these two files (the only test
that imports them, `test_main_transient_dispatch.py`, mocks the solver classes
out entirely) -- this file closes that gap and doubles as the Stage A
bit-identical gate for the transient rewire.

Two verification mechanisms, not one, because the raw call sites split across
two different code paths per file:

1. `quasi_steady` fluid_model (the legacy `solve_transient()` path) is
   exercised end-to-end via `main_transient.run_transient`. This is the ONLY
   code path `main_solve_shellntube_transient.py`'s 9 call sites are ever
   reached through (its `solve_transient_core` delegates entirely to
   `transient_core.adapters_shelltube`, a separate, already dispatch-routed
   module -- confirmed by reading the code, not assumed). It also reaches most
   of `main_solve_transient.py`'s call sites (`_march_fluids`/
   `_relax_counter_flow`), which both `solve_transient` and
   `solve_transient_core` share.

2. A handful of `main_solve_transient.py` (helical) call sites live in helper
   methods reached ONLY from `solve_transient_core`
   (`_helical_nominal_coolant_dp`, `_helical_lumped_resistance_over_density`,
   `_helical_face_resistance_over_density` -- the low-Mach coolant-momentum
   closures). `solve_transient_core` itself could not be used for an
   end-to-end smoke gate here: at the time this file was written it raised a
   `FloatingPointError` ("coolant internal energy left the configured CoolProp
   temperature range") on every short/coarse configuration tried, for both
   helical and shell-and-tube. Root-caused 2026-08-18 (see
   `tests/test_transient_core_shelltube.py`): a genuine CFL instability in the
   explicit forward-Euler coolant mass/energy advection
   (`transient_core/compressible_coolant.py`), not a bug in the guard itself --
   any macro step exceeding roughly one cell's coolant residence time blows up.
   The **shell-and-tube** path (`transient_core/adapters_shelltube.py`) now
   self-subdivides such steps (`_cfl_stable_substep_count`) and is fixed and
   tested there. The **helical** path (`main_solve_transient.py`'s own inline
   mass/energy loop -- separate code, not the shared `adapters_shelltube.py`
   kernel) has the same root cause but was NOT fixed here (out of scope; only
   the shell-and-tube transient is a Stage D migration target). So this file
   still calls those three helper methods directly rather than through a live
   `solve_transient_core` run, the same way `test_core_thermo.py` pins
   getter-vs-raw-call equivalence.

Runtimes recorded while capturing these baselines: helical quasi_steady ~9-10
s/case (dominated by FPV manifold + chemistry setup), shell-and-tube
quasi_steady ~1-2 s/case, helper-method checks <1 s (construction only, no
march).
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
    transientProp,
)
from hps_combustor.main_solve_transient import transient_solver
from hps_combustor.main_transient import build_inputs, run_transient

REL_TOL = 1e-6

# Recorded 2026-08-18 on the pre-rewire (raw CP.PropsSI) transient solvers,
# short smoke configuration: t_end=1.0 s, max_step=0.25 s, n_save=5, pinned
# Helium coolant (T_in=303.15 K, p_in=82e5 Pa) -- see `_helium_coolant()`.
HELICAL_QUASI_STEADY_BASELINE = {
    "co": {
        "T_c_out_final": 355.27149418373165,
        "T_g_out_final": 2415.759160793546,
        "Q_hot_kW_final": 282.1702150712807,
        "Q_cold_kW_final": 16.406086590408105,
        "dT_wall_max_final": 33.6587881696837,
    },
    "counter": {
        # These 4 values shift by ~1e-5 relative between a cold and a warm
        # `cache/fpv_manifolds` cache (the FPV manifold build vs. disk-cache
        # round-trip introduces a tiny floating-point difference) -- pinned to
        # the warm-cache values, reproducible across repeated pytest and
        # standalone-script runs alike once the cache is warm, which is the
        # steady state this suite normally runs in.
        "T_c_out_final": 431.9752150223438,
        "T_g_out_final": 2450.8022705531757,
        "Q_hot_kW_final": 272.9130714790373,
        "Q_cold_kW_final": 51.05866191534962,
        "dT_wall_max_final": 23.593681405729228,
    },
}

SHELLTUBE_QUASI_STEADY_BASELINE = {
    "co": {
        "T_c_out_final": 348.28587434045835,
        "T_g_out_final": 1336.962142972538,
        "Q_hot_kW_final": 637.1449039416051,
        "Q_cold_kW_final": 14.248566057063874,
    },
    "counter": {
        "T_c_out_final": 484.71800538362197,
        "T_g_out_final": 1336.962142972538,
        "Q_hot_kW_final": 638.9758178828876,
        "Q_cold_kW_final": 153.63083512126568,
    },
}

# Helical solve_transient_core-exclusive helper methods (low-Mach coolant
# momentum closures), evaluated directly at representative synthetic
# arguments -- see module docstring for why an end-to-end solve_transient_core
# run isn't used to gate these.
HELICAL_CORE_HELPER_BASELINE = {
    "nominal_dp": 104866.3907628565,
    "lumped_resistance": 77670654.34323278,
    "face_resistance_sum": 73863331.98895773,
    "face_resistance_first": 461645.8249309858,
    "face_resistance_last": 461645.8249309858,
}


def _helium_coolant():
    return coolantProp(coolant="Helium", coolant_model="single_phase_coolprop", T_in=303.15, p_in=82e5)


def _run_quasi_steady(hx_config, flow_config, n_axial=None):
    inputs = build_inputs()
    inputs["coolant"] = _helium_coolant()
    inputs["combustor"].HX_config = hx_config
    inputs["combustor"].flow_config = flow_config
    inputs["transient"].fluid_model = "quasi_steady"
    inputs["transient"].t_end = 1.0
    inputs["transient"].max_step = 0.25
    inputs["transient"].n_save = 5
    if n_axial is not None:
        inputs["run"].shelltube_transient_nodes = n_axial
    _solver, summary = run_transient(inputs)
    return summary


@pytest.mark.parametrize("flow_config", ["co", "counter"])
def test_helical_transient_quasi_steady_baseline_unchanged(flow_config):
    summary = _run_quasi_steady("shellnHelicalTube", flow_config)
    expected = HELICAL_QUASI_STEADY_BASELINE[flow_config]
    for key, value in expected.items():
        assert summary[key] == pytest.approx(value, rel=REL_TOL), key


@pytest.mark.parametrize("flow_config", ["co", "counter"])
def test_shelltube_transient_quasi_steady_baseline_unchanged(flow_config):
    summary = _run_quasi_steady("shellntube", flow_config, n_axial=20)
    expected = SHELLTUBE_QUASI_STEADY_BASELINE[flow_config]
    for key, value in expected.items():
        assert summary[key] == pytest.approx(value, rel=REL_TOL), key


def test_helical_transient_core_helper_methods_unchanged():
    """Directly gates the raw-PropsSI call sites inside the three low-Mach
    coolant-momentum helpers that only `solve_transient_core` reaches (see
    module docstring)."""
    solver = transient_solver(
        coolantProp=_helium_coolant(),
        hotgasProp=hotgasProp(),
        combustorProp=combustorProp(HX_config="shellnHelicalTube", flow_config="counter"),
        numericalProp=numericalProp(),
        system_requirements=system_requirements(),
        transientProp=transientProp(),
        corrCoeffs=CorrelationCoefficients(),
    )

    nominal_dp = solver._helical_nominal_coolant_dp(T=350.0, p=80e5, mdot_total=0.06)
    assert nominal_dp == pytest.approx(HELICAL_CORE_HELPER_BASELINE["nominal_dp"], rel=REL_TOL)

    lumped_resistance = solver._helical_lumped_resistance_over_density(
        temperature=np.full(solver.N, 350.0),
        pressure_profile=np.full(solver.N, 80e5),
        density=np.full(solver.N, 4.0),
        mdot=0.06,
    )
    assert lumped_resistance == pytest.approx(
        HELICAL_CORE_HELPER_BASELINE["lumped_resistance"], rel=REL_TOL
    )

    face_resistance = solver._helical_face_resistance_over_density(
        temperature=np.full(solver.N, 350.0),
        pressure_profile=np.full(solver.N, 80e5),
        density_face=np.full(solver.N + 1, 4.0),
        mdot_reference=0.06,
        mdot_floor=1e-4,
    )
    assert float(np.sum(face_resistance)) == pytest.approx(
        HELICAL_CORE_HELPER_BASELINE["face_resistance_sum"], rel=REL_TOL
    )
    assert float(face_resistance[0]) == pytest.approx(
        HELICAL_CORE_HELPER_BASELINE["face_resistance_first"], rel=REL_TOL
    )
    assert float(face_resistance[-1]) == pytest.approx(
        HELICAL_CORE_HELPER_BASELINE["face_resistance_last"], rel=REL_TOL
    )

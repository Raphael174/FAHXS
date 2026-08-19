"""Regression test + fix gate for the helical `solve_transient_core` CFL
instability found and fixed 2026-08-19.

Companion to `tests/test_transient_core_shelltube.py` (2026-08-18), which
fixed the SAME root cause in the shell-and-tube kernel
(`transient_core/adapters_shelltube.py`). At the time that fix landed, the
helical transient core (`main_solve_transient.py`'s own inline mass/energy
loop -- separate code, not the shared `adapters_shelltube.py` kernel) had
the identical `FloatingPointError: coolant internal energy left the
configured CoolProp temperature range` crash signature but was deliberately
left unfixed (out of scope for that session). Fixed here using the exact
same mechanism: `_cfl_stable_substep_count` (defined once in
`transient_core/adapters_shelltube.py`, imported into `main_solve_transient.py`
alongside the other private helpers it already borrows from that module)
subdivides a macro step into CFL-safe substeps for the coolant mass/energy
advection only -- momentum, the hot-gas march (`_march_fluids`), and the
wall/coolant conductance stay frozen across substeps, matching this
method's own documented "hot gas remains quasi-steady through the fixed
step" assumption.

Verified against the exact config that crashed earlier in this session
(`shellnHelicalTube`, counter-flow, `max_step=0.25`, the config `tests/
test_transient_baseline_regression.py`'s helical cases already exercise via
the LEGACY `quasi_steady` path): now completes in ~56s with energy residual
at true machine precision (~7e-10 J), vs. an immediate crash before the fix.
Also verified for `coolant_momentum_model="low_mach"` (same shared kernel).

This file mirrors `test_transient_core_shelltube.py`'s structure and intent
-- see that file for the fuller root-cause writeup (identical mechanism,
found and explained there first).
"""
from __future__ import annotations

import CoolProp.CoolProp as CP
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

REL_TOL = 1e-6

# Recorded 2026-08-19 with the CFL subcycling fix in place. Small/fast
# fixture (N nodes from default helical geometry, t_end=0.1s) chosen
# specifically so `_cfl_stable_substep_count` returns > 1 every macro step
# (verified: 38 substeps/step at max_step=0.025s) -- deliberately NOT the
# full-scale config that crashed earlier in this session (that one takes
# ~56s with the fix in place; fine for a one-off verification, too slow for
# a routine regression test).
EXPECTED = {
    "T_c_out_final": 316.06849368680116,
    "T_wall_max_final": 322.0310283718078,
    "Q_hot_kW_final": 249.28786873531595,
    "Q_cold_kW_final": 11.3099348986242,
}


def _build_solver(*, momentum_model="quasi_steady"):
    tp = transientProp()
    tp.fluid_model = "transient_coolant"
    tp.coolant_momentum_model = momentum_model
    tp.t_end = 0.1
    tp.max_step = 0.025
    tp.n_save = 3
    tp.schedule_mass_flow_c = ((0.0, 0.06),)
    tp.schedule_mass_flow_g = ((0.0, 0.08),)
    return transient_solver(
        coolantProp=coolantProp(
            coolant="Helium", coolant_model="single_phase_coolprop", T_in=303.15, p_in=82e5
        ),
        hotgasProp=hotgasProp(),
        combustorProp=combustorProp(HX_config="shellnHelicalTube", flow_config="counter"),
        numericalProp=numericalProp(),
        system_requirements=system_requirements(),
        transientProp=tp,
        corrCoeffs=CorrelationCoefficients(),
    )


def test_solve_transient_core_completes_without_cfl_blowup():
    """The regression itself: this used to raise FloatingPointError at the
    full-scale config that crashed earlier in this session. Must now
    complete cleanly on this fast fixture too."""
    solver = _build_solver()
    solver.solve_transient_core(verbose=False)  # raises on failure; that IS the test


def test_low_mach_momentum_variant_also_fixed():
    """Same shared kernel, different momentum closure -- must be fixed too."""
    solver = _build_solver(momentum_model="low_mach")
    solver.solve_transient_core(verbose=False)


def test_energy_closure_near_machine_precision():
    solver = _build_solver()
    solver.solve_transient_core(verbose=False)
    residual = solver.time_series["scalars"]["energy_residual_J"]
    assert float(max(abs(residual))) < 1e-6


def test_coolant_temperature_stays_in_coolprop_valid_range():
    solver = _build_solver()
    solver.solve_transient_core(verbose=False)
    t_min = CP.PropsSI("Tmin", "Helium")
    t_max = CP.PropsSI("Tmax", "Helium")
    T_field = solver.time_series["fields"]["T_c"]
    assert float(T_field.min()) > t_min
    assert float(T_field.max()) < t_max


def test_headline_scalars_pinned():
    solver = _build_solver()
    solver.solve_transient_core(verbose=False)
    summary = {f"{k}_final": float(v[-1]) for k, v in solver.time_series["scalars"].items()}
    for key, value in EXPECTED.items():
        assert summary[key] == pytest.approx(value, rel=REL_TOL), key


def test_subcycling_actually_engages_for_this_fixture():
    """Guard against this test silently becoming a no-op check."""
    import hps_combustor.main_solve_transient as mst

    original = mst._cfl_stable_substep_count
    seen = []

    def spy(mass, faces, dt, **kwargs):
        n = original(mass, faces, dt, **kwargs)
        seen.append(n)
        return n

    mst._cfl_stable_substep_count = spy
    try:
        _build_solver().solve_transient_core(verbose=False)
    finally:
        mst._cfl_stable_substep_count = original

    assert seen, "subcycling helper was never called"
    assert all(n > 1 for n in seen), f"expected subcycling every step, got {seen}"

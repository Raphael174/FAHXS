"""Regression test + fix gate for the shell-and-tube `solve_transient_core`
CFL instability found and fixed 2026-08-18.

Background: this path (`transient_core/adapters_shelltube.py::
run_shelltube_transient_core`) is the Stage D migration target of
`docs/solver_design/FV_CORE_REWORK_PLAN.md` ("migrate the shell-and-tube
transient first"). Its own documented validation matrix
(`1Dmodel/validation/transient_core_short_runs.py`, results stored in
`docs/validation/transient_core_short_run_results.json` from 2026-07-10) had
silently regressed to a crash: `enforce_internal_energy_bounds` raised
`FloatingPointError` ("coolant internal energy left the configured CoolProp
temperature range") on every documented case.

Root cause, established by direct probing (not assumed): the coolant
mass/energy update (`transient_core/compressible_coolant.py::
conservative_mass_energy_step`) is explicit forward Euler in the conserved
variables. It is unconditionally unstable once a macro step exceeds roughly
one cell's residence time (`mass / mdot`). Confirmed empirically:
`max_step`/`tau` ~0.21 -> stable and matches the July reference to <0.6%;
~0.43 -> a fast-growing single-cell spike (traced cell-by-cell: a `T` bump at
the inlet cell amplifies and travels downstream over successive steps,
the textbook signature of violating an explicit-upwind CFL limit); the
documented cases run at `max_step`/`tau` ~50, which explains the crash. This
is NOT a CFL-sweep coincidence -- it reproduces at every `max_step` from 0.25
down to 0.005 (50x), ruling out "just take smaller documented steps" as a fix,
and is unrelated to Stage A (`solve_transient_core` calls none of the methods
touched there).

Fix: `adapters_shelltube.py::_cfl_stable_substep_count` subdivides a macro
step into CFL-safe substeps for the mass/energy advection only (momentum,
hot-gas march, and wall/coolant conductance stay frozen across substeps --
the same quasi-steady-per-macro-step assumption `solve_transient_core`
already documents for the hot side). The outer `t` grid callers see is
unchanged; diagnostics are summed over substeps.

This file is the gate: run a small, fast (~1s) config chosen so subcycling is
exercised (n_sub > 1, verified while writing this test) but the case stays
cheap, and assert (a) it completes without the guard tripping, (b) energy
closure is near machine precision (not merely "small"), (c) the reconstructed
coolant temperature never left the fluid's real CoolProp range, (d) pinned
headline scalars for future-regression detection.
"""
from __future__ import annotations

import CoolProp.CoolProp as CP
import pytest

from hps_combustor.input_data import (
    CorrelationCoefficients,
    coolantProp,
    hotgasProp,
    numericalProp,
    shellTubeProp,
    system_requirements,
    transientProp,
)
from hps_combustor.main_solve_shellntube_transient import shellntube_transient_solver

REL_TOL = 1e-6

# Recorded 2026-08-18 with the CFL subcycling fix in place. Small/fast fixture
# (N_axial=8, t_end=0.1s, constant schedules) chosen specifically so
# `_cfl_stable_substep_count` returns > 1 every macro step (verified: 9
# substeps/step at max_step=0.025s) -- this is deliberately NOT one of the
# four documented validation-matrix cases (those take 35-140s each with the
# fix in place, since max_step there is ~50-200x the coolant residence time;
# fine for a one-off validation run, too slow for a routine regression test).
EXPECTED = {
    "T_c_out_final": 334.5000446839529,
    "T_wall_max_final": 363.2578918561147,
    "Q_hot_kW_final": 479.9223867444702,
    "Q_cold_kW_final": 30.732995968519887,
}


def _build_solver():
    tp = transientProp()
    tp.fluid_model = "transient_coolant"
    tp.t_end = 0.1
    tp.max_step = 0.025
    tp.n_save = 3
    tp.schedule_mass_flow_c = ((0.0, 0.15),)
    tp.schedule_mass_flow_g = ((0.0, 0.08),)
    return shellntube_transient_solver(
        coolantProp=coolantProp(
            coolant="Helium", coolant_model="single_phase_coolprop", T_in=303.15, p_in=80e5
        ),
        hotgasProp=hotgasProp(),
        shellTubeProp=shellTubeProp(),
        numericalProp=numericalProp(),
        system_requirements=system_requirements(),
        transientProp=tp,
        corrCoeffs=CorrelationCoefficients(),
        N_axial=8,
        flow_config="counter",
    )


def test_solve_transient_core_completes_without_cfl_blowup():
    """The regression itself: this used to raise FloatingPointError on every
    documented configuration. Must now complete cleanly."""
    solver = _build_solver()
    solver.solve_transient_core(verbose=False)  # raises on failure; that IS the test


def test_energy_closure_near_machine_precision():
    """Distinguishes a real fix from a wider guard band: the OLD (broken)
    code's own validation JSON reported energy_residual_abs_max ~30-93 J as
    "passing" -- that was itself the unresolved CFL error, not closure. A
    correct fix should land near roundoff, not merely under some threshold
    smaller than 30 J."""
    solver = _build_solver()
    solver.solve_transient_core(verbose=False)
    residual = solver.time_series["scalars"]["energy_residual_J"]
    assert float(max(abs(residual))) < 1e-6


def test_coolant_temperature_stays_in_coolprop_valid_range():
    """The guard this regression tripped was checking (60, 2500) K -- both
    ends wrong for Helium's real CoolProp range. Assert against the REAL
    range (CP.PropsSI('Tmin'/'Tmax', 'Helium')), not the guard's constants."""
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
    """Guard against this test silently becoming a no-op check (n_sub == 1
    every step would mean the fix's actual mechanism isn't exercised here)."""
    import hps_combustor.transient_core.adapters_shelltube as ad

    original = ad._cfl_stable_substep_count
    seen = []

    def spy(mass, faces, dt, **kwargs):
        n = original(mass, faces, dt, **kwargs)
        seen.append(n)
        return n

    ad._cfl_stable_substep_count = spy
    try:
        _build_solver().solve_transient_core(verbose=False)
    finally:
        ad._cfl_stable_substep_count = original

    assert seen, "subcycling helper was never called"
    assert all(n > 1 for n in seen), f"expected subcycling every step, got {seen}"

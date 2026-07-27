"""Accuracy and speed validation for the opt-in TTSE/BICUBIC CoolProp
property backend (coolantProp.liquid_property_backend).

Both tabulated backends are compared against the exact "HEOS" evaluation
across a (p, h) grid spanning the project's operating pressure range and the
full subcooled -> two-phase -> superheated span, with a DENSE grid right at
the saturation boundaries specifically - interpolation error from a bicubic/
Taylor-expansion table is largest exactly there, which is also exactly where
CHF margin and quality are evaluated, so accuracy at the dome boundary is the
one number that actually matters for trusting these backends in the liquid
march (see docs/solver_design/STEADY_LIQUID_BOILING_REQUIREMENTS.md item 7).

Run from the repository root:

    python -m hps_combustor.validation.liquid_ttse_backend_validation
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from hps_combustor.physics.liquid_flow.correlations import equilibrium_state_ph, saturation_state
from hps_combustor.physics.liquid_flow.dispatch import coolprop_fluid_string

FLUID = "Water"
PRESSURES_PA = [30e5, 80e5, 150e5]  # spans this project's ~80 bar operating point
N_ACROSS_DOME = 60       # coarse sweep, subcooled -> superheated
N_NEAR_BOUNDARY = 200    # dense sweep within a thin band of h_l/h_v
NEAR_BOUNDARY_WINDOW = 0.01  # fraction of h_fg on each side of the boundary
N_TIMING_CALLS = 3000    # representative of one march's worth of node evaluations


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _grid_for_pressure(p_Pa: float) -> np.ndarray:
    sat = saturation_state(FLUID, p_Pa)
    h_l, h_v, h_fg = sat.h_l_J_kg, sat.h_v_J_kg, sat.h_fg_J_kg
    coarse = np.linspace(h_l - 0.3 * h_fg, h_v + 0.3 * h_fg, N_ACROSS_DOME)
    near_l = np.linspace(h_l - NEAR_BOUNDARY_WINDOW * h_fg, h_l + NEAR_BOUNDARY_WINDOW * h_fg, N_NEAR_BOUNDARY)
    near_v = np.linspace(h_v - NEAR_BOUNDARY_WINDOW * h_fg, h_v + NEAR_BOUNDARY_WINDOW * h_fg, N_NEAR_BOUNDARY)
    return np.unique(np.concatenate([coarse, near_l, near_v]))


def _accuracy_for_backend(backend: str) -> dict:
    fluid_tagged = coolprop_fluid_string(FLUID, backend)
    per_pressure = []
    for p_Pa in PRESSURES_PA:
        h_grid = _grid_for_pressure(p_Pa)
        max_rel_T = max_rel_rho = 0.0
        max_abs_quality = 0.0
        n_points = 0
        for h in h_grid:
            ref = equilibrium_state_ph(p_Pa, float(h), FLUID)  # HEOS, exact
            test = equilibrium_state_ph(p_Pa, float(h), fluid_tagged)
            n_points += 1
            max_rel_T = max(max_rel_T, abs(test.T_K - ref.T_K) / ref.T_K)
            if ref.rho_kg_m3 > 0:
                max_rel_rho = max(max_rel_rho, abs(test.rho_kg_m3 - ref.rho_kg_m3) / ref.rho_kg_m3)
            if ref.phase == "two_phase" or test.phase == "two_phase":
                max_abs_quality = max(max_abs_quality, abs(test.quality - ref.quality))
        per_pressure.append(dict(
            p_bar=p_Pa / 1e5, n_points=n_points,
            max_rel_error_T=max_rel_T, max_rel_error_rho=max_rel_rho,
            max_abs_error_quality=max_abs_quality,
        ))
    return dict(backend=backend, per_pressure=per_pressure)


def _timing_for_backend(backend: str) -> dict:
    fluid_tagged = coolprop_fluid_string(FLUID, backend)
    p_Pa = 80e5
    sat = saturation_state(FLUID, p_Pa)
    h_points = np.linspace(sat.h_l_J_kg - 0.2 * sat.h_fg_J_kg, sat.h_v_J_kg + 0.2 * sat.h_fg_J_kg, 50)
    # warm up (builds/caches the table on first use for TTSE/BICUBIC - excluded from the timed loop)
    for h in h_points:
        equilibrium_state_ph(p_Pa, float(h), fluid_tagged)
    t0 = time.perf_counter()
    for i in range(N_TIMING_CALLS):
        h = h_points[i % len(h_points)]
        equilibrium_state_ph(p_Pa, float(h), fluid_tagged)
    elapsed = time.perf_counter() - t0
    return dict(backend=backend, n_calls=N_TIMING_CALLS, total_s=elapsed,
                per_call_us=elapsed / N_TIMING_CALLS * 1e6)


def generate_report() -> dict:
    accuracy = [_accuracy_for_backend(b) for b in ("BICUBIC", "TTSE")]
    timing = [_timing_for_backend(b) for b in ("HEOS", "BICUBIC", "TTSE")]
    heos_us = next(t["per_call_us"] for t in timing if t["backend"] == "HEOS")
    for t in timing:
        t["speedup_vs_heos"] = heos_us / t["per_call_us"]
    return dict(fluid=FLUID, pressures_bar=[p / 1e5 for p in PRESSURES_PA],
                accuracy=accuracy, timing=timing)


def main() -> None:
    report = generate_report()
    out = _repo_root() / "docs" / "validation" / "liquid_ttse_backend_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))

    print("=" * 70)
    print("TTSE/BICUBIC accuracy vs HEOS (equilibrium_state_ph, Water)")
    print("=" * 70)
    for entry in report["accuracy"]:
        print(f"  {entry['backend']}:")
        for pp in entry["per_pressure"]:
            print(f"    p={pp['p_bar']:.0f} bar  n={pp['n_points']:4d}  "
                  f"max|dT|/T={pp['max_rel_error_T']:.2e}  "
                  f"max|drho|/rho={pp['max_rel_error_rho']:.2e}  "
                  f"max|dx|={pp['max_abs_error_quality']:.2e}")
    print("-" * 70)
    print("Measured per-call speed (equilibrium_state_ph, 80 bar, water):")
    for t in report["timing"]:
        print(f"  {t['backend']:8s} {t['per_call_us']:8.3f} us/call  "
              f"({t['speedup_vs_heos']:.1f}x vs HEOS)")
    print("=" * 70)
    print(f"SAVED {out}")


if __name__ == "__main__":
    main()

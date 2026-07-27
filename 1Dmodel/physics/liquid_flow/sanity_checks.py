"""Thermodynamic sanity gates for the coupled liquid coolant march.

Run at the end of a steady liquid-coolant solve
(``main_solver``/``shellntube_solver`` with
``coolantProp.coolant_model == "equilibrium_liquid"``). See
docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md, Phase 2, Design
Decision 6. These gates are engineering-decision gates, not just diagnostics:
CHF margin and dryout are hard failures because the model has no post-CHF/
mist-flow closure, so results past that point are not physical.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from hps_combustor.physics.liquid_flow.chf import chf_regime
from hps_combustor.physics.liquid_flow.correlations import saturation_state


@dataclass
class LiquidMarchSanityReport:
    passed: bool
    energy_balance_ok: bool
    energy_balance_rel_error: float
    temperature_ordering_ok: bool
    temperature_ordering_violations: int
    saturation_consistency_ok: bool
    max_saturation_deviation_K: float
    pressure_monotonic_ok: bool
    quality_bounds_ok: bool
    void_fraction_bounds_ok: bool
    chf_margin_ok: bool
    min_chf_margin: float
    chf_regime_at_min_margin: str | None
    dryout_risk: bool
    mach_c_max: float
    mach_choking_ok: bool
    messages: list[str] = field(default_factory=list)


def check_liquid_march(
    data_master: dict,
    *,
    fluid: str,
    mass_flow_c: float,
    energy_balance_tol: float = 0.02,
    chf_margin_limit: float = 1.0,
    saturation_tol_K: float = 0.5,
    mach_c_warn_limit: float = 0.5,
    mach_c_choking_limit: float = 1.0,
) -> LiquidMarchSanityReport:
    """Run all thermodynamic sanity gates on a completed liquid-coolant march.

    ``data_master`` is the solver's own dict of per-node lists. Index order
    follows the arc-length/axial march direction, not necessarily the
    physical coolant flow direction (counter-flow marches with a negative
    flow sign) — checks that care about direction (energy closure) compare
    magnitudes rather than assuming a specific inlet/outlet index.
    """
    messages: list[str] = []

    enthalpy_c = np.asarray(data_master["enthalpy_c"], dtype=float)
    p_c = np.asarray(data_master["p_c"], dtype=float)
    T_c = np.asarray(data_master["T_c"], dtype=float)
    T_wc = np.asarray(data_master["T_wc"], dtype=float)
    T_wg = np.asarray(data_master["T_wg"], dtype=float)
    T_g = np.asarray(data_master["T_g"], dtype=float)
    quality_c = np.asarray(data_master["quality_c"], dtype=float)
    void_c = np.asarray(data_master["void_c"], dtype=float)
    chf_margin_c = np.asarray(data_master["chf_margin_c"], dtype=float)
    dQ = np.asarray(data_master["dQ"], dtype=float)

    # --- energy closure: total wall duty vs coolant enthalpy rise ---
    # data_master records each node's state BEFORE that node's step advances
    # it (main_solve.py: append, then _advance_state()). So dQ[i] is the duty
    # that carries the march from the recorded state i to recorded state i+1;
    # the LAST recorded node's dQ advances to an unrecorded final state. Sum
    # only dQ[:-1] to match the recorded enthalpy span exactly (to float
    # precision), not ~1/N (the last node's share) off from it.
    Q_nodes = float(np.sum(dQ[:-1])) if dQ.size > 1 else float(np.sum(dQ))
    expected_dh = Q_nodes / mass_flow_c
    actual_dh = float(enthalpy_c[-1] - enthalpy_c[0])
    energy_balance_rel_error = abs(abs(actual_dh) - abs(expected_dh)) / max(abs(expected_dh), 1.0)
    energy_balance_ok = energy_balance_rel_error < energy_balance_tol
    if not energy_balance_ok:
        messages.append(
            f"energy balance off by {energy_balance_rel_error * 100:.2f}% "
            f"(> {energy_balance_tol * 100:.2f}% tolerance)"
        )

    # --- temperature ordering (heating: T_g > T_wg > T_wc > T_c) ---
    ordering_violations = int(
        np.sum(T_c > T_wc + 0.5) + np.sum(T_wc > T_wg + 0.5) + np.sum(T_wg > T_g + 0.5)
    )
    temperature_ordering_ok = ordering_violations == 0
    if not temperature_ordering_ok:
        messages.append(f"temperature ordering violated at {ordering_violations} node(s)")

    # --- saturation consistency where 0 < x < 1: T must track Tsat(p) ---
    two_phase = (quality_c > 0.0) & (quality_c < 1.0) & np.isfinite(quality_c)
    max_saturation_deviation_K = 0.0
    if np.any(two_phase):
        deviations = [
            abs(float(T) - saturation_state(fluid, float(p)).T_sat_K)
            for p, T in zip(p_c[two_phase], T_c[two_phase])
        ]
        max_saturation_deviation_K = float(max(deviations))
    saturation_consistency_ok = max_saturation_deviation_K <= saturation_tol_K
    if not saturation_consistency_ok:
        messages.append(
            f"T vs Tsat(p) deviation {max_saturation_deviation_K:.2f} K "
            f"(> {saturation_tol_K:.2f} K tolerance) in the two-phase region"
        )

    # --- pressure strictly monotonic along the march direction ---
    dp = np.diff(p_c)
    pressure_monotonic_ok = bool(np.all(dp <= 1.0e-6)) or bool(np.all(dp >= -1.0e-6))
    if not pressure_monotonic_ok:
        messages.append("coolant pressure is not monotonic along the march")

    # --- bounds ---
    # Quality has no fixed physical lower bound: deeply subcooled liquid can
    # have arbitrarily negative equilibrium quality (it is an enthalpy-based
    # subcooling measure, not a fraction, outside [0, 1]). The meaningful
    # upper-end check is the dryout flag (x >= 1) below. NaN quality is NOT a
    # failure: it is the supercritical-pressure signal (no dome, so quality is
    # undefined) -- only reject non-finite INFINITIES. Subcritically quality is
    # always finite, so this is bit-identical to the former isfinite() check.
    finite_quality = quality_c[np.isfinite(quality_c)]
    quality_bounds_ok = bool(quality_c.size > 0 and not np.any(np.isinf(quality_c)))
    if not quality_bounds_ok:
        messages.append("equilibrium quality is infinite at one or more nodes")
    finite_void = void_c[np.isfinite(void_c)]
    void_fraction_bounds_ok = bool(
        finite_void.size == 0 or np.all((finite_void >= 0.0) & (finite_void <= 1.0))
    )
    if not void_fraction_bounds_ok:
        messages.append("void fraction outside [0, 1]")

    # --- CHF margin: hard gate, not just a diagnostic (no post-CHF model exists) ---
    finite_margin_mask = np.isfinite(chf_margin_c)
    finite_margin = chf_margin_c[finite_margin_mask]
    min_chf_margin = float(np.min(finite_margin)) if finite_margin.size else float("nan")
    chf_margin_ok = bool(np.isnan(min_chf_margin) or min_chf_margin >= chf_margin_limit)
    chf_regime_at_min_margin = None
    if finite_margin.size:
        worst_idx = np.flatnonzero(finite_margin_mask)[np.argmin(finite_margin)]
        chf_regime_at_min_margin = chf_regime(float(quality_c[worst_idx]))
    if not chf_margin_ok:
        messages.append(
            f"minimum CHF margin {min_chf_margin:.2f} < {chf_margin_limit:.2f} "
            f"-- dryout risk ({chf_regime_at_min_margin} regime)"
        )

    dryout_risk = bool(np.any(finite_quality >= 1.0)) or not chf_margin_ok
    if dryout_risk and chf_margin_ok:
        messages.append("quality reached 1.0 (complete vaporization) -- past model validity")

    # --- two-phase compressibility / choking: mixture sound speed collapses
    # (Wood's equation) well before either pure-phase sound speed, so a
    # modest coolant velocity can still approach choking. See "Assumptions
    # That Break" in docs/solver_design/water_coolant_conversion_plan.md.
    mach_c = np.asarray(data_master.get("Mach_c", []), dtype=float)
    finite_mach_c = mach_c[np.isfinite(mach_c)]
    mach_c_max = float(np.max(finite_mach_c)) if finite_mach_c.size else float("nan")
    mach_choking_ok = bool(np.isnan(mach_c_max) or mach_c_max < mach_c_choking_limit)
    if not mach_choking_ok:
        messages.append(
            f"Mach_c reached {mach_c_max:.3f} >= {mach_c_choking_limit:.2f} -- "
            f"two-phase choking risk, results past this point are not physical"
        )
    elif not np.isnan(mach_c_max) and mach_c_max > mach_c_warn_limit:
        messages.append(
            f"Mach_c_max = {mach_c_max:.3f} > {mach_c_warn_limit:.2f} -- "
            f"coolant compressibility significant"
        )

    passed = (
        energy_balance_ok
        and temperature_ordering_ok
        and saturation_consistency_ok
        and pressure_monotonic_ok
        and quality_bounds_ok
        and void_fraction_bounds_ok
        and chf_margin_ok
        and not dryout_risk
        and mach_choking_ok
    )

    return LiquidMarchSanityReport(
        passed=passed,
        energy_balance_ok=energy_balance_ok,
        energy_balance_rel_error=float(energy_balance_rel_error),
        temperature_ordering_ok=temperature_ordering_ok,
        temperature_ordering_violations=ordering_violations,
        saturation_consistency_ok=saturation_consistency_ok,
        max_saturation_deviation_K=max_saturation_deviation_K,
        pressure_monotonic_ok=pressure_monotonic_ok,
        quality_bounds_ok=quality_bounds_ok,
        void_fraction_bounds_ok=void_fraction_bounds_ok,
        chf_margin_ok=chf_margin_ok,
        min_chf_margin=min_chf_margin,
        chf_regime_at_min_margin=chf_regime_at_min_margin,
        dryout_risk=dryout_risk,
        mach_c_max=mach_c_max,
        mach_choking_ok=mach_choking_ok,
        messages=messages,
    )

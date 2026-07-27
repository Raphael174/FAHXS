"""Audit bang-bang helium schedule timescales for coolant momentum relevance."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def load_time_mdot(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        try:
            t = float(parts[0].strip().replace(",", "."))
            mdot = float(parts[1].strip().replace(",", "."))
        except ValueError:
            continue
        rows.append((t, max(0.0, mdot)))
    if not rows:
        raise ValueError(f"No time/mass-flow rows found in {path}")
    data = np.asarray(rows, dtype=float)
    order = np.argsort(data[:, 0])
    data = data[order]
    unique_t, unique_idx = np.unique(data[:, 0], return_index=True)
    return unique_t, data[unique_idx, 1]


def audit_schedule(
    t: np.ndarray,
    mdot: np.ndarray,
    *,
    helical_length_m: float = 6.0,
    helical_inner_diameter_m: float = 3.5e-3,
    shelltube_length_m: float = 0.235,
    shelltube_flow_area_m2: float = 7.5e-3,
    sound_speed_min_m_s: float = 550.0,
    sound_speed_max_m_s: float = 1000.0,
) -> dict:
    dt = np.diff(t)
    dmdot = np.diff(mdot)
    valid = dt > 0.0
    slopes = np.abs(dmdot[valid] / dt[valid])
    max_slope = float(np.max(slopes)) if slopes.size else 0.0
    mdot_max = float(np.max(mdot))
    mdot_mean = float(np.mean(mdot))
    tau_10_90 = _ten_to_ninety_time(t, mdot)

    A_helical = np.pi * helical_inner_diameter_m**2 / 4.0
    A_shell = shelltube_flow_area_m2
    dp_inertia_helical = helical_length_m / A_helical * max_slope
    dp_inertia_shell = shelltube_length_m / A_shell * max_slope

    return {
        "n_points": int(len(t)),
        "t_start_s": float(t[0]),
        "t_end_s": float(t[-1]),
        "dt_min_positive_s": float(np.min(dt[valid])) if np.any(valid) else None,
        "dt_max_s": float(np.max(dt)) if dt.size else None,
        "mdot_max_kg_s": mdot_max,
        "mdot_mean_kg_s": mdot_mean,
        "max_abs_dmdot_dt_kg_s2": max_slope,
        "ten_to_ninety_ramp_s": tau_10_90,
        "helical": {
            "length_m": helical_length_m,
            "inner_diameter_m": helical_inner_diameter_m,
            "area_m2": A_helical,
            "acoustic_time_min_s": helical_length_m / sound_speed_max_m_s,
            "acoustic_time_max_s": helical_length_m / sound_speed_min_m_s,
            "inertial_dp_scale_Pa": dp_inertia_helical,
            "inertial_dp_scale_bar": dp_inertia_helical / 1.0e5,
        },
        "shelltube": {
            "length_m": shelltube_length_m,
            "flow_area_m2": shelltube_flow_area_m2,
            "acoustic_time_min_s": shelltube_length_m / sound_speed_max_m_s,
            "acoustic_time_max_s": shelltube_length_m / sound_speed_min_m_s,
            "inertial_dp_scale_Pa": dp_inertia_shell,
            "inertial_dp_scale_bar": dp_inertia_shell / 1.0e5,
        },
        "interpretation": {
            "helical_momentum": (
                "transient momentum likely significant if the supplied dmdot/dt "
                "survives upstream smoothing"
            ),
            "shelltube_momentum": (
                "lower inertial pressure scale; transient mass/energy with "
                "quasi-steady momentum is a reasonable first implementation"
            ),
        },
    }


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("schedule_file")
    parser.add_argument(
        "--output",
        default="docs/validation/bangbang_momentum_audit.json",
    )
    args = parser.parse_args(argv)

    t, mdot = load_time_mdot(args.schedule_file)
    payload = audit_schedule(t, mdot)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {out}")
    print(f"max mdot: {payload['mdot_max_kg_s']:.6g} kg/s")
    print(f"max |dmdot/dt|: {payload['max_abs_dmdot_dt_kg_s2']:.6g} kg/s2")
    print(f"helical inertial dp scale: {payload['helical']['inertial_dp_scale_bar']:.3g} bar")
    print(f"shelltube inertial dp scale: {payload['shelltube']['inertial_dp_scale_bar']:.3g} bar")
    return payload


def _ten_to_ninety_time(t: np.ndarray, mdot: np.ndarray) -> float | None:
    peak = float(np.max(mdot))
    if peak <= 0.0:
        return None
    lo = 0.1 * peak
    hi = 0.9 * peak
    above_lo = np.where(mdot >= lo)[0]
    above_hi = np.where(mdot >= hi)[0]
    if above_lo.size == 0 or above_hi.size == 0:
        return None
    return float(t[above_hi[0]] - t[above_lo[0]])


if __name__ == "__main__":
    main()

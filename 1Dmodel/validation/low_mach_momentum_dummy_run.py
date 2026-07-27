"""Run the HX dummy schedule with shell-side low-Mach coolant momentum.

This is a short validation case for the pressure-driven transient momentum
closure:

    python -m hps_combustor.validation.low_mach_momentum_dummy_run

The case uses `inputs/HX_dummy_inlet.xlsx` for helium inlet pressure and
temperature, fixes downstream shell-side pressure to 70 bar, and runs only the
first 20 simulated seconds.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from ..main_transient import build_inputs, run_transient
from ..schedule_inputs import apply_schedule_file


def run_case(
    *,
    schedule_file: str = "inputs/HX_dummy_inlet.xlsx",
    output: str = "docs/validation/low_mach_momentum_dummy_run.json",
    n_nodes: int = 20,
    t_end: float = 20.0,
    max_step: float = 0.05,
) -> dict:
    inputs = build_inputs()
    apply_schedule_file(
        inputs["transient"],
        inputs["hotgas"],
        schedule_file,
        coolant=inputs["coolant"],
    )

    inputs["combustor"].HX_config = "shellntube"
    inputs["combustor"].flow_config = "counter"
    inputs["transient"].fluid_model = "transient_coolant"
    inputs["transient"].coolant_momentum_model = "low_mach"
    inputs["transient"].transient_coolant_outlet_pressure = 70.0e5
    inputs["transient"].t_end = float(t_end)
    inputs["transient"].max_step = float(max_step)
    inputs["transient"].insert_schedule_breakpoints = False
    inputs["transient"].n_save = 81
    inputs["transient"].T_wall_initial = 300.0
    inputs["transient"].T_coolant_initial = None
    inputs["transient"].schedule_mass_flow_lox = ((0.0, 0.090), (float(t_end), 0.090))
    inputs["transient"].schedule_mass_flow_diesel = ((0.0, 0.030), (float(t_end), 0.030))
    inputs["transient"].schedule_mass_flow_g = ((0.0, 0.120), (float(t_end), 0.120))
    inputs["transient"].schedule_OF = ((0.0, 3.0), (float(t_end), 3.0))
    inputs["transient"].schedule_T_lox_in = ((0.0, 100.0), (float(t_end), 100.0))
    inputs["transient"].schedule_ignition_state = ((0.0, 1.0), (float(t_end), 1.0))
    inputs["transient"].ignition_time = 0.0
    inputs["transient"].progress_interval_time_s = 2.0
    inputs["run"].shelltube_transient_nodes = int(n_nodes)

    start = time.perf_counter()
    solver, summary = run_transient(inputs)
    runtime = time.perf_counter() - start
    payload = _metrics(solver, summary, runtime, schedule_file)

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def _metrics(solver, summary: dict, runtime_s: float, schedule_file: str) -> dict:
    ts = solver.time_series
    scalars = ts["scalars"]
    fields = ts["fields"]
    t = np.asarray(ts["t"], dtype=float)
    mdot_in = np.asarray(scalars.get("mdot_c_inlet_face", np.nan), dtype=float)
    mdot_out = np.asarray(scalars.get("mdot_c_outlet_face", np.nan), dtype=float)
    p_field = np.asarray(fields.get("p_c", np.nan), dtype=float)

    positive = mdot_in > 0.0
    return {
        "case": {
            "schedule_file": schedule_file,
            "t_end_s": float(t[-1]) if t.size else None,
            "n_time": int(t.size),
            "n_nodes": int(len(ts["x"])),
            "momentum_model": "low_mach",
            "downstream_pressure_Pa": 70.0e5,
            "lox_kg_s": 0.090,
            "diesel_kg_s": 0.030,
            "ignition": 1,
            "T_wall_initial_K": 300.0,
        },
        "runtime_s": float(runtime_s),
        "summary": summary,
        "diagnostics": {
            "mdot_in_mean_kg_s": _finite_mean(mdot_in),
            "mdot_in_positive_mean_kg_s": _finite_mean(mdot_in[positive]),
            "mdot_in_min_kg_s": _finite_min(mdot_in),
            "mdot_in_max_kg_s": _finite_max(mdot_in),
            "mdot_out_mean_kg_s": _finite_mean(mdot_out),
            "mdot_out_min_kg_s": _finite_min(mdot_out),
            "mdot_out_max_kg_s": _finite_max(mdot_out),
            "T_c_out_final_K": _finite_last(scalars.get("T_c_out", np.nan)),
            "T_c_out_max_K": _finite_max(scalars.get("T_c_out", np.nan)),
            "T_wall_max_peak_K": _finite_max(scalars.get("T_wall_max", np.nan)),
            "T_g_out_final_K": _finite_last(scalars.get("T_g_out", np.nan)),
            "Q_hot_peak_kW": _finite_max(scalars.get("Q_hot_kW", np.nan)),
            "coolant_pressure_min_bar": _finite_min(p_field) / 1.0e5,
            "coolant_pressure_max_bar": _finite_max(p_field) / 1.0e5,
            "energy_residual_abs_max_J": _finite_abs_max(scalars.get("energy_residual_J", np.nan)),
            "coolant_mass_residual_abs_max_kg": _finite_abs_max(
                scalars.get("coolant_mass_residual_kg", np.nan)
            ),
        },
    }


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--schedule-file", default="inputs/HX_dummy_inlet.xlsx")
    parser.add_argument("--output", default="docs/validation/low_mach_momentum_dummy_run.json")
    parser.add_argument("--n-nodes", type=int, default=20)
    parser.add_argument("--t-end", type=float, default=20.0)
    parser.add_argument("--max-step", type=float, default=0.05)
    args = parser.parse_args(argv)

    payload = run_case(
        schedule_file=args.schedule_file,
        output=args.output,
        n_nodes=args.n_nodes,
        t_end=args.t_end,
        max_step=args.max_step,
    )
    print(f"Wrote {args.output}")
    diag = payload["diagnostics"]
    print(
        "low_mach_dummy: "
        f"runtime={payload['runtime_s']:.2f}s, "
        f"mdot_in_mean={diag['mdot_in_mean_kg_s']:.5f} kg/s, "
        f"T_c_out_final={diag['T_c_out_final_K']:.2f} K, "
        f"T_wall_peak={diag['T_wall_max_peak_K']:.2f} K"
    )
    return payload


def _finite_values(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def _finite_mean(values) -> float:
    arr = _finite_values(values)
    return float(np.mean(arr)) if arr.size else float("nan")


def _finite_min(values) -> float:
    arr = _finite_values(values)
    return float(np.min(arr)) if arr.size else float("nan")


def _finite_max(values) -> float:
    arr = _finite_values(values)
    return float(np.max(arr)) if arr.size else float("nan")


def _finite_abs_max(values) -> float:
    arr = _finite_values(values)
    return float(np.max(np.abs(arr))) if arr.size else float("nan")


def _finite_last(values) -> float:
    arr = _finite_values(values)
    return float(arr[-1]) if arr.size else float("nan")


if __name__ == "__main__":
    main()

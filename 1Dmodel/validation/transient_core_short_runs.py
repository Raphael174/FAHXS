"""Short validation matrix for the shell-and-tube transient-coolant core.

Run from the repository root:

    python -m hps_combustor.validation.transient_core_short_runs

The matrix is intentionally short. It is meant to catch numerical/pathology
issues and track runtime before committing to longer 100 s studies.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from ..main_transient import build_inputs, run_transient


@dataclass(frozen=True)
class ValidationCase:
    name: str
    flow_config: str
    t_end: float
    max_step: float
    n_nodes: int
    schedule_mass_flow_c: tuple[tuple[float, float], ...]
    schedule_mass_flow_g: tuple[tuple[float, float], ...]
    schedule_ignition_state: tuple[tuple[float, float], ...] | None = None
    schedule_mass_flow_lox: tuple[tuple[float, float], ...] | None = None
    schedule_T_lox_in: tuple[tuple[float, float], ...] | None = None


def build_case_matrix() -> list[ValidationCase]:
    """Return the default short validation matrix."""

    bangbang_he = (
        (0.00, 0.15),
        (0.25, 0.00),
        (0.50, 0.15),
        (1.00, 0.00),
        (1.25, 0.15),
    )
    hot_ramp = ((0.00, 0.065), (1.00, 0.120))

    gox_he = (
        (0.00, 0.15),
        (0.75, 0.15),
        (1.00, 0.02),
        (1.25, 0.15),
    )
    ignition = ((0.00, 0.0), (0.70, 0.0), (0.75, 1.0), (2.00, 1.0))
    gox_flow = ((0.00, 0.03), (0.70, 0.03), (0.75, 0.00))
    gox_T = ((0.00, 120.0), (0.70, 120.0))
    torch_to_full = (
        (0.00, 0.000),
        (0.75, 0.005),
        (1.25, 0.080),
        (2.00, 0.120),
    )

    cases: list[ValidationCase] = []
    for max_step in (0.25, 0.10):
        cases.append(
            ValidationCase(
                name=f"finite_counter_bangbang_dt{_tag(max_step)}_n20",
                flow_config="counter",
                t_end=2.0,
                max_step=max_step,
                n_nodes=20,
                schedule_mass_flow_c=bangbang_he,
                schedule_mass_flow_g=hot_ramp,
            )
        )

    for n_nodes in (20, 40):
        cases.append(
            ValidationCase(
                name=f"gox_ignition_ramp_dt0p25_n{n_nodes}",
                flow_config="counter",
                t_end=2.0,
                max_step=0.25,
                n_nodes=n_nodes,
                schedule_mass_flow_c=gox_he,
                schedule_mass_flow_g=torch_to_full,
                schedule_ignition_state=ignition,
                schedule_mass_flow_lox=gox_flow,
                schedule_T_lox_in=gox_T,
            )
        )

    return cases


def run_case(case: ValidationCase) -> dict[str, Any]:
    """Run one validation case and return serializable metrics."""

    inputs = build_inputs()
    inputs["combustor"].HX_config = "shellntube"
    inputs["combustor"].flow_config = case.flow_config
    inputs["transient"].fluid_model = "transient_coolant"
    inputs["transient"].t_end = case.t_end
    inputs["transient"].max_step = case.max_step
    inputs["transient"].n_save = max(2, int(round(case.t_end / case.max_step)) + 1)
    inputs["transient"].schedule_mass_flow_c = case.schedule_mass_flow_c
    inputs["transient"].schedule_mass_flow_g = case.schedule_mass_flow_g
    inputs["transient"].schedule_ignition_state = case.schedule_ignition_state
    inputs["transient"].schedule_mass_flow_lox = case.schedule_mass_flow_lox
    inputs["transient"].schedule_T_lox_in = case.schedule_T_lox_in
    inputs["run"].shelltube_transient_nodes = case.n_nodes

    start = time.perf_counter()
    solver, summary = run_transient(inputs)
    runtime = time.perf_counter() - start

    ts = solver.time_series
    scalars = ts["scalars"]
    fields = ts["fields"]
    return {
        "case": asdict(case),
        "runtime_s": runtime,
        "runtime_per_simulated_second": runtime / case.t_end if case.t_end > 0.0 else np.nan,
        "n_time": int(len(ts["t"])),
        "n_nodes": int(len(ts["x"])),
        "n_steps": int(len(solver.core_result.integration.t)),
        "summary": summary,
        "diagnostics": {
            "T_c_out_peak_K": _finite_max(scalars["T_c_out"]),
            "T_wall_max_peak_K": _finite_max(scalars["T_wall_max"]),
            "Q_hot_peak_kW": _finite_max(scalars["Q_hot_kW"]),
            "Q_cold_peak_kW": _finite_max(scalars["Q_cold_kW"]),
            "energy_residual_abs_max_J": _finite_abs_max(scalars["energy_residual_J"]),
            "Re_g_max": _finite_max(scalars.get("Re_g_max", np.nan)),
            "Re_shell_max": _finite_max(scalars.get("Re_shell_max", np.nan)),
            "dp_g_total_peak_Pa": _finite_max(scalars.get("dp_g_total_Pa", np.nan)),
            "dp_shell_total_peak_Pa": _finite_max(scalars.get("dp_shell_total_Pa", np.nan)),
            "progress_g_out_final": _finite_last(scalars.get("progress_g_out", np.nan)),
            "h_removed_g_out_final_J_kg": _finite_last(
                scalars.get("h_removed_g_out_J_kg", np.nan)
            ),
            "T_c_field_min_K": _finite_min(fields["T_c"]),
            "T_c_field_max_K": _finite_max(fields["T_c"]),
            "T_wall_field_min_K": _finite_min(fields["Tbar"]),
            "T_wall_field_max_K": _finite_max(fields["Tbar"]),
        },
    }


def run_matrix(cases: list[ValidationCase] | None = None) -> dict[str, Any]:
    cases = build_case_matrix() if cases is None else cases
    started = datetime.now().isoformat(timespec="seconds")
    results = [run_case(case) for case in cases]
    return {
        "created_at": started,
        "description": "Shell-and-tube transient_coolant short-run validation matrix.",
        "results": results,
    }


def main(argv: list[str] | None = None) -> dict[str, Any]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="docs/validation/transient_core_short_run_results.json",
        help="JSON output path.",
    )
    args = parser.parse_args(argv)

    payload = run_matrix()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"Wrote {output}")
    for result in payload["results"]:
        summary = result["summary"]
        print(
            f"{result['case']['name']}: "
            f"runtime={result['runtime_s']:.3f}s, "
            f"T_c_out={summary['T_c_out_final']:.2f}K, "
            f"T_wall_max={summary['T_wall_max_final']:.2f}K"
        )
    return payload


def _tag(value: float) -> str:
    return str(value).replace(".", "p")


def _finite_min(values) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.min(arr)) if arr.size else float("nan")


def _finite_max(values) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.max(arr)) if arr.size else float("nan")


def _finite_abs_max(values) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.max(np.abs(arr))) if arr.size else float("nan")


def _finite_last(values) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(arr[-1]) if arr.size else float("nan")


if __name__ == "__main__":
    main()

"""Short audit for shell-and-tube bang-bang helium coolant behavior.

Run from the repository root:

    python -m hps_combustor.validation.shelltube_bangbang_coolant_audit

The audit checks the production `transient_coolant` shell-and-tube path. It is
not a physics validation against test data; it verifies the intended numerical
behavior for a valve-like helium schedule:

- inlet face follows the open-flow command before shutoff;
- inlet face closes after shutoff;
- outlet face can continue discharging residual helium from the exchanger;
- total coolant inventory decreases after shutoff.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from ..main_transient import build_inputs, run_transient


def run_audit(
    *,
    output: str | Path = "docs/validation/shelltube_bangbang_coolant_audit.json",
    flow_config: str = "counter",
) -> dict[str, Any]:
    inputs = build_inputs()
    inputs["combustor"].HX_config = "shellntube"
    inputs["combustor"].flow_config = flow_config
    inputs["transient"].fluid_model = "transient_coolant"
    inputs["transient"].t_end = 0.02
    inputs["transient"].max_step = 0.25
    inputs["transient"].n_save = 2
    inputs["transient"].schedule_mass_flow_c = ((0.0, 0.15), (0.01, 0.0))
    inputs["run"].shelltube_transient_nodes = 3

    started = time.perf_counter()
    solver, summary = run_transient(inputs)
    runtime_s = time.perf_counter() - started

    fields = solver.time_series["fields"]
    t = np.asarray(solver.core_result.integration.t, dtype=float)
    face_mdot = np.asarray(fields["face_mdot_c"], dtype=float)
    mass = np.asarray(fields["coolant_mass_kg"], dtype=float)

    closure_index = int(np.searchsorted(t, 0.01, side="left"))
    after_start = int(np.searchsorted(t, 0.01, side="right"))
    after = slice(after_start, None)
    if flow_config == "counter":
        inlet_face = face_mdot[:, -1]
        outlet_face = face_mdot[:, 0]
        inlet_command_sign = -1.0
        residual_outlet_positive = np.abs(outlet_face[after])
    else:
        inlet_face = face_mdot[:, 0]
        outlet_face = face_mdot[:, -1]
        inlet_command_sign = 1.0
        residual_outlet_positive = np.maximum(outlet_face[after], 0.0)

    total_mass = np.sum(mass, axis=1)
    dt = np.diff(t)

    metrics = {
        "flow_config": flow_config,
        "runtime_s": runtime_s,
        "n_internal_steps": int(t.size),
        "max_internal_dt_s": float(np.max(dt)) if dt.size else 0.0,
        "closure_index": closure_index,
        "after_closure_index": after_start,
        "summary": summary,
        "pre_shutoff_inlet_mdot_kg_s": float(inlet_face[1]) if t.size > 1 else float("nan"),
        "post_shutoff_inlet_mdot_abs_max_kg_s": float(np.max(np.abs(inlet_face[after]))),
        "post_shutoff_outlet_mdot_abs_max_kg_s": float(np.max(np.abs(outlet_face[after]))),
        "post_shutoff_residual_outlet_mdot_max_kg_s": float(np.max(residual_outlet_positive)),
        "mass_at_shutoff_kg": float(total_mass[closure_index]),
        "mass_final_kg": float(total_mass[-1]),
        "mass_loss_after_shutoff_kg": float(total_mass[closure_index] - total_mass[-1]),
        "energy_residual_final_J": float(summary.get("energy_residual_J_final", np.nan)),
        "mass_residual_final_kg": float(summary.get("coolant_mass_residual_kg_final", np.nan)),
        "checks": {
            "inlet_open_before_shutoff": bool(
                inlet_command_sign * inlet_face[1] > 0.05
            ) if t.size > 1 else False,
            "inlet_closed_after_shutoff": bool(np.max(np.abs(inlet_face[after])) < 1.0e-8),
            "residual_outlet_flow_after_shutoff": bool(np.max(np.abs(outlet_face[after])) > 1.0e-5),
            "inventory_decreases_after_shutoff": bool(total_mass[-1] < total_mass[closure_index]),
        },
    }
    metrics["checks"]["all_passed"] = bool(all(metrics["checks"].values()))

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    return metrics


def main(argv: list[str] | None = None) -> dict[str, Any]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="docs/validation/shelltube_bangbang_coolant_audit.json",
        help="JSON output path.",
    )
    parser.add_argument("--flow-config", choices=("co", "counter"), default="counter")
    args = parser.parse_args(argv)

    metrics = run_audit(output=args.output, flow_config=args.flow_config)
    print(f"Wrote {args.output}")
    print(json.dumps(metrics["checks"], indent=2))
    return metrics


if __name__ == "__main__":
    main()

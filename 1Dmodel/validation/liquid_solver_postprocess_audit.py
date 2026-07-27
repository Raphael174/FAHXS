"""Audit opt-in liquid postprocess hooks on maintained steady solver classes.

Run from the repository root:

    python -m hps_combustor.validation.liquid_solver_postprocess_audit

This is not a coupled HX solve. It verifies that the solver-class postprocess
hooks consume solver-output-shaped duty fields and return liquid p-h diagnostics
without changing the existing steady helium marches.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from hps_combustor.input_data import combustorProp, coolantProp, numericalProp, shellTubeProp
from hps_combustor.main_solve import main_solver
from hps_combustor.main_solve_shellntube import shellntube_solver


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_lut_path() -> Path:
    return _repo_root() / "docs" / "reference" / "external" / "2006LUTdata.txt"


def _water_coolant() -> coolantProp:
    return coolantProp(
        coolant="Water",
        coolant_model="equilibrium_liquid",
        mass_flow_c=0.12,
        T_in=420.0,
        p_in=1.0e6,
    )


def _diagnostics_dict(result) -> dict[str, object]:
    d = result.diagnostics
    return {
        "heat_rate_W": d.heat_rate_W,
        "pressure_drop_Pa": d.pressure_drop_Pa,
        "inlet_T_K": d.inlet_T_K,
        "outlet_T_K": d.outlet_T_K,
        "min_quality": d.min_quality,
        "max_quality": d.max_quality,
        "outlet_quality": d.outlet_quality,
        "max_void_fraction": d.max_void_fraction,
        "min_chf_margin": d.min_chf_margin,
        "boiling_reached": d.boiling_reached,
        "dryout_or_vapor_reached": d.dryout_or_vapor_reached,
        "chf_margin_below_limit": d.chf_margin_below_limit,
        "energy_residual_abs_J_kg": d.energy_residual_abs_J_kg,
        "energy_residual_ok": d.energy_residual_ok,
        "coolant_enters_at": result.coolant_enters_at,
        "n_cells": int(result.cell_fields_hx_order["z_m"].size),
    }


def run_audit(
    *,
    output: str | Path = "docs/validation/liquid_solver_postprocess_audit.json",
    lut_path: str | Path | None = None,
) -> dict[str, object]:
    lut = Path(lut_path) if lut_path is not None else _default_lut_path()
    coolant = _water_coolant()

    dQ_helical = np.linspace(80.0, 120.0, 8)
    helical_solver_like = SimpleNamespace(
        coolantProp=coolant,
        combustorProp=combustorProp(flow_config="counter", N_coils=1, Dh_coil=0.006),
        numericalProp=numericalProp(),
        data_master={"dQ": dQ_helical, "L_ch": np.linspace(0.0, 0.7, dQ_helical.size)},
    )
    helical = main_solver.liquid_coolant_postprocess(
        helical_solver_like,
        lut_path=lut,
    )

    stp = shellTubeProp()
    dQ_shell = np.linspace(4.0, 7.0, 12)
    shell_solver_like = SimpleNamespace(
        coolantProp=coolant,
        stp=stp,
        tube={"dQ": dQ_shell},
        flow_config="counter",
    )
    shell = shellntube_solver.liquid_coolant_postprocess(
        shell_solver_like,
        lut_path=lut,
    )

    report = {
        "source": "synthetic solver-output-shaped duty fixtures",
        "purpose": "verify maintained steady solver liquid_coolant_postprocess hooks",
        "helical": _diagnostics_dict(helical),
        "shelltube": _diagnostics_dict(shell),
        "checks": {
            "helical_stored_result": helical_solver_like.liquid_coolant is helical,
            "shelltube_stored_result": shell_solver_like.liquid_coolant is shell,
            "helical_energy_ok": helical.diagnostics.energy_residual_ok,
            "shelltube_energy_ok": shell.diagnostics.energy_residual_ok,
            "helical_counterflow_mapping": helical.coolant_enters_at == "z_max",
            "shelltube_counterflow_mapping": shell.coolant_enters_at == "z_max",
        },
    }
    report["checks"]["all_passed"] = bool(all(report["checks"].values()))

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> dict[str, object]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="docs/validation/liquid_solver_postprocess_audit.json",
        help="JSON output path.",
    )
    args = parser.parse_args(argv)
    report = run_audit(output=args.output)
    print(f"Wrote {args.output}")
    print(json.dumps(report["checks"], indent=2))
    return report


if __name__ == "__main__":
    main()

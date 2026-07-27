"""Unified validation matrix for the liquid heated-channel solver.

Run from the repository root:

    python -m hps_combustor.validation.liquid_validation_matrix

The matrix intentionally combines three different gates:

- straight-pipe water boiling and literature-correlation checks;
- HX-style imposed-duty handoff checks;
- maintained steady-solver postprocess hook checks.

It is a readiness report, not a replacement for the individual validation
scripts. The output JSON records both metrics and pass/fail criteria so the
liquid coolant integration path has one machine-readable audit artifact.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from hps_combustor.validation.liquid_boiling_straight_pipe import (
    generate_validation_report,
    generate_yu2002_validation_report,
)
from hps_combustor.validation.liquid_hx_imposed_duty import generate_imposed_duty_report
from hps_combustor.validation.liquid_solver_postprocess_audit import run_audit


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_output() -> Path:
    return _repo_root() / "docs" / "validation" / "liquid_validation_matrix.json"


def _default_artifact_root() -> Path:
    return _repo_root() / "docs" / "validation"


def _checks(
    *,
    straight_pipe: dict[str, Any],
    yu2002: dict[str, Any],
    imposed_duty: dict[str, Any],
    postprocess: dict[str, Any],
) -> dict[str, bool]:
    checks = {
        "straight_pipe_energy_residual_ok": abs(float(straight_pipe["energy_residual_J_kg"])) <= 1.0e-6,
        "straight_pipe_reaches_boiling": float(straight_pipe["outlet_quality"]) > 0.0,
        "straight_pipe_chf_margin_ok": float(straight_pipe["min_chf_margin"]) > 1.0,
        "groeneveld_lut_page9_exact": float(straight_pipe["groeneveld_page9_max_abs_error_kW_m2"]) == 0.0,
        "yu2002_pressure_fit_reasonable": float(
            yu2002["pressure_multiplier_mean_abs_rel_error_yu2002_fit"]
        ) <= 0.10,
        "yu2002_htc_digitized_reasonable": float(
            yu2002["htc_fig10_digitized_mean_abs_rel_error"]
        ) <= 0.10,
        "yu2002_chf_trend_ok": bool(yu2002["chf_digitized_monotonic_decrease_with_quality"]),
        "hx_imposed_duty_energy_ok": bool(imposed_duty["energy_residual_ok"]),
        "hx_imposed_duty_reaches_boiling": bool(imposed_duty["boiling_reached"]),
        "hx_imposed_duty_chf_margin_ok": float(imposed_duty["min_chf_margin"]) > 1.0,
        "hx_imposed_duty_no_dryout": not bool(imposed_duty["dryout_or_vapor_reached"]),
        "solver_postprocess_hooks_ok": bool(postprocess["checks"]["all_passed"]),
    }
    checks["all_passed"] = bool(all(checks.values()))
    return checks


def run_validation_matrix(
    *,
    output: str | Path = _default_output(),
    artifact_root: str | Path = _default_artifact_root(),
) -> dict[str, Any]:
    """Generate all liquid validation artifacts and a consolidated JSON report."""

    output = Path(output)
    artifact_root = Path(artifact_root)
    artifact_root.mkdir(parents=True, exist_ok=True)

    straight_pipe = generate_validation_report(artifact_root / "liquid_boiling_poc")
    yu2002 = generate_yu2002_validation_report(artifact_root / "liquid_boiling_yu2002")
    imposed_duty = generate_imposed_duty_report(artifact_root / "liquid_hx_imposed_duty")
    postprocess = run_audit(output=artifact_root / "liquid_solver_postprocess_audit.json")

    report = {
        "purpose": "single readiness gate for 1D heated liquid-flow boiling solver integration",
        "scope": {
            "steady_governing_state": "p,h with CoolProp/equilibrium-liquid closure",
            "validated_physics": [
                "single-phase liquid convection/friction",
                "saturated two-phase equilibrium state",
                "Gungor-Winterton boiling HTC",
                "Muller-Steinhagen-Heck two-phase friction pressure drop",
                "homogeneous acceleration pressure gradient",
                "Groeneveld 2006 CHF LUT ingestion and margin",
                "HX-grid duty handoff and counterflow postprocess mapping",
            ],
            "not_yet_validated": [
                "fully coupled production liquid wall march",
                "transient boiling/liquid finite-volume model",
                "geometry-specific shell-side or helical-coil boiling correlations",
                "dryout/post-CHF model validation",
            ],
        },
        "artifacts": {
            "straight_pipe_dir": str(artifact_root / "liquid_boiling_poc"),
            "yu2002_dir": str(artifact_root / "liquid_boiling_yu2002"),
            "hx_imposed_duty_dir": str(artifact_root / "liquid_hx_imposed_duty"),
            "postprocess_audit": str(artifact_root / "liquid_solver_postprocess_audit.json"),
        },
        "straight_pipe": straight_pipe,
        "yu2002": yu2002,
        "hx_imposed_duty": imposed_duty,
        "solver_postprocess": postprocess,
    }
    report["checks"] = _checks(
        straight_pipe=straight_pipe,
        yu2002=yu2002,
        imposed_duty=imposed_duty,
        postprocess=postprocess,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> dict[str, Any]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default=str(_default_output()),
        help="Consolidated JSON output path.",
    )
    parser.add_argument(
        "--artifact-root",
        default=str(_default_artifact_root()),
        help="Directory where validation sub-artifacts are written.",
    )
    args = parser.parse_args(argv)

    report = run_validation_matrix(output=args.output, artifact_root=args.artifact_root)
    print(f"Wrote {args.output}")
    print(json.dumps(report["checks"], indent=2))
    return report


if __name__ == "__main__":
    main()

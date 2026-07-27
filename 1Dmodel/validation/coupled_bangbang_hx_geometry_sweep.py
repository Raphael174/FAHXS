"""Geometry sweep for the coupled bang-bang HX validation module.

Run from the repository root:

    python -m hps_combustor.validation.coupled_bangbang_hx_geometry_sweep

The sweep reuses the same 0D tank/feed/valve schedule and the same coupled HX
runner as `coupled_bangbang_hx.py`. It is intentionally small: the purpose is
to screen shell-and-tube and helical candidates before committing to long
100-second runs.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from .coupled_bangbang_hx import (
    CoupledBangBangHxConfig,
    DEFAULT_SETTINGS_FILE,
    default_feed_design,
    default_pressurant_system_config,
    load_settings,
    run_coupled_case,
)


def default_candidates(include_shell: bool = True) -> list[dict[str, Any]]:
    """Return a compact shell baseline + helical envelope sweep."""

    cases: list[dict[str, Any]] = []
    if include_shell:
        cases.append(
            {
                "name": "shellntube_echtherm_default",
                "hx_config": "shellntube",
            }
        )
    cases.extend(
        [
            {
                "name": "helical_id13p5_wall1p0_D60_L12",
                "hx_config": "shellnHelicalTube",
                "helical_inner_diameter_m": 13.5e-3,
                "helical_wall_thickness_m": 1.0e-3,
                "helical_centerline_diameter_m": 60.0e-3,
                "helical_pipe_length_m": 12.0,
            },
            {
                "name": "helical_id14p0_wall1p25_D65_L12",
                "hx_config": "shellnHelicalTube",
                "helical_inner_diameter_m": 14.0e-3,
                "helical_wall_thickness_m": 1.25e-3,
                "helical_centerline_diameter_m": 65.0e-3,
                "helical_pipe_length_m": 12.0,
            },
            {
                "name": "helical_id14p5_wall1p5_D70_L12",
                "hx_config": "shellnHelicalTube",
                "helical_inner_diameter_m": 14.5e-3,
                "helical_wall_thickness_m": 1.5e-3,
                "helical_centerline_diameter_m": 70.0e-3,
                "helical_pipe_length_m": 12.0,
            },
        ]
    )
    return cases


def run_sweep(
    *,
    settings: str = DEFAULT_SETTINGS_FILE,
    output_dir: str = "docs/validation/coupled_bangbang_hx_geometry_sweep",
    t_end_s: float = 1.0,
    hx_max_step_s: float = 0.003,
    hx_nodes: int = 5,
    hx_save_points: int = 101,
    include_shell: bool = True,
) -> dict[str, Any]:
    """Run the coupled HX geometry sweep and write summary artifacts."""

    if Path(settings).exists():
        base_cfg, system_config, feed_design = load_settings(settings)
    else:
        base_cfg = CoupledBangBangHxConfig()
        system_config = default_pressurant_system_config()
        feed_design = default_feed_design()

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    payloads: list[dict[str, Any]] = []
    for candidate in default_candidates(include_shell=include_shell):
        name = str(candidate["name"])
        updates = {
            key: value
            for key, value in candidate.items()
            if key != "name"
        }
        cfg = replace(
            base_cfg,
            output_dir=str(root / name),
            t_end_s=float(t_end_s),
            hx_max_step_s=float(hx_max_step_s),
            hx_nodes=int(hx_nodes),
            hx_save_points=int(hx_save_points),
            coolant_momentum_model="low_mach",
            **updates,
        )
        payload = run_coupled_case(cfg, system_config=system_config, feed_design=feed_design)
        payloads.append(payload)
        rows.append(_summary_row(name, cfg, payload))

    rows_sorted = sorted(rows, key=lambda row: row["score"])
    summary = {
        "settings": settings,
        "output_dir": str(root),
        "t_end_s": float(t_end_s),
        "hx_max_step_s": float(hx_max_step_s),
        "hx_nodes": int(hx_nodes),
        "hx_save_points": int(hx_save_points),
        "rows": rows_sorted,
    }
    (root / "sweep_summary.json").write_text(json.dumps(_jsonable(summary), indent=2), encoding="utf-8")
    _write_csv(root / "sweep_summary.csv", rows_sorted)
    _write_markdown(root / "SWEEP_SUMMARY.md", summary)
    return {"summary": summary, "payloads": payloads}


def summarize_existing_sweep(output_dir: str | Path) -> dict[str, Any]:
    """Rebuild aggregate sweep summary files from existing case `summary.json` files."""

    root = Path(output_dir)
    rows: list[dict[str, Any]] = []
    for summary_file in sorted(root.glob("*/summary.json")):
        payload = json.loads(summary_file.read_text(encoding="utf-8"))
        cfg = CoupledBangBangHxConfig(**payload["config"])
        rows.append(_summary_row(summary_file.parent.name, cfg, payload))
    rows_sorted = sorted(rows, key=lambda row: row["score"])
    summary = {
        "settings": "",
        "output_dir": str(root),
        "t_end_s": rows_sorted[0].get("t_end_s", "") if rows_sorted else "",
        "hx_max_step_s": "",
        "hx_nodes": "",
        "hx_save_points": "",
        "rows": rows_sorted,
    }
    if rows_sorted:
        first_payload = json.loads((Path(rows_sorted[0]["output_dir"]) / "summary.json").read_text(encoding="utf-8"))
        first_cfg = first_payload["config"]
        summary.update(
            {
                "t_end_s": first_cfg.get("t_end_s", ""),
                "hx_max_step_s": first_cfg.get("hx_max_step_s", ""),
                "hx_nodes": first_cfg.get("hx_nodes", ""),
                "hx_save_points": first_cfg.get("hx_save_points", ""),
            }
        )
    (root / "sweep_summary.json").write_text(json.dumps(_jsonable(summary), indent=2), encoding="utf-8")
    _write_csv(root / "sweep_summary.csv", rows_sorted)
    _write_markdown(root / "SWEEP_SUMMARY.md", summary)
    return summary


def _summary_row(name: str, config: CoupledBangBangHxConfig, payload: dict[str, Any]) -> dict[str, Any]:
    diag = payload.get("coupled_diagnostics", {})
    flags = _normalize_flags(payload.get("validation_flags", []), config)
    target_mdot = 0.5 * (config.min_target_helium_mdot_kg_s + config.max_target_helium_mdot_kg_s)
    mdot = _value(diag, "hx_mean_outlet_face_mdot_kg_s", _value(diag, "system_mean_helium_mdot_kg_s"))
    pressure_mean = _value(
        diag,
        "hx_mean_flowing_coolant_pressure_span_bar",
        _value(diag, "hx_mean_coolant_pressure_span_bar"),
    )
    pressure_p95 = _value(
        diag,
        "hx_p95_flowing_coolant_pressure_span_bar",
        _value(
            diag,
            "hx_max_coolant_pressure_span_bar",
            _value(diag, "hx_max_shell_pressure_drop_estimate_bar"),
        ),
    )
    pressure_span = _value(
        diag,
        "hx_max_flowing_coolant_pressure_span_bar",
        _value(
            diag,
            "hx_max_coolant_pressure_span_bar",
            _value(diag, "hx_max_shell_pressure_drop_estimate_bar"),
        ),
    )
    if not np.isfinite(pressure_mean):
        pressure_mean = pressure_span
    if not np.isfinite(pressure_p95):
        pressure_p95 = pressure_span
    wall_peak = _value(diag, "hx_peak_T_wall_K")
    flowing_T = _value(diag, "hx_max_flowing_T_c_out_K")
    mdot_score_value = _finite_or(mdot, target_mdot)
    pressure_score_value = _finite_or(pressure_p95, 0.0)
    pressure_mean_score_value = _finite_or(pressure_mean, 7.5)
    wall_score_value = _finite_or(wall_peak, 0.0)
    flowing_T_score_value = _finite_or(flowing_T, config.max_target_helium_outlet_temperature_K)
    score = (
        100.0 * abs(mdot_score_value - target_mdot)
        + 2.0 * max(pressure_score_value - config.max_hx_pressure_loss_bar, 0.0)
        + 1.0 * max(5.0 - pressure_mean_score_value, 0.0)
        + 0.02 * max(wall_score_value - config.max_wall_temperature_K, 0.0)
        + 0.01 * abs(flowing_T_score_value - config.max_target_helium_outlet_temperature_K)
        + 10.0 * sum("No coupled sanity-check" not in flag for flag in flags)
    )
    return {
        "name": name,
        "hx_config": config.hx_config,
        "score": float(score),
        "runtime_s": payload.get("hx_runtime_s"),
        "mean_helium_mdot_kg_s": mdot,
        "max_flowing_helium_outlet_K": flowing_T,
        "peak_wall_K": wall_peak,
        "mean_flowing_coolant_pressure_span_bar": pressure_mean,
        "p95_flowing_coolant_pressure_span_bar": pressure_p95,
        "max_coolant_pressure_span_bar": pressure_span,
        "mean_water_flow_L_s": _value(payload.get("system_summary", {}), "mean_water_flow_L_s"),
        "mean_tank_pressure_bar": _value(payload.get("system_summary", {}), "mean_pressure_bar"),
        "flags": " | ".join(flags),
        "output_dir": str(Path(config.output_dir)),
        "helical_id_mm": config.helical_inner_diameter_m * 1e3 if config.hx_config == "shellnHelicalTube" else "",
        "helical_wall_mm": config.helical_wall_thickness_m * 1e3 if config.hx_config == "shellnHelicalTube" else "",
        "helical_centerline_diameter_mm": config.helical_centerline_diameter_m * 1e3 if config.hx_config == "shellnHelicalTube" else "",
        "helical_pipe_length_m": config.helical_pipe_length_m if config.hx_config == "shellnHelicalTube" else "",
}


def _normalize_flags(flags: list[str], config: CoupledBangBangHxConfig) -> list[str]:
    if config.coolant_momentum_model != "low_mach":
        return list(flags)
    normalized = []
    for flag in flags:
        text = str(flag).replace(
            "HX shell-side pressure-drop estimate exceeds",
            "HX p95 flowing coolant pressure span exceeds",
        )
        text = text.replace(
            "HX coolant pressure span exceeds",
            "HX p95 flowing coolant pressure span exceeds",
        )
        normalized.append(text)
    return normalized


def _write_markdown(path: Path, summary: dict[str, Any]) -> None:
    rows = summary["rows"]
    lines = [
        "# Coupled Bang-Bang HX Geometry Sweep",
        "",
        f"- Simulated duration per case: {summary['t_end_s']} s",
        f"- HX max step: {summary['hx_max_step_s']} s",
        f"- HX nodes: {summary['hx_nodes']}",
        "",
        "| Rank | Case | HX | Score | He mdot kg/s | He out max K | Wall peak K | dp mean/p95/max bar | Runtime s |",
        "|---:|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for rank, row in enumerate(rows, start=1):
        lines.append(
            "| "
            f"{rank} | {row['name']} | {row['hx_config']} | {row['score']:.3g} | "
            f"{_fmt(row['mean_helium_mdot_kg_s'])} | {_fmt(row['max_flowing_helium_outlet_K'])} | "
            f"{_fmt(row['peak_wall_K'])} | "
            f"{_fmt(row['mean_flowing_coolant_pressure_span_bar'])}/"
            f"{_fmt(row['p95_flowing_coolant_pressure_span_bar'])}/"
            f"{_fmt(row['max_coolant_pressure_span_bar'])} | "
            f"{_fmt(row['runtime_s'])} |"
        )
    lines.extend(["", "## Notes", ""])
    lines.append(
        "- For low-Mach cases, the pressure column is the flowing-period coolant "
        "pressure span as `mean/p95/max`; target flags are based on the p95 value. "
        "Legacy summaries that only saved max span repeat the max value in all "
        "three pressure slots."
    )
    for row in rows:
        lines.append(f"- `{row['name']}`: {row['flags']}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    columns = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _fmt(value: Any, spec: str = ".4g") -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return ""
    if not np.isfinite(numeric):
        return ""
    return format(numeric, spec)


def _value(mapping: dict[str, Any], key: str, default: float = float("nan")) -> float:
    try:
        value = mapping.get(key, default)
    except AttributeError:
        return float(default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def _finite_or(value: float, fallback: float) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return float(fallback)
    return numeric if np.isfinite(numeric) else float(fallback)


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_jsonable(val) for val in value]
    if isinstance(value, tuple):
        return [_jsonable(val) for val in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, float):
        return value if np.isfinite(value) else str(value)
    return value


def main(argv: list[str] | None = None) -> dict[str, Any]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", default=DEFAULT_SETTINGS_FILE)
    parser.add_argument("--output-dir", default="docs/validation/coupled_bangbang_hx_geometry_sweep")
    parser.add_argument("--t-end", type=float, default=1.0)
    parser.add_argument("--hx-max-step", type=float, default=0.003)
    parser.add_argument("--hx-nodes", type=int, default=5)
    parser.add_argument("--hx-save-points", type=int, default=101)
    parser.add_argument("--no-shell-baseline", action="store_true")
    parser.add_argument("--summarize-existing", action="store_true")
    args = parser.parse_args(argv)

    if args.summarize_existing:
        summary = summarize_existing_sweep(args.output_dir)
        output = Path(args.output_dir)
        print(f"Wrote {output / 'sweep_summary.json'}")
        print(f"Wrote {output / 'SWEEP_SUMMARY.md'}")
        return {"summary": summary, "payloads": []}

    result = run_sweep(
        settings=args.settings,
        output_dir=args.output_dir,
        t_end_s=args.t_end,
        hx_max_step_s=args.hx_max_step,
        hx_nodes=args.hx_nodes,
        hx_save_points=args.hx_save_points,
        include_shell=not args.no_shell_baseline,
    )
    output = Path(args.output_dir)
    print(f"Wrote {output / 'sweep_summary.json'}")
    print(f"Wrote {output / 'SWEEP_SUMMARY.md'}")
    for row in result["summary"]["rows"]:
        print(
            f"{row['name']}: score={row['score']:.3g}, "
            f"mdot={row['mean_helium_mdot_kg_s']:.4g} kg/s, "
            f"THe,max={row['max_flowing_helium_outlet_K']:.1f} K, "
            f"Twall={row['peak_wall_K']:.1f} K"
        )
    return result


if __name__ == "__main__":
    main()

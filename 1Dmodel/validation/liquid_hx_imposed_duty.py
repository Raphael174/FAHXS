"""HX-style imposed-duty validation for the liquid heated-channel adapter.

This module does not model a specific helical or shell-and-tube geometry. It
mimics the interface those solvers will provide after their wall solve:
segment edges, local hydraulic geometry, and per-segment heat duty ``dQ``.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from hps_combustor.input_data import coolantProp
from hps_combustor.physics.liquid_flow.governing_equations import (
    HeatedChannelProfileCase,
    heated_channel_cell_fields,
    solve_steady_heated_channel_profile,
    summarize_heated_channel_result,
)
from hps_combustor.physics.liquid_flow.correlations import saturation_state


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_lut_path() -> Path:
    return _repo_root() / "docs" / "reference" / "external" / "2006LUTdata.txt"


def imposed_duty_reference_case() -> HeatedChannelProfileCase:
    """Return a nonuniform-grid, nonuniform-duty water case for HX adapter tests."""
    n_cells = 48
    xi = np.linspace(0.0, 1.0, n_cells + 1)
    z_edges = 1.2 * xi**1.15
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

    diameter = 0.006 * (1.0 + 0.08 * np.sin(np.pi * z_centers / z_edges[-1]))
    area = np.pi * diameter**2 / 4.0
    perimeter = np.pi * diameter

    p_in = 1.0e6
    sat = saturation_state("Water", p_in)
    h_in = sat.h_l_J_kg - 5.0e4

    base_heat_flux = 1.8e5 * (1.0 + 0.35 * np.sin(np.pi * z_centers / z_edges[-1]))
    heat_per_segment = base_heat_flux * perimeter * np.diff(z_edges)

    return HeatedChannelProfileCase(
        coolant_prop=coolantProp(coolant="Water", coolant_model="equilibrium_liquid"),
        z_edges_m=z_edges,
        hydraulic_diameter_m=diameter,
        flow_area_m2=area,
        heated_perimeter_m=perimeter,
        mass_flow_kg_s=0.012,
        p_in_Pa=p_in,
        h_in_J_kg=h_in,
        heat_per_segment_W=heat_per_segment,
        lut_path=_default_lut_path(),
    )


def _write_fields_csv(path: Path, fields: dict[str, np.ndarray]) -> None:
    keys = list(fields.keys())
    rows = zip(*(fields[key] for key in keys))
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(keys)
        writer.writerows(rows)


def generate_imposed_duty_report(output_dir: str | Path | None = None) -> dict[str, object]:
    """Run the HX-style imposed-duty case and write comparison artifacts."""
    output = (
        Path(output_dir)
        if output_dir is not None
        else _repo_root() / "docs" / "validation" / "liquid_hx_imposed_duty"
    )
    output.mkdir(parents=True, exist_ok=True)

    case = imposed_duty_reference_case()
    result = solve_steady_heated_channel_profile(case)
    fields = heated_channel_cell_fields(result)
    diagnostics = summarize_heated_channel_result(result, min_pressure_Pa=case.min_pressure_Pa)

    fields_csv = output / "liquid_hx_imposed_duty_fields.csv"
    _write_fields_csv(fields_csv, fields)

    z = result.z_m
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), sharex=True)
    axes[0, 0].plot(z, result.T_K, color="#1f6f78")
    axes[0, 0].set_ylabel("coolant T [K]")
    axes[0, 1].plot(z, result.quality, color="#8f2d56")
    axes[0, 1].axhline(0.0, color="#777777", linewidth=0.8, linestyle=":")
    axes[0, 1].axhline(1.0, color="#777777", linewidth=0.8, linestyle=":")
    axes[0, 1].set_ylabel("quality [-]")
    axes[1, 0].plot(z, (result.p_Pa[0] - result.p_Pa) / 1000.0, color="#264653")
    axes[1, 0].set_ylabel("pressure drop [kPa]")
    axes[1, 0].set_xlabel("z [m]")
    axes[1, 1].plot(z, result.htc_W_m2_K / 1000.0, color="#e76f51")
    axes[1, 1].set_ylabel("HTC [kW/m2/K]")
    axes[1, 1].set_xlabel("z [m]")
    for ax in axes.ravel():
        ax.grid(True, alpha=0.3)
    fig.suptitle("HX-style imposed-duty liquid coolant adapter")
    fig.tight_layout()
    profiles_png = output / "liquid_hx_imposed_duty_profiles.png"
    fig.savefig(profiles_png, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(fields["z_m"], fields["heat_flux_W_m2"] / 1000.0, color="#2a9d8f", label="imposed heat flux")
    if np.any(np.isfinite(fields["chf_W_m2"])):
        ax.plot(fields["z_m"], fields["chf_W_m2"] / 1000.0, color="#8f2d56", label="Groeneveld CHF")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("heat flux [kW/m2]")
    ax.set_title("Imposed duty and CHF margin")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    chf_png = output / "liquid_hx_imposed_duty_chf.png"
    fig.savefig(chf_png, dpi=180)
    plt.close(fig)

    summary = {
        "source": "synthetic HX-style imposed-duty adapter case",
        "purpose": "validate integration contract from wall-solver dQ segments to p-h liquid coolant march",
        "n_cells": int(len(case.z_edges_m) - 1),
        "heat_rate_W": diagnostics.heat_rate_W,
        "pressure_drop_Pa": diagnostics.pressure_drop_Pa,
        "inlet_T_K": diagnostics.inlet_T_K,
        "outlet_T_K": diagnostics.outlet_T_K,
        "min_quality": diagnostics.min_quality,
        "max_quality": diagnostics.max_quality,
        "outlet_quality": diagnostics.outlet_quality,
        "max_void_fraction": diagnostics.max_void_fraction,
        "min_chf_margin": diagnostics.min_chf_margin,
        "boiling_reached": diagnostics.boiling_reached,
        "dryout_or_vapor_reached": diagnostics.dryout_or_vapor_reached,
        "chf_margin_below_limit": diagnostics.chf_margin_below_limit,
        "energy_residual_abs_J_kg": diagnostics.energy_residual_abs_J_kg,
        "energy_residual_ok": diagnostics.energy_residual_ok,
        "outputs": [
            fields_csv.name,
            profiles_png.name,
            chf_png.name,
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    report = generate_imposed_duty_report()
    for key, value in report.items():
        if isinstance(value, float):
            print(f"{key}: {value:.6g}")
        else:
            print(f"{key}: {value}")

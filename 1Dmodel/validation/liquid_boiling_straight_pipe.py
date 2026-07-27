"""Straight-pipe liquid/boiling proof-of-concept validation case.

This module deliberately avoids the helical and shell-and-tube architectures.
It exercises the governing variables and closures needed before integrating
liquid coolants into HX-specific solvers:

    state = p, h
    heat input -> dh/dz
    EOS -> phase, density, quality, void fraction
    pressure drop -> single-phase or equilibrium two-phase
    boiling HTC -> Gungor-Winterton in saturated two-phase flow
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import math
from pathlib import Path

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np

from hps_combustor.physics.liquid_flow.governing_equations import (
    HeatedChannelCase,
    solve_steady_heated_channel,
)
from hps_combustor.physics.liquid_flow.chf import (
    groeneveld_2006_chf,
    load_groeneveld_2006_lut,
)
from hps_combustor.physics.liquid_flow.correlations import (
    chisholm_two_phase_multiplier,
    darcy_friction_smooth_pipe,
    equilibrium_state_ph,
    gungor_winterton_boiling_htc,
    homogeneous_acceleration_pressure_gradient,
    liquid_single_phase_nusselt,
    muller_steinhagen_heck_friction_gradient,
    saturation_state,
    yu2002_modified_anl_boiling_htc,
    yu2002_small_channel_pressure_multiplier,
)


@dataclass(frozen=True)
class StraightPipeBoilingCase:
    fluid: str = "Water"
    length_m: float = 1.0
    diameter_m: float = 0.010
    mass_flow_kg_s: float = 0.05
    p_in_Pa: float = 1.0e5
    h_in_J_kg: float | None = None
    T_in_K: float = 370.0
    heat_flux_W_m2: float = 1.0e5
    n_cells: int = 80


@dataclass(frozen=True)
class StraightPipeBoilingResult:
    z_m: np.ndarray
    p_Pa: np.ndarray
    h_J_kg: np.ndarray
    T_K: np.ndarray
    quality: np.ndarray
    void_fraction: np.ndarray
    rho_kg_m3: np.ndarray
    htc_W_m2_K: np.ndarray
    dpdz_friction_Pa_m: np.ndarray
    dpdz_acceleration_Pa_m: np.ndarray
    heat_rate_W: float
    energy_residual_J_kg: float

    @property
    def outlet_quality(self) -> float:
        return float(self.quality[-1])

    @property
    def pressure_drop_Pa(self) -> float:
        return float(self.p_Pa[0] - self.p_Pa[-1])


GROENEVELD_PAGE9_REFERENCE_POINTS = (
    # Values transcribed from Groeneveld et al. 2007 PDF page 9, first pressure block.
    # Units: p [MPa], G [kg/m2/s], x [-], CHF [kW/m2].
    (0.10, 0.0, -0.50, 8111.0),
    (0.10, 0.0, -0.20, 4802.0),
    (0.10, 0.0, 0.00, 1142.0),
    (0.10, 0.0, 0.25, 188.0),
    (0.10, 0.0, 0.50, 123.0),
    (0.10, 0.0, 0.90, 55.0),
    (0.10, 50.0, -0.50, 8317.0),
    (0.10, 50.0, -0.20, 5035.0),
    (0.10, 50.0, 0.00, 1570.0),
    (0.10, 50.0, 0.25, 553.0),
    (0.10, 50.0, 0.50, 387.0),
    (0.10, 50.0, 0.90, 204.0),
    (0.10, 100.0, -0.50, 8390.0),
    (0.10, 100.0, -0.20, 5322.0),
    (0.10, 100.0, 0.00, 2103.0),
    (0.10, 100.0, 0.25, 847.0),
    (0.10, 100.0, 0.50, 715.0),
    (0.10, 100.0, 0.90, 359.0),
)


YU2002_PRESSURE_MULTIPLIER_DIGITIZED = (
    # Representative points digitized visually from Yu et al. 2002 Figs. 7-8.
    # The paper equations are for phi_l^2; the plotted ordinate is treated as
    # phi_l magnitude here. These points are suitable for trend/regression
    # validation, not as a replacement for original tabular data.
    # X, phi_l_exp, mass_flux_kg_m2_s
    (0.0065, 118.0, 151.0),
    (0.0080, 96.0, 129.0),
    (0.0100, 78.0, 151.0),
    (0.0150, 54.0, 129.0),
    (0.0200, 40.0, 103.0),
    (0.0300, 27.0, 103.0),
    (0.0500, 17.0, 76.0),
    (0.0750, 11.5, 76.0),
    (0.1000, 8.3, 50.0),
    (0.1500, 5.8, 50.0),
)


YU2002_HTC_COMPARISON_DIGITIZED = (
    # Representative points digitized visually from Yu et al. 2002 Fig. 10
    # comparing experimental local HTC to the modified ANL prediction.
    # h_exp_W_m2_K, h_pred_paper_W_m2_K
    (9000.0, 9800.0),
    (12000.0, 13200.0),
    (15000.0, 16200.0),
    (18000.0, 18800.0),
    (22000.0, 21800.0),
    (26000.0, 25200.0),
    (30000.0, 29200.0),
    (34000.0, 32900.0),
    (38000.0, 36500.0),
    (42000.0, 43800.0),
)


YU2002_CHF_TREND_DIGITIZED = (
    # Representative points digitized visually from Yu et al. 2002 Fig. 12.
    # x_exit, q_chf_kW_m2. Used only for monotonic trend checks.
    (0.52, 185.0),
    (0.60, 165.0),
    (0.70, 135.0),
    (0.82, 105.0),
    (0.93, 78.0),
)


def inlet_enthalpy(case: StraightPipeBoilingCase) -> float:
    if case.h_in_J_kg is not None:
        return case.h_in_J_kg
    return CP.PropsSI("H", "P", case.p_in_Pa, "T", case.T_in_K, case.fluid)


def _single_phase_pressure_gradient(
    *,
    p_Pa: float,
    h_J_kg: float,
    mass_flux_kg_m2_s: float,
    diameter_m: float,
    fluid: str,
) -> float:
    rho = CP.PropsSI("D", "P", p_Pa, "H", h_J_kg, fluid)
    mu = CP.PropsSI("V", "P", p_Pa, "H", h_J_kg, fluid)
    Re = mass_flux_kg_m2_s * diameter_m / mu
    f = darcy_friction_smooth_pipe(Re)
    return f * mass_flux_kg_m2_s**2 / (2.0 * diameter_m * rho)


def _single_phase_htc(
    *,
    p_Pa: float,
    h_J_kg: float,
    mass_flux_kg_m2_s: float,
    diameter_m: float,
    fluid: str,
) -> float:
    mu = CP.PropsSI("V", "P", p_Pa, "H", h_J_kg, fluid)
    k = CP.PropsSI("L", "P", p_Pa, "H", h_J_kg, fluid)
    cp = CP.PropsSI("C", "P", p_Pa, "H", h_J_kg, fluid)
    Re = mass_flux_kg_m2_s * diameter_m / mu
    Pr = cp * mu / k
    Nu = liquid_single_phase_nusselt(Re, Pr)
    return Nu * k / diameter_m


def solve_steady_straight_pipe(case: StraightPipeBoilingCase) -> StraightPipeBoilingResult:
    if case.length_m <= 0.0 or case.diameter_m <= 0.0:
        raise ValueError("length and diameter must be positive")
    if case.mass_flow_kg_s <= 0.0:
        raise ValueError("mass flow must be positive")
    if case.n_cells < 1:
        raise ValueError("n_cells must be at least 1")

    channel = HeatedChannelCase(
        fluid=case.fluid,
        length_m=case.length_m,
        hydraulic_diameter_m=case.diameter_m,
        flow_area_m2=math.pi * case.diameter_m**2 / 4.0,
        heated_perimeter_m=math.pi * case.diameter_m,
        mass_flow_kg_s=case.mass_flow_kg_s,
        p_in_Pa=case.p_in_Pa,
        h_in_J_kg=case.h_in_J_kg,
        T_in_K=case.T_in_K,
        heat_flux_W_m2=case.heat_flux_W_m2,
        n_cells=case.n_cells,
        coolant_model="equilibrium_liquid",
        lut_path=_default_lut_path(),
    )
    result = solve_steady_heated_channel(channel)
    return StraightPipeBoilingResult(
        z_m=result.z_m,
        p_Pa=result.p_Pa,
        h_J_kg=result.h_J_kg,
        T_K=result.T_K,
        quality=result.quality,
        void_fraction=result.void_fraction,
        rho_kg_m3=result.rho_kg_m3,
        htc_W_m2_K=result.htc_W_m2_K,
        dpdz_friction_Pa_m=result.dpdz_friction_Pa_m,
        dpdz_acceleration_Pa_m=result.dpdz_acceleration_Pa_m,
        heat_rate_W=result.heat_rate_W,
        energy_residual_J_kg=result.energy_residual_J_kg,
    )


def saturated_water_reference_case() -> StraightPipeBoilingCase:
    """Return a simple saturated-inlet case for regression and manual checks."""
    p_in = 5.0e6
    sat = saturation_state("Water", p_in)
    h_in = sat.h_l_J_kg + 0.05 * sat.h_fg_J_kg
    return StraightPipeBoilingCase(
        fluid="Water",
        length_m=1.0,
        diameter_m=0.020,
        mass_flow_kg_s=0.20,
        p_in_Pa=p_in,
        h_in_J_kg=h_in,
        T_in_K=sat.T_sat_K,
        heat_flux_W_m2=1.0e5,
        n_cells=80,
    )


def run_reference_case() -> dict[str, float]:
    case = saturated_water_reference_case()
    result = solve_steady_straight_pipe(case)
    return {
        "heat_rate_W": result.heat_rate_W,
        "outlet_quality": result.outlet_quality,
        "pressure_drop_Pa": result.pressure_drop_Pa,
        "max_void_fraction": float(np.max(result.void_fraction)),
        "mean_boiling_htc_W_m2_K": float(np.mean(result.htc_W_m2_K)),
        "energy_residual_J_kg": result.energy_residual_J_kg,
    }


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_lut_path() -> Path:
    return _repo_root() / "docs" / "reference" / "external" / "2006LUTdata.txt"


def _write_pipe_csv(path: Path, result: StraightPipeBoilingResult) -> None:
    rows = zip(
        result.z_m,
        result.p_Pa,
        result.h_J_kg,
        result.T_K,
        result.quality,
        result.void_fraction,
        result.rho_kg_m3,
        result.htc_W_m2_K,
        result.dpdz_friction_Pa_m,
        result.dpdz_acceleration_Pa_m,
    )
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "z_m",
                "p_Pa",
                "h_J_kg",
                "T_K",
                "quality",
                "void_fraction",
                "rho_kg_m3",
                "htc_W_m2_K",
                "dpdz_friction_Pa_m",
                "dpdz_acceleration_Pa_m",
            ]
        )
        writer.writerows(rows)


def _lookup_table_value(table: np.ndarray, p_axis: np.ndarray, g_axis: np.ndarray, x_axis: np.ndarray, p: float, g: float, x: float) -> float:
    i_p = int(np.where(np.isclose(p_axis, p))[0][0])
    i_g = int(np.where(np.isclose(g_axis, g))[0][0])
    i_x = int(np.where(np.isclose(x_axis, x))[0][0])
    return float(table[i_p, i_g, i_x])


def _write_groeneveld_page9_comparison(output_dir: Path, lut_path: Path) -> list[dict[str, float]]:
    p_axis, g_axis, x_axis, table = load_groeneveld_2006_lut(lut_path)
    rows = []
    for p, g, x, paper_chf in GROENEVELD_PAGE9_REFERENCE_POINTS:
        ingested = _lookup_table_value(table, p_axis, g_axis, x_axis, p, g, x)
        rows.append(
            {
                "p_MPa": p,
                "mass_flux_kg_m2_s": g,
                "quality": x,
                "paper_page9_chf_kW_m2": paper_chf,
                "ingested_lut_chf_kW_m2": ingested,
                "relative_error": 0.0 if paper_chf == 0.0 else (ingested - paper_chf) / paper_chf,
            }
        )
    csv_path = output_dir / "groeneveld_page9_comparison.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    paper = np.array([r["paper_page9_chf_kW_m2"] for r in rows])
    ingested = np.array([r["ingested_lut_chf_kW_m2"] for r in rows])
    ax.scatter(paper, ingested, color="#1f6f78", edgecolor="white", linewidth=0.6)
    lim = [0.0, max(float(np.max(paper)), float(np.max(ingested))) * 1.05]
    ax.plot(lim, lim, color="#aa3a2a", linewidth=1.2, label="exact match")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Groeneveld 2007 page 9 CHF [kW/m2]")
    ax.set_ylabel("ingested GitHub LUT CHF [kW/m2]")
    ax.set_title("2006 CHF LUT ingestion check")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "groeneveld_page9_comparison.png", dpi=180)
    plt.close(fig)
    return rows


def _mean_abs_relative_error(observed: np.ndarray, predicted: np.ndarray) -> float:
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    return float(np.mean(np.abs((predicted - observed) / observed)))


def generate_yu2002_validation_report(output_dir: str | Path | None = None) -> dict[str, object]:
    """Generate an executable validation report based on Yu et al. 2002.

    The local paper does not include machine-readable original tables. This
    report therefore uses explicit representative digitized points from the
    published plots and validates the implemented equations/trends against
    those points. The output files label the digitization source and should be
    replaced with author/source data if it becomes available.
    """
    output = (
        Path(output_dir)
        if output_dir is not None
        else _repo_root() / "docs" / "validation" / "liquid_boiling_yu2002"
    )
    output.mkdir(parents=True, exist_ok=True)

    pressure_rows = []
    for X, phi_exp, G in YU2002_PRESSURE_MULTIPLIER_DIGITIZED:
        phi_yu = math.sqrt(yu2002_small_channel_pressure_multiplier(X))
        phi_chisholm = math.sqrt(chisholm_two_phase_multiplier(X, C=12.0))
        pressure_rows.append(
            {
                "source": "Yu2002_Fig7_Fig8_digitized_representative",
                "martinelli_X": X,
                "mass_flux_kg_m2_s": G,
                "phi_l_exp_digitized": phi_exp,
                "phi_l_yu2002_fit": phi_yu,
                "phi_l_chisholm_C12": phi_chisholm,
                "rel_error_yu2002_fit": (phi_yu - phi_exp) / phi_exp,
                "rel_error_chisholm": (phi_chisholm - phi_exp) / phi_exp,
            }
        )
    pressure_csv = output / "yu2002_pressure_multiplier_digitized.csv"
    with pressure_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(pressure_rows[0].keys()))
        writer.writeheader()
        writer.writerows(pressure_rows)

    pressure_x = np.array([r["martinelli_X"] for r in pressure_rows])
    pressure_exp = np.array([r["phi_l_exp_digitized"] for r in pressure_rows])
    pressure_yu = np.array([r["phi_l_yu2002_fit"] for r in pressure_rows])
    pressure_ch = np.array([r["phi_l_chisholm_C12"] for r in pressure_rows])

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    x_line = np.logspace(math.log10(0.005), math.log10(0.25), 200)
    ax.scatter(pressure_x, pressure_exp, color="#1f6f78", label="digitized Yu 2002 points")
    ax.plot(x_line, np.sqrt([yu2002_small_channel_pressure_multiplier(x) for x in x_line]), color="#e76f51", label="Yu 2002 fit")
    ax.plot(x_line, np.sqrt([chisholm_two_phase_multiplier(x, C=12.0) for x in x_line]), color="#264653", linestyle="--", label="Chisholm C=12")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Lockhart-Martinelli parameter X [-]")
    ax.set_ylabel("pressure multiplier magnitude [-]")
    ax.set_title("Yu 2002 pressure-multiplier validation")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "yu2002_pressure_multiplier_comparison.png", dpi=180)
    plt.close(fig)

    htc_rows = []
    for h_exp, h_pred_paper in YU2002_HTC_COMPARISON_DIGITIZED:
        htc_rows.append(
            {
                "source": "Yu2002_Fig10_digitized_representative",
                "h_exp_digitized_W_m2_K": h_exp,
                "h_pred_paper_digitized_W_m2_K": h_pred_paper,
                "rel_error_digitized_prediction": (h_pred_paper - h_exp) / h_exp,
            }
        )
    htc_csv = output / "yu2002_htc_fig10_digitized.csv"
    with htc_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(htc_rows[0].keys()))
        writer.writeheader()
        writer.writerows(htc_rows)

    h_exp = np.array([r["h_exp_digitized_W_m2_K"] for r in htc_rows])
    h_pred = np.array([r["h_pred_paper_digitized_W_m2_K"] for r in htc_rows])
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    ax.scatter(h_exp / 1000.0, h_pred / 1000.0, color="#2a9d8f", edgecolor="white", linewidth=0.6)
    lim = [0.0, max(float(np.max(h_exp)), float(np.max(h_pred))) / 1000.0 * 1.08]
    ax.plot(lim, lim, color="#222222", linewidth=1.1, label="perfect agreement")
    ax.plot(lim, [1.3 * v for v in lim], color="#777777", linestyle=":", linewidth=1.0, label="+/-30%")
    ax.plot(lim, [0.7 * v for v in lim], color="#777777", linestyle=":", linewidth=1.0)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("digitized experimental HTC [kW/m2/K]")
    ax.set_ylabel("digitized predicted HTC [kW/m2/K]")
    ax.set_title("Yu 2002 modified ANL HTC comparison")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "yu2002_htc_fig10_comparison.png", dpi=180)
    plt.close(fig)

    htc_param_rows = []
    for q in (40_000.0, 80_000.0, 120_000.0, 160_000.0):
        for G in (50.0, 103.0, 151.0):
            h_model = yu2002_modified_anl_boiling_htc(
                p_Pa=200_000.0,
                mass_flux_kg_m2_s=G,
                diameter_m=0.00298,
                heat_flux_W_m2=q,
                fluid="Water",
            )
            htc_param_rows.append(
                {
                    "p_Pa": 200_000.0,
                    "diameter_m": 0.00298,
                    "mass_flux_kg_m2_s": G,
                    "heat_flux_W_m2": q,
                    "yu2002_modified_anl_htc_W_m2_K": h_model,
                }
            )
    htc_param_csv = output / "yu2002_modified_anl_parametric_htc.csv"
    with htc_param_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(htc_param_rows[0].keys()))
        writer.writeheader()
        writer.writerows(htc_param_rows)

    chf_rows = []
    for x, q_chf in YU2002_CHF_TREND_DIGITIZED:
        chf_rows.append(
            {
                "source": "Yu2002_Fig12_digitized_representative",
                "exit_quality": x,
                "chf_kW_m2_digitized": q_chf,
            }
        )
    chf_csv = output / "yu2002_chf_trend_digitized.csv"
    with chf_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(chf_rows[0].keys()))
        writer.writeheader()
        writer.writerows(chf_rows)

    chf_x = np.array([r["exit_quality"] for r in chf_rows])
    chf_q = np.array([r["chf_kW_m2_digitized"] for r in chf_rows])
    fig, ax = plt.subplots(figsize=(6.6, 4.5))
    ax.plot(chf_x, chf_q, marker="o", color="#8f2d56")
    ax.set_xlabel("exit quality [-]")
    ax.set_ylabel("CHF [kW/m2]")
    ax.set_title("Yu 2002 CHF trend, digitized representative points")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output / "yu2002_chf_trend_digitized.png", dpi=180)
    plt.close(fig)

    summary = {
        "source": "Yu et al. 2002, International Journal of Multiphase Flow 28, 927-941",
        "data_status": "representative points digitized from local PDF plots; replace with source tables if acquired",
        "pressure_multiplier_mean_abs_rel_error_yu2002_fit": _mean_abs_relative_error(
            pressure_exp, pressure_yu
        ),
        "pressure_multiplier_mean_abs_rel_error_chisholm": _mean_abs_relative_error(
            pressure_exp, pressure_ch
        ),
        "htc_fig10_digitized_mean_abs_rel_error": _mean_abs_relative_error(h_exp, h_pred),
        "chf_digitized_monotonic_decrease_with_quality": bool(np.all(np.diff(chf_q) < 0.0)),
        "outputs": [
            pressure_csv.name,
            "yu2002_pressure_multiplier_comparison.png",
            htc_csv.name,
            "yu2002_htc_fig10_comparison.png",
            htc_param_csv.name,
            chf_csv.name,
            "yu2002_chf_trend_digitized.png",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def generate_validation_report(
    output_dir: str | Path | None = None,
    lut_path: str | Path | None = None,
) -> dict[str, object]:
    """Generate figures and comparison data for the liquid boiling PoC."""
    output = Path(output_dir) if output_dir is not None else _repo_root() / "docs" / "validation" / "liquid_boiling_poc"
    output.mkdir(parents=True, exist_ok=True)
    lut = Path(lut_path) if lut_path is not None else _default_lut_path()

    case = saturated_water_reference_case()
    result = solve_steady_straight_pipe(case)
    _write_pipe_csv(output / "straight_pipe_profiles.csv", result)

    area_flow = math.pi * case.diameter_m**2 / 4.0
    mass_flux = case.mass_flow_kg_s / area_flow
    chf = np.array(
        [
            groeneveld_2006_chf(
                p_Pa=p,
                mass_flux_kg_m2_s=mass_flux,
                quality=float(np.clip(x, -0.5, 1.0)),
                diameter_m=case.diameter_m,
                lut_path=lut,
            )
            for p, x in zip(result.p_Pa, result.quality)
        ]
    )
    chf_margin = chf / case.heat_flux_W_m2

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), sharex=True)
    axes[0, 0].plot(result.z_m, result.quality, color="#26547c")
    axes[0, 0].set_ylabel("thermodynamic quality [-]")
    axes[0, 1].plot(result.z_m, result.void_fraction, color="#2a9d8f")
    axes[0, 1].set_ylabel("HEM void fraction [-]")
    axes[1, 0].plot(result.z_m, (result.p_Pa[0] - result.p_Pa) / 1000.0, color="#8f2d56")
    axes[1, 0].set_ylabel("pressure drop [kPa]")
    axes[1, 0].set_xlabel("z [m]")
    axes[1, 1].plot(result.z_m, result.htc_W_m2_K / 1000.0, color="#e76f51", label="Gungor-Winterton HTC")
    axes[1, 1].set_ylabel("HTC [kW/m2/K]")
    axes[1, 1].set_xlabel("z [m]")
    for ax in axes.ravel():
        ax.grid(True, alpha=0.3)
    fig.suptitle("Straight-pipe saturated-water boiling PoC")
    fig.tight_layout()
    fig.savefig(output / "straight_pipe_profiles.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(result.z_m, np.full_like(result.z_m, case.heat_flux_W_m2 / 1000.0), label="applied heat flux")
    ax.plot(result.z_m, chf / 1000.0, label="Groeneveld 2006 CHF, diameter-corrected")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("heat flux [kW/m2]")
    ax.set_title("CHF margin along straight-pipe PoC")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "straight_pipe_chf_margin.png", dpi=180)
    plt.close(fig)

    qualities = np.linspace(0.01, 0.90, 120)
    htc = np.array(
        [
            gungor_winterton_boiling_htc(
                p_Pa=case.p_in_Pa,
                mass_flux_kg_m2_s=mass_flux,
                diameter_m=case.diameter_m,
                quality=float(x),
                heat_flux_W_m2=case.heat_flux_W_m2,
                fluid=case.fluid,
            )
            for x in qualities
        ]
    )
    dpdz = np.array(
        [
            muller_steinhagen_heck_friction_gradient(
                p_Pa=case.p_in_Pa,
                mass_flux_kg_m2_s=mass_flux,
                diameter_m=case.diameter_m,
                quality=float(x),
                fluid=case.fluid,
            )
            for x in qualities
        ]
    )
    fig, ax1 = plt.subplots(figsize=(7.0, 4.8))
    ax2 = ax1.twinx()
    ax1.plot(qualities, htc / 1000.0, color="#e76f51", label="Gungor-Winterton HTC")
    ax2.plot(qualities, dpdz / 1000.0, color="#264653", label="MSH friction gradient")
    ax1.set_xlabel("quality [-]")
    ax1.set_ylabel("HTC [kW/m2/K]", color="#e76f51")
    ax2.set_ylabel("friction gradient [kPa/m]", color="#264653")
    ax1.grid(True, alpha=0.3)
    ax1.set_title("Literature correlations used by the PoC")
    fig.tight_layout()
    fig.savefig(output / "correlation_quality_sweep.png", dpi=180)
    plt.close(fig)

    page9_rows = _write_groeneveld_page9_comparison(output, lut)
    summary = {
        "case": case.__dict__,
        "heat_rate_W": result.heat_rate_W,
        "outlet_quality": result.outlet_quality,
        "pressure_drop_Pa": result.pressure_drop_Pa,
        "energy_residual_J_kg": result.energy_residual_J_kg,
        "min_chf_margin": float(np.min(chf_margin)),
        "mean_boiling_htc_W_m2_K": float(np.mean(result.htc_W_m2_K)),
        "groeneveld_page9_max_abs_error_kW_m2": float(
            max(abs(r["ingested_lut_chf_kW_m2"] - r["paper_page9_chf_kW_m2"]) for r in page9_rows)
        ),
        "outputs": [
            "straight_pipe_profiles.csv",
            "straight_pipe_profiles.png",
            "straight_pipe_chf_margin.png",
            "correlation_quality_sweep.png",
            "groeneveld_page9_comparison.csv",
            "groeneveld_page9_comparison.png",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    summary = generate_validation_report()
    yu_summary = generate_yu2002_validation_report()
    print("[generic_straight_pipe_poc]")
    for key, value in summary.items():
        if isinstance(value, float):
            print(f"{key}: {value:.6g}")
        else:
            print(f"{key}: {value}")
    print("[yu2002_validation]")
    for key, value in yu_summary.items():
        if isinstance(value, float):
            print(f"{key}: {value:.6g}")
        else:
            print(f"{key}: {value}")

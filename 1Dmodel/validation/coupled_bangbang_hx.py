"""Coupled bang-bang pressurant + HX transient validation module.

Run from the repository root:

    python -m hps_combustor.validation.coupled_bangbang_hx

This is intentionally a separate validation/coupling sandbox. It keeps the
external feed/tank system as a cheap 0D model, then drives the detailed 1D
shell-and-tube transient HX with the resulting helium inlet boundary histories.
The output is meant for sanity checking coolant/material transients before the
HX model is embedded in a larger Simulink or ESPSS-style multiphysics loop.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import asdict, dataclass, fields as dataclass_fields, replace
from pathlib import Path
from typing import Any

import numpy as np

from ..main_transient import build_inputs, run_transient
from ..result_package import get_transient_time_series
from .pressurant_bangbang_sizing import (
    FeedDesign,
    PressurantSystemConfig,
    run_design,
    write_history_csv,
)


@dataclass(frozen=True)
class CoupledBangBangHxConfig:
    """Runtime controls for the coupled validation module."""

    output_dir: str = "docs/validation/coupled_bangbang_hx"
    t_end_s: float = 20.0
    system_dt_s: float = 0.005
    hx_max_step_s: float = 0.05
    hx_nodes: int = 20
    hx_save_points: int = 101
    hx_boundary_average_window_s: float = 0.05
    hx_config: str = "shellntube"
    flow_config: str = "counter"
    chemistry: str = "finite_rate"
    coolant_momentum_model: str = "quasi_steady"
    lox_mdot_kg_s: float = 0.04833333333333333
    diesel_mdot_kg_s: float = 0.024166666666666666
    hot_gas_control: str = "track_helium_mdot"
    hot_gas_nominal_helium_mdot_kg_s: float = 0.145
    hot_gas_min_mdot_kg_s: float = 0.0
    hot_gas_max_mdot_kg_s: float = 0.0725
    hot_gas_off_helium_mdot_kg_s: float = 1.0e-5
    low_mach_mdot_cap_kg_s: float = 0.35
    lox_temperature_K: float = 100.0
    wall_initial_temperature_K: float = 300.0
    max_wall_temperature_K: float = 1473.15
    max_hx_pressure_loss_bar: float = 10.0
    min_target_helium_mdot_kg_s: float = 0.120
    max_target_helium_mdot_kg_s: float = 0.150
    max_target_helium_outlet_temperature_K: float = 700.0
    target_water_flow_L_s: float = 30.0
    target_water_tank_pressure_bar: float = 70.0
    run_hx: bool = True
    helical_inner_diameter_m: float = 14.0e-3
    helical_wall_thickness_m: float = 1.25e-3
    helical_pipe_length_m: float = 12.0
    helical_centerline_diameter_m: float = 65.0e-3
    helical_coil_gap_m: float = 4.0e-3
    helical_material: str = "ST316L"
    helical_shell_nusselt_correction: float = 0.28


DEFAULT_SETTINGS_FILE = "inputs/coupled_bangbang_hx_dummy_validation.json"


def default_feed_design() -> FeedDesign:
    """Return the dummy hot-pressurant feed design for coupled validation."""

    return FeedDesign(
        n_branches=3,
        orifice_diameter_m=1.75e-3,
        valve_equivalent_diameter_m=1.75e-3,
        control_frequency_Hz=40.0,
    )


def default_pressurant_system_config() -> PressurantSystemConfig:
    """Return a 0D system config consistent with hot helium leaving the HX."""

    return PressurantSystemConfig(
        helium_tank_temperature_K=90.0,
        pressurant_temperature_K=700.0,
        target_line_pressure_Pa=77.5e5,
        hx_pressure_loss_nominal_Pa=7.5e5,
        hx_nominal_helium_mdot_kg_s=0.15,
    )


def load_settings(path: str | Path) -> tuple[CoupledBangBangHxConfig, PressurantSystemConfig, FeedDesign]:
    """Load a coupled validation settings JSON file."""

    source = Path(path)
    payload = json.loads(source.read_text(encoding="utf-8"))
    cfg = _dataclass_from_mapping(CoupledBangBangHxConfig, payload.get("coupled", {}))
    system = replace(
        default_pressurant_system_config(),
        **_select_dataclass_fields(PressurantSystemConfig, payload.get("pressurant", {})),
    )
    design = replace(
        default_feed_design(),
        **_select_dataclass_fields(FeedDesign, payload.get("feed_design", {})),
    )
    return cfg, system, design


def write_settings(
    path: str | Path,
    config: CoupledBangBangHxConfig,
    system_config: PressurantSystemConfig,
    feed_design: FeedDesign,
) -> None:
    """Write the exact coupled validation settings used by a run."""

    payload = {
        "description": (
            "Dummy coupled bang-bang validation settings. The external "
            "pressurant/feed/water-tank system is 0D; the HX is the 1D "
            "transient coolant/material HX solver."
        ),
        "coupled": asdict(config),
        "pressurant": asdict(system_config),
        "feed_design": asdict(feed_design),
    }
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(_jsonable(payload), indent=2), encoding="utf-8")


def run_coupled_case(
    config: CoupledBangBangHxConfig | None = None,
    *,
    system_config: PressurantSystemConfig | None = None,
    feed_design: FeedDesign | None = None,
) -> dict[str, Any]:
    """Run the 0D pressurant system and, optionally, the detailed HX transient."""

    cfg = config or CoupledBangBangHxConfig()
    sys_cfg = replace(
        system_config or default_pressurant_system_config(),
        t_end_s=float(cfg.t_end_s),
        dt_s=float(cfg.system_dt_s),
    )
    design = feed_design or default_feed_design()
    output = Path(cfg.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    write_settings(output / "settings_used.json", cfg, sys_cfg, design)

    system_summary, system_history = run_design(sys_cfg, design)
    outlet_sweep = water_outlet_orifice_sweep(sys_cfg, design)
    hx_boundary_history = averaged_hx_boundary_history(cfg, system_history)
    hx_boundary_history.update(
        hot_gas_schedule_from_helium(
            cfg,
            hx_boundary_history["time_s"],
            hx_boundary_history["helium_mdot_kg_s"],
        )
    )
    write_history_csv(output / "system_timeseries.csv", system_history)
    write_history_csv(output / "hx_boundary_schedule.csv", hx_boundary_history)
    (output / "water_outlet_orifice_sweep.json").write_text(
        json.dumps(_jsonable(outlet_sweep), indent=2),
        encoding="utf-8",
    )

    hx_runtime_s = None
    hx_summary: dict[str, Any] = {}
    hx_time_series = None
    if cfg.run_hx:
        inputs = build_hx_inputs_from_system_history(cfg, hx_boundary_history)
        start = time.perf_counter()
        solver, hx_summary = run_transient(inputs)
        hx_runtime_s = time.perf_counter() - start
        hx_time_series = get_transient_time_series(solver)
        _write_hx_npz(output / "hx_transient_timeseries.npz", hx_time_series)
        write_coupled_scalar_csv(output / "coupled_timeseries.csv", system_history, hx_time_series)

    diagnostics = coupled_diagnostics(system_history, hx_time_series, coolant_momentum_model=cfg.coolant_momentum_model)
    flags = validation_flags(cfg, system_history, hx_time_series, diagnostics)
    payload = {
        "module": "coupled_bangbang_hx",
        "coupling_level": "0D feed/tanks -> 1D transient HX, one-way boundary drive",
        "config": asdict(cfg),
        "pressurant_config": asdict(sys_cfg),
        "feed_design": asdict(design),
        "system_summary": asdict(system_summary),
        "water_outlet_orifice_sweep": outlet_sweep,
        "hx_runtime_s": hx_runtime_s,
        "hx_summary": hx_summary,
        "coupled_diagnostics": diagnostics,
        "validation_flags": flags,
        "outputs": {
            "system_timeseries_csv": "system_timeseries.csv",
            "hx_boundary_schedule_csv": "hx_boundary_schedule.csv",
            "water_outlet_orifice_sweep_json": "water_outlet_orifice_sweep.json",
            "coupled_timeseries_csv": "coupled_timeseries.csv" if hx_time_series is not None else None,
            "hx_timeseries_npz": "hx_transient_timeseries.npz" if hx_time_series is not None else None,
            "dashboard_html": "coupled_dashboard.html",
        },
    }
    (output / "summary.json").write_text(json.dumps(_jsonable(payload), indent=2), encoding="utf-8")
    write_dashboard(output / "coupled_dashboard.html", system_history, hx_time_series, payload)
    write_readme(output / "README.md", payload)
    return payload


def build_hx_inputs_from_system_history(
    config: CoupledBangBangHxConfig,
    system_history: dict[str, np.ndarray],
) -> dict[str, Any]:
    """Build `main_transient` inputs driven by the 0D system history."""

    inputs = build_inputs()
    t = np.asarray(system_history["time_s"], dtype=float)
    p_line = np.asarray(system_history["line_pressure_before_hx_bar"], dtype=float) * 1.0e5
    p_out = np.asarray(system_history["water_tank_pressure_bar"], dtype=float) * 1.0e5
    mdot_he = np.asarray(system_history["helium_mdot_kg_s"], dtype=float)
    T_he = np.asarray(system_history["supply_temperature_K"], dtype=float)

    transient = inputs["transient"]
    hotgas = inputs["hotgas"]

    hx_config = str(config.hx_config)
    if hx_config not in ("shellntube", "shellnHelicalTube"):
        raise ValueError("hx_config must be 'shellntube' or 'shellnHelicalTube'")
    inputs["combustor"].HX_config = hx_config
    inputs["combustor"].flow_config = config.flow_config
    if hx_config == "shellntube":
        inputs["run"].shelltube_transient_nodes = int(config.hx_nodes)
    else:
        _apply_helical_geometry(inputs, config)
        transient.n_axial = int(config.hx_nodes)
        transient.skip_steady_reference_probe = True

    transient.fluid_model = "transient_coolant"
    transient.coolant_momentum_model = config.coolant_momentum_model
    transient.chemistry_transient = config.chemistry
    transient.t_end = float(config.t_end_s)
    transient.max_step = float(config.hx_max_step_s)
    transient.n_save = int(config.hx_save_points)
    transient.insert_schedule_breakpoints = False
    transient.T_wall_initial = float(config.wall_initial_temperature_K)
    transient.T_coolant_initial = None
    if config.coolant_momentum_model == "low_mach":
        transient.transient_coolant_outlet_pressure = float(p_out[0])
    transient.progress_interval_time_s = max(1.0, float(config.t_end_s) / 10.0)

    if config.coolant_momentum_model == "low_mach":
        cap = max(float(config.low_mach_mdot_cap_kg_s), 1.0e-6)
        mdot_available = np.minimum(np.maximum(mdot_he, 0.0), cap)
        transient.schedule_mass_flow_c = _schedule_from_arrays(t, mdot_available)
        inputs["coolant"].mass_flow_c = float(np.nanmax(mdot_available))
    else:
        transient.schedule_mass_flow_c = _schedule_from_arrays(t, mdot_he)
        inputs["coolant"].mass_flow_c = float(np.nanmax(mdot_he))
    transient.schedule_p_c_in = _schedule_from_arrays(t, p_line)
    transient.schedule_p_c_out = _schedule_from_arrays(t, p_out)
    transient.schedule_T_c_in = _schedule_from_arrays(t, T_he)
    hot_schedule = hot_gas_schedule_from_helium(config, t, mdot_he)
    transient.schedule_mass_flow_lox = _schedule_from_arrays(t, hot_schedule["lox_mdot_kg_s"])
    transient.schedule_mass_flow_diesel = _schedule_from_arrays(t, hot_schedule["diesel_mdot_kg_s"])
    transient.schedule_mass_flow_g = _schedule_from_arrays(t, hot_schedule["hot_gas_mdot_kg_s"])
    transient.schedule_OF = (
        (0.0, config.lox_mdot_kg_s / max(config.diesel_mdot_kg_s, 1.0e-12)),
        (float(config.t_end_s), config.lox_mdot_kg_s / max(config.diesel_mdot_kg_s, 1.0e-12)),
    )
    transient.schedule_T_lox_in = ((0.0, config.lox_temperature_K), (float(config.t_end_s), config.lox_temperature_K))
    transient.schedule_ignition_state = ((0.0, 1.0), (float(config.t_end_s), 1.0))
    transient.ignition_time = 0.0

    inputs["coolant"].p_in = float(p_line[0])
    inputs["coolant"].T_in = float(T_he[0])
    hotgas.mass_flow_g = float(np.nanmax(hot_schedule["hot_gas_mdot_kg_s"]))
    hotgas.mixing_ratio = config.lox_mdot_kg_s / max(config.diesel_mdot_kg_s, 1.0e-12)
    hotgas.T_inj_LOX = config.lox_temperature_K
    return inputs


def _apply_helical_geometry(inputs: dict[str, Any], config: CoupledBangBangHxConfig) -> None:
    """Apply user-level helical geometry targets to the legacy input dataclasses."""

    combustor = inputs["combustor"]
    numerical = inputs["numerical"]
    tube_id = float(config.helical_inner_diameter_m)
    wall = float(config.helical_wall_thickness_m)
    tube_od = tube_id + 2.0 * wall
    centerline_diameter = float(config.helical_centerline_diameter_m)
    pitch = tube_od + float(config.helical_coil_gap_m)
    if tube_id <= 0.0 or wall <= 0.0 or centerline_diameter <= 0.0 or pitch <= 0.0:
        raise ValueError("helical geometry dimensions must be positive")

    combustor.Dh_coil = tube_id
    combustor.thickness_coil_wall = wall
    combustor.coil_gap = float(config.helical_coil_gap_m)
    combustor.material_HX = str(config.helical_material)
    combustor.Nusselt_correction = float(config.helical_shell_nusselt_correction)
    combustor.gap_shell2coil = 0.5 * (
        float(combustor.inner_diameter) - centerline_diameter - tube_od
    )
    if combustor.gap_shell2coil <= 0.0:
        raise ValueError(
            "helical_centerline_diameter_m and tube OD do not fit inside the combustor inner diameter"
        )

    radius = 0.5 * centerline_diameter
    h = pitch / (2.0 * np.pi)
    arc_factor = np.sqrt(1.0 + (radius / h) ** 2)
    axial_coil_length = float(config.helical_pipe_length_m) / arc_factor
    numerical.L_HX_max = (
        axial_coil_length
        + float(combustor.mixing_length)
        + 2.0 * float(combustor.length_2_coil)
        + tube_od
    )


def hot_gas_schedule_from_helium(
    config: CoupledBangBangHxConfig,
    t: np.ndarray,
    mdot_he: np.ndarray,
) -> dict[str, np.ndarray]:
    """Return OF-preserving hot-side schedules for the coupled validation case."""

    t_arr = np.asarray(t, dtype=float)
    mdot_he_arr = np.maximum(np.asarray(mdot_he, dtype=float), 0.0)
    nominal_hot = float(config.lox_mdot_kg_s + config.diesel_mdot_kg_s)
    of_ratio = float(config.lox_mdot_kg_s / max(config.diesel_mdot_kg_s, 1.0e-12))

    control = str(config.hot_gas_control).lower()
    if control == "constant":
        hot_total = np.full_like(t_arr, nominal_hot, dtype=float)
    elif control == "track_helium_mdot":
        nominal_he = max(float(config.hot_gas_nominal_helium_mdot_kg_s), 1.0e-12)
        hot_total = nominal_hot * mdot_he_arr / nominal_he
        on = mdot_he_arr > float(config.hot_gas_off_helium_mdot_kg_s)
        hot_total = np.where(on, hot_total, 0.0)
        if config.hot_gas_min_mdot_kg_s > 0.0:
            hot_total = np.where(on, np.maximum(hot_total, float(config.hot_gas_min_mdot_kg_s)), 0.0)
        hot_total = np.minimum(hot_total, float(config.hot_gas_max_mdot_kg_s))
    else:
        raise ValueError("hot_gas_control must be 'track_helium_mdot' or 'constant'")

    diesel = hot_total / (1.0 + of_ratio)
    lox = hot_total - diesel
    return {
        "hot_gas_mdot_kg_s": hot_total,
        "lox_mdot_kg_s": lox,
        "diesel_mdot_kg_s": diesel,
    }


def water_outlet_orifice_sweep(
    system_config: PressurantSystemConfig,
    feed_design: FeedDesign,
    *,
    diameters_mm: tuple[float, ...] = (20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0),
) -> list[dict[str, float]]:
    """Evaluate fixed water outlet diameters around the 30 L/s target."""

    rows = []
    for diameter_mm in diameters_mm:
        cfg = replace(
            system_config,
            water_exit_orifice_diameter_m=float(diameter_mm) * 1.0e-3,
        )
        summary, _history = run_design(cfg, feed_design)
        rows.append(
            {
                "water_exit_orifice_diameter_mm": float(diameter_mm),
                "mean_water_flow_L_s": summary.mean_water_flow_L_s,
                "mean_water_tank_pressure_bar": summary.mean_pressure_bar,
                "mean_line_pressure_bar": summary.mean_line_pressure_bar,
                "helium_used_kg": summary.helium_used_kg,
                "final_supply_pressure_bar": summary.final_supply_pressure_bar,
                "score": summary.score,
            }
        )
    return rows


def coupled_diagnostics(
    system_history: dict[str, np.ndarray],
    hx_time_series: dict[str, Any] | None,
    *,
    coolant_momentum_model: str = "quasi_steady",
) -> dict[str, float]:
    """Return scalar checks spanning the feed system and HX response."""

    out = {
        "system_mean_line_pressure_bar": _finite_mean(system_history["line_pressure_before_hx_bar"]),
        "system_min_line_pressure_bar": _finite_min(system_history["line_pressure_before_hx_bar"]),
        "system_max_line_pressure_bar": _finite_max(system_history["line_pressure_before_hx_bar"]),
        "system_min_line_minus_tank_pressure_bar": _finite_min(
            np.asarray(system_history["line_pressure_before_hx_bar"], dtype=float)
            - np.asarray(system_history["water_tank_pressure_bar"], dtype=float)
        ),
        "system_mean_helium_mdot_kg_s": _finite_mean(system_history["helium_mdot_kg_s"]),
        "system_max_helium_mdot_kg_s": _finite_max(system_history["helium_mdot_kg_s"]),
        "system_min_supply_temperature_K": _finite_min(system_history["supply_temperature_K"]),
    }
    if hx_time_series is None:
        return out

    scalars = hx_time_series.get("scalars", {})
    fields = hx_time_series.get("fields", {})
    t_hx = np.asarray(hx_time_series.get("t", []), dtype=float)
    x_hx = np.asarray(hx_time_series.get("x", []), dtype=float)
    mdot_out = np.asarray(scalars.get("mdot_c_outlet_face", []), dtype=float)
    flow_mask = np.abs(mdot_out) >= 0.01
    T_c_out = np.asarray(scalars.get("T_c_out", []), dtype=float)
    peak_hot_face = _field_peak(fields.get("T_wg"), t_hx, x_hx)
    peak_mean_wall = _field_peak(fields.get("Tbar"), t_hx, x_hx)
    peak_cold_face = _field_peak(fields.get("T_wc"), t_hx, x_hx)
    out.update(
        {
            "hx_final_T_c_out_K": _finite_last(scalars.get("T_c_out")),
            "hx_max_flowing_T_c_out_K": _masked_finite_max(T_c_out, flow_mask),
            "hx_mean_flowing_T_c_out_K": _masked_finite_mean(T_c_out, flow_mask),
            "hx_final_mdot_c_outlet_face_kg_s": _finite_last(scalars.get("mdot_c_outlet_face")),
            "hx_peak_T_wall_K": _finite_max(scalars.get("T_wall_max")),
            "hx_peak_T_wall_hot_face_K": peak_hot_face["value"],
            "hx_peak_T_wall_hot_face_time_s": peak_hot_face["time_s"],
            "hx_peak_T_wall_hot_face_x_m": peak_hot_face["x_m"],
            "hx_peak_T_wall_mean_K": peak_mean_wall["value"],
            "hx_peak_T_wall_mean_time_s": peak_mean_wall["time_s"],
            "hx_peak_T_wall_mean_x_m": peak_mean_wall["x_m"],
            "hx_peak_T_wall_cold_face_K": peak_cold_face["value"],
            "hx_peak_T_wall_cold_face_time_s": peak_cold_face["time_s"],
            "hx_peak_T_wall_cold_face_x_m": peak_cold_face["x_m"],
            "hx_peak_Q_hot_kW": _finite_max(scalars.get("Q_hot_kW")),
            "hx_max_hot_gas_pressure_drop_bar": _finite_max(scalars.get("dp_g_total_Pa")) / 1.0e5,
            "hx_max_shell_pressure_drop_estimate_bar": _finite_max(scalars.get("dp_shell_total_Pa")) / 1.0e5,
            "hx_mean_inlet_face_mdot_kg_s": _finite_mean(scalars.get("mdot_c_inlet_face")),
            "hx_mean_outlet_face_mdot_kg_s": _finite_mean(scalars.get("mdot_c_outlet_face")),
            "hx_max_inlet_face_mdot_kg_s": _finite_max(scalars.get("mdot_c_inlet_face")),
            "hx_max_outlet_face_mdot_kg_s": _finite_max(scalars.get("mdot_c_outlet_face")),
        }
    )
    p_c_min = _finite_min(fields.get("p_c")) / 1.0e5
    p_c_max = _finite_max(fields.get("p_c")) / 1.0e5
    p_c_span_stats = _field_span_stats(fields.get("p_c"))
    p_c_flowing_span_stats = _field_span_stats(fields.get("p_c"), flow_mask)
    if str(coolant_momentum_model).lower() == "low_mach":
        out.update(
            {
                "hx_min_coolant_pressure_bar": p_c_min,
                "hx_max_coolant_pressure_bar": p_c_max,
                "hx_max_coolant_pressure_span_bar": p_c_span_stats["max"] / 1.0e5,
                "hx_mean_coolant_pressure_span_bar": p_c_span_stats["mean"] / 1.0e5,
                "hx_p95_coolant_pressure_span_bar": p_c_span_stats["p95"] / 1.0e5,
                "hx_max_flowing_coolant_pressure_span_bar": p_c_flowing_span_stats["max"] / 1.0e5,
                "hx_mean_flowing_coolant_pressure_span_bar": p_c_flowing_span_stats["mean"] / 1.0e5,
                "hx_p95_flowing_coolant_pressure_span_bar": p_c_flowing_span_stats["p95"] / 1.0e5,
            }
        )
    else:
        out.update(
            {
                "hx_min_coolant_thermodynamic_pressure_bar": p_c_min,
                "hx_max_coolant_thermodynamic_pressure_bar": p_c_max,
                "hx_max_coolant_thermodynamic_pressure_span_bar": p_c_span_stats["max"] / 1.0e5,
            }
        )
    return out


def validation_flags(
    config: CoupledBangBangHxConfig,
    system_history: dict[str, np.ndarray],
    hx_time_series: dict[str, Any] | None,
    diagnostics: dict[str, float],
) -> list[str]:
    """Return human-readable sanity-check flags for the coupled run."""

    flags: list[str] = []
    water_flow = _finite_mean(system_history.get("water_flow_L_s"))
    tank_pressure = _finite_mean(system_history.get("water_tank_pressure_bar"))
    helium_mdot = _finite_mean(system_history.get("helium_mdot_kg_s"))
    if hx_time_series is not None and config.coolant_momentum_model == "low_mach":
        helium_mdot = diagnostics.get("hx_mean_outlet_face_mdot_kg_s", helium_mdot)
    if np.isfinite(water_flow) and abs(water_flow - config.target_water_flow_L_s) > 1.0:
        flags.append(
            "Mean water outflow is outside the configured target band "
            f"({water_flow:.2f} L/s vs {config.target_water_flow_L_s:.2f} L/s)."
        )
    if np.isfinite(tank_pressure) and abs(tank_pressure - config.target_water_tank_pressure_bar) > 2.0:
        flags.append(
            "Mean water-tank pressure is outside the configured target band "
            f"({tank_pressure:.2f} bar vs {config.target_water_tank_pressure_bar:.2f} bar)."
        )
    if (
        np.isfinite(helium_mdot)
        and not (config.min_target_helium_mdot_kg_s <= helium_mdot <= config.max_target_helium_mdot_kg_s)
    ):
        flags.append(
            "Mean helium flow is outside the expected target range "
            f"({helium_mdot:.4f} kg/s, target "
            f"{config.min_target_helium_mdot_kg_s:.3f}-"
            f"{config.max_target_helium_mdot_kg_s:.3f} kg/s)."
        )

    line_margin = diagnostics.get("system_min_line_minus_tank_pressure_bar", float("nan"))
    if np.isfinite(line_margin) and line_margin < -0.5:
        flags.append(
            "System line pressure falls below water-tank/backpressure; reverse or weak HX flow is possible."
        )

    if hx_time_series is None:
        flags.append("HX transient was not run; only the 0D system schedule was checked.")
        return flags

    if diagnostics.get("hx_peak_T_wall_K", 0.0) > config.max_wall_temperature_K:
        flags.append(
            "HX wall temperature exceeds the configured material limit "
            f"({config.max_wall_temperature_K:.1f} K)."
        )
    flowing_T = diagnostics.get("hx_max_flowing_T_c_out_K", float("nan"))
    if np.isfinite(flowing_T) and flowing_T > config.max_target_helium_outlet_temperature_K:
        flags.append(
            "HX flowing helium outlet temperature exceeds the configured target "
            f"({flowing_T:.1f} K vs {config.max_target_helium_outlet_temperature_K:.1f} K)."
        )
    scalars = hx_time_series.get("scalars", {})
    if _scalar_hits_limit(scalars.get("T_c_min"), 60.0) or _scalar_hits_limit(scalars.get("T_c_max"), 2500.0):
        flags.append(
            "HX coolant temperature reached a property-validity limiter; treat this run as nonphysical."
        )

    p_drop = diagnostics.get(
        "hx_p95_flowing_coolant_pressure_span_bar",
        diagnostics.get(
            "hx_max_coolant_pressure_span_bar",
            diagnostics.get("hx_max_shell_pressure_drop_estimate_bar", float("nan")),
        ),
    )
    if np.isfinite(p_drop) and p_drop > config.max_hx_pressure_loss_bar:
        pressure_label = (
            "p95 flowing coolant pressure span"
            if config.coolant_momentum_model == "low_mach" else
            "shell-side pressure-drop estimate"
        )
        flags.append(
            f"HX {pressure_label} exceeds the configured target "
            f"({config.max_hx_pressure_loss_bar:.1f} bar)."
        )

    if config.coolant_momentum_model == "low_mach":
        p_hx_max = diagnostics.get("hx_max_coolant_pressure_bar", float("nan"))
        p_system_max = diagnostics.get("system_max_line_pressure_bar", float("nan"))
        if np.isfinite(p_hx_max) and np.isfinite(p_system_max) and p_hx_max > 1.25 * p_system_max:
            flags.append(
                "HX coolant pressure exceeds the driving system pressure by more than 25%; "
                "treat this run as a numerical/boundary-condition sanity failure."
            )
        hx_mean = diagnostics.get("hx_mean_outlet_face_mdot_kg_s", float("nan"))
        system_mean = diagnostics.get("system_mean_helium_mdot_kg_s", float("nan"))
        if np.isfinite(hx_mean) and np.isfinite(system_mean) and abs(hx_mean - system_mean) > 0.03:
            flags.append(
                "Pressure-driven HX helium flow differs from the one-way 0D system "
                f"prediction by more than 30 g/s ({hx_mean:.4f} vs {system_mean:.4f} kg/s); "
                "iterate the feed/HX loss surrogate before a long run."
            )

    mdot_in = diagnostics.get("hx_mean_inlet_face_mdot_kg_s", float("nan"))
    mdot_system = diagnostics.get("system_mean_helium_mdot_kg_s", float("nan"))
    if np.isfinite(mdot_in) and np.isfinite(mdot_system) and mdot_system > 1e-9 and mdot_in < -0.05 * mdot_system:
        flags.append(
            "Mean HX inlet face flow is reversed relative to the positive system helium command."
        )
    if not flags:
        flags.append("No coupled sanity-check flags triggered.")
    return flags


def write_coupled_scalar_csv(
    path: Path,
    system_history: dict[str, np.ndarray],
    hx_time_series: dict[str, Any],
) -> None:
    """Write system and HX scalar histories interpolated on the HX save times."""

    t_hx = np.asarray(hx_time_series["t"], dtype=float)
    rows = {
        "time_s": t_hx,
        "system_supply_pressure_bar": _interp_system(system_history, "supply_pressure_bar", t_hx),
        "system_supply_temperature_K": _interp_system(system_history, "supply_temperature_K", t_hx),
        "system_line_pressure_before_hx_bar": _interp_system(system_history, "line_pressure_before_hx_bar", t_hx),
        "system_water_tank_pressure_bar": _interp_system(system_history, "water_tank_pressure_bar", t_hx),
        "system_helium_mdot_kg_s": _interp_system(system_history, "helium_mdot_kg_s", t_hx),
        "system_water_flow_L_s": _interp_system(system_history, "water_flow_L_s", t_hx),
        "system_open_branches": _interp_system(system_history, "open_branches", t_hx),
    }
    for name, values in hx_time_series.get("scalars", {}).items():
        arr = np.asarray(values, dtype=float)
        if arr.ndim == 1 and len(arr) == len(t_hx):
            rows[f"hx_{name}"] = arr
    _write_csv(path, rows)


def write_dashboard(
    path: Path,
    system_history: dict[str, np.ndarray],
    hx_time_series: dict[str, Any] | None,
    payload: dict[str, Any],
) -> None:
    """Write a standalone HTML dashboard for the coupled validation run."""

    data = {
        "system": {k: _jsonable(v) for k, v in system_history.items()},
        "hx_boundary": None,
        "hx": _jsonable(hx_time_series) if hx_time_series is not None else None,
        "summary": _jsonable(payload),
    }
    html = _COUPLED_DASHBOARD_TEMPLATE.replace("/*__DATA__*/", json.dumps(data))
    path.write_text(html, encoding="utf-8")


def write_readme(path: Path, payload: dict[str, Any]) -> None:
    cfg = payload["config"]
    diag = payload["coupled_diagnostics"]
    lines = [
        "# Coupled Bang-Bang HX Validation",
        "",
        "This folder is generated by:",
        "",
        "```powershell",
        f"python -m hps_combustor.validation.coupled_bangbang_hx --settings {DEFAULT_SETTINGS_FILE}",
        "```",
        "",
        "The module is intentionally separate from `main_transient.py`. The external",
        "system is a 0D pressurant/feed/tank surrogate; the HX remains the detailed",
        f"1D transient coolant/material {cfg.get('hx_config', 'shellntube')} solver.",
        "",
        "## Current Coupling",
        "",
        "- 0D system outputs helium `mdot`, pre-HX pressure, supply temperature, and",
        "  water-tank/backpressure histories.",
        "- Those histories are passed to the HX as inlet mass flow, inlet pressure,",
        "  inlet temperature, and outlet pressure schedules.",
        "- By default, raw valve ripple is averaged before it is passed to the HX.",
        f"  Current averaging window: {cfg['hx_boundary_average_window_s']} s.",
        "- The coupling is one-way in this module. It is for sanity checking and",
        "  Simulink integration preparation, not a replacement for a closed-loop",
        "  multiphysics solver.",
        "- The exact settings used for this run are saved as `settings_used.json`.",
        "- `water_outlet_orifice_sweep.json` checks fixed water outlet diameters",
        "  around the 30 L/s target.",
        "- Open `coupled_dashboard.html` in a browser to view the dedicated",
        "  bang-bang validation dashboard.",
        "",
        "## Case",
        "",
        f"- Simulated duration: {cfg['t_end_s']} s",
        f"- HX nodes: {cfg['hx_nodes']}",
        f"- HX config: {cfg.get('hx_config', 'shellntube')}",
        f"- HX max step: {cfg['hx_max_step_s']} s",
        f"- Chemistry: {cfg['chemistry']}",
        f"- Flow configuration: {cfg['flow_config']}",
        f"- Hot-gas control: {cfg['hot_gas_control']}",
        f"- Low-Mach helium mass-flow cap: {cfg.get('low_mach_mdot_cap_kg_s', float('nan')):.6g} kg/s",
        f"- Nominal hot-gas flow: {cfg['lox_mdot_kg_s'] + cfg['diesel_mdot_kg_s']:.6g} kg/s",
        f"- Maximum hot-gas flow: {cfg['hot_gas_max_mdot_kg_s']:.6g} kg/s",
        f"- Wall temperature limit: {cfg['max_wall_temperature_K']:.2f} K "
        f"({cfg['max_wall_temperature_K'] - 273.15:.2f} degC)",
        f"- Helium-side hydraulic pressure-drop target: <= {cfg['max_hx_pressure_loss_bar']} bar",
    ]
    if cfg.get("hx_config") == "shellnHelicalTube":
        lines.extend(
            [
                f"- Helical tube ID: {cfg['helical_inner_diameter_m'] * 1e3:.3g} mm",
                f"- Helical wall thickness: {cfg['helical_wall_thickness_m'] * 1e3:.3g} mm",
                f"- Helical target pipe length: {cfg['helical_pipe_length_m']:.3g} m",
                f"- Helical centerline diameter: {cfg['helical_centerline_diameter_m'] * 1e3:.3g} mm",
                f"- Helical material: {cfg['helical_material']}",
                f"- Helical shell Nu correction: {cfg['helical_shell_nusselt_correction']:.3g}",
            ]
        )
    lines.extend(
        [
            "",
            "## Pressure-Drop Interpretation",
            "",
            "In `quasi_steady` coolant momentum mode, the hydraulic HX pressure loss is",
            "`hx_max_shell_pressure_drop_estimate_bar`, computed from the shell-side",
            "correlation. The saved `p_c` field is a reconstructed thermodynamic coolant",
            "state used for mass/energy bookkeeping; its axial span is not a hydraulic",
            "pressure-drop prediction in this mode. In `low_mach` mode, `p_c` is the",
            "resolved pressure field and the extrema become the momentum-model diagnostic.",
            "For `low_mach`, `schedule_mass_flow_c` is only a numerical cap; the",
            "physical helium mass flow is the solved HX face flow.",
            "",
            "## Diagnostics",
            "",
        ]
    )
    for key, value in diag.items():
        lines.append(f"- {key}: {value:.6g}")
    lines.extend(["", "## Validation Flags", ""])
    for flag in payload.get("validation_flags", []):
        lines.append(f"- {flag}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _schedule_from_arrays(t: np.ndarray, values: np.ndarray) -> tuple[tuple[float, float], ...]:
    return tuple((float(tt), float(vv)) for tt, vv in zip(np.asarray(t, dtype=float), np.asarray(values, dtype=float)))


def averaged_hx_boundary_history(
    config: CoupledBangBangHxConfig,
    system_history: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Return HX boundary histories, optionally averaged over valve cycles."""

    window = float(config.hx_boundary_average_window_s)
    if window <= 0.0:
        return {key: np.asarray(value, dtype=float).copy() for key, value in system_history.items()}

    t = np.asarray(system_history["time_s"], dtype=float)
    dt = float(np.median(np.diff(t))) if len(t) > 1 else window
    n = max(int(round(window / max(dt, 1.0e-12))), 1)
    out = {key: np.asarray(value, dtype=float).copy() for key, value in system_history.items()}
    for key in (
        "water_tank_pressure_bar",
        "line_pressure_before_hx_bar",
        "supply_pressure_bar",
        "supply_temperature_K",
        "water_flow_L_s",
        "helium_mdot_kg_s",
        "open_branches",
    ):
        out[key] = _centered_moving_average(out[key], n)
    return out


def _centered_moving_average(values: np.ndarray, n: int) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if n <= 1 or arr.size <= 1:
        return arr.copy()
    left = n // 2
    right = n - 1 - left
    padded = np.pad(arr, (left, right), mode="edge")
    kernel = np.ones(n, dtype=float) / float(n)
    return np.convolve(padded, kernel, mode="valid")


def _interp_system(history: dict[str, np.ndarray], key: str, t_target: np.ndarray) -> np.ndarray:
    t = np.asarray(history["time_s"], dtype=float)
    values = np.asarray(history[key], dtype=float)
    return np.interp(np.asarray(t_target, dtype=float), t, values)


def _write_hx_npz(path: Path, time_series: dict[str, Any]) -> None:
    arrays = {
        "t": np.asarray(time_series["t"], dtype=float),
        "x": np.asarray(time_series["x"], dtype=float),
    }
    arrays.update({f"field_{k}": np.asarray(v, dtype=float) for k, v in time_series.get("fields", {}).items()})
    arrays.update({f"scalar_{k}": np.asarray(v, dtype=float) for k, v in time_series.get("scalars", {}).items()})
    np.savez_compressed(path, **arrays)


def _write_csv(path: Path, columns: dict[str, np.ndarray]) -> None:
    names = list(columns)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for row in zip(*(columns[name] for name in names)):
            writer.writerow([f"{float(value):.10g}" for value in row])


def _finite_array(values) -> np.ndarray:
    if values is None:
        return np.array([], dtype=float)
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def _finite_mean(values) -> float:
    arr = _finite_array(values)
    return float(np.mean(arr)) if arr.size else float("nan")


def _finite_min(values) -> float:
    arr = _finite_array(values)
    return float(np.min(arr)) if arr.size else float("nan")


def _finite_max(values) -> float:
    arr = _finite_array(values)
    return float(np.max(arr)) if arr.size else float("nan")


def _finite_last(values) -> float:
    arr = _finite_array(values)
    return float(arr[-1]) if arr.size else float("nan")


def _masked_finite_max(values, mask) -> float:
    arr = np.asarray(values, dtype=float)
    mask_arr = np.asarray(mask, dtype=bool)
    if arr.shape != mask_arr.shape:
        return float("nan")
    selected = arr[mask_arr & np.isfinite(arr)]
    return float(np.max(selected)) if selected.size else float("nan")


def _masked_finite_mean(values, mask) -> float:
    arr = np.asarray(values, dtype=float)
    mask_arr = np.asarray(mask, dtype=bool)
    if arr.shape != mask_arr.shape:
        return float("nan")
    selected = arr[mask_arr & np.isfinite(arr)]
    return float(np.mean(selected)) if selected.size else float("nan")


def _field_peak(values, time_s: np.ndarray, x_m: np.ndarray) -> dict[str, float]:
    if values is None:
        return {"value": float("nan"), "time_s": float("nan"), "x_m": float("nan")}
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2 or arr.size == 0:
        return {"value": float("nan"), "time_s": float("nan"), "x_m": float("nan")}
    finite = np.where(np.isfinite(arr), arr, -np.inf)
    if not np.any(np.isfinite(arr)):
        return {"value": float("nan"), "time_s": float("nan"), "x_m": float("nan")}
    i, j = np.unravel_index(int(np.argmax(finite)), arr.shape)
    t_val = float(time_s[i]) if i < len(time_s) else float("nan")
    x_val = float(x_m[j]) if j < len(x_m) else float("nan")
    return {"value": float(arr[i, j]), "time_s": t_val, "x_m": x_val}


def _max_field_span(values) -> float:
    if values is None:
        return float("nan")
    arr = np.asarray(values, dtype=float)
    if arr.ndim < 2:
        return float("nan")
    spans = np.nanmax(arr, axis=1) - np.nanmin(arr, axis=1)
    spans = spans[np.isfinite(spans)]
    return float(np.max(spans)) if spans.size else float("nan")


def _field_span_stats(values, mask=None) -> dict[str, float]:
    if values is None:
        return {"max": float("nan"), "mean": float("nan"), "p95": float("nan")}
    arr = np.asarray(values, dtype=float)
    if arr.ndim < 2:
        return {"max": float("nan"), "mean": float("nan"), "p95": float("nan")}
    spans = []
    for row in arr:
        finite = row[np.isfinite(row)]
        if finite.size:
            spans.append(float(np.max(finite) - np.min(finite)))
        else:
            spans.append(float("nan"))
    span_arr = np.asarray(spans, dtype=float)
    if mask is not None:
        mask_arr = np.asarray(mask, dtype=bool)
        if mask_arr.shape == span_arr.shape:
            span_arr = span_arr[mask_arr]
    span_arr = span_arr[np.isfinite(span_arr)]
    if not span_arr.size:
        return {"max": float("nan"), "mean": float("nan"), "p95": float("nan")}
    return {
        "max": float(np.max(span_arr)),
        "mean": float(np.mean(span_arr)),
        "p95": float(np.percentile(span_arr, 95.0)),
    }


def _scalar_hits_limit(values, limit: float, *, atol: float = 1.0e-9) -> bool:
    arr = _finite_array(values)
    if arr.size == 0:
        return False
    return bool(np.any(np.isclose(arr, float(limit), atol=atol, rtol=0.0)))


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, float):
        return value if np.isfinite(value) else str(value)
    return value


def _dataclass_from_mapping(cls, values: dict[str, Any]):
    return cls(**_select_dataclass_fields(cls, values))


def _select_dataclass_fields(cls, values: dict[str, Any]) -> dict[str, Any]:
    allowed = {field.name for field in dataclass_fields(cls)}
    return {key: value for key, value in values.items() if key in allowed}


def main(argv: list[str] | None = None) -> dict[str, Any]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", default=DEFAULT_SETTINGS_FILE, help="JSON settings file for this coupled validation.")
    parser.add_argument("--write-default-settings", action="store_true", help="Write the default settings JSON and exit.")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--t-end", type=float, default=None)
    parser.add_argument("--system-dt", type=float, default=None)
    parser.add_argument("--hx-max-step", type=float, default=None)
    parser.add_argument("--hx-nodes", type=int, default=None)
    parser.add_argument("--hx-save-points", type=int, default=None)
    parser.add_argument("--hx-config", choices=("shellntube", "shellnHelicalTube"), default=None)
    parser.add_argument("--flow-config", choices=("co", "counter"), default=None)
    parser.add_argument("--chemistry", choices=("finite_rate", "equilibrium", "frozen"), default=None)
    parser.add_argument("--coolant-momentum-model", choices=("quasi_steady", "low_mach"), default=None)
    parser.add_argument("--system-only", action="store_true", help="Write the 0D system dashboard without running the HX.")
    args = parser.parse_args(argv)

    if args.write_default_settings:
        write_settings(
            args.settings,
            CoupledBangBangHxConfig(),
            default_pressurant_system_config(),
            default_feed_design(),
        )
        print(f"Wrote settings: {args.settings}")
        return {}

    if Path(args.settings).exists():
        cfg, system_config, feed_design = load_settings(args.settings)
    else:
        cfg = CoupledBangBangHxConfig()
        system_config = default_pressurant_system_config()
        feed_design = default_feed_design()

    cfg = _apply_cli_overrides(cfg, args)
    if args.system_only:
        cfg = replace(cfg, run_hx=False)

    payload = run_coupled_case(cfg, system_config=system_config, feed_design=feed_design)
    diag = payload["coupled_diagnostics"]
    output = Path(cfg.output_dir)
    print(f"Wrote {output}")
    print(f"settings:  {output / 'settings_used.json'}")
    print(f"dashboard: {output / 'coupled_dashboard.html'}")
    print(
        "coupled bang-bang HX: "
        f"line_p_mean={diag['system_mean_line_pressure_bar']:.2f} bar, "
        f"He_mdot_mean={diag['system_mean_helium_mdot_kg_s']:.4f} kg/s"
    )
    if not args.system_only:
        print(
            f"HX runtime={payload['hx_runtime_s']:.2f}s, "
            f"T_c_out_final={diag['hx_final_T_c_out_K']:.1f} K, "
            f"T_wall_peak={diag['hx_peak_T_wall_K']:.1f} K"
        )
    return payload


def _apply_cli_overrides(config: CoupledBangBangHxConfig, args) -> CoupledBangBangHxConfig:
    updates = {}
    mapping = {
        "output_dir": args.output_dir,
        "t_end_s": args.t_end,
        "system_dt_s": args.system_dt,
        "hx_max_step_s": args.hx_max_step,
        "hx_nodes": args.hx_nodes,
        "hx_save_points": args.hx_save_points,
        "hx_config": args.hx_config,
        "flow_config": args.flow_config,
        "chemistry": args.chemistry,
        "coolant_momentum_model": args.coolant_momentum_model,
    }
    for key, value in mapping.items():
        if value is not None:
            updates[key] = value
    return replace(config, **updates) if updates else config


_COUPLED_DASHBOARD_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Combustor-HX Coupled Bang-Bang Dashboard</title>
<style>
  :root{--bg:#f8f8f6;--panel:#ffffff;--ink:#111;--muted:#666;--grid:#ddd;--blue:#276fd1;--red:#d74738;--green:#14996b;--orange:#d97828;--violet:#7658c9}
  body{margin:0;background:var(--bg);color:var(--ink);font-family:system-ui,-apple-system,Segoe UI,sans-serif;font-size:14px}
  header{padding:16px 22px;border-bottom:1px solid #ddd}
  h1{font-size:18px;margin:0 0 4px}
  .sub{color:var(--muted)}
  .tiles{display:flex;gap:10px;flex-wrap:wrap;padding:14px 22px}
  .tile{background:var(--panel);border:1px solid #ddd;border-radius:8px;padding:10px 12px;min-width:150px}
  .lbl{font-size:11px;color:var(--muted);text-transform:uppercase}
  .val{font-size:19px;font-weight:650;font-variant-numeric:tabular-nums}
  .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:14px;padding:0 22px 24px}
  .panel{background:var(--panel);border:1px solid #ddd;border-radius:8px;padding:12px}
  .flags{margin:0 22px 14px;background:var(--panel);border:1px solid #ddd;border-radius:8px;padding:12px}
  .flags ul{margin:6px 0 0;padding-left:20px}
  h2{font-size:13px;margin:0 0 6px}
  canvas{width:100%;height:230px;display:block}
  .note{color:var(--muted);font-size:12px;margin-top:6px}
</style>
</head>
<body>
<header>
  <h1>Combustor-HX Coupled Bang-Bang Dashboard</h1>
  <div class="sub" id="subtitle"></div>
</header>
<div class="tiles" id="tiles"></div>
<div class="flags">
  <h2>Validation flags</h2>
  <ul id="flags"></ul>
</div>
<div class="grid">
  <div class="panel"><h2>System pressures</h2><canvas id="pSystem"></canvas><div class="note">Supply, pre-HX line, and water-tank pressure.</div></div>
  <div class="panel"><h2>System flows and valves</h2><canvas id="flowSystem"></canvas><div class="note">Helium flow, water drain flow, and staged open branches.</div></div>
  <div class="panel"><h2>HX outlet and wall temperatures</h2><canvas id="tempHx"></canvas><div class="note">Available after the detailed HX run.</div></div>
  <div class="panel"><h2>HX heat duty and coolant flow</h2><canvas id="powerHx"></canvas><div class="note">HX heat duty with inlet/outlet shell-side face flows.</div></div>
  <div class="panel"><h2>HX pressure diagnostics</h2><canvas id="pHx"></canvas><div class="note">Quasi-steady mode shows shell-side pressure-drop estimate; low-Mach mode shows reconstructed pressure extrema.</div></div>
  <div class="panel"><h2>Axial temperatures at final HX snapshot</h2><canvas id="profileHx"></canvas><div class="note">Final gas, coolant, and wall profiles.</div></div>
</div>
<script>
const DATA = /*__DATA__*/;
const S = DATA.system;
const HX = DATA.hx;
function arr(x){return x || []}
function finite(a){return arr(a).filter(Number.isFinite)}
function last(a){a=finite(a); return a.length?a[a.length-1]:NaN}
function min(a){a=finite(a); return a.length?Math.min(...a):NaN}
function max(a){a=finite(a); return a.length?Math.max(...a):NaN}
function mean(a){a=finite(a); return a.length?a.reduce((x,y)=>x+y,0)/a.length:NaN}
function fmt(v,u){return Number.isFinite(v)?v.toFixed(Math.abs(v)>100?0:2)+' '+u:'n/a'}
document.getElementById('subtitle').textContent =
  DATA.summary.coupling_level + ' | t_end=' + DATA.summary.config.t_end_s + ' s';
const tiles = [
  ['Mean line pressure', fmt(mean(S.line_pressure_before_hx_bar),'bar')],
  ['Mean He flow', fmt(mean(S.helium_mdot_kg_s),'kg/s')],
  ['Final supply', fmt(last(S.supply_pressure_bar),'bar')],
  ['Final supply T', fmt(last(S.supply_temperature_K),'K')],
  ['HX He outlet', HX ? fmt(last(HX.scalars.T_c_out),'K') : 'not run'],
  ['HX wall peak', HX ? fmt(max(HX.scalars.T_wall_max),'K') : 'not run'],
];
document.getElementById('tiles').innerHTML = tiles.map(t=>`<div class="tile"><div class="lbl">${t[0]}</div><div class="val">${t[1]}</div></div>`).join('');
document.getElementById('flags').innerHTML =
  (DATA.summary.validation_flags || []).map(v=>`<li>${v}</li>`).join('');
function color(name){return getComputedStyle(document.documentElement).getPropertyValue(name).trim()}
function setup(id){const c=document.getElementById(id); const r=Math.min(devicePixelRatio||1,2); const w=c.clientWidth||420,h=230; c.width=w*r;c.height=h*r;c.style.height=h+'px'; const ctx=c.getContext('2d'); ctx.setTransform(r,0,0,r,0,0); return {ctx,w,h}}
function extent(series){let lo=Infinity,hi=-Infinity; for(const y of series){for(const v of finite(y)){lo=Math.min(lo,v);hi=Math.max(hi,v)}} if(!Number.isFinite(lo)){lo=0;hi=1} if(lo===hi){lo-=1;hi+=1} const p=(hi-lo)*0.08; return [lo-p,hi+p]}
function plot(id, x, series, ylab){
  const {ctx,w,h}=setup(id); ctx.clearRect(0,0,w,h); const pad={l:54,r:12,t:12,b:28};
  const xd=[x[0]||0,x[x.length-1]||1], yd=extent(series.map(s=>s.y));
  const sx=v=>pad.l+(v-xd[0])/(xd[1]-xd[0]||1)*(w-pad.l-pad.r);
  const sy=v=>h-pad.b-(v-yd[0])/(yd[1]-yd[0]||1)*(h-pad.t-pad.b);
  ctx.strokeStyle=color('--grid'); ctx.fillStyle=color('--muted'); ctx.font='11px system-ui';
  for(let i=0;i<=4;i++){const yy=pad.t+i*(h-pad.t-pad.b)/4; ctx.beginPath();ctx.moveTo(pad.l,yy);ctx.lineTo(w-pad.r,yy);ctx.stroke(); const val=yd[1]-(yd[1]-yd[0])*i/4; ctx.textAlign='right';ctx.fillText(val.toFixed(Math.abs(val)>100?0:2),pad.l-5,yy+4)}
  ctx.textAlign='left'; ctx.fillText(ylab,8,12); ctx.textAlign='right'; ctx.fillText('t [s]',w-pad.r,h-8);
  for(const s of series){ctx.strokeStyle=s.c;ctx.lineWidth=2;ctx.beginPath();let started=false; for(let i=0;i<x.length;i++){const v=s.y[i]; if(!Number.isFinite(v))continue; const px=sx(x[i]),py=sy(v); if(!started){ctx.moveTo(px,py);started=true}else ctx.lineTo(px,py)} ctx.stroke(); ctx.fillStyle=s.c; ctx.fillText(s.name,w-pad.r,pad.t+14*series.indexOf(s));}
}
const tS = S.time_s;
plot('pSystem', tS, [
  {name:'supply', y:S.supply_pressure_bar, c:color('--blue')},
  {name:'line', y:S.line_pressure_before_hx_bar, c:color('--orange')},
  {name:'water tank', y:S.water_tank_pressure_bar, c:color('--green')},
], 'pressure [bar]');
plot('flowSystem', tS, [
  {name:'He kg/s', y:S.helium_mdot_kg_s, c:color('--blue')},
  {name:'water L/s /100', y:S.water_flow_L_s.map(v=>v/100), c:color('--green')},
  {name:'branches /10', y:S.open_branches.map(v=>v/10), c:color('--red')},
], 'mixed scale');
if(HX){
  const t=HX.t, sc=HX.scalars, f=HX.fields;
  plot('tempHx', t, [
    {name:'T He out', y:sc.T_c_out, c:color('--blue')},
    {name:'T gas out', y:sc.T_g_out, c:color('--red')},
    {name:'T wall max', y:sc.T_wall_max, c:color('--orange')},
    {name:'T wall min', y:sc.T_wall_min, c:color('--green')},
  ], 'T [K]');
  plot('powerHx', t, [
    {name:'Q hot kW', y:sc.Q_hot_kW, c:color('--red')},
    {name:'mdot in kg/s x1000', y:sc.mdot_c_inlet_face.map(v=>v*1000), c:color('--blue')},
    {name:'mdot out kg/s x1000', y:sc.mdot_c_outlet_face.map(v=>v*1000), c:color('--green')},
  ], 'mixed scale');
  const pLine=t.map(tt=>interp(tS,S.line_pressure_before_hx_bar,tt));
  const pOut=t.map(tt=>interp(tS,S.water_tank_pressure_bar,tt));
  if(DATA.summary.config.coolant_momentum_model === 'low_mach'){
    const pmin=f.p_c.map(row=>Math.min(...row)/1e5), pmax=f.p_c.map(row=>Math.max(...row)/1e5);
    plot('pHx', t, [
      {name:'HX p min', y:pmin, c:color('--blue')},
      {name:'HX p max', y:pmax, c:color('--orange')},
      {name:'system line', y:pLine, c:color('--red')},
      {name:'system outlet', y:pOut, c:color('--green')},
    ], 'pressure [bar]');
  } else {
    plot('pHx', t, [
      {name:'shell dp estimate', y:sc.dp_shell_total_Pa.map(v=>v/1e5), c:color('--blue')},
      {name:'line - tank', y:pLine.map((v,i)=>v-pOut[i]), c:color('--green')},
    ], 'pressure drop [bar]');
  }
  const x=HX.x, k=t.length-1;
  plot('profileHx', x, [
    {name:'T gas', y:f.T_g[k], c:color('--red')},
    {name:'T He', y:f.T_c[k], c:color('--blue')},
    {name:'T wall hot', y:f.T_wg[k], c:color('--orange')},
    {name:'T wall cold', y:f.T_wc[k], c:color('--green')},
  ], 'T [K]');
}
function interp(x,y,xx){if(xx<=x[0])return y[0]; for(let i=1;i<x.length;i++){if(xx<=x[i]){const a=(xx-x[i-1])/(x[i]-x[i-1]||1); return y[i-1]+a*(y[i]-y[i-1])}} return y[y.length-1]}
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()

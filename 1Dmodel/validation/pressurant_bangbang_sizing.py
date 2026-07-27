"""Simple helium bang-bang pressurization sizing model.

Run from the repository root:

    python -m hps_combustor.validation.pressurant_bangbang_sizing

This is a system-level surrogate, not a detailed ESPSS replacement. It sizes a
parallel valve/orifice helium feed from a 400 bar, 100 K, 265 L helium tank into
a nominal 80 bar pre-HX line, then through a simplified HX/feed pressure loss
into a 3000 L water tank ullage at about 70 bar. Water exits through a
calibrated orifice at a target 30 L/s.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


GAMMA_HE = 1.667
R_HE = 2077.1  # J/kg/K
RHO_WATER = 1000.0


@dataclass(frozen=True)
class PressurantSystemConfig:
    helium_tank_pressure_Pa: float = 400.0e5
    helium_tank_temperature_K: float = 100.0
    helium_tank_volume_m3: float = 0.265
    water_tank_volume_m3: float = 3.0
    initial_ullage_fraction: float = 0.02
    pressurant_temperature_K: float = 300.0
    target_line_pressure_Pa: float = 80.0e5
    target_water_tank_pressure_Pa: float = 70.0e5
    target_water_flow_m3_s: float = 0.030
    downstream_water_pressure_Pa: float = 1.0e5
    water_exit_Cd: float = 0.70
    water_exit_orifice_diameter_m: float | None = None
    feed_Cd: float = 0.80
    hx_pressure_loss_nominal_Pa: float = 10.0e5
    hx_nominal_helium_mdot_kg_s: float = 0.33
    hysteresis_Pa: float = 0.25e5
    t_end_s: float = 100.0
    dt_s: float = 0.005


@dataclass(frozen=True)
class FeedDesign:
    n_branches: int
    orifice_diameter_m: float
    valve_equivalent_diameter_m: float
    control_frequency_Hz: float


@dataclass(frozen=True)
class SimulationSummary:
    score: float
    mean_pressure_bar: float
    mean_line_pressure_bar: float
    pressure_ripple_bar: float
    line_pressure_ripple_bar: float
    mean_water_flow_L_s: float
    delivered_water_L: float
    helium_used_kg: float
    helium_remaining_kg: float
    duty_cycle: float
    switch_count: int
    min_supply_pressure_bar: float
    final_supply_pressure_bar: float
    min_supply_temperature_K: float
    final_supply_temperature_K: float
    max_branch_mdot_kg_s: float
    water_exit_orifice_diameter_mm: float
    branch_orifice_diameter_mm: float
    valve_equivalent_diameter_mm: float
    valve_equivalent_Kv_m3_h: float
    n_branches: int
    control_frequency_Hz: float


def run_design(config: PressurantSystemConfig, design: FeedDesign):
    """Simulate one bang-bang design and return `(summary, history)`."""

    if design.n_branches <= 0:
        raise ValueError("n_branches must be positive")
    if design.control_frequency_Hz <= 0.0:
        raise ValueError("control_frequency_Hz must be positive")

    exit_area = water_exit_orifice_area(config)
    branch_cda = branch_effective_cda(config, design)
    initial_ullage = config.water_tank_volume_m3 * config.initial_ullage_fraction
    initial_water_volume = config.water_tank_volume_m3 - initial_ullage
    helium_supply_mass = ideal_gas_mass(
        config.helium_tank_pressure_Pa,
        config.helium_tank_volume_m3,
        config.helium_tank_temperature_K,
    )
    helium_supply_mass_initial = helium_supply_mass
    helium_supply_temperature = config.helium_tank_temperature_K
    tank_helium_mass = ideal_gas_mass(
        config.target_water_tank_pressure_Pa,
        initial_ullage,
        config.pressurant_temperature_K,
    )

    n_steps = int(np.ceil(config.t_end_s / config.dt_s)) + 1
    history = {
        "time_s": np.zeros(n_steps),
        "water_tank_pressure_bar": np.zeros(n_steps),
        "line_pressure_before_hx_bar": np.zeros(n_steps),
        "supply_pressure_bar": np.zeros(n_steps),
        "supply_temperature_K": np.zeros(n_steps),
        "water_flow_L_s": np.zeros(n_steps),
        "helium_mdot_kg_s": np.zeros(n_steps),
        "open_branches": np.zeros(n_steps),
        "water_volume_m3": np.zeros(n_steps),
    }

    water_volume = initial_water_volume
    open_branches = 0
    next_control_time = 0.0
    control_period = 1.0 / design.control_frequency_Hz
    switch_count = 0

    for k in range(n_steps):
        t = min(k * config.dt_s, config.t_end_s)
        ullage = max(config.water_tank_volume_m3 - water_volume, 1.0e-6)
        tank_pressure = ideal_gas_pressure(
            tank_helium_mass,
            ullage,
            config.pressurant_temperature_K,
        )
        supply_pressure = ideal_gas_pressure(
            helium_supply_mass,
            config.helium_tank_volume_m3,
            helium_supply_temperature,
        )
        nominal_line_pressure = tank_pressure + config.hx_pressure_loss_nominal_Pa

        if t >= next_control_time - 0.5 * config.dt_s:
            old = open_branches
            open_branches = commanded_open_branches(
                nominal_line_pressure,
                target_pressure=config.target_line_pressure_Pa,
                hysteresis=config.hysteresis_Pa,
                n_branches=design.n_branches,
            )
            if old != open_branches:
                switch_count += 1
            next_control_time += control_period

        q_water = water_exit_flow(config, tank_pressure, exit_area)
        mdot_he = merged_feed_mdot(
            config,
            open_branches=open_branches,
            supply_pressure=supply_pressure,
            supply_temperature=helium_supply_temperature,
            tank_pressure=tank_pressure,
            branch_cda=branch_cda,
        )
        mdot_branch = mdot_he / max(open_branches, 1)
        mdot_he = min(mdot_he, helium_supply_mass / max(config.dt_s, 1.0e-12))

        history["time_s"][k] = t
        history["water_tank_pressure_bar"][k] = tank_pressure / 1.0e5
        history["line_pressure_before_hx_bar"][k] = (
            tank_pressure + hx_pressure_loss(config, mdot_he)
        ) / 1.0e5
        history["supply_pressure_bar"][k] = supply_pressure / 1.0e5
        history["supply_temperature_K"][k] = helium_supply_temperature
        history["water_flow_L_s"][k] = q_water * 1000.0
        history["helium_mdot_kg_s"][k] = mdot_he
        history["open_branches"][k] = open_branches
        history["water_volume_m3"][k] = water_volume

        if k == n_steps - 1:
            break

        dt = min(config.dt_s, config.t_end_s - t)
        water_out = min(q_water * dt, water_volume)
        water_volume -= water_out
        helium_supply_mass -= mdot_he * dt
        helium_supply_temperature = adiabatic_tank_temperature(
            mass=helium_supply_mass,
            initial_mass=helium_supply_mass_initial,
            initial_temperature=config.helium_tank_temperature_K,
        )
        tank_helium_mass += mdot_he * dt

    summary = summarize(config, design, history, helium_supply_mass, switch_count, branch_cda)
    return summary, history


def sweep_designs(config: PressurantSystemConfig):
    """Sweep a practical grid and return sorted `(summary, design)` pairs."""

    candidates: list[tuple[SimulationSummary, FeedDesign]] = []
    for n_branches in (2, 3):
        for orifice_mm in np.linspace(2.5, 5.5, 7):
            for valve_mm in np.linspace(2.5, 5.5, 7):
                for freq in (10.0, 20.0, 40.0):
                    design = FeedDesign(
                        n_branches=n_branches,
                        orifice_diameter_m=float(orifice_mm) * 1.0e-3,
                        valve_equivalent_diameter_m=float(valve_mm) * 1.0e-3,
                        control_frequency_Hz=freq,
                    )
                    summary, _history = run_design(config, design)
                    candidates.append((summary, design))
    candidates.sort(key=lambda item: item[0].score)
    return candidates


def commanded_open_branches(
    pressure: float,
    *,
    target_pressure: float,
    hysteresis: float,
    n_branches: int,
) -> int:
    """Staged bang-bang command: 0..n open branches from pressure error."""

    if n_branches <= 0:
        return 0
    error = float(target_pressure) - float(pressure)
    if error <= -float(hysteresis):
        return 0
    if error >= float(hysteresis):
        return int(n_branches)
    fraction = (error + float(hysteresis)) / max(2.0 * float(hysteresis), 1.0e-12)
    return int(np.clip(np.rint(fraction * int(n_branches)), 0, int(n_branches)))


def water_exit_orifice_area(config: PressurantSystemConfig) -> float:
    if config.water_exit_orifice_diameter_m is not None:
        diameter = float(config.water_exit_orifice_diameter_m)
        if diameter <= 0.0:
            raise ValueError("water_exit_orifice_diameter_m must be positive when provided")
        return float(np.pi * diameter**2 / 4.0)
    dp = max(config.target_water_tank_pressure_Pa - config.downstream_water_pressure_Pa, 1.0)
    velocity_factor = np.sqrt(2.0 * dp / RHO_WATER)
    return config.target_water_flow_m3_s / (config.water_exit_Cd * velocity_factor)


def branch_effective_cda(config: PressurantSystemConfig, design: FeedDesign) -> float:
    area_orifice = np.pi * design.orifice_diameter_m**2 / 4.0
    area_valve = np.pi * design.valve_equivalent_diameter_m**2 / 4.0
    cda_orifice = config.feed_Cd * area_orifice
    cda_valve = config.feed_Cd * area_valve
    return 1.0 / np.sqrt(1.0 / cda_orifice**2 + 1.0 / cda_valve**2)


def helium_orifice_mdot(
    *,
    upstream_pressure: float,
    downstream_pressure: float,
    upstream_temperature: float,
    cda: float,
) -> float:
    """Ideal-gas compressible orifice flow from upstream to downstream."""

    if cda <= 0.0 or upstream_pressure <= downstream_pressure:
        return 0.0
    gamma = GAMMA_HE
    pressure_ratio = max(downstream_pressure / upstream_pressure, 0.0)
    critical_ratio = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    scale = cda * upstream_pressure / np.sqrt(R_HE * upstream_temperature)
    if pressure_ratio <= critical_ratio:
        return float(
            scale
            * np.sqrt(gamma)
            * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
        )
    term = pressure_ratio ** (2.0 / gamma) - pressure_ratio ** ((gamma + 1.0) / gamma)
    return float(scale * np.sqrt((2.0 * gamma / (gamma - 1.0)) * max(term, 0.0)))


def merged_feed_mdot(
    config: PressurantSystemConfig,
    *,
    open_branches: int,
    supply_pressure: float,
    supply_temperature: float,
    tank_pressure: float,
    branch_cda: float,
) -> float:
    """Return total helium flow through open branches and common HX/feed loss."""

    if open_branches <= 0 or supply_pressure <= tank_pressure:
        return 0.0
    if config.hx_pressure_loss_nominal_Pa <= 0.0:
        return float(open_branches) * helium_orifice_mdot(
            upstream_pressure=supply_pressure,
            downstream_pressure=tank_pressure,
            upstream_temperature=supply_temperature,
            cda=branch_cda,
        )

    mdot_nom = max(config.hx_nominal_helium_mdot_kg_s, 1.0e-9)
    k_hx = config.hx_pressure_loss_nominal_Pa / mdot_nom**2
    high = float(open_branches) * helium_orifice_mdot(
        upstream_pressure=supply_pressure,
        downstream_pressure=tank_pressure,
        upstream_temperature=supply_temperature,
        cda=branch_cda,
    )
    if high <= 0.0:
        return 0.0
    low = 0.0
    for _ in range(32):
        mid = 0.5 * (low + high)
        downstream = tank_pressure + k_hx * mid * mid
        available = float(open_branches) * helium_orifice_mdot(
            upstream_pressure=supply_pressure,
            downstream_pressure=downstream,
            upstream_temperature=supply_temperature,
            cda=branch_cda,
        )
        if available >= mid:
            low = mid
        else:
            high = mid
    return low


def hx_pressure_loss(config: PressurantSystemConfig, mdot_he: float) -> float:
    mdot_nom = max(config.hx_nominal_helium_mdot_kg_s, 1.0e-9)
    return float(config.hx_pressure_loss_nominal_Pa * (max(float(mdot_he), 0.0) / mdot_nom) ** 2)


def water_exit_flow(config: PressurantSystemConfig, pressure: float, area: float) -> float:
    dp = max(float(pressure) - config.downstream_water_pressure_Pa, 0.0)
    if dp <= 0.0:
        return 0.0
    return float(config.water_exit_Cd * area * np.sqrt(2.0 * dp / RHO_WATER))


def ideal_gas_mass(pressure: float, volume: float, temperature: float) -> float:
    return float(pressure * volume / (R_HE * temperature))


def ideal_gas_pressure(mass: float, volume: float, temperature: float) -> float:
    return float(max(mass, 0.0) * R_HE * temperature / max(volume, 1.0e-12))


def adiabatic_tank_temperature(*, mass: float, initial_mass: float, initial_temperature: float) -> float:
    mass_ratio = max(float(mass), 1.0e-12) / max(float(initial_mass), 1.0e-12)
    return float(initial_temperature * mass_ratio ** (GAMMA_HE - 1.0))


def summarize(
    config: PressurantSystemConfig,
    design: FeedDesign,
    history: dict[str, np.ndarray],
    helium_supply_mass_final: float,
    switch_count: int,
    branch_cda: float,
) -> SimulationSummary:
    t = history["time_s"]
    mask = t >= min(5.0, 0.1 * config.t_end_s)
    p = history["water_tank_pressure_bar"][mask]
    line_p = history["line_pressure_before_hx_bar"][mask]
    q = history["water_flow_L_s"][mask]
    mdot = history["helium_mdot_kg_s"]
    open_branches = history["open_branches"]
    initial_supply_mass = ideal_gas_mass(
        config.helium_tank_pressure_Pa,
        config.helium_tank_volume_m3,
        config.helium_tank_temperature_K,
    )
    helium_used = initial_supply_mass - helium_supply_mass_final
    delivered_water_L = float(
        np.trapezoid(history["water_flow_L_s"], history["time_s"])
    )
    pressure_error = abs(float(np.mean(p)) - config.target_water_tank_pressure_Pa / 1.0e5)
    line_pressure_error = abs(float(np.mean(line_p)) - config.target_line_pressure_Pa / 1.0e5)
    flow_error = abs(float(np.mean(q)) - config.target_water_flow_m3_s * 1000.0)
    pressure_ripple = float(np.max(p) - np.min(p)) if p.size else float("nan")
    line_pressure_ripple = float(np.max(line_p) - np.min(line_p)) if line_p.size else float("nan")
    duty = float(np.mean(open_branches > 0.0))
    score = (
        pressure_error
        + 1.5 * line_pressure_error
        + 0.5 * pressure_ripple
        + 0.35 * line_pressure_ripple
        + 2.0 * flow_error
    )

    valve_area = np.pi * design.valve_equivalent_diameter_m**2 / 4.0
    valve_kv = equivalent_water_kv(config.feed_Cd * valve_area)
    exit_diameter = np.sqrt(4.0 * water_exit_orifice_area(config) / np.pi)
    return SimulationSummary(
        score=float(score),
        mean_pressure_bar=float(np.mean(p)),
        mean_line_pressure_bar=float(np.mean(line_p)),
        pressure_ripple_bar=pressure_ripple,
        line_pressure_ripple_bar=line_pressure_ripple,
        mean_water_flow_L_s=float(np.mean(q)),
        delivered_water_L=delivered_water_L,
        helium_used_kg=float(helium_used),
        helium_remaining_kg=float(helium_supply_mass_final),
        duty_cycle=duty,
        switch_count=int(switch_count),
        min_supply_pressure_bar=float(np.min(history["supply_pressure_bar"])),
        final_supply_pressure_bar=float(history["supply_pressure_bar"][-1]),
        min_supply_temperature_K=float(np.min(history["supply_temperature_K"])),
        final_supply_temperature_K=float(history["supply_temperature_K"][-1]),
        max_branch_mdot_kg_s=float(np.max(mdot) / max(design.n_branches, 1)),
        water_exit_orifice_diameter_mm=float(exit_diameter * 1000.0),
        branch_orifice_diameter_mm=float(design.orifice_diameter_m * 1000.0),
        valve_equivalent_diameter_mm=float(design.valve_equivalent_diameter_m * 1000.0),
        valve_equivalent_Kv_m3_h=float(valve_kv),
        n_branches=int(design.n_branches),
        control_frequency_Hz=float(design.control_frequency_Hz),
    )


def equivalent_water_kv(cda: float) -> float:
    """Return Kv [m3/h] for an equivalent water valve at 1 bar pressure drop."""

    q_m3_s = float(cda) * np.sqrt(2.0 * 1.0e5 / RHO_WATER)
    return float(q_m3_s * 3600.0)


def write_history_csv(path: Path, history: dict[str, np.ndarray]) -> None:
    keys = list(history)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(keys)
        for row in zip(*(history[key] for key in keys)):
            writer.writerow([f"{float(value):.9g}" for value in row])


def run_sizing(output_dir: str | Path = "docs/validation/pressurant_bangbang") -> dict:
    config = PressurantSystemConfig()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    ranked = sweep_designs(config)
    best_summary, best_design = ranked[0]
    best_summary, best_history = run_design(config, best_design)

    payload = {
        "config": asdict(config),
        "best_design": asdict(best_design),
        "best_summary": asdict(best_summary),
        "top_candidates": [
            {"summary": asdict(summary), "design": asdict(design)}
            for summary, design in ranked[:20]
        ],
    }
    (output / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_history_csv(output / "best_timeseries.csv", best_history)
    return payload


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="docs/validation/pressurant_bangbang",
        help="Folder for summary.json and best_timeseries.csv.",
    )
    args = parser.parse_args(argv)
    payload = run_sizing(args.output_dir)
    best = payload["best_summary"]
    design = payload["best_design"]
    print(f"Wrote {args.output_dir}")
    print(
        "best: "
        f"{design['n_branches']} branches, "
        f"orifice={best['branch_orifice_diameter_mm']:.2f} mm, "
        f"valve_d_eq={best['valve_equivalent_diameter_mm']:.2f} mm, "
        f"Kv={best['valve_equivalent_Kv_m3_h']:.3f} m3/h, "
        f"f={design['control_frequency_Hz']:.1f} Hz"
    )
    print(
        f"water_p={best['mean_pressure_bar']:.2f} bar "
        f"line_p={best['mean_line_pressure_bar']:.2f} bar, "
        f"water={best['mean_water_flow_L_s']:.2f} L/s, "
        f"He used={best['helium_used_kg']:.2f} kg, "
        f"supply_final={best['final_supply_pressure_bar']:.1f} bar/"
        f"{best['final_supply_temperature_K']:.1f} K"
    )
    return payload


if __name__ == "__main__":
    main()

"""Reusable 1D heated-channel liquid coolant governing equations.

State variables are pressure and enthalpy ``(p, h)`` throughout, never
temperature: temperature is not a valid marching/convergence state inside the
two-phase dome (it plateaus at ``Tsat(p)`` while enthalpy keeps rising). See
docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md, Design Decision 1.

The solver is intentionally geometry-agnostic: callers provide flow area,
heated perimeter, hydraulic diameter, and a heat-flux boundary condition. That
is the integration surface needed by straight pipes, coils, and future HX
wrappers while keeping the liquid/two-phase state closure in ``p,h`` form.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Callable

import CoolProp.CoolProp as CP
import numpy as np

from hps_combustor.input_data import coolantProp
from hps_combustor.physics.liquid_flow.correlations import (
    homogeneous_acceleration_pressure_gradient,
    saturation_state,
)
from hps_combustor.physics.liquid_flow.dispatch import evaluate_coolant_closure

HeatFluxBC = float | Callable[[float], float] | np.ndarray


@dataclass(frozen=True)
class HeatedChannelCase:
    fluid: str = "Water"
    length_m: float = 1.0
    hydraulic_diameter_m: float = 0.010
    flow_area_m2: float | None = None
    heated_perimeter_m: float | None = None
    mass_flow_kg_s: float = 0.05
    p_in_Pa: float = 1.0e5
    h_in_J_kg: float | None = None
    T_in_K: float = 370.0
    heat_flux_W_m2: HeatFluxBC = 1.0e5
    n_cells: int = 80
    coolant_model: str = "equilibrium_liquid"
    lut_path: str | Path | None = None
    min_pressure_Pa: float = 1.0

    @property
    def effective_flow_area_m2(self) -> float:
        if self.flow_area_m2 is not None:
            return float(self.flow_area_m2)
        return math.pi * self.hydraulic_diameter_m**2 / 4.0

    @property
    def effective_heated_perimeter_m(self) -> float:
        if self.heated_perimeter_m is not None:
            return float(self.heated_perimeter_m)
        return math.pi * self.hydraulic_diameter_m


@dataclass(frozen=True)
class HeatedChannelResult:
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
    chf_W_m2: np.ndarray
    chf_margin: np.ndarray
    heat_flux_W_m2: np.ndarray
    heat_rate_W: float
    energy_residual_J_kg: float

    @property
    def outlet_quality(self) -> float:
        return float(self.quality[-1])

    @property
    def pressure_drop_Pa(self) -> float:
        return float(self.p_Pa[0] - self.p_Pa[-1])

    @property
    def min_chf_margin(self) -> float:
        finite = self.chf_margin[np.isfinite(self.chf_margin)]
        return float(np.min(finite)) if finite.size else float("nan")


@dataclass(frozen=True)
class HeatedChannelDiagnostics:
    """Scalar diagnostics intended for HX integration gates and warnings."""

    heat_rate_W: float
    pressure_drop_Pa: float
    inlet_T_K: float
    outlet_T_K: float
    inlet_p_Pa: float
    outlet_p_Pa: float
    inlet_h_J_kg: float
    outlet_h_J_kg: float
    min_quality: float
    max_quality: float
    outlet_quality: float
    max_void_fraction: float
    min_chf_margin: float
    boiling_reached: bool
    dryout_or_vapor_reached: bool
    chf_margin_below_limit: bool
    pressure_floor_reached: bool
    energy_residual_abs_J_kg: float
    energy_residual_ok: bool


@dataclass(frozen=True)
class HXGridHeatedChannelResult:
    """Liquid-channel result mapped back to an HX solver's axial grid order."""

    flow_result: HeatedChannelResult
    cell_fields_hx_order: dict[str, np.ndarray]
    node_fields_hx_order: dict[str, np.ndarray]
    diagnostics: HeatedChannelDiagnostics
    coolant_enters_at: str


@dataclass(frozen=True)
class HeatedChannelProfileCase:
    """HX-facing 1D coolant march inputs on a segment grid.

    ``z_edges_m`` has length ``n_cells + 1``. Geometry and heating arrays are
    cell-centered with length ``n_cells``. Scalars are broadcast to all cells.
    Supply exactly one heat input form: ``heat_flux_W_m2``,
    ``heat_per_length_W_m``, or ``heat_per_segment_W``.
    """

    coolant_prop: coolantProp
    z_edges_m: np.ndarray
    hydraulic_diameter_m: float | np.ndarray
    flow_area_m2: float | np.ndarray
    heated_perimeter_m: float | np.ndarray
    mass_flow_kg_s: float
    p_in_Pa: float
    h_in_J_kg: float | None = None
    T_in_K: float | None = None
    heat_flux_W_m2: float | np.ndarray | None = None
    heat_per_length_W_m: float | np.ndarray | None = None
    heat_per_segment_W: float | np.ndarray | None = None
    lut_path: str | Path | None = None
    min_pressure_Pa: float = 1.0


def heated_channel_cell_fields(result: HeatedChannelResult) -> dict[str, np.ndarray]:
    """Return cell-centered liquid fields for HX wall/plotting integration.

    Nodal thermodynamic fields are averaged onto cells. Quantities already
    defined by the marching segment, such as heat flux and pressure gradients,
    use the upstream cell value.
    """
    return {
        "z_m": 0.5 * (result.z_m[:-1] + result.z_m[1:]),
        "p_Pa": 0.5 * (result.p_Pa[:-1] + result.p_Pa[1:]),
        "h_J_kg": 0.5 * (result.h_J_kg[:-1] + result.h_J_kg[1:]),
        "T_K": 0.5 * (result.T_K[:-1] + result.T_K[1:]),
        "quality": 0.5 * (result.quality[:-1] + result.quality[1:]),
        "void_fraction": 0.5 * (result.void_fraction[:-1] + result.void_fraction[1:]),
        "rho_kg_m3": 0.5 * (result.rho_kg_m3[:-1] + result.rho_kg_m3[1:]),
        "htc_W_m2_K": 0.5 * (result.htc_W_m2_K[:-1] + result.htc_W_m2_K[1:]),
        "dpdz_friction_Pa_m": result.dpdz_friction_Pa_m[:-1],
        "dpdz_acceleration_Pa_m": result.dpdz_acceleration_Pa_m[:-1],
        "chf_W_m2": result.chf_W_m2[:-1],
        "chf_margin": result.chf_margin[:-1],
        "heat_flux_W_m2": result.heat_flux_W_m2[:-1],
    }


def heated_channel_node_fields(result: HeatedChannelResult) -> dict[str, np.ndarray]:
    """Return nodal liquid fields for HX integration and plotting."""
    return {
        "z_m": result.z_m,
        "p_Pa": result.p_Pa,
        "h_J_kg": result.h_J_kg,
        "T_K": result.T_K,
        "quality": result.quality,
        "void_fraction": result.void_fraction,
        "rho_kg_m3": result.rho_kg_m3,
        "htc_W_m2_K": result.htc_W_m2_K,
        "chf_W_m2": result.chf_W_m2,
        "chf_margin": result.chf_margin,
        "heat_flux_W_m2": result.heat_flux_W_m2,
    }


def summarize_heated_channel_result(
    result: HeatedChannelResult,
    *,
    min_pressure_Pa: float = 1.0,
    chf_margin_limit: float = 1.0,
    energy_residual_tol_J_kg: float = 1.0e-6,
) -> HeatedChannelDiagnostics:
    """Summarize liquid-channel health for HX integration checks."""
    finite_quality = result.quality[np.isfinite(result.quality)]
    min_quality = float(np.min(finite_quality)) if finite_quality.size else float("nan")
    max_quality = float(np.max(finite_quality)) if finite_quality.size else float("nan")
    min_margin = result.min_chf_margin
    pressure_floor_reached = bool(np.any(result.p_Pa <= min_pressure_Pa * (1.0 + 1.0e-12)))
    return HeatedChannelDiagnostics(
        heat_rate_W=float(result.heat_rate_W),
        pressure_drop_Pa=result.pressure_drop_Pa,
        inlet_T_K=float(result.T_K[0]),
        outlet_T_K=float(result.T_K[-1]),
        inlet_p_Pa=float(result.p_Pa[0]),
        outlet_p_Pa=float(result.p_Pa[-1]),
        inlet_h_J_kg=float(result.h_J_kg[0]),
        outlet_h_J_kg=float(result.h_J_kg[-1]),
        min_quality=min_quality,
        max_quality=max_quality,
        outlet_quality=result.outlet_quality,
        max_void_fraction=float(np.nanmax(result.void_fraction)),
        min_chf_margin=min_margin,
        boiling_reached=bool(np.any((result.quality >= 0.0) & (result.quality <= 1.0))),
        dryout_or_vapor_reached=bool(np.any(result.quality >= 1.0)),
        chf_margin_below_limit=bool(np.isfinite(min_margin) and min_margin < chf_margin_limit),
        pressure_floor_reached=pressure_floor_reached,
        energy_residual_abs_J_kg=float(abs(result.energy_residual_J_kg)),
        energy_residual_ok=bool(abs(result.energy_residual_J_kg) <= energy_residual_tol_J_kg),
    )


def _reverse_cell_input(value: float | np.ndarray | None) -> float | np.ndarray | None:
    if value is None or np.isscalar(value):
        return value
    return np.asarray(value, dtype=float)[::-1]


def _fields_to_hx_order(
    fields: dict[str, np.ndarray],
    *,
    z_values: np.ndarray,
    reverse: bool,
) -> dict[str, np.ndarray]:
    mapped = {}
    for key, value in fields.items():
        if key == "z_m":
            mapped[key] = np.asarray(z_values, dtype=float)
        else:
            mapped[key] = np.asarray(value, dtype=float)[::-1] if reverse else np.asarray(value, dtype=float)
    return mapped


def inlet_enthalpy(case: HeatedChannelCase) -> float:
    if case.h_in_J_kg is not None:
        return float(case.h_in_J_kg)
    return float(CP.PropsSI("H", "P", case.p_in_Pa, "T", case.T_in_K, case.fluid))


def _profile_inlet_enthalpy(case: HeatedChannelProfileCase) -> float:
    if case.h_in_J_kg is not None:
        return float(case.h_in_J_kg)
    if case.T_in_K is None:
        raise ValueError("T_in_K is required when h_in_J_kg is not supplied")
    return float(
        CP.PropsSI(
            "H",
            "P",
            case.p_in_Pa,
            "T",
            float(case.T_in_K),
            case.coolant_prop.coolant,
        )
    )


def heat_flux_profile(case: HeatedChannelCase) -> np.ndarray:
    """Return cell-centered heat fluxes for a scalar, callable, or array BC."""
    if case.n_cells < 1:
        raise ValueError("n_cells must be at least 1")
    z_centers = (np.arange(case.n_cells, dtype=float) + 0.5) * case.length_m / case.n_cells
    bc = case.heat_flux_W_m2
    if callable(bc):
        q = np.array([float(bc(float(z))) for z in z_centers], dtype=float)
    elif isinstance(bc, np.ndarray):
        q = np.asarray(bc, dtype=float)
        if q.shape != (case.n_cells,):
            raise ValueError("array heat_flux_W_m2 must have shape (n_cells,)")
    else:
        q = np.full(case.n_cells, float(bc), dtype=float)
    if not np.all(np.isfinite(q)):
        raise ValueError("heat flux profile must be finite")
    return q


def validate_heated_channel_case(case: HeatedChannelCase) -> None:
    if case.length_m <= 0.0:
        raise ValueError("length must be positive")
    if case.hydraulic_diameter_m <= 0.0:
        raise ValueError("hydraulic diameter must be positive")
    if case.effective_flow_area_m2 <= 0.0:
        raise ValueError("flow area must be positive")
    if case.effective_heated_perimeter_m <= 0.0:
        raise ValueError("heated perimeter must be positive")
    if case.mass_flow_kg_s <= 0.0:
        raise ValueError("mass flow must be positive")
    if case.p_in_Pa <= 0.0:
        raise ValueError("inlet pressure must be positive")
    if case.min_pressure_Pa <= 0.0:
        raise ValueError("minimum pressure must be positive")
    if case.coolant_model not in ("single_phase_coolprop", "equilibrium_liquid"):
        raise ValueError(f"unsupported coolant_model: {case.coolant_model!r}")


def _cell_array(value: float | np.ndarray, n_cells: int, name: str) -> np.ndarray:
    if np.isscalar(value):
        arr = np.full(n_cells, float(value), dtype=float)
    else:
        arr = np.asarray(value, dtype=float)
        if arr.shape != (n_cells,):
            raise ValueError(f"{name} must be scalar or have shape (n_cells,)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _validate_profile_case(case: HeatedChannelProfileCase) -> tuple[np.ndarray, np.ndarray]:
    z_edges = np.asarray(case.z_edges_m, dtype=float)
    if z_edges.ndim != 1 or z_edges.size < 2:
        raise ValueError("z_edges_m must be a 1D array with at least two entries")
    dz = np.diff(z_edges)
    if not np.all(np.isfinite(z_edges)) or not np.all(dz > 0.0):
        raise ValueError("z_edges_m must be finite and strictly increasing")
    n_cells = dz.size
    if case.mass_flow_kg_s <= 0.0:
        raise ValueError("mass flow must be positive")
    if case.p_in_Pa <= 0.0:
        raise ValueError("inlet pressure must be positive")
    if case.min_pressure_Pa <= 0.0:
        raise ValueError("minimum pressure must be positive")
    heat_input_count = sum(
        value is not None
        for value in (
            case.heat_flux_W_m2,
            case.heat_per_length_W_m,
            case.heat_per_segment_W,
        )
    )
    if heat_input_count != 1:
        raise ValueError(
            "supply exactly one of heat_flux_W_m2, heat_per_length_W_m, or heat_per_segment_W"
        )
    model = getattr(case.coolant_prop, "coolant_model", "single_phase_coolprop")
    if model not in ("single_phase_coolprop", "equilibrium_liquid"):
        raise ValueError(f"unsupported coolant_model: {model!r}")
    return z_edges, dz


def solve_steady_heated_channel_profile(
    case: HeatedChannelProfileCase,
) -> HeatedChannelResult:
    """March a heated coolant profile on an HX-provided segment grid."""
    z_edges, dz = _validate_profile_case(case)
    n_cells = dz.size
    D_h = _cell_array(case.hydraulic_diameter_m, n_cells, "hydraulic_diameter_m")
    flow_area = _cell_array(case.flow_area_m2, n_cells, "flow_area_m2")
    heated_perimeter = _cell_array(case.heated_perimeter_m, n_cells, "heated_perimeter_m")
    if np.any(D_h <= 0.0) or np.any(flow_area <= 0.0) or np.any(heated_perimeter <= 0.0):
        raise ValueError("hydraulic diameter, flow area, and heated perimeter must be positive")
    if case.heat_flux_W_m2 is not None:
        q_profile = _cell_array(case.heat_flux_W_m2, n_cells, "heat_flux_W_m2")
        heat_per_length = q_profile * heated_perimeter
    elif case.heat_per_length_W_m is not None:
        heat_per_length = _cell_array(
            case.heat_per_length_W_m,
            n_cells,
            "heat_per_length_W_m",
        )
        q_profile = heat_per_length / heated_perimeter
    else:
        heat_per_segment = _cell_array(
            case.heat_per_segment_W,
            n_cells,
            "heat_per_segment_W",
        )
        heat_per_length = heat_per_segment / dz
        q_profile = heat_per_length / heated_perimeter
    heat_rate = float(np.sum(heat_per_length * dz))

    z = np.array(z_edges, dtype=float)
    p = np.empty(n_cells + 1)
    h = np.empty(n_cells + 1)
    T = np.empty(n_cells + 1)
    quality = np.empty(n_cells + 1)
    alpha = np.empty(n_cells + 1)
    rho = np.empty(n_cells + 1)
    htc = np.empty(n_cells + 1)
    chf = np.full(n_cells + 1, np.nan)
    margin = np.full(n_cells + 1, np.nan)
    q_nodes = np.empty(n_cells + 1)
    dpdz_f = np.zeros(n_cells + 1)
    dpdz_a = np.zeros(n_cells + 1)

    p[0] = float(case.p_in_Pa)
    h[0] = _profile_inlet_enthalpy(case)
    fluid = case.coolant_prop.coolant
    model = getattr(case.coolant_prop, "coolant_model", "single_phase_coolprop")

    for i in range(n_cells + 1):
        j = min(i, n_cells - 1)
        q_local = float(q_profile[j])
        mass_flux = case.mass_flow_kg_s / flow_area[j]
        closure = evaluate_coolant_closure(
            coolant_prop=case.coolant_prop,
            p_Pa=p[i],
            h_J_kg=h[i],
            mass_flux_kg_m2_s=mass_flux,
            hydraulic_diameter_m=D_h[j],
            heat_flux_W_m2=max(q_local, 0.0),
            lut_path=case.lut_path,
        )
        state = closure.state
        T[i] = state.T_K
        quality[i] = state.quality
        alpha[i] = state.void_fraction
        rho[i] = state.rho_kg_m3
        htc[i] = closure.htc_W_m2_K
        dpdz_f[i] = closure.dpdz_friction_Pa_m
        chf[i] = np.nan if closure.chf_W_m2 is None else closure.chf_W_m2
        margin[i] = np.nan if closure.chf_margin is None else closure.chf_margin
        q_nodes[i] = q_local

        if i == n_cells:
            break

        dhdz = heat_per_length[i] / case.mass_flow_kg_s
        if model == "equilibrium_liquid" and 0.0 <= state.quality <= 1.0:
            sat = saturation_state(fluid, p[i])
            quality_gradient = dhdz / sat.h_fg_J_kg
            dpdz_a[i] = homogeneous_acceleration_pressure_gradient(
                mass_flux_kg_m2_s=mass_flux,
                p_Pa=p[i],
                quality_gradient_1_m=quality_gradient,
                fluid=fluid,
            )

        h[i + 1] = h[i] + dhdz * dz[i]
        p[i + 1] = max(p[i] - (dpdz_f[i] + dpdz_a[i]) * dz[i], case.min_pressure_Pa)

    expected_dh = heat_rate / case.mass_flow_kg_s
    energy_residual = (h[-1] - h[0]) - expected_dh
    return HeatedChannelResult(
        z_m=z,
        p_Pa=p,
        h_J_kg=h,
        T_K=T,
        quality=quality,
        void_fraction=alpha,
        rho_kg_m3=rho,
        htc_W_m2_K=htc,
        dpdz_friction_Pa_m=dpdz_f,
        dpdz_acceleration_Pa_m=dpdz_a,
        chf_W_m2=chf,
        chf_margin=margin,
        heat_flux_W_m2=q_nodes,
        heat_rate_W=heat_rate,
        energy_residual_J_kg=float(energy_residual),
    )


def solve_steady_heated_channel_on_hx_grid(
    *,
    coolant_prop: coolantProp,
    z_edges_m: np.ndarray,
    hydraulic_diameter_m: float | np.ndarray,
    flow_area_m2: float | np.ndarray,
    heated_perimeter_m: float | np.ndarray,
    mass_flow_kg_s: float,
    p_in_Pa: float,
    h_in_J_kg: float | None = None,
    T_in_K: float | None = None,
    heat_flux_W_m2: float | np.ndarray | None = None,
    heat_per_length_W_m: float | np.ndarray | None = None,
    heat_per_segment_W: float | np.ndarray | None = None,
    coolant_enters_at: str = "z_min",
    lut_path: str | Path | None = None,
    min_pressure_Pa: float = 1.0,
) -> HXGridHeatedChannelResult:
    """Solve liquid coolant on an HX axial grid and map fields back to HX order.

    ``z_edges_m`` and all non-scalar cell arrays are supplied in the HX solver's
    native axial order. For counterflow, pass ``coolant_enters_at="z_max"``;
    this reverses cell inputs for the march and reverses outputs back to HX
    order.
    """
    z_hx = np.asarray(z_edges_m, dtype=float)
    if z_hx.ndim != 1 or z_hx.size < 2:
        raise ValueError("z_edges_m must be a 1D array with at least two entries")
    if not np.all(np.isfinite(z_hx)) or not np.all(np.diff(z_hx) > 0.0):
        raise ValueError("z_edges_m must be finite and strictly increasing")
    if coolant_enters_at not in ("z_min", "z_max"):
        raise ValueError("coolant_enters_at must be 'z_min' or 'z_max'")

    reverse = coolant_enters_at == "z_max"
    if reverse:
        z_flow = z_hx[-1] - z_hx[::-1]
    else:
        z_flow = z_hx - z_hx[0]

    profile_case = HeatedChannelProfileCase(
        coolant_prop=coolant_prop,
        z_edges_m=z_flow,
        hydraulic_diameter_m=_reverse_cell_input(hydraulic_diameter_m) if reverse else hydraulic_diameter_m,
        flow_area_m2=_reverse_cell_input(flow_area_m2) if reverse else flow_area_m2,
        heated_perimeter_m=_reverse_cell_input(heated_perimeter_m) if reverse else heated_perimeter_m,
        mass_flow_kg_s=mass_flow_kg_s,
        p_in_Pa=p_in_Pa,
        h_in_J_kg=h_in_J_kg,
        T_in_K=T_in_K,
        heat_flux_W_m2=_reverse_cell_input(heat_flux_W_m2) if reverse else heat_flux_W_m2,
        heat_per_length_W_m=_reverse_cell_input(heat_per_length_W_m) if reverse else heat_per_length_W_m,
        heat_per_segment_W=_reverse_cell_input(heat_per_segment_W) if reverse else heat_per_segment_W,
        lut_path=lut_path,
        min_pressure_Pa=min_pressure_Pa,
    )
    result = solve_steady_heated_channel_profile(profile_case)
    cell_centers_hx = 0.5 * (z_hx[:-1] + z_hx[1:])
    cell_fields = _fields_to_hx_order(
        heated_channel_cell_fields(result),
        z_values=cell_centers_hx,
        reverse=reverse,
    )
    node_fields = _fields_to_hx_order(
        heated_channel_node_fields(result),
        z_values=z_hx,
        reverse=reverse,
    )
    diagnostics = summarize_heated_channel_result(
        result,
        min_pressure_Pa=min_pressure_Pa,
    )
    return HXGridHeatedChannelResult(
        flow_result=result,
        cell_fields_hx_order=cell_fields,
        node_fields_hx_order=node_fields,
        diagnostics=diagnostics,
        coolant_enters_at=coolant_enters_at,
    )


def solve_steady_heated_channel(case: HeatedChannelCase) -> HeatedChannelResult:
    """March a 1D heated coolant channel using ``p,h`` as the governing state."""
    validate_heated_channel_case(case)
    q_profile = heat_flux_profile(case)

    flow_area = case.effective_flow_area_m2
    heated_perimeter = case.effective_heated_perimeter_m
    props = coolantProp(coolant=case.fluid, coolant_model=case.coolant_model)
    return solve_steady_heated_channel_profile(
        HeatedChannelProfileCase(
            coolant_prop=props,
            z_edges_m=np.linspace(0.0, case.length_m, case.n_cells + 1),
            hydraulic_diameter_m=case.hydraulic_diameter_m,
            flow_area_m2=flow_area,
            heated_perimeter_m=heated_perimeter,
            mass_flow_kg_s=case.mass_flow_kg_s,
            p_in_Pa=case.p_in_Pa,
            h_in_J_kg=case.h_in_J_kg,
            T_in_K=case.T_in_K,
            heat_flux_W_m2=q_profile,
            lut_path=case.lut_path,
            min_pressure_Pa=case.min_pressure_Pa,
        )
    )

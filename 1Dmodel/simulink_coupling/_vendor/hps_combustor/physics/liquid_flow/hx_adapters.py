"""Adapters from HX configuration objects to liquid heated-channel inputs.

These helpers keep the validated liquid coolant march separate from production
HX solvers while making the handoff explicit and testable. They are not a
replacement for geometry-specific boiling correlations; they convert existing
solver geometry and wall duty into the generic 1D liquid profile interface.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from hps_combustor.input_data import coolantProp, combustorProp, numericalProp, shellTubeProp
from hps_combustor.mechanical.geometry.helix_geometry import HelixGeometryRadiusCST
from hps_combustor.mechanical.geometry.shelltube_geometry import compute_bell_delaware_geometry
from hps_combustor.physics.liquid_flow.governing_equations import (
    HXGridHeatedChannelResult,
    solve_steady_heated_channel_on_hx_grid,
)


def _shelltube_geometry(stp: shellTubeProp) -> dict:
    return compute_bell_delaware_geometry(
        D_shell_inner=stp.D_shell_inner,
        D_tube_outer=stp.D_tube_outer,
        pitch_ratio=stp.pitch_ratio,
        layout=stp.layout,
        N_tubes=stp.N_tubes,
        N_baffles=stp.N_baffles,
        baffle_cut=stp.baffle_cut,
        L_tube=stp.L_tube,
        clearance_tube_baffle=stp.clearance_tube_baffle,
        clearance_baffle_shell=stp.clearance_baffle_shell,
        clearance_bundle_shell=stp.clearance_bundle_shell,
        N_sealing_strip_pairs=stp.N_sealing_strip_pairs,
        baffle_spacing=stp.baffle_spacing,
        L_inlet_spacing=stp.L_inlet_spacing,
        L_outlet_spacing=stp.L_outlet_spacing,
    )


def _helical_pipe_length(cp: combustorProp, npv: numericalProp) -> float:
    Dh = float(cp.Dh_coil)
    coil_pitch = Dh + 2.0 * float(cp.thickness_coil_wall) + float(cp.coil_gap)
    D_coil = (
        float(cp.inner_diameter)
        - 2.0 * float(cp.gap_shell2coil)
        - Dh
        - 2.0 * float(cp.thickness_coil_wall)
    )
    if D_coil <= 0.0:
        raise ValueError("computed helical coil centerline diameter must be positive")
    L_coil = (
        (float(npv.L_HX_max) - float(cp.mixing_length))
        - 2.0 * float(cp.length_2_coil)
        - (Dh + 2.0 * float(cp.thickness_coil_wall))
    )
    if L_coil <= 0.0:
        raise ValueError("computed helical coil axial length must be positive")
    _func_s_to_x, _func_s_to_theta, L_ch_max = HelixGeometryRadiusCST(
        coil_pitch=coil_pitch,
        D_coil=D_coil,
        L_coil=L_coil,
    )
    return float(L_ch_max)


def solve_shelltube_shellside_liquid_from_duty(
    *,
    coolant_prop: coolantProp,
    shelltube_prop: shellTubeProp,
    dQ_profile_per_tube_W: np.ndarray,
    coolant_enters_at: str,
    n_axial: int | None = None,
    lut_path: str | Path | None = None,
    min_pressure_Pa: float = 1.0,
) -> HXGridHeatedChannelResult:
    """Run the liquid profile solver from shell-and-tube wall duty.

    Parameters
    ----------
    dQ_profile_per_tube_W:
        Heat transferred from one representative hot-gas tube to the shell-side
        coolant in each axial segment, matching ``shellntube_solver._tube_side_march``.
        This adapter multiplies by ``shelltube_prop.N_tubes`` to get total
        coolant heat input per segment.

    Notes
    -----
    This is a pseudo-1D shell-side adapter. It uses Bell-Delaware ``S_m`` as the
    effective flow area, tube outside diameter as the shell-side characteristic
    length, and total outside tube perimeter as heated perimeter. It is an
    integration bridge for the validated liquid solver, not a final shell-side
    boiling correlation package.
    """
    dQ = np.asarray(dQ_profile_per_tube_W, dtype=float)
    if dQ.ndim != 1 or dQ.size < 1:
        raise ValueError("dQ_profile_per_tube_W must be a non-empty 1D array")
    if not np.all(np.isfinite(dQ)):
        raise ValueError("dQ_profile_per_tube_W must be finite")
    n = int(n_axial) if n_axial is not None else int(dQ.size)
    if n != dQ.size:
        raise ValueError("n_axial must match dQ_profile_per_tube_W length")
    if coolant_enters_at not in ("z_min", "z_max"):
        raise ValueError("coolant_enters_at must be 'z_min' or 'z_max'")

    geom = _shelltube_geometry(shelltube_prop)
    z_edges = np.linspace(0.0, float(shelltube_prop.L_tube), n + 1)
    hydraulic_diameter = float(shelltube_prop.D_tube_outer)
    flow_area = float(geom["S_m"])
    heated_perimeter = float(shelltube_prop.N_tubes * np.pi * shelltube_prop.D_tube_outer)
    heat_per_segment = dQ * float(shelltube_prop.N_tubes)

    return solve_steady_heated_channel_on_hx_grid(
        coolant_prop=coolant_prop,
        z_edges_m=z_edges,
        hydraulic_diameter_m=hydraulic_diameter,
        flow_area_m2=flow_area,
        heated_perimeter_m=heated_perimeter,
        mass_flow_kg_s=coolant_prop.mass_flow_c,
        p_in_Pa=coolant_prop.p_in,
        T_in_K=coolant_prop.T_in,
        heat_per_segment_W=heat_per_segment,
        coolant_enters_at=coolant_enters_at,
        lut_path=lut_path,
        min_pressure_Pa=min_pressure_Pa,
    )


def solve_shelltube_shellside_liquid_from_tube_result(
    *,
    coolant_prop: coolantProp,
    shelltube_prop: shellTubeProp,
    tube_result: dict,
    coolant_enters_at: str,
    lut_path: str | Path | None = None,
    min_pressure_Pa: float = 1.0,
) -> HXGridHeatedChannelResult:
    """Run the shell-side liquid adapter from ``shellntube_solver.tube`` output."""
    if "dQ" not in tube_result:
        raise KeyError("tube_result must contain a 'dQ' array")
    return solve_shelltube_shellside_liquid_from_duty(
        coolant_prop=coolant_prop,
        shelltube_prop=shelltube_prop,
        dQ_profile_per_tube_W=np.asarray(tube_result["dQ"], dtype=float),
        coolant_enters_at=coolant_enters_at,
        lut_path=lut_path,
        min_pressure_Pa=min_pressure_Pa,
    )


def solve_helical_coil_liquid_from_duty(
    *,
    coolant_prop: coolantProp,
    combustor_prop: combustorProp,
    numerical_prop: numericalProp,
    dQ_profile_W: np.ndarray,
    coolant_enters_at: str | None = None,
    z_edges_m: np.ndarray | None = None,
    lut_path: str | Path | None = None,
    min_pressure_Pa: float = 1.0,
) -> HXGridHeatedChannelResult:
    """Run the liquid profile solver from helical-coil wall duty.

    ``dQ_profile_W`` is the total heat transferred to the coolant in each coil
    path segment, matching the maintained helical solver's ``data_master["dQ"]``
    convention. Geometry is the coil inner flow path: ``Dh_coil`` as hydraulic
    diameter, ``N_coils * pi*Dh^2/4`` as total flow area, and
    ``N_coils * pi*Dh`` as total heated perimeter.
    """
    dQ = np.asarray(dQ_profile_W, dtype=float)
    if dQ.ndim != 1 or dQ.size < 1:
        raise ValueError("dQ_profile_W must be a non-empty 1D array")
    if not np.all(np.isfinite(dQ)):
        raise ValueError("dQ_profile_W must be finite")
    if coolant_enters_at is None:
        coolant_enters_at = "z_min" if combustor_prop.flow_config == "co" else "z_max"
    if coolant_enters_at not in ("z_min", "z_max"):
        raise ValueError("coolant_enters_at must be 'z_min' or 'z_max'")

    if z_edges_m is None:
        z_edges = np.linspace(0.0, _helical_pipe_length(combustor_prop, numerical_prop), dQ.size + 1)
    else:
        z_edges = np.asarray(z_edges_m, dtype=float)
        if z_edges.shape != (dQ.size + 1,):
            raise ValueError("z_edges_m must have length len(dQ_profile_W) + 1")

    n_channels = int(combustor_prop.N_coils)
    if n_channels < 1:
        raise ValueError("N_coils must be at least 1")
    Dh = float(combustor_prop.Dh_coil)
    if Dh <= 0.0:
        raise ValueError("Dh_coil must be positive")
    total_area = n_channels * np.pi * Dh**2 / 4.0
    total_perimeter = n_channels * np.pi * Dh

    return solve_steady_heated_channel_on_hx_grid(
        coolant_prop=coolant_prop,
        z_edges_m=z_edges,
        hydraulic_diameter_m=Dh,
        flow_area_m2=total_area,
        heated_perimeter_m=total_perimeter,
        mass_flow_kg_s=coolant_prop.mass_flow_c,
        p_in_Pa=coolant_prop.p_in,
        T_in_K=coolant_prop.T_in,
        heat_per_segment_W=dQ,
        coolant_enters_at=coolant_enters_at,
        lut_path=lut_path,
        min_pressure_Pa=min_pressure_Pa,
    )


def solve_helical_coil_liquid_from_data_master(
    *,
    coolant_prop: coolantProp,
    combustor_prop: combustorProp,
    numerical_prop: numericalProp,
    data_master: dict,
    lut_path: str | Path | None = None,
    min_pressure_Pa: float = 1.0,
) -> HXGridHeatedChannelResult:
    """Run the helical liquid adapter from ``main_solver.data_master`` output."""
    if "dQ" not in data_master:
        raise KeyError("data_master must contain a 'dQ' array")
    dQ = np.asarray(data_master["dQ"], dtype=float)
    z_edges = None
    if "L_ch" in data_master and len(data_master["L_ch"]) == dQ.size:
        L_ch_nodes = np.asarray(data_master["L_ch"], dtype=float)
        if L_ch_nodes.size > 1 and np.all(np.isfinite(L_ch_nodes)):
            step = float(np.median(np.diff(L_ch_nodes)))
            if step > 0.0:
                z_edges = np.concatenate([L_ch_nodes, [L_ch_nodes[-1] + step]])
    return solve_helical_coil_liquid_from_duty(
        coolant_prop=coolant_prop,
        combustor_prop=combustor_prop,
        numerical_prop=numerical_prop,
        dQ_profile_W=dQ,
        z_edges_m=z_edges,
        lut_path=lut_path,
        min_pressure_Pa=min_pressure_Pa,
    )

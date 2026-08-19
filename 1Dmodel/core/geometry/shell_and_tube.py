"""``shellntube`` geometry -> FlowPath pair (Stage B).

Hot combustion gas is INSIDE the tubes; coolant is on the baffled shell side.
That orientation is load-bearing downstream: wall conduction must be solved
with ``hot_side="inner"`` (CLAUDE.md -- do not regress the corrected hot/cold
perimeter mapping in ``physics/heat_conduction.py``), and tube-side gas
quantities are PER TUBE. Here the per-tube convention is structural:
``FlowPath.n_parallel = N_tubes`` and ``FlowPath.mass_flux()`` divides by it
exactly once.

Bell-Delaware shell-side flow areas stay in
``mechanical/geometry/shelltube_geometry.py`` -- this builder reuses that
helper rather than duplicating it, and stores the resulting cross-flow area
as the shell FlowPath's ``A_flow``.
"""

from __future__ import annotations

import numpy as np

from hps_combustor.core.mesh import FlowPath, HXAssembly
from hps_combustor.mechanical.geometry.shelltube_geometry import (
    compute_bell_delaware_geometry,
)


def build_shelltube_assembly(shell_tube_prop, combustor_prop, *, n_cells: int = 200):
    """Build the (hot tube-side, cold shell-side) FlowPath pair.

    ``n_cells`` defaults to the production steady grid
    (``runProp.shelltube_steady_nodes``). Both streams share one uniform
    axial partition here, so the coupling operator is the identity -- but it
    is still built through the same conservative machinery so a future
    non-uniform or refined tube-side grid needs no special-casing.
    """
    D_o = shell_tube_prop.D_tube_outer
    t_w = shell_tube_prop.thickness_tube_wall
    D_i = D_o - 2.0 * t_w
    if D_i <= 0.0:
        raise ValueError("tube inner diameter is non-positive - check thickness_tube_wall")
    L = shell_tube_prop.L_tube

    z_edges = np.linspace(0.0, L, n_cells + 1)

    tube = FlowPath(
        name="tube_gas",
        s_edges=z_edges.copy(),
        z_of_s_edges=z_edges.copy(),
        A_flow=np.pi * D_i**2 / 4.0,
        Dh=D_i,
        P_wetted=np.pi * D_i,
        P_heated=np.pi * D_i,
        n_parallel=shell_tube_prop.N_tubes,
        geometry_tag="straight_tube",
        roughness=shell_tube_prop.tube_roughness,
        flow_direction=1,
    )

    geom = compute_bell_delaware_geometry(
        D_shell_inner=shell_tube_prop.D_shell_inner,
        D_tube_outer=D_o,
        pitch_ratio=shell_tube_prop.pitch_ratio,
        layout=shell_tube_prop.layout,
        N_tubes=shell_tube_prop.N_tubes,
        N_baffles=shell_tube_prop.N_baffles,
        baffle_cut=shell_tube_prop.baffle_cut,
        L_tube=L,
        clearance_tube_baffle=shell_tube_prop.clearance_tube_baffle,
        clearance_baffle_shell=shell_tube_prop.clearance_baffle_shell,
        clearance_bundle_shell=shell_tube_prop.clearance_bundle_shell,
        N_sealing_strip_pairs=shell_tube_prop.N_sealing_strip_pairs,
        baffle_spacing=shell_tube_prop.baffle_spacing,
        L_inlet_spacing=shell_tube_prop.L_inlet_spacing,
        L_outlet_spacing=shell_tube_prop.L_outlet_spacing,
    )

    shell = FlowPath(
        name="shell_coolant",
        s_edges=z_edges.copy(),
        z_of_s_edges=z_edges.copy(),
        A_flow=geom["S_m"],           # Bell-Delaware centreline cross-flow area
        Dh=D_o,                        # cross-flow reference length = tube OD
        P_wetted=np.pi * D_o,
        P_heated=np.pi * D_o,          # heat leaves via the tube OD
        n_parallel=1,
        geometry_tag="shell_crossflow",
        flow_direction=1 if combustor_prop.flow_config == "co" else -1,
    )
    return HXAssembly(hot=tube, cold=shell)

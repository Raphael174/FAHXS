"""``shellnHelicalTube`` geometry -> FlowPath pair (Stage B).

Reproduces the geometry derivations currently inlined in
``main_solve.py.__init__`` (lines ~110-148) so the FV core consumes geometry
as DATA rather than recomputing it mid-march. In particular the coil's
arc-length-to-axial map, which the legacy ``_advance_state()`` approximated
linearly and got wrong for any config other than ``shellnHelicalTube``
(CLAUDE.md, 2026-07-13) -- here it is an explicit ``z_of_s_edges`` array
built once from the same ``HelixGeometryRadiusCST`` helper the legacy solver
uses, so there is one source of truth.

Not yet consumed by ``main_solve.py`` -- Stage B builds and validates the
builders; the legacy solver is repointed at them during its own migration
(Stage E), per the staged-migration decision.
"""

from __future__ import annotations

import numpy as np

from hps_combustor.core.mesh import FlowPath, HXAssembly
from hps_combustor.mechanical.geometry.helix_geometry import (
    HelixGeometryRadiusCST,
    compute_Dh_shell,
)


def build_helical_assembly(combustor_prop, numerical_prop, *, n_cells_shell: int | None = None):
    """Build the (hot shell-side, cold coil-side) FlowPath pair.

    Mirrors ``main_solve.main_solver.__init__``'s geometry block. Returns an
    ``HXAssembly``; ``flow_config`` on ``combustor_prop`` sets the coil's
    ``flow_direction`` (hot gas is always +z by convention here).
    """
    Dh_ch = combustor_prop.Dh_coil
    t_wall = combustor_prop.thickness_coil_wall
    coil_pitch = Dh_ch + 2.0 * t_wall + combustor_prop.coil_gap
    D_coil = (
        combustor_prop.inner_diameter
        - 2.0 * combustor_prop.gap_shell2coil
        - Dh_ch
        - 2.0 * t_wall
    )
    if D_coil <= 0.0:
        raise ValueError("Coil center-to-center diameter is negative - check your geometry")

    L_coil = (
        (numerical_prop.L_HX_max - combustor_prop.mixing_length)
        - 2.0 * combustor_prop.length_2_coil
        - (Dh_ch + 2.0 * t_wall)
    )
    func_s_to_x, _func_s_to_theta, L_ch_max = HelixGeometryRadiusCST(
        coil_pitch=coil_pitch, D_coil=D_coil, L_coil=L_coil
    )

    # Arc-length grid, same step convention as numericalProp.dx in the legacy
    # solver: one full turn split into N_arc_steps_per_turn sub-steps.
    ds_nominal = np.pi * D_coil / numerical_prop.N_arc_steps_per_turn
    n_cells_coil = max(int(round(L_ch_max / ds_nominal)), 1)
    s_edges = np.linspace(0.0, L_ch_max, n_cells_coil + 1)
    # THE map that matters: arc length -> axial station, from the same helper
    # the legacy solver uses (not a linear-dx approximation).
    z_edges_coil = np.asarray(func_s_to_x(s_edges), dtype=float)

    # Curvature radius including pitch (Rc in main_solve.py line ~118).
    Rc = D_coil / 2.0 * (1.0 + (coil_pitch / (np.pi * D_coil)) ** 2)

    coil = FlowPath(
        name="coil",
        s_edges=s_edges,
        z_of_s_edges=z_edges_coil,
        A_flow=np.pi * Dh_ch**2 / 4.0,
        Dh=Dh_ch,
        P_wetted=np.pi * Dh_ch,
        P_heated=np.pi * Dh_ch,
        n_parallel=combustor_prop.N_coils,
        geometry_tag="helical_coil",
        R_curv=np.full(n_cells_coil, Rc),
        roughness=combustor_prop.channel_roughness,
        flow_direction=1 if combustor_prop.flow_config == "co" else -1,
    )

    # Shell side: annular gas passage around the coil, uniform axial grid over
    # the same z span the coil covers.
    Dh_cc = compute_Dh_shell(
        D_coil=D_coil,
        d_coil_outer=Dh_ch + 2.0 * t_wall,
        shell_diameter=combustor_prop.inner_diameter,
        coil_pitch=coil_pitch,
    )
    n_shell = int(n_cells_shell) if n_cells_shell else max(n_cells_coil // 10, 2)
    z_edges_shell = np.linspace(z_edges_coil[0], z_edges_coil[-1], n_shell + 1)
    shell = FlowPath(
        name="shell",
        s_edges=z_edges_shell.copy(),
        z_of_s_edges=z_edges_shell.copy(),
        A_flow=np.pi * Dh_cc**2 / 4.0,
        Dh=Dh_cc,
        P_wetted=np.pi * Dh_cc,
        P_heated=np.pi * (Dh_ch + 2.0 * t_wall),  # heat enters via the coil OD
        n_parallel=1,
        geometry_tag="shell_crossflow",
        roughness=combustor_prop.combustor_roughness,
        flow_direction=1,
    )
    return HXAssembly(hot=shell, cold=coil)

"""Stage B acceptance: geometry builders reproduce the legacy solvers' own
geometry derivations to machine precision.

The builders in ``core/geometry/`` re-derive what ``main_solve.py`` and
``main_solve_shellntube.py`` currently compute inline. Until those solvers
are repointed at the builders (their own migration, Stage E), these tests are
what guarantees the two derivations have not silently diverged -- they
recompute the legacy expressions here, from the same shared helpers the
legacy files import, and demand exact agreement.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from hps_combustor.core.geometry.shell_and_helical_tube import build_helical_assembly
from hps_combustor.core.geometry.shell_and_tube import build_shelltube_assembly
from hps_combustor.input_data import combustorProp, numericalProp, shellTubeProp
from hps_combustor.mechanical.geometry.helix_geometry import (
    HelixGeometryRadiusCST,
    compute_Dh_shell,
)


@pytest.fixture(scope="module")
def helical_reference():
    """Legacy inline derivation from main_solve.main_solver.__init__."""
    cp = combustorProp(HX_config="shellnHelicalTube", flow_config="co")
    npr = numericalProp()
    Dh_ch = cp.Dh_coil
    t = cp.thickness_coil_wall
    coil_pitch = Dh_ch + 2.0 * t + cp.coil_gap
    D_coil = cp.inner_diameter - 2.0 * cp.gap_shell2coil - cp.Dh_coil - 2.0 * t
    L_coil = (
        (npr.L_HX_max - cp.mixing_length) - 2.0 * cp.length_2_coil - (Dh_ch + 2.0 * t)
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, _, L_ch_max = HelixGeometryRadiusCST(
            coil_pitch=coil_pitch, D_coil=D_coil, L_coil=L_coil
        )
    return dict(
        cp=cp,
        npr=npr,
        L_ch_max=L_ch_max,
        Dh_shell=compute_Dh_shell(
            D_coil=D_coil,
            d_coil_outer=Dh_ch + 2.0 * t,
            shell_diameter=cp.inner_diameter,
            coil_pitch=coil_pitch,
        ),
        Rc=D_coil / 2.0 * (1.0 + (coil_pitch / (np.pi * D_coil)) ** 2),
        A_ch=np.pi * cp.Dh_coil**2 / 4.0,
    )


@pytest.fixture(scope="module")
def helical_assembly(helical_reference):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_helical_assembly(helical_reference["cp"], helical_reference["npr"])


def test_helical_arc_length_matches_legacy(helical_reference, helical_assembly):
    assert helical_assembly.cold.length_s == pytest.approx(
        helical_reference["L_ch_max"], abs=1e-12
    )


def test_helical_shell_hydraulic_diameter_matches_legacy(helical_reference, helical_assembly):
    assert helical_assembly.hot.Dh[0] == pytest.approx(
        helical_reference["Dh_shell"], abs=1e-15
    )


def test_helical_curvature_radius_matches_legacy(helical_reference, helical_assembly):
    assert helical_assembly.cold.R_curv[0] == pytest.approx(helical_reference["Rc"], abs=1e-15)


def test_helical_flow_area_matches_legacy(helical_reference, helical_assembly):
    assert helical_assembly.cold.A_flow[0] == pytest.approx(helical_reference["A_ch"], rel=1e-15)


def test_helical_node_count_matches_documented_real_geometry(helical_assembly):
    """CLAUDE.md records that this combustor's real coil is ~1378 arc-length
    nodes, NOT ~100 -- the discrepancy that silently broke every Phase 0-2
    liquid-coolant test before the HX_config guard existed. Pin it."""
    assert helical_assembly.cold.n_cells == 1378


def test_helical_arc_length_far_exceeds_axial_length(helical_assembly):
    """s != z is the whole reason FlowPath separates them. ~14.6x here."""
    coil = helical_assembly.cold
    assert coil.length_s / coil.length_z > 10.0


def test_helical_coupling_is_conservative(helical_assembly):
    assert helical_assembly.hot_to_cold.conservation_defect < 1e-12
    assert helical_assembly.cold_to_hot.conservation_defect < 1e-12
    q = np.linspace(1.0, 100.0, helical_assembly.cold.n_cells)
    assert helical_assembly.cold_to_hot.apply(q).sum() == pytest.approx(q.sum(), rel=1e-12)


def test_helical_flow_direction_follows_flow_config():
    npr = numericalProp()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        co = build_helical_assembly(combustorProp(flow_config="co"), npr)
        ct = build_helical_assembly(combustorProp(flow_config="counter"), npr)
    assert not co.is_counterflow
    assert ct.is_counterflow


# --- shell-and-tube ------------------------------------------------------


def test_shelltube_tube_inner_diameter_and_per_tube_mass_flux():
    stp = shellTubeProp()
    asm = build_shelltube_assembly(stp, combustorProp(flow_config="co"), n_cells=200)
    tube = asm.hot
    D_i = stp.D_tube_outer - 2.0 * stp.thickness_tube_wall
    assert tube.Dh[0] == pytest.approx(D_i, rel=1e-15)
    assert tube.n_parallel == stp.N_tubes
    # The "per tube" sharp edge (CLAUDE.md): total mdot divided by N_tubes once.
    mdot_total = 0.1
    expected = (mdot_total / stp.N_tubes) / (np.pi / 4.0 * D_i**2)
    assert tube.mass_flux(mdot_total)[0] == pytest.approx(expected, rel=1e-15)


def test_shelltube_shell_area_is_bell_delaware_crossflow_area():
    from hps_combustor.mechanical.geometry.shelltube_geometry import (
        compute_bell_delaware_geometry,
    )

    stp = shellTubeProp()
    asm = build_shelltube_assembly(stp, combustorProp(flow_config="co"), n_cells=50)
    geom = compute_bell_delaware_geometry(
        D_shell_inner=stp.D_shell_inner, D_tube_outer=stp.D_tube_outer,
        pitch_ratio=stp.pitch_ratio, layout=stp.layout, N_tubes=stp.N_tubes,
        N_baffles=stp.N_baffles, baffle_cut=stp.baffle_cut, L_tube=stp.L_tube,
        clearance_tube_baffle=stp.clearance_tube_baffle,
        clearance_baffle_shell=stp.clearance_baffle_shell,
        clearance_bundle_shell=stp.clearance_bundle_shell,
        N_sealing_strip_pairs=stp.N_sealing_strip_pairs,
        baffle_spacing=stp.baffle_spacing, L_inlet_spacing=stp.L_inlet_spacing,
        L_outlet_spacing=stp.L_outlet_spacing,
    )
    assert asm.cold.A_flow[0] == pytest.approx(geom["S_m"], rel=1e-15)


def test_shelltube_matched_grids_give_identity_coupling():
    stp = shellTubeProp()
    asm = build_shelltube_assembly(stp, combustorProp(flow_config="co"), n_cells=25)
    np.testing.assert_allclose(asm.hot_to_cold.weights, np.eye(25), atol=1e-14)


def test_shelltube_rejects_impossible_wall_thickness():
    stp = shellTubeProp(D_tube_outer=5e-3, thickness_tube_wall=3e-3)
    with pytest.raises(ValueError, match="inner diameter"):
        build_shelltube_assembly(stp, combustorProp(flow_config="co"))

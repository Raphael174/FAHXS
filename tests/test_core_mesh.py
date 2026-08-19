"""Stage B of docs/solver_design/FV_CORE_REWORK_PLAN.md.

The acceptance gate for this stage is explicitly: "a unit test proving the
overlap operator conserves energy exactly on a non-uniform helical<->shell
mapping". That is ``test_coupling_conserves_energy_helical_to_shell`` and
``test_coupling_conserves_energy_nonuniform_both_sides`` below — everything
else here is supporting validation of FlowPath's invariants.
"""
from __future__ import annotations

import numpy as np
import pytest

from hps_combustor.core.mesh import (
    FlowPath,
    HXAssembly,
    build_coupling,
)


def _straight_path(name, length, n_cells, *, flow_direction=1, n_parallel=1, D=0.01):
    z = np.linspace(0.0, length, n_cells + 1)
    return FlowPath(
        name=name,
        s_edges=z.copy(),
        z_of_s_edges=z.copy(),
        A_flow=np.pi / 4.0 * D**2,
        Dh=D,
        P_wetted=np.pi * D,
        P_heated=np.pi * D,
        n_parallel=n_parallel,
        flow_direction=flow_direction,
    )


def _helical_path(name, z_length, n_turns, n_per_turn, D_coil, D_tube, *, flow_direction=1):
    """Helical coil: arc length s advances ~pi*D_coil per turn while z
    advances only one pitch — the real geometry that makes s != z matter."""
    n_cells = n_turns * n_per_turn
    z = np.linspace(0.0, z_length, n_cells + 1)
    turn_arc = np.pi * D_coil
    s = np.linspace(0.0, n_turns * turn_arc, n_cells + 1)
    return FlowPath(
        name=name,
        s_edges=s,
        z_of_s_edges=z,
        A_flow=np.pi / 4.0 * D_tube**2,
        Dh=D_tube,
        P_wetted=np.pi * D_tube,
        P_heated=np.pi * D_tube,
        geometry_tag="helical_coil",
        R_curv=np.full(n_cells, D_coil / 2.0),
        flow_direction=flow_direction,
    )


# --- FlowPath invariants -------------------------------------------------


def test_flowpath_basic_measures():
    p = _straight_path("hot", 1.0, 10)
    assert p.n_cells == 10
    assert p.length_s == pytest.approx(1.0)
    assert p.length_z == pytest.approx(1.0)
    np.testing.assert_allclose(p.ds, 0.1)
    np.testing.assert_allclose(p.dz, 0.1)


def test_flowpath_helical_s_exceeds_z():
    """The defining helical property: arc length is much longer than axial
    extent. This is what the legacy linear-dx approximation got wrong."""
    p = _helical_path("coil", z_length=0.2, n_turns=10, n_per_turn=20,
                      D_coil=0.136, D_tube=0.0135)
    assert p.length_z == pytest.approx(0.2)
    assert p.length_s == pytest.approx(10 * np.pi * 0.136)
    assert p.length_s > 20 * p.length_z


def test_flowpath_rejects_nonmonotonic_z():
    z = np.array([0.0, 0.5, 0.4, 1.0])
    with pytest.raises(ValueError, match="strictly increasing"):
        FlowPath(name="bad", s_edges=np.array([0.0, 1.0, 2.0, 3.0]), z_of_s_edges=z,
                 A_flow=1e-4, Dh=0.01, P_wetted=0.03, P_heated=0.03)


def test_flowpath_rejects_nonpositive_area():
    z = np.linspace(0.0, 1.0, 4)
    with pytest.raises(ValueError, match="strictly positive"):
        FlowPath(name="bad", s_edges=z.copy(), z_of_s_edges=z.copy(),
                 A_flow=np.array([1e-4, 0.0, 1e-4]), Dh=0.01,
                 P_wetted=0.03, P_heated=0.03)


def test_mass_flux_divides_by_n_parallel_once():
    """The shell-and-tube 'per tube' sharp edge, handled structurally."""
    p = _straight_path("tubes", 0.235, 10, n_parallel=235, D=5e-3)
    G = p.mass_flux(1.0)  # 1 kg/s total across 235 tubes
    expected = (1.0 / 235) / (np.pi / 4.0 * 5e-3**2)
    np.testing.assert_allclose(G, expected)


def test_inlet_outlet_indices_follow_flow_direction():
    fwd = _straight_path("a", 1.0, 5, flow_direction=1)
    rev = _straight_path("b", 1.0, 5, flow_direction=-1)
    assert (fwd.inlet_index, fwd.outlet_index) == (0, 4)
    assert (rev.inlet_index, rev.outlet_index) == (4, 0)


# --- Conservation: the Stage B acceptance gate ---------------------------


def test_coupling_conserves_energy_helical_to_shell():
    """THE Stage B gate: non-uniform helical<->shell mapping, exact
    conservation. 200 coil cells over the same 0.2 m axial span as 37 shell
    cells — deliberately non-commensurate (200/37 is not an integer) so no
    cell boundary coincidence can hide a defect."""
    coil = _helical_path("coil", z_length=0.2, n_turns=10, n_per_turn=20,
                         D_coil=0.136, D_tube=0.0135)
    shell = _straight_path("shell", 0.2, 37, D=0.05)

    coupling = build_coupling(coil, shell)
    assert coupling.conservation_defect == pytest.approx(0.0, abs=1e-14)

    rng = np.random.default_rng(20260731)
    q_coil = rng.uniform(1.0, 500.0, size=coil.n_cells)  # W per cell
    q_shell = coupling.apply(q_coil)

    assert q_shell.shape == (shell.n_cells,)
    assert q_shell.sum() == pytest.approx(q_coil.sum(), rel=1e-13)


def test_coupling_conserves_energy_nonuniform_both_sides():
    """Both partitions non-uniform (clustered cells, as a real grid-refined
    solver would produce) — still exact."""
    rng = np.random.default_rng(7)
    z_a = np.sort(np.concatenate([[0.0, 1.0], rng.uniform(0.0, 1.0, 40)]))
    z_b = np.sort(np.concatenate([[0.0, 1.0], rng.uniform(0.0, 1.0, 25)]))
    z_a, z_b = np.unique(z_a), np.unique(z_b)

    a = FlowPath(name="a", s_edges=z_a.copy(), z_of_s_edges=z_a.copy(),
                 A_flow=1e-4, Dh=0.01, P_wetted=0.03, P_heated=0.03)
    b = FlowPath(name="b", s_edges=z_b.copy(), z_of_s_edges=z_b.copy(),
                 A_flow=1e-4, Dh=0.01, P_wetted=0.03, P_heated=0.03)

    coupling = build_coupling(a, b)
    assert coupling.conservation_defect == pytest.approx(0.0, abs=1e-13)

    q_a = rng.uniform(-100.0, 900.0, size=a.n_cells)
    assert coupling.apply(q_a).sum() == pytest.approx(q_a.sum(), rel=1e-12)


def test_coupling_round_trip_preserves_total():
    """Mapping A->B then B->A conserves the grand total (not the per-cell
    distribution — that is genuinely lossy, and intentionally so)."""
    coil = _helical_path("coil", 0.2, 8, 15, 0.136, 0.0135)
    shell = _straight_path("shell", 0.2, 23, D=0.05)
    fwd = build_coupling(coil, shell)
    back = build_coupling(shell, coil)

    q = np.linspace(10.0, 300.0, coil.n_cells)
    assert back.apply(fwd.apply(q)).sum() == pytest.approx(q.sum(), rel=1e-12)


def test_coupling_identity_when_partitions_match():
    a = _straight_path("a", 1.0, 12)
    b = _straight_path("b", 1.0, 12)
    coupling = build_coupling(a, b)
    np.testing.assert_allclose(coupling.weights, np.eye(12), atol=1e-14)


def test_coupling_rejects_mismatched_spans():
    a = _straight_path("a", 1.0, 10)
    b = _straight_path("b", 0.8, 10)
    with pytest.raises(ValueError, match="axial spans differ"):
        build_coupling(a, b)


def test_coupling_apply_rejects_wrong_shape():
    a = _straight_path("a", 1.0, 10)
    b = _straight_path("b", 1.0, 6)
    coupling = build_coupling(a, b)
    with pytest.raises(ValueError, match="expected source array"):
        coupling.apply(np.ones(5))


def test_uniform_source_maps_to_proportional_target():
    """A uniform per-LENGTH source must land proportional to target cell
    length — the physical sanity check behind the weights."""
    a = _straight_path("a", 1.0, 100)
    z_b = np.array([0.0, 0.1, 0.5, 1.0])
    b = FlowPath(name="b", s_edges=z_b.copy(), z_of_s_edges=z_b.copy(),
                 A_flow=1e-4, Dh=0.01, P_wetted=0.03, P_heated=0.03)
    coupling = build_coupling(a, b)
    q_a = np.full(a.n_cells, 1.0)  # 1 W per cell, uniform cells => uniform per length
    q_b = coupling.apply(q_a)
    np.testing.assert_allclose(q_b, np.array([10.0, 40.0, 50.0]), rtol=1e-12)


# --- Assembly ------------------------------------------------------------


def test_assembly_detects_counterflow():
    hot = _straight_path("hot", 1.0, 10, flow_direction=1)
    cold_co = _straight_path("cold", 1.0, 10, flow_direction=1)
    cold_ct = _straight_path("cold", 1.0, 10, flow_direction=-1)
    assert not HXAssembly(hot=hot, cold=cold_co).is_counterflow
    assert HXAssembly(hot=hot, cold=cold_ct).is_counterflow


def test_assembly_couplings_both_directions_conserve():
    hot = _straight_path("hot", 0.235, 40, n_parallel=235, D=5e-3)
    cold = _straight_path("cold", 0.235, 17, D=0.05, flow_direction=-1)
    asm = HXAssembly(hot=hot, cold=cold)
    q_hot = np.linspace(5.0, 50.0, hot.n_cells)
    assert asm.hot_to_cold.apply(q_hot).sum() == pytest.approx(q_hot.sum(), rel=1e-12)
    q_cold = np.linspace(1.0, 9.0, cold.n_cells)
    assert asm.cold_to_hot.apply(q_cold).sum() == pytest.approx(q_cold.sum(), rel=1e-12)

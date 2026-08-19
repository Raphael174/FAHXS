"""Stage D, Slice 2 of docs/solver_design/FV_CORE_REWORK_PLAN.md.

`core/state.py` generalizes `transient_core/state.py::TransientStateLayout`
(which only packs `[Tbar_wall, T_coolant]`, the OLD temperature-only coolant
model) to the conservative `(mass, internal_energy)` pair the project
actually runs today. See that module's docstring for the full rationale.
"""
from __future__ import annotations

import numpy as np
import pytest

from hps_combustor.core.state import CoolantState, WallCoolantStateLayout, WallState


def test_coolant_state_rejects_mismatched_lengths():
    with pytest.raises(ValueError, match="equal length"):
        CoolantState(mass_kg=np.array([1.0, 2.0]), internal_energy_J=np.array([1.0]))


def test_coolant_state_rejects_nonpositive_mass():
    with pytest.raises(ValueError, match="strictly positive"):
        CoolantState(mass_kg=np.array([1.0, 0.0]), internal_energy_J=np.array([1.0, 1.0]))


def test_coolant_state_specific_internal_energy():
    state = CoolantState(mass_kg=np.array([2.0, 4.0]), internal_energy_J=np.array([10.0, 40.0]))
    np.testing.assert_allclose(state.specific_internal_energy_J_kg, [5.0, 10.0])
    assert state.n_cells == 2


def test_wall_state_rejects_nonfinite():
    with pytest.raises(ValueError, match="non-finite"):
        WallState(Tbar_K=np.array([300.0, np.nan]))


def test_layout_pack_unpack_round_trip():
    layout = WallCoolantStateLayout(n_cells=3)
    wall = WallState(Tbar_K=np.array([300.0, 310.0, 320.0]))
    coolant = CoolantState(
        mass_kg=np.array([1e-3, 2e-3, 3e-3]),
        internal_energy_J=np.array([100.0, 200.0, 300.0]),
    )
    y = layout.pack(wall, coolant)
    assert y.shape == (9,)

    wall_out, coolant_out = layout.unpack(y)
    np.testing.assert_allclose(wall_out.Tbar_K, wall.Tbar_K)
    np.testing.assert_allclose(coolant_out.mass_kg, coolant.mass_kg)
    np.testing.assert_allclose(coolant_out.internal_energy_J, coolant.internal_energy_J)


def test_layout_pack_rejects_wrong_shapes():
    layout = WallCoolantStateLayout(n_cells=3)
    wall = WallState(Tbar_K=np.array([300.0, 310.0]))  # wrong length
    coolant = CoolantState(
        mass_kg=np.array([1e-3, 2e-3, 3e-3]),
        internal_energy_J=np.array([100.0, 200.0, 300.0]),
    )
    with pytest.raises(ValueError, match="wall.Tbar_K"):
        layout.pack(wall, coolant)


def test_layout_unpack_rejects_wrong_size():
    layout = WallCoolantStateLayout(n_cells=3)
    with pytest.raises(ValueError, match="shape"):
        layout.unpack(np.zeros(7))


def test_layout_size_and_slices():
    layout = WallCoolantStateLayout(n_cells=4)
    assert layout.size == 12
    assert layout.wall_slice == slice(0, 4)
    assert layout.mass_slice == slice(4, 8)
    assert layout.energy_slice == slice(8, 12)

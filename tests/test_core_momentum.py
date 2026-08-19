"""Stage D, Slice 2 of docs/solver_design/FV_CORE_REWORK_PLAN.md.

`core/momentum.py::quasi_steady_face_mdot` generalizes
`transient_core/adapters_shelltube.py::_shelltube_quasi_steady_faces` and
`main_solve_transient.py::_helical_quasi_steady_faces`, parameterized by
`FlowPath` instead of a raw `flow_direction` int. Proves bit-identical
equivalence to the shell-and-tube version specifically (the two legacy
functions diverge in the closed-valve branch -- helical additionally
computes non-zero interior-face flows there; that divergence is
deliberately not reproduced, see core/momentum.py's module docstring and
`docs/solver_design/FV_CORE_REWORK_PLAN.md`).
"""
from __future__ import annotations

import numpy as np
import pytest

from hps_combustor.core.mesh import FlowPath
from hps_combustor.core.momentum import quasi_steady_face_mdot
from hps_combustor.transient_core.adapters_shelltube import _shelltube_quasi_steady_faces

P = np.array([5.0e5, 4.9e5, 4.8e5, 4.7e5])
RHO = np.array([10.0, 10.0, 10.0, 10.0])
RESISTANCE = np.array([1.0e5, 1.0e5, 1.0e5, 1.0e5])


def _path(flow_direction: int) -> FlowPath:
    s = np.linspace(0.0, 1.0, 5)
    return FlowPath(
        name="test",
        s_edges=s,
        z_of_s_edges=s.copy(),
        A_flow=1e-4,
        Dh=0.01,
        P_wetted=0.03,
        P_heated=0.03,
        flow_direction=flow_direction,
    )


@pytest.mark.parametrize("flow_direction", [1, -1])
def test_throughput_branch_matches_legacy_shelltube(flow_direction):
    """mdot_inlet > mdot_floor: both legacy functions agree here (uniform
    prescribed flow), so this is the least ambiguous equivalence check."""
    path = _path(flow_direction)
    got = quasi_steady_face_mdot(
        path, P, RHO, RESISTANCE, mdot_inlet=0.05, outlet_pressure=4.5e5, mdot_floor=1e-9
    )
    ref = _shelltube_quasi_steady_faces(
        P, RHO, RESISTANCE, mdot_inlet=0.05, outlet_pressure=4.5e5,
        flow_direction=flow_direction, mdot_floor=1e-9,
    )
    np.testing.assert_array_equal(got, ref)
    expected_value = 0.05 if flow_direction == 1 else -0.05
    np.testing.assert_allclose(got, expected_value)


@pytest.mark.parametrize(
    "flow_direction,outlet_pressure",
    [(1, 4.5e5), (-1, 4.5e5)],
)
def test_closed_valve_branch_matches_legacy_shelltube(flow_direction, outlet_pressure):
    """mdot_inlet <= mdot_floor: only the downstream boundary face carries
    pressure-driven discharge; interior faces stay zero. Pressures chosen so
    the discharge is genuinely non-zero (not a trivial both-zero match)."""
    path = _path(flow_direction)
    got = quasi_steady_face_mdot(
        path, P, RHO, RESISTANCE, mdot_inlet=0.0, outlet_pressure=outlet_pressure, mdot_floor=1e-9
    )
    ref = _shelltube_quasi_steady_faces(
        P, RHO, RESISTANCE, mdot_inlet=0.0, outlet_pressure=outlet_pressure,
        flow_direction=flow_direction, mdot_floor=1e-9,
    )
    np.testing.assert_array_equal(got, ref)
    # interior faces (all but the one active boundary face) must be exactly zero
    interior = got[1:-1] if flow_direction == 1 else got[1:-1]
    assert np.count_nonzero(got) == 1, f"expected exactly one nonzero face, got {got}"


def test_rejects_wrong_length_inputs():
    path = _path(1)
    with pytest.raises(ValueError, match="pressure and density"):
        quasi_steady_face_mdot(
            path, P[:-1], RHO, RESISTANCE, mdot_inlet=0.0, outlet_pressure=4.5e5, mdot_floor=1e-9
        )
    with pytest.raises(ValueError, match="resistance"):
        quasi_steady_face_mdot(
            path, P, RHO, RESISTANCE[:-1], mdot_inlet=0.0, outlet_pressure=4.5e5, mdot_floor=1e-9
        )

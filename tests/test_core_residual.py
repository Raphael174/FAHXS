"""Stage D, Slice 3 of docs/solver_design/FV_CORE_REWORK_PLAN.md.

`core/residual.py::hot_gas_march` reproduces `transient_core/
adapters_shelltube.py::shelltube_hot_gas_march`'s sequential per-cell
algorithm on the new FlowPath/core.closures/core.wall abstractions. This is
genuinely new code (not a relocation like the other Slice 3 pieces), so it
is gated by exact-match reproduction against the legacy march on a
hand-built fixture, for both `inside_tube_choice` values and with nonzero
wall-side radiation and calibration coefficients exercised, rather than an
identical-object proof.
"""
from __future__ import annotations

import numpy as np
import pytest

from hps_combustor.core.hotgas.combustor import equilibrium_gas_state_provider
from hps_combustor.core.mesh import FlowPath
from hps_combustor.core.residual import hot_gas_march
from hps_combustor.core.wall import CylindricalWall
from hps_combustor.input_data import CorrelationCoefficients
from hps_combustor.transient_core.adapters_shelltube import (
    ShellTubeCoreGeometry,
    shelltube_hot_gas_march,
)
from hps_combustor.transient_core.grid import AxialGrid

N = 4
L = 0.5
D_I = 0.01
T_WALL_THICKNESS = 0.0008
N_TUBES = 10
ROUGHNESS = 1.5e-6


class _FakeEquilibriumManifold:
    """Deterministic, cheap stand-in for a real equilibrium manifold -- gas
    cools linearly with enthalpy removed, floored at 400 K."""

    def at(self, h_removed):
        T = max(2500.0 - h_removed / 2000.0, 400.0)
        return (T, 3.0, 5.0e-5, 0.15, 1400.0, 0.1, 0.2)


def _k_of_T(T):
    return 15.0 + 0.002 * T


def _new_path_and_wall():
    z = np.linspace(0.0, L, N + 1)
    path = FlowPath(
        name="tube",
        s_edges=z,
        z_of_s_edges=z.copy(),
        A_flow=np.pi * D_I**2 / 4.0,
        Dh=D_I,
        P_wetted=np.pi * D_I,
        P_heated=np.pi * D_I,
        n_parallel=N_TUBES,
        roughness=ROUGHNESS,
        flow_direction=1,
    )
    wall = CylindricalWall(D_inner=D_I, thickness=T_WALL_THICKNESS, k_of_T=_k_of_T, hot_side="inner")
    return path, wall


def _legacy_geometry():
    grid = AxialGrid.uniform(
        length=L, n_cells=N, coolant_area=1.0, wall_area=1.0,
        hot_perimeter=1.0, coolant_perimeter=1.0, flow_direction=1,
    )
    return ShellTubeCoreGeometry(
        grid=grid,
        tube_inner_diameter=D_I,
        tube_outer_diameter=D_I + 2 * T_WALL_THICKNESS,
        wall_thickness=T_WALL_THICKNESS,
        n_tubes=N_TUBES,
        shell_inner_diameter=0.1,
    )


@pytest.mark.parametrize("inside_tube_choice", ["smooth", "grooved"])
def test_hot_gas_march_matches_legacy_exactly(inside_tube_choice):
    corrCoeffs = CorrelationCoefficients()
    Tbar = np.full(N, 700.0)
    Tc = np.full(N, 350.0)
    hs = np.full(N, 800.0)
    hrad = np.array([100.0, 90.0, 80.0, 70.0])
    mdot = 0.1

    path, wall = _new_path_and_wall()
    provider_new, _ = equilibrium_gas_state_provider(_FakeEquilibriumManifold())
    result = hot_gas_march(
        path, wall,
        Tbar_wall=Tbar, T_coolant=Tc, h_shell=hs, mdot_hot_total=mdot,
        gas_state_at=provider_new, inside_tube_choice=inside_tube_choice,
        corrCoeffs=corrCoeffs, corrugation_thickness_m=3e-4, corrugation_pitch_m=4e-3,
        h_g_rad=hrad,
    )

    geometry = _legacy_geometry()
    provider_old, _ = equilibrium_gas_state_provider(_FakeEquilibriumManifold())
    legacy = shelltube_hot_gas_march(
        geometry,
        Tbar_wall=Tbar, T_coolant=Tc, h_shell=hs, mdot_hot_total=mdot,
        gas_state_at=provider_old, wall_conductivity_at_T=_k_of_T,
        inside_tube_choice=inside_tube_choice, nusselt_selector="gnielinski_blended",
        roughness=ROUGHNESS, corrCoeffs=corrCoeffs,
        corrugation_thickness=3e-4, corrugation_pitch=4e-3, h_g_rad=hrad,
    )

    np.testing.assert_array_equal(result.T_gas, legacy.T_gas)
    np.testing.assert_array_equal(result.h_gas_W_m2K, legacy.h_gas_W_m2K)
    np.testing.assert_array_equal(result.gas_velocity_m_s, legacy.gas_velocity_m_s)
    np.testing.assert_array_equal(result.reynolds, legacy.reynolds)
    np.testing.assert_array_equal(result.prandtl, legacy.prandtl)
    np.testing.assert_array_equal(result.friction_factor, legacy.friction_factor)
    # nusselt is a DERIVED diagnostic here (Nu = h*D/k, reversing legacy's
    # h = Nu*k/D) since core/closures.py's registered closures return h
    # directly, not Nu -- floating point division/multiplication isn't
    # perfectly invertible, so this one field carries ~1e-15 relative
    # reassociation noise even though h_gas_W_m2K (the physically
    # consequential quantity that actually drives the wall flux) matches
    # exactly. Not worth redundantly recomputing Nu via a second call into
    # the underlying dispatch function just for bit-perfect diagnostic
    # reporting on a field nothing downstream consumes.
    np.testing.assert_allclose(result.nusselt, legacy.nusselt, rtol=1e-12)
    np.testing.assert_array_equal(result.dp_per_length_Pa_m, legacy.dp_per_length_Pa_m)
    np.testing.assert_array_equal(result.T_wg, legacy.wall_flux.T_wg)
    np.testing.assert_array_equal(result.T_wc, legacy.wall_flux.T_wc)
    np.testing.assert_array_equal(result.dq_hot_per_length_W_m, legacy.wall_flux.dq_hot_per_length_W_m)
    np.testing.assert_array_equal(result.dq_cold_per_length_W_m, legacy.wall_flux.dq_cold_per_length_W_m)
    np.testing.assert_array_equal(result.hot_heat_W, legacy.wall_flux.hot_heat_W)
    np.testing.assert_array_equal(result.cold_heat_W, legacy.wall_flux.cold_heat_W)
    np.testing.assert_array_equal(result.k_wall, legacy.wall_flux.k_wall)
    np.testing.assert_array_equal(result.enthalpy_removed_J_kg, legacy.enthalpy_removed_J_kg)
    np.testing.assert_array_equal(result.progress_variable, legacy.progress_variable)
    assert result.T_gas_outlet == legacy.T_gas_outlet
    assert result.enthalpy_removed_outlet_J_kg == legacy.enthalpy_removed_outlet_J_kg
    assert result.progress_outlet == legacy.progress_outlet


def test_hot_gas_march_rejects_nonpositive_mdot():
    path, wall = _new_path_and_wall()
    provider, _ = equilibrium_gas_state_provider(_FakeEquilibriumManifold())
    with pytest.raises(ValueError, match="mdot_hot_total"):
        hot_gas_march(
            path, wall,
            Tbar_wall=np.full(N, 700.0), T_coolant=np.full(N, 350.0), h_shell=np.full(N, 800.0),
            mdot_hot_total=0.0, gas_state_at=provider, inside_tube_choice="smooth",
        )


def test_hot_gas_march_h_removed_trajectory_is_monotonic_for_cooling_gas():
    """Sanity check independent of the legacy comparison: a gas that's
    losing heat to the wall should show a monotonically increasing
    h_removed trajectory."""
    path, wall = _new_path_and_wall()
    provider, _ = equilibrium_gas_state_provider(_FakeEquilibriumManifold())
    result = hot_gas_march(
        path, wall,
        Tbar_wall=np.full(N, 500.0), T_coolant=np.full(N, 350.0), h_shell=np.full(N, 800.0),
        mdot_hot_total=0.1, gas_state_at=provider, inside_tube_choice="smooth",
    )
    assert np.all(np.diff(result.enthalpy_removed_J_kg) > 0.0)
    assert np.all(np.diff(result.T_gas) < 0.0)

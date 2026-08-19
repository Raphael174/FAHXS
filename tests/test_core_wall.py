"""Stage C of docs/solver_design/FV_CORE_REWORK_PLAN.md.

Gate: existing-config behavior preserved (the quadratic reconstruction and
hot_side orientation must match ``physics/heat_conduction.py``), and the NEW
rib/fin path verified independently -- it is dead code for the two existing
configs, so it gets its own hand-calculation and limiting-case checks rather
than riding on a solver regression.
"""
from __future__ import annotations

import numpy as np
import pytest

from hps_combustor.core.wall import (
    CylindricalWall,
    channel_effective_perimeter,
    rectangular_fin_efficiency,
)
from hps_combustor.physics.heat_conduction import (
    compute_eta_fin_rectangular,
    OneDimensionalSteadyConduction_ShellnHelicalTube,
)

K_STEEL = 20.0


def _k_const(T):
    return K_STEEL


# --- fin / rib model -----------------------------------------------------


def test_fin_efficiency_matches_legacy_factor_two_helper():
    """The rib between two channels is wetted on both faces -> factor 2.
    Must agree with physics/heat_conduction.py::compute_eta_fin_rectangular
    (the factor-2 variant), NOT compute_fin_efficiency_ch (factor 1)."""
    h_c, t_rib, h_ch = 15000.0, 1.2e-3, 4e-3
    eta, m = rectangular_fin_efficiency(h_c=h_c, k_w=K_STEEL, t_rib=t_rib, h_channel=h_ch)
    eta_legacy, m_legacy = compute_eta_fin_rectangular(h_c, K_STEEL, t_rib, h_ch)
    assert eta == pytest.approx(eta_legacy, rel=1e-12)
    assert m == pytest.approx(m_legacy, rel=1e-12)


def test_fin_efficiency_hand_calculation():
    h_c, t_rib, h_ch, k = 10000.0, 1e-3, 3e-3, 25.0
    eta, m = rectangular_fin_efficiency(h_c=h_c, k_w=k, t_rib=t_rib, h_channel=h_ch)
    m_expected = np.sqrt(2.0 * h_c / (k * t_rib))
    assert m == pytest.approx(m_expected, rel=1e-12)
    assert eta == pytest.approx(np.tanh(m_expected * h_ch) / (m_expected * h_ch), rel=1e-12)


def test_fin_efficiency_bounded_zero_to_one():
    for h_c in (1e2, 1e3, 1e4, 1e5, 1e6):
        eta, _ = rectangular_fin_efficiency(h_c=h_c, k_w=K_STEEL, t_rib=1e-3, h_channel=5e-3)
        assert 0.0 < eta <= 1.0


def test_fin_efficiency_approaches_unity_for_short_stubby_fin():
    """Short, thick, highly conductive fin -> nearly isothermal -> eta -> 1."""
    eta, _ = rectangular_fin_efficiency(h_c=100.0, k_w=400.0, t_rib=5e-3, h_channel=1e-4)
    assert eta > 0.999


def test_fin_efficiency_degrades_for_tall_thin_fin():
    """Tall, thin, low-conductivity fin under high h -> strongly derated."""
    eta, _ = rectangular_fin_efficiency(h_c=5e4, k_w=15.0, t_rib=0.5e-3, h_channel=1.5e-2)
    assert eta < 0.2


def test_effective_perimeter_between_plain_and_full_geometric():
    """P_eff must sit between the no-fin-credit case (root only) and the
    full geometric perimeter (eta=1)."""
    w, h, t = 2e-3, 5e-3, 1e-3
    fin = channel_effective_perimeter(
        w_channel=w, h_channel=h, t_rib=t, h_c=2e4, k_w=K_STEEL, include_closeout=True
    )
    no_credit = 2.0 * w                 # root + closeout, ribs contribute nothing
    full_geometric = 2.0 * w + 2.0 * h  # eta = 1
    assert no_credit < fin.P_effective_m < full_geometric


def test_effective_perimeter_closeout_toggle():
    kw = dict(w_channel=2e-3, h_channel=5e-3, t_rib=1e-3, h_c=2e4, k_w=K_STEEL)
    with_co = channel_effective_perimeter(include_closeout=True, **kw).P_effective_m
    without_co = channel_effective_perimeter(include_closeout=False, **kw).P_effective_m
    assert with_co == pytest.approx(without_co + 2e-3, rel=1e-12)


def test_ribbed_wall_raises_without_film_coefficient():
    wall = CylindricalWall(
        D_inner=0.1, thickness=1e-3, k_of_T=_k_const,
        rib=dict(w_channel=2e-3, h_channel=4e-3, t_rib=1e-3, n_channels=60),
    )
    with pytest.raises(ValueError, match="requires h_c and k_w"):
        wall.perimeters()


def test_ribbed_wall_increases_cold_side_perimeter():
    """The point of ribs: far more cold-side area than a plain bore."""
    plain = CylindricalWall(D_inner=0.1, thickness=1e-3, k_of_T=_k_const, hot_side="inner")
    ribbed = CylindricalWall(
        D_inner=0.1, thickness=1e-3, k_of_T=_k_const, hot_side="inner",
        rib=dict(w_channel=2e-3, h_channel=4e-3, t_rib=1e-3, n_channels=60),
    )
    _, P_cold_plain = plain.perimeters(h_c=2e4, k_w=K_STEEL)
    _, P_cold_ribbed = ribbed.perimeters(h_c=2e4, k_w=K_STEEL)
    assert P_cold_ribbed > P_cold_plain


# --- orientation ---------------------------------------------------------


def test_hot_side_orientation_swaps_perimeters():
    """CLAUDE.md do-not-regress item: shell-and-tube needs hot_side='inner'."""
    D, s = 5e-3, 0.75e-3
    inner = CylindricalWall(D_inner=D, thickness=s, k_of_T=_k_const, hot_side="inner")
    outer = CylindricalWall(D_inner=D, thickness=s, k_of_T=_k_const, hot_side="outer")
    P_hot_i, P_cold_i = inner.perimeters()
    P_hot_o, P_cold_o = outer.perimeters()
    assert P_hot_i == pytest.approx(np.pi * D)
    assert P_cold_i == pytest.approx(np.pi * (D + 2 * s))
    assert (P_hot_o, P_cold_o) == (P_cold_i, P_hot_i)


def test_invalid_hot_side_rejected():
    with pytest.raises(ValueError, match="hot_side"):
        CylindricalWall(D_inner=0.01, thickness=1e-3, k_of_T=_k_const, hot_side="sideways")


# --- steady / transient consistency --------------------------------------


def test_steady_solution_conserves_flux_across_wall():
    wall = CylindricalWall(D_inner=0.0135, thickness=0.85e-3, k_of_T=_k_const)
    faces = wall.solve_steady(T_g=2000.0, T_c=400.0, h_g=800.0, h_c=5000.0)
    assert faces.q_hot_W_m == pytest.approx(faces.q_cold_W_m, rel=1e-12)


def test_steady_temperature_ordering():
    wall = CylindricalWall(D_inner=0.0135, thickness=0.85e-3, k_of_T=_k_const)
    f = wall.solve_steady(T_g=2000.0, T_c=400.0, h_g=800.0, h_c=5000.0)
    assert 400.0 < f.T_wc < f.T_wg < 2000.0


def test_fluxes_at_Tbar_reduces_to_steady_at_the_steady_mean_wall_temp():
    """Self-consistency guardrail the legacy module documents: evaluated at
    the steady mean wall temperature, the transient reconstruction must
    reproduce the steady faces and have equal in/out fluxes."""
    wall = CylindricalWall(D_inner=0.0135, thickness=0.85e-3, k_of_T=_k_const)
    T_g, T_c, h_g, h_c = 2000.0, 400.0, 800.0, 5000.0
    steady = wall.solve_steady(T_g=T_g, T_c=T_c, h_g=h_g, h_c=h_c)
    T_bar = 0.5 * (steady.T_wg + steady.T_wc)
    tr = wall.fluxes_at_Tbar(T_bar=T_bar, T_g=T_g, T_c=T_c, h_g=h_g, h_c=h_c)
    assert tr.T_wg == pytest.approx(steady.T_wg, rel=2e-3)
    assert tr.T_wc == pytest.approx(steady.T_wc, rel=2e-3)
    assert tr.q_hot_W_m == pytest.approx(tr.q_cold_W_m, rel=5e-3)


def test_fluxes_at_Tbar_matches_legacy_heat_conduction_module():
    """The generalized wall must reproduce the validated legacy
    reconstruction (physics/heat_conduction.py::fluxes_at_Tbar) for the
    plain-tube helical case it was written for."""
    D, s = 0.0135, 0.85e-3
    T_g, T_c, h_g, h_c, T_bar = 1800.0, 350.0, 750.0, 4200.0, 900.0

    legacy = OneDimensionalSteadyConduction_ShellnHelicalTube(
        h_g=h_g, h_c=h_c, T_c=T_c, T_g=T_g, s_w=s, Dh_ch=D,
        f_kw_at_T=_k_const, T_wg_0=T_g, T_wc_0=T_c, T_c_check_0=T_c,
        hot_side="outer",
    )
    ref = legacy.fluxes_at_Tbar(T_bar, h_g_rad=0.0)

    wall = CylindricalWall(D_inner=D, thickness=s, k_of_T=_k_const, hot_side="outer")
    got = wall.fluxes_at_Tbar(T_bar=T_bar, T_g=T_g, T_c=T_c, h_g=h_g, h_c=h_c)

    assert got.T_wg == pytest.approx(ref["T_wg"], rel=1e-12)
    assert got.T_wc == pytest.approx(ref["T_wc"], rel=1e-12)
    assert got.q_hot_W_m == pytest.approx(ref["dq_hot__dx"], rel=1e-12)
    assert got.q_cold_W_m == pytest.approx(ref["dq_cold__dx"], rel=1e-12)


def test_fluxes_at_Tbar_matches_legacy_hot_side_inner():
    """Same equivalence with the shell-and-tube orientation."""
    D, s = 3.5e-3, 0.75e-3
    T_g, T_c, h_g, h_c, T_bar = 2200.0, 500.0, 1200.0, 3000.0, 1100.0
    legacy = OneDimensionalSteadyConduction_ShellnHelicalTube(
        h_g=h_g, h_c=h_c, T_c=T_c, T_g=T_g, s_w=s, Dh_ch=D,
        f_kw_at_T=_k_const, T_wg_0=T_g, T_wc_0=T_c, T_c_check_0=T_c,
        hot_side="inner",
    )
    ref = legacy.fluxes_at_Tbar(T_bar, h_g_rad=0.0)
    wall = CylindricalWall(D_inner=D, thickness=s, k_of_T=_k_const, hot_side="inner")
    got = wall.fluxes_at_Tbar(T_bar=T_bar, T_g=T_g, T_c=T_c, h_g=h_g, h_c=h_c)
    assert got.T_wg == pytest.approx(ref["T_wg"], rel=1e-12)
    assert got.q_hot_W_m == pytest.approx(ref["dq_hot__dx"], rel=1e-12)


def test_temperature_dependent_conductivity_is_used():
    """k_of_T must be evaluated, not ignored."""
    calls = []

    def k_of_T(T):
        calls.append(T)
        return 10.0 + 0.01 * T

    wall = CylindricalWall(D_inner=0.01, thickness=1e-3, k_of_T=k_of_T)
    faces = wall.fluxes_at_Tbar(T_bar=800.0, T_g=1500.0, T_c=300.0, h_g=500.0, h_c=3000.0)
    assert calls and faces.k_w == pytest.approx(10.0 + 0.01 * 800.0)


def test_radiation_coefficient_adds_in_parallel_on_hot_side():
    wall = CylindricalWall(D_inner=0.0135, thickness=0.85e-3, k_of_T=_k_const)
    base = wall.fluxes_at_Tbar(T_bar=900.0, T_g=1800.0, T_c=350.0, h_g=750.0, h_c=4200.0)
    with_rad = wall.fluxes_at_Tbar(
        T_bar=900.0, T_g=1800.0, T_c=350.0, h_g=750.0, h_c=4200.0, h_g_rad=200.0
    )
    assert with_rad.q_hot_W_m > base.q_hot_W_m
    assert with_rad.T_wg > base.T_wg


def test_heat_capacity_per_length_matches_annulus_hand_calc():
    D, s, rho = 0.0135, 0.85e-3, 7900.0
    wall = CylindricalWall(D_inner=D, thickness=s, k_of_T=_k_const)
    C = wall.heat_capacity_per_length(rho, lambda T: 500.0, 700.0)
    r_i, r_o = D / 2.0, D / 2.0 + s
    assert C == pytest.approx(rho * 500.0 * np.pi * (r_o**2 - r_i**2), rel=1e-12)

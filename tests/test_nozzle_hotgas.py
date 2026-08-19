"""Stage F groundwork tests (see docs/solver_design/FV_CORE_REWORK_PLAN.md
section 5, "Staged implementation" Stage F item 2): area-Mach and T_aw
against closed-form isentropic relations, tight tolerance since these are
exact analytic relations, not fitted correlations.
"""
from __future__ import annotations

import numpy as np
import pytest

from hps_combustor.core.geometry.nozzle_contour import (
    NozzleContour,
    build_conical_contour,
)
from hps_combustor.core.hotgas.nozzle_gas import (
    ChamberState,
    adiabatic_wall_temperature,
    area_mach_ratio,
    bartz_cornelisse_htc,
    choked_mass_flux,
    isentropic_static_state,
    mach_from_area_ratio,
    throat_diameter_for_mass_flow,
)

GAMMA = 1.2268  # C2H4/O2 chamber value, O/F=2.3, 50 bar (2026-07-31 session)


def test_area_mach_ratio_is_unity_at_M1():
    assert area_mach_ratio(1.0, GAMMA) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("M_true", [0.1, 0.3, 0.7, 0.95])
def test_mach_from_area_ratio_subsonic_round_trip(M_true):
    AR = area_mach_ratio(M_true, GAMMA)
    M_back = mach_from_area_ratio(AR, GAMMA, supersonic=False)
    assert M_back == pytest.approx(M_true, rel=1e-8)


@pytest.mark.parametrize("M_true", [1.05, 1.5, 2.0, 3.356, 5.0])
def test_mach_from_area_ratio_supersonic_round_trip(M_true):
    AR = area_mach_ratio(M_true, GAMMA)
    M_back = mach_from_area_ratio(AR, GAMMA, supersonic=True)
    assert M_back == pytest.approx(M_true, rel=1e-8)


def test_mach_at_area_ratio_one_is_sonic():
    assert mach_from_area_ratio(1.0, GAMMA, supersonic=False) == pytest.approx(1.0, abs=1e-6)
    assert mach_from_area_ratio(1.0, GAMMA, supersonic=True) == pytest.approx(1.0, abs=1e-6)


def test_area_ratio_below_one_raises():
    with pytest.raises(ValueError):
        mach_from_area_ratio(0.5, GAMMA, supersonic=False)


def test_isentropic_static_state_stagnation_limit():
    """M->0 recovers stagnation T, p, rho (to numerical tolerance)."""
    chamber = ChamberState(T0_K=3760.0, p0_Pa=50e5, gamma=GAMMA, R_specific_J_kgK=391.1)
    T, p, rho, V = isentropic_static_state(1e-6, chamber)
    assert T == pytest.approx(chamber.T0_K, rel=1e-5)
    assert p == pytest.approx(chamber.p0_Pa, rel=1e-5)
    assert rho == pytest.approx(chamber.rho0_kg_m3, rel=1e-5)
    assert V == pytest.approx(0.0, abs=1.0)


def test_isentropic_static_state_matches_closed_form_at_M1():
    """Closed-form critical-point ratios: T*/T0=2/(g+1),
    p*/p0=(2/(g+1))^(g/(g-1)), rho*/rho0=(2/(g+1))^(1/(g-1))."""
    chamber = ChamberState(T0_K=3760.0, p0_Pa=50e5, gamma=GAMMA, R_specific_J_kgK=391.1)
    T, p, rho, V = isentropic_static_state(1.0, chamber)
    g = GAMMA
    ratio = 2.0 / (g + 1.0)
    assert T / chamber.T0_K == pytest.approx(ratio, rel=1e-10)
    assert p / chamber.p0_Pa == pytest.approx(ratio ** (g / (g - 1.0)), rel=1e-10)
    assert rho / chamber.rho0_kg_m3 == pytest.approx(ratio ** (1.0 / (g - 1.0)), rel=1e-10)
    a = np.sqrt(g * chamber.R_specific_J_kgK * T)
    assert V == pytest.approx(a, rel=1e-10)  # M=1 => V equals local sound speed


def test_adiabatic_wall_temperature_bounds():
    """T_aw must lie between static T (Pr->0 limit, r->0) and T0 (Pr=1,
    recovery factor=1 exactly recovers stagnation temperature)."""
    chamber = ChamberState(T0_K=3760.0, p0_Pa=50e5, gamma=GAMMA, R_specific_J_kgK=391.1)
    T, p, rho, V = isentropic_static_state(1.0, chamber)
    T_aw_Pr1 = adiabatic_wall_temperature(T, chamber, M=1.0, Pr=1.0, turbulent=True)
    assert T_aw_Pr1 == pytest.approx(chamber.T0_K, rel=1e-10)
    T_aw_low_Pr = adiabatic_wall_temperature(T, chamber, M=1.0, Pr=1e-9, turbulent=True)
    assert T_aw_low_Pr == pytest.approx(T, rel=1e-3)


def test_bartz_cornelisse_htc_hand_calc():
    """Direct hand-computed value from the formula, independent of the
    implementation's expression grouping."""
    mu, cp, Pr, rho, V, D = 1.0e-4, 2000.0, 0.6, 3.0, 1500.0, 0.12
    rho_f, mu_f = 2.5, 1.1e-4
    h = bartz_cornelisse_htc(
        mu_Pa_s=mu, cp_J_kgK=cp, Pr=Pr, rho_kg_m3=rho, V_m_s=V, D_m=D,
        rho_film_kg_m3=rho_f, mu_film_Pa_s=mu_f,
    )
    expected = (
        0.026 * (mu**0.2 * cp / Pr**0.6) * (rho * V) ** 0.8 / D**0.2
        * (rho_f / rho) * (mu_f / mu)
    )
    assert h == pytest.approx(expected, rel=1e-12)


def test_bartz_htc_positive_and_finite():
    h = bartz_cornelisse_htc(
        mu_Pa_s=1e-4, cp_J_kgK=2000.0, Pr=0.6, rho_kg_m3=3.0, V_m_s=1500.0,
        D_m=0.12, rho_film_kg_m3=2.5, mu_film_Pa_s=1.1e-4,
    )
    assert np.isfinite(h) and h > 0.0


def test_throat_diameter_for_mass_flow_round_trips_choked_mass_flux():
    """throat_diameter_for_mass_flow must be the exact inverse of
    choked_mass_flux -- feeding the resulting diameter's implied mdot back
    in must reproduce the target mdot (the 2026-07-31 geometry-consistency
    finding depends on this round-trip being exact)."""
    chamber = ChamberState(T0_K=3760.0, p0_Pa=50e5, gamma=GAMMA, R_specific_J_kgK=390.9)
    mdot_target = 45.0
    D_t = throat_diameter_for_mass_flow(mdot_target, chamber)
    G_t = choked_mass_flux(chamber)
    A_t = np.pi / 4.0 * D_t**2
    assert G_t * A_t == pytest.approx(mdot_target, rel=1e-9)


def test_choked_mass_flux_matches_120mm_throat_mismatch_finding():
    """Regression pin for the specific inconsistency found 2026-07-31: a
    120mm throat at 50 bar (this chamber chemistry) chokes ~30.5 kg/s, not
    the originally-stated ~45 kg/s."""
    chamber = ChamberState(T0_K=3760.0, p0_Pa=50e5, gamma=GAMMA, R_specific_J_kgK=390.9)
    G_t = choked_mass_flux(chamber)
    A_t_120mm = np.pi / 4.0 * 0.120**2
    mdot_120mm = G_t * A_t_120mm
    assert mdot_120mm == pytest.approx(30.5, abs=0.5)


def test_conical_contour_throat_and_expansion_ratio():
    contour = build_conical_contour(D_throat_m=0.120, expansion_ratio=10.0)
    assert contour.D_t_m == pytest.approx(0.120, rel=1e-9)
    assert contour.A_m2[-1] / contour.A_t_m2 == pytest.approx(10.0, rel=1e-9)
    assert contour.z_m[contour.throat_index] == pytest.approx(0.0, abs=1e-9)
    assert np.argmin(contour.r_m) == contour.throat_index


def test_conical_contour_rejects_bad_throat_index():
    with pytest.raises(ValueError):
        NozzleContour(z_m=np.array([0.0, 1.0, 2.0]), r_m=np.array([1.0, 0.5, 0.8]), throat_index=0)


def test_conical_contour_rejects_nonincreasing_z():
    with pytest.raises(ValueError):
        NozzleContour(z_m=np.array([0.0, 0.0, 1.0]), r_m=np.array([1.0, 0.5, 0.8]), throat_index=1)

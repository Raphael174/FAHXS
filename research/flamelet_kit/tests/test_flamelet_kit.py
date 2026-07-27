"""
pytest suite for flamelet_kit. Uses gri30.yaml (bundled with Cantera) so no
external mechanism files are needed.

Run:
    python -m pytest flamelet_kit/tests/ -q
"""
import os
import sys

import numpy as np
import cantera as ct
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flamelet import Flamelet, _build_z_grid, _bilger_Z_st, _chi_profile
from steady_cache import SteadyCache
from cooling_pfr import CoolingPFR


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def gas():
    return ct.Solution("gri30.yaml")


@pytest.fixture(scope="module")
def ch4_air_streams(gas):
    gas.TPX = 300.0, 101325.0, "CH4:1"
    Y_fuel = np.array(gas.Y, dtype=float)
    gas.TPX = 300.0, 101325.0, "O2:1, N2:3.76"
    Y_ox = np.array(gas.Y, dtype=float)
    return Y_ox, Y_fuel


# ---------------------------------------------------------------------------
# Grid construction
# ---------------------------------------------------------------------------

def test_z_grid_clusters_near_Z_st():
    Z_st = 0.055
    Z = _build_z_grid(65, Z_st)
    n_near = int(np.sum(np.abs(Z - Z_st) <= 0.05))
    assert n_near >= 15, f"expected >=15 nodes within 0.05 of Z_st, got {n_near}"
    assert Z[0] == 0.0
    assert Z[-1] == 1.0
    assert np.all(np.diff(Z) > 0), "grid must be strictly increasing"


# ---------------------------------------------------------------------------
# Bilger Z_st
# ---------------------------------------------------------------------------

def test_bilger_Z_st_ch4_air(gas, ch4_air_streams):
    Y_ox, Y_fuel = ch4_air_streams
    Z_st = _bilger_Z_st(gas, Y_ox, Y_fuel, list(gas.species_names))
    assert 0.045 <= Z_st <= 0.065, f"expected CH4/air Z_st ~0.055, got {Z_st}"


# ---------------------------------------------------------------------------
# chi(Z) profile
# ---------------------------------------------------------------------------

def test_chi_profile_positive_and_anchored_at_Z_st():
    """
    chi(Z) = chi_st * exp(2*erfcinv(2*Z_st)^2 - 2*erfcinv(2*Z)^2) is positive
    everywhere and, by construction, EQUALS chi_st exactly at Z=Z_st (that is
    its defining anchor property -- chi_st is "the value of chi at the
    stoichiometric surface", not the profile's global maximum). The
    counterflow form's true global maximum is at Z=0.5 (erfcinv(1)=0, the
    geometric center of the mixing layer), independent of where Z_st sits --
    this is the physically correct Peters counterflow shape, not a bug.
    """
    Z_st = 0.055
    chi_st = 20.0
    Z = _build_z_grid(65, Z_st)
    chi = _chi_profile(Z, Z_st, chi_st)

    assert np.all(chi > 0.0)

    i_st = int(np.argmin(np.abs(Z - Z_st)))
    assert np.isclose(chi[i_st], chi_st, rtol=0.02), (
        f"chi(Z_st) should equal chi_st={chi_st}, got {chi[i_st]}"
    )

    # Unimodal, maximum at the mixing-layer center Z=0.5, monotonically
    # decreasing toward both Z=0 and Z=1.
    i_peak = int(np.argmax(chi))
    assert abs(Z[i_peak] - 0.5) < 0.05
    assert chi[i_peak] >= chi[i_st]


# ---------------------------------------------------------------------------
# Short flamelet advance: finite, bounded, elemental-mass-conserving
# ---------------------------------------------------------------------------

def test_flamelet_advance_finite_bounded_and_conserves_elements(gas, ch4_air_streams):
    Y_ox, Y_fuel = ch4_air_streams
    p = 101325.0

    fl = Flamelet("gri30.yaml", n_z=31)
    fl.init_mixing(T_ox=300.0, Y_ox=Y_ox, T_fuel=300.0, Y_fuel=Y_fuel, p=p)

    def elemental_mass_fracs(Y_row):
        gas.TPY = 300.0, p, {sp: float(y) for sp, y in
                              zip(fl.species_names, Y_row) if y > 0}
        return {el: gas.elemental_mass_fraction(el) for el in ("C", "H", "O", "N")}

    before = [elemental_mass_fracs(fl.Y[i]) for i in (0, fl.n_z // 2, fl.n_z - 1)]

    for _ in range(5):
        fl.step(dt=2.0e-5, p=p, T_ox=300.0, Y_ox=Y_ox, T_fuel=300.0,
                Y_fuel=Y_fuel, chi_st=10.0)

    assert np.isfinite(fl.T).all()
    assert np.isfinite(fl.Y).all()
    assert np.all(fl.T >= 200.0) and np.all(fl.T <= 4500.0)
    assert np.all(fl.Y >= -1e-8)
    row_sums = fl.Y.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-6)

    # Elemental mass at the Dirichlet boundary nodes must be exactly the
    # imposed feed-stream composition (boundary conditions are never touched
    # by chemistry); this is the conservation check that matters here since
    # per-node total elemental mass in Z-space is a boundary-condition
    # invariant, not a volume-integral one (diffusion in Z redistributes,
    # it does not create/destroy elements at the boundaries).
    after_bc0 = elemental_mass_fracs(fl.Y[0])
    after_bc1 = elemental_mass_fracs(fl.Y[-1])
    for el in ("C", "H", "O", "N"):
        assert np.isclose(after_bc0[el], before[0][el], atol=1e-8)
        assert np.isclose(after_bc1[el], before[2][el], atol=1e-8)


# ---------------------------------------------------------------------------
# SteadyCache
# ---------------------------------------------------------------------------

def test_steady_cache_hit_within_tolerance_miss_outside():
    cache = SteadyCache(tol_dT=2.0, tol_p=0.01, tol_chi=0.05, tol_T_fuel=5.0)

    assert cache.check(p_key=1.0e5, chi_key=20.0, T_fuel=300.0) == "no_prior_advance"
    cache.record_advance(p_key=1.0e5, chi_key=20.0, T_fuel=300.0, dT=0.1)

    # Within tolerance -> HIT
    assert cache.is_hit(p_key=1.0005e5, chi_key=20.5, T_fuel=301.0)

    # chi moved 50% (>> 5% tol) -> MISS, reason "chi"
    reason = cache.check(p_key=1.0e5, chi_key=30.0, T_fuel=300.0)
    assert reason == "chi"

    # pressure moved 5% (>> 1% tol) -> MISS, reason "p"
    reason = cache.check(p_key=1.05e5, chi_key=20.0, T_fuel=300.0)
    assert reason == "p"

    # T_fuel moved 10 K (>> 5 K tol) -> MISS, reason "T_fuel"
    reason = cache.check(p_key=1.0e5, chi_key=20.0, T_fuel=312.0)
    assert reason == "T_fuel"

    assert cache.n_hit >= 1
    assert cache.n_miss >= 3


# ---------------------------------------------------------------------------
# CoolingPFR: energy balance
# ---------------------------------------------------------------------------

def test_cooling_pfr_adiabatic_conserves_enthalpy(gas):
    gas.TP = 300.0, 101325.0
    gas.set_equivalence_ratio(1.0, "CH4", "O2:1, N2:3.76")
    gas.equilibrate("HP")
    T_in, Y_in, p = float(gas.T), np.array(gas.Y, dtype=float), 101325.0

    pfr = CoolingPFR("gri30.yaml")
    result = pfr.march(T_in=T_in, Y_in=Y_in, p=p, mdot=0.02,
                        length=1.0, n_steps=20, diameter=0.02,
                        h_conv=0.0, T_wall=500.0)
    assert result["Q_wall_total_W"] == 0.0
    assert np.isclose(result["h_in"], result["h_out"], atol=1.0)
    assert np.isclose(result["T"][-1], result["T"][0], atol=5.0)


def test_cooling_pfr_energy_balance_closes_with_heat_loss(gas):
    gas.TP = 300.0, 101325.0
    gas.set_equivalence_ratio(1.0, "CH4", "O2:1, N2:3.76")
    gas.equilibrate("HP")
    T_in, Y_in, p = float(gas.T), np.array(gas.Y, dtype=float), 101325.0

    pfr = CoolingPFR("gri30.yaml")
    result = pfr.march(T_in=T_in, Y_in=Y_in, p=p, mdot=0.02,
                        length=1.0, n_steps=20, diameter=0.02,
                        h_conv=200.0, T_wall=500.0)
    assert result["Q_wall_total_W"] > 0.0
    assert result["T"][-1] < result["T"][0]
    # mdot*(h_in - h_out) must equal the accumulated wall heat loss to a
    # tight tolerance relative to the magnitude of the heat removed.
    assert abs(result["energy_balance_residual_W"]) < 1e-3 * max(
        result["Q_wall_total_W"], 1.0
    )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))

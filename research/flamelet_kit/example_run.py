"""
End-to-end demo of flamelet_kit using gri30.yaml (bundled with Cantera --
no extra mechanism files needed). Methane/air two-stream non-premixed
flamelet at a fixed chi_st, plus a SteadyCache hit/miss demonstration.

Run:
    python flamelet_kit/example_run.py
"""
from __future__ import annotations

import time

import numpy as np
import cantera as ct

from flamelet import Flamelet
from steady_cache import SteadyCache


def build_streams(gas: ct.Solution):
    """Pure-CH4 fuel stream and air (O2:1, N2:3.76 mole) oxidizer stream."""
    gas.TPX = 300.0, 101325.0, "CH4:1"
    Y_fuel = np.array(gas.Y, dtype=float)

    gas.TPX = 300.0, 101325.0, "O2:1, N2:3.76"
    Y_ox = np.array(gas.Y, dtype=float)
    return Y_ox, Y_fuel


def seed_flame_kernel(fl: Flamelet, gas: ct.Solution, p: float,
                       half_width: float = 0.25) -> None:
    """
    Seed nodes near Z_st with the adiabatic-flame (HP-equilibrium) state at
    their local mixture composition, instead of the cold linear-mixing
    profile.

    This is a DEMO-ONLY convenience, deliberately kept out of the library
    (flamelet.py has no ignition/spark logic per the no-ignition/extinction
    scope of this kit -- see ADAPTATION_GUIDE.md). It stands in for what a
    real steady heat-exchanger deployment does naturally: you initialize a
    new representative flamelet from a neighboring already-converged
    condition rather than from cold mixing, because the whole premise of the
    steady/steadily-changing regime is that you are always close to an
    existing burning solution. Starting this demo from pure cold mixing
    would instead simulate a CH4/air autoignition delay at 300 K, which is
    minutes-to-hours of chemical time -- not a meaningful few-step demo.
    """
    for i in range(fl.n_z):
        if abs(fl.Z[i] - fl.Z_st) > half_width:
            continue
        Y_mix = fl.Y[i].copy()
        try:
            gas.TPY = float(fl.T[i]), p, {sp: float(y) for sp, y in
                                           zip(fl.species_names, Y_mix) if y > 0}
            gas.equilibrate("HP")
            fl.T[i] = float(np.clip(gas.T, 250.0, 4500.0))
            fl.Y[i] = np.maximum(np.array(gas.Y, dtype=float), 0.0)
        except Exception:
            pass
    fl.sync_reactors(p)


def main() -> None:
    p = 101325.0
    chi_st = 3.0  # s^-1 -- moderate strain; see note below on why not higher
    n_z = 65

    gas = ct.Solution("gri30.yaml")
    Y_ox, Y_fuel = build_streams(gas)

    fl = Flamelet("gri30.yaml", n_z=n_z)
    fl.init_mixing(T_ox=300.0, Y_ox=Y_ox, T_fuel=300.0, Y_fuel=Y_fuel, p=p)
    print(f"Z_st (Bilger) = {fl.Z_st:.4f}")
    print(f"nodes within |Z-Z_st|<0.05: {fl.n_nodes_near_Z_st(0.05)}")

    seed_flame_kernel(fl, gas, p, half_width=0.3)

    # Note on chi_st=3/s: GRI-30 CH4/air non-premixed extinction strain is
    # numerically close to what this coarse-ish demo grid/dt can resolve
    # stably; chi_st=3 was found (by direct trial) to relax to a genuine
    # quasi-steady plateau rather than slowly quenching, so it is the
    # honest choice for a short, fast-running demo -- NOT a claim about the
    # physical CH4/air extinction strain rate. Real deployments should
    # dt/grid-converge before trusting a specific chi_st's burn/quench outcome.
    dt = 2.0e-5
    n_steps = 60
    t0 = time.perf_counter()
    for _ in range(n_steps):
        fl.step(dt, p, T_ox=300.0, Y_ox=Y_ox, T_fuel=300.0, Y_fuel=Y_fuel,
                chi_st=chi_st)
    wall = time.perf_counter() - t0

    Z_peak = float(fl.Z[np.argmax(fl.T)])
    print(f"After {n_steps} steps (dt={dt:.1e}s, {wall:.1f}s wall):")
    print(f"  T_max        = {fl.T_max:.1f} K")
    print(f"  Z at T_max   = {Z_peak:.4f}")
    print(f"  t_flamelet   = {fl.t_flamelet:.2e} s")

    # --- Regime-1 -> regime-2 hand-off preview (see cooling_pfr.py / -----
    #     example_cooling.py for the full downstream march) ---------------
    T_burnt = float(fl.T_at_Z(fl.Z_st))
    print(f"\nBurnt-gas state at Z_st (hand-off to CoolingPFR regime):")
    print(f"  T(Z_st) = {T_burnt:.1f} K")

    # --- SteadyCache hit/miss demonstration -------------------------------
    cache = SteadyCache()
    T_fuel_bc = 300.0

    reason1 = cache.check(p, chi_st, T_fuel_bc)
    print(f"\ncache.check() before any record   -> {reason1!r} (expect 'no_prior_advance')")
    cache.record_advance(p_key=p, chi_key=chi_st, T_fuel=T_fuel_bc, dT=0.5)

    reason2 = cache.check(p, chi_st, T_fuel_bc)
    print(f"cache.check() same p/chi/T_fuel   -> {reason2!r} (expect '' = HIT)")

    chi_moved = chi_st * 1.5  # 50% jump, outside default 5% tol_chi
    reason3 = cache.check(p, chi_moved, T_fuel_bc)
    print(f"cache.check() chi_st moved 50%     -> {reason3!r} (expect 'chi' = MISS)")

    print(f"\ncache stats: hits={cache.n_hit} misses={cache.n_miss} "
          f"hit_rate={cache.hit_rate:.2f}")


if __name__ == "__main__":
    main()

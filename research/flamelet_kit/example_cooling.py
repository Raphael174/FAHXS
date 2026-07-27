"""
End-to-end demo of the regime-2 cooling PFR, using GRI-30 CH4/air adiabatic-
flame products as a stand-in "burnt gas" inlet (in a real diesel/O2 heat
exchanger this inlet would instead come from `Flamelet.T_at_Z(Z_st)` /
`Flamelet.Y_at_Z(Z_st, ...)` -- see example_run.py's hand-off preview).

Marches the burnt gas down a channel with a fixed wall temperature, cooler
than the gas, and shows T(x) dropping while accumulating wall heat loss.
Also demonstrates the SteadyCache-gated `march_cached` reuse.

Run:
    python flamelet_kit/example_cooling.py
"""
from __future__ import annotations

import numpy as np
import cantera as ct

from cooling_pfr import CoolingPFR
from steady_cache import SteadyCache


def adiabatic_flame_state(mechanism: str, T_in: float, p: float):
    """Stoichiometric CH4/air adiabatic-flame (HP-equilibrium) state --
    a cheap stand-in for a flamelet's burnt-gas hand-off state."""
    gas = ct.Solution(mechanism)
    gas.TP = T_in, p
    gas.set_equivalence_ratio(1.0, "CH4", "O2:1, N2:3.76")
    gas.equilibrate("HP")
    return float(gas.T), np.array(gas.Y, dtype=float), list(gas.species_names)


def main() -> None:
    mechanism = "gri30.yaml"
    p = 101325.0

    T_burnt, Y_burnt, species_names = adiabatic_flame_state(mechanism, 300.0, p)
    print(f"Inlet (stand-in burnt-gas) state: T = {T_burnt:.1f} K")

    pfr = CoolingPFR(mechanism)

    mdot = 0.02       # kg/s
    diameter = 0.02    # m
    length = 1.0       # m
    n_steps = 40
    h_conv = 200.0     # W/m^2/K
    T_wall = 500.0     # K, constant

    result = pfr.march(
        T_in=T_burnt, Y_in=Y_burnt, p=p, mdot=mdot,
        length=length, n_steps=n_steps, diameter=diameter,
        h_conv=h_conv, T_wall=T_wall,
    )
    print(f"\nMarched {length} m channel in {n_steps} segments "
          f"(mdot={mdot} kg/s, D={diameter} m, h_conv={h_conv} W/m2K, "
          f"T_wall={T_wall} K):")
    print(f"  T(0)         = {result['T'][0]:.1f} K")
    print(f"  T(L)         = {result['T'][-1]:.1f} K")
    print(f"  Q_wall_total = {result['Q_wall_total_W']:.1f} W")
    print(f"  energy balance residual = {result['energy_balance_residual_W']:.3e} W "
          f"(should be ~0)")

    # --- Adiabatic sanity check: h_conv=0 -> no cooling, h conserved -------
    result_adiabatic = pfr.march(
        T_in=T_burnt, Y_in=Y_burnt, p=p, mdot=mdot,
        length=length, n_steps=n_steps, diameter=diameter,
        h_conv=0.0, T_wall=T_wall,
    )
    print(f"\nAdiabatic check (h_conv=0): T(L) = {result_adiabatic['T'][-1]:.1f} K "
          f"(vs T(0) = {result_adiabatic['T'][0]:.1f} K), "
          f"Q_wall_total = {result_adiabatic['Q_wall_total_W']:.3e} W")

    # --- SteadyCache-gated reuse demonstration ------------------------------
    cache = SteadyCache()
    march_kwargs = dict(length=length, n_steps=n_steps, diameter=diameter,
                        h_conv=h_conv, T_wall=T_wall)

    r1 = pfr.march_cached(cache, T_burnt, Y_burnt, p, mdot, **march_kwargs)
    print(f"\nmarch_cached() 1st call -> real march "
          f"(hits={cache.n_hit}, misses={cache.n_miss})")

    r2 = pfr.march_cached(cache, T_burnt, Y_burnt, p, mdot, **march_kwargs)
    same = np.allclose(r1["T"], r2["T"])
    print(f"march_cached() same (T_in, mdot, p) -> reused cached result "
          f"(hits={cache.n_hit}, misses={cache.n_miss}), profiles identical: {same}")

    r3 = pfr.march_cached(cache, T_burnt, Y_burnt, p, mdot * 2.0, **march_kwargs)
    print(f"march_cached() mdot doubled -> re-marched "
          f"(hits={cache.n_hit}, misses={cache.n_miss})")


if __name__ == "__main__":
    main()

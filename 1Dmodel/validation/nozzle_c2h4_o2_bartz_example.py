"""First-pass C2H4/O2 nozzle hot-gas sizing example.

Design point given by the user, 2026-07-31: conical nozzle, D_throat=120 mm,
expansion ratio=10, C2H4/O2 at O/F=2.3 (mass), p0=50 bar, total mdot~45 kg/s.
"Keep it simple for now" -- this is a standalone script exercising the Stage F
groundwork (``core/geometry/nozzle_contour.py``, ``core/hotgas/nozzle_gas.py``),
not yet wired into the coupled FV core (Stages B-E aren't done -- see
docs/solver_design/FV_CORE_REWORK_PLAN.md). No wall/coolant coupling: heat
flux uses an ASSUMED uniform hot-wall temperature for Bartz's film
properties, not a solved one.

**Geometry consistency finding (2026-07-31):** D_throat=120mm, p0=50 bar, and
mdot~45 kg/s are NOT mutually consistent for this chamber chemistry --
choking ties throat area, chamber pressure, and mass flow together; a 120mm
throat at 50 bar only passes ~30.5 kg/s here, not 45. Mass flow is normally
the design-driving quantity (tied to thrust/Isp), so this script treats
``MDOT_TOTAL_KG_S`` + ``P0_PA`` as authoritative and DERIVES throat diameter
from choked flow (~146mm) rather than using the originally-stated 120mm
directly -- both are reported below so the discrepancy stays visible instead
of being silently resolved one way.

Chemistry caveat: chamber equilibrium uses Cantera's bundled GRI-Mech 3.0
(``gri30.yaml``), which includes C2H4/O2/CO2/H2O/CO/H2/OH/O/H with valid NASA
polynomial thermo data, but is tuned/validated for natural-gas combustion and
NOx chemistry, not purpose-built for ethylene/oxygen rocket combustion.
Equilibrium composition/temperature depends on thermodynamics (present here)
more than on the reaction mechanism's kinetics (irrelevant to an equilibrium
calculation), so this should be a reasonable first-pass chamber state, but a
dedicated combustion mechanism should replace it before trusting numbers past
a first sizing pass.

Run directly:
    python -m hps_combustor.validation.nozzle_c2h4_o2_bartz_example
"""

from __future__ import annotations

import cantera as ct
from scipy.constants import gas_constant

from hps_combustor.core.geometry.nozzle_contour import build_conical_contour
from hps_combustor.core.hotgas.nozzle_gas import (
    ChamberState,
    choked_mass_flux,
    solve_frozen_expansion,
    throat_diameter_for_mass_flow,
)

D_THROAT_STATED_M = 0.120   # as originally given -- inconsistent with mdot below, see finding
EXPANSION_RATIO = 10.0
OF_MASS = 2.3
P0_PA = 50e5
MDOT_TOTAL_KG_S = 45.0      # as given -- treated as authoritative, throat is derived from it
T_REACTANTS_K = 298.15
T_WALL_GUESS_K = 800.0  # placeholder only -- see module docstring; swept below


def chamber_state() -> tuple[ChamberState, ct.Solution]:
    """Equilibrium C2H4/O2 chamber state; returns the frozen-composition gas
    object too (reused, TP re-set per station, for Bartz property lookups)."""
    gas = ct.Solution("gri30.yaml")
    gas.TPY = T_REACTANTS_K, P0_PA, f"O2:{OF_MASS}, C2H4:1"
    gas.equilibrate("HP")
    R_specific = gas_constant * 1000.0 / gas.mean_molecular_weight  # J/kg/K
    chamber = ChamberState(
        T0_K=float(gas.T),
        p0_Pa=float(P0_PA),
        gamma=float(gas.cp_mass / gas.cv_mass),
        R_specific_J_kgK=float(R_specific),
    )
    # Freeze this composition for the expansion -- gas.TP is re-set per
    # station downstream, gas.Y is never touched again after this point.
    return chamber, gas


def run() -> None:
    chamber, gas = chamber_state()
    print("=== Chamber state (C2H4/O2, O/F=2.3 mass, p0=50 bar) ===")
    print(f"  T0 = {chamber.T0_K:.1f} K")
    print(f"  gamma = {chamber.gamma:.4f}")
    print(f"  MW = {gas.mean_molecular_weight:.3f} g/mol")
    print(f"  rho0 = {chamber.rho0_kg_m3:.4f} kg/m3")

    G_t = choked_mass_flux(chamber)
    A_t_stated = 3.141592653589793 / 4.0 * D_THROAT_STATED_M**2
    mdot_at_stated_throat = G_t * A_t_stated
    D_t_for_mdot = throat_diameter_for_mass_flow(MDOT_TOTAL_KG_S, chamber)
    print("\n=== Geometry consistency check (choking ties D_throat, p0, mdot"
          " together -- they were over-specified) ===")
    print(f"  stated D_throat={D_THROAT_STATED_M*1e3:.1f} mm at p0=50 bar chokes"
          f" {mdot_at_stated_throat:.2f} kg/s (stated total mdot was"
          f" {MDOT_TOTAL_KG_S:.0f} kg/s, {mdot_at_stated_throat/MDOT_TOTAL_KG_S*100:.0f}%)")
    print(f"  D_throat consistent with {MDOT_TOTAL_KG_S:.0f} kg/s at p0=50 bar:"
          f" {D_t_for_mdot*1e3:.1f} mm  <- used below")

    contour = build_conical_contour(
        D_throat_m=D_t_for_mdot,
        expansion_ratio=EXPANSION_RATIO,
    )
    print("\n=== Contour (conical, default 15deg divergent / 30deg convergent"
          " half-angles, contraction ratio 6 -- all ASSUMED, see"
          " nozzle_contour.py docstring) ===")
    print(f"  D_throat = {contour.D_t_m*1e3:.1f} mm")
    print(f"  D_exit = {contour.D_m2[-1]*1e3:.1f} mm")
    print(f"  length (convergent start to exit) = "
          f"{(contour.z_m[-1]-contour.z_m[0])*1e3:.1f} mm")

    stations = solve_frozen_expansion(
        contour, chamber, gas, T_wall_guess_K=T_WALL_GUESS_K
    )

    throat = stations[contour.throat_index]
    peak = max(stations, key=lambda s: s.q_w_W_m2)
    exit_station = stations[-1]

    print(f"\n=== Hot-gas expansion (frozen, constant gamma={chamber.gamma:.4f};"
          f" Bartz Cornelisse film-property form; T_wall assumed"
          f" {T_WALL_GUESS_K:.0f} K uniform -- placeholder, not solved) ===")
    print(f"  throat: M={throat.M:.3f}  T={throat.T_K:.1f} K  p={throat.p_Pa/1e5:.2f} bar"
          f"  h_g={throat.h_g_W_m2K:.0f} W/m2K  T_aw={throat.T_aw_K:.1f} K"
          f"  q_w={throat.q_w_W_m2/1e6:.2f} MW/m2")
    print(f"  peak q_w: z={peak.z_m*1e3:+.1f} mm (relative to throat)  "
          f"M={peak.M:.3f}  h_g={peak.h_g_W_m2K:.0f} W/m2K  "
          f"q_w={peak.q_w_W_m2/1e6:.2f} MW/m2")
    print(f"  exit: M={exit_station.M:.3f}  T={exit_station.T_K:.1f} K  "
          f"p={exit_station.p_Pa/1e5:.3f} bar  h_g={exit_station.h_g_W_m2K:.0f} W/m2K  "
          f"q_w={exit_station.q_w_W_m2/1e6:.3f} MW/m2")

    print("\n=== Throat q_w sensitivity to the ASSUMED wall temperature"
          " (no wall/coolant coupling exists yet -- this is the single"
          " biggest current source of uncertainty in the number above) ===")
    for T_wall in (600.0, 800.0, 1000.0, 1200.0, 1400.0):
        s = solve_frozen_expansion(contour, chamber, gas, T_wall_guess_K=T_wall)[contour.throat_index]
        print(f"  T_wall={T_wall:6.0f} K   h_g={s.h_g_W_m2K:7.0f} W/m2K"
              f"   q_w={s.q_w_W_m2/1e6:6.2f} MW/m2")


if __name__ == "__main__":
    run()

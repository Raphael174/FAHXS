"""First supercritical-N2 (LN2 coolant) run on the helical steady solver.

LN2 at 80 bar is fully SUPERCRITICAL (N2 Pc = 33.96 bar): no boiling dome, a
pseudo-critical transition near 145.7 K. Inlet 100 K (supercritical liquid-like),
co-flow, frozen chemistry (coolant-side supercritical physics is independent of
hot-gas chemistry; frozen is fast and sufficient for this first check).

Run:  python -m studies.ln2_supercritical.run_ln2_helical    (from repo root)
   or python studies/ln2_supercritical/run_ln2_helical.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "water_vs_helium"))

import numpy as np
import CoolProp.CoolProp as CP

from hps_combustor.input_data import (
    CorrelationCoefficients, combustorProp, coolantProp, hotgasProp,
    numericalProp, system_requirements,
)
from hps_combustor.main_solve import main_solver
from run_case import L_HX_MAX_6M

T_IN = 100.0
P_IN = 80e5
MDOT_C = 0.10


def main():
    coolant = coolantProp(
        coolant="Nitrogen", coolant_model="equilibrium_liquid",
        mass_flow_c=MDOT_C, T_in=T_IN, p_in=P_IN,
    )
    solver = main_solver(
        coolantProp=coolant, hotgasProp=hotgasProp(),
        combustorProp=combustorProp(HX_config="shellnHelicalTube", flow_config="co"),
        numericalProp=numericalProp(L_HX_max=L_HX_MAX_6M, chemistry_model="frozen"),
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
    )
    solver.solver()
    solver._check_global()

    d = solver.data_master
    x = np.asarray(d["L_HX"], dtype=float)
    T_c = np.asarray(d["T_c"], dtype=float)
    Re_c = np.asarray(d["Re_c"], dtype=float)
    T_wc = np.asarray(d["T_wc"], dtype=float)
    T_wg = np.asarray(d["T_wg"], dtype=float)
    h_c = np.asarray(d["h_c"], dtype=float)
    p_c = np.asarray(d["p_c"], dtype=float)

    Q_tot = float(np.sum(np.asarray(d["dQ"], dtype=float)))  # compute_performance() not needed
    dh = CP.PropsSI("H", "T", T_c[-1], "P", p_c[-1], "Nitrogen") - CP.PropsSI(
        "H", "T", T_c[0], "P", p_c[0], "Nitrogen")
    Q_coolant = MDOT_C * dh

    print("=" * 64)
    print("LN2 SUPERCRITICAL HELICAL RUN")
    print("=" * 64)
    print(f"  N_nodes = {len(x)}   coil arc length = {x[-1]:.3f} m")
    print(f"  T_c: {T_c[0]:.1f} -> {T_c[-1]:.1f} K   (pseudo-critical ~145.7 K)")
    print(f"  Re_c: {Re_c.min():.0f} -> {Re_c.max():.0f}  (Cheng2020 window 7000-27000)")
    print(f"  h_c (coolant HTC): {h_c.min():.0f} -> {h_c.max():.0f} W/m2K")
    print(f"  T_wg max = {T_wg.max():.1f} K   T_wc max = {T_wc.max():.1f} K")
    print(f"  p_c: {p_c[0]/1e5:.2f} -> {p_c[-1]/1e5:.2f} bar")
    print(f"  Q_tot = {Q_tot/1e3:.2f} kW")
    print(f"  energy check: Q_tot={Q_tot/1e3:.2f} kW vs mdot*dh={Q_coolant/1e3:.2f} kW "
          f"(rel err {abs(Q_tot-Q_coolant)/max(abs(Q_coolant),1)*100:.3f}%)")
    print("=" * 64)


if __name__ == "__main__":
    main()

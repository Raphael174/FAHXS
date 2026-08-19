"""Shell-and-tube (Inconel 718) steady design-point search — LN2/supercritical
N2 coolant.

Target (given 2026-08-19, updated 2026-08-19, due 2026-08-21): from 100 K
inlet, find the operating point on the CURRENT shell-and-straight-tube
Inconel 718 design delivering ~0.9 MW at ~75 bar exit pressure, with the
coolant flow rate constrained by a target ~30 L/s volumetric flow rate
MEASURED AT THE OUTLET (delivered-gas basis, same convention as the water
case in `friday_shelltube_water.py` -- confirmed via AskUserQuestion).

N2's critical pressure is ~34 bar; a ~75 bar exit-pressure target is
inherently supercritical for N2 at any inlet pressure that can plausibly
reach it with a physically small pressure drop, so `quality` reports NaN
throughout (confirmed, not a gap) -- this reads as "cryogenic-N2-supplied,
operating supercritically" rather than "boiling LN2", so there is no
separate two-phase-outlet constraint to satisfy here the way there is for
water (N2 never has an "is it real steam" question -- it's a dense
supercritical fluid at the outlet either way).

Found by direct search: at the ORIGINAL point (`mass_flow_c=9.7 kg/s`,
`p_in=80 bar`) the outlet density is high enough (449.9 kg/m3, T_c_out=139K)
that the outlet volumetric flow is only 21.6 L/s -- short of the 30 L/s
target. Increasing `mass_flow_c` raises `Vdot_out` (lower-density,
higher-T_c_out state) but also increases the shell-side pressure drop
sharply (Bell-Delaware dP scales close to mdot_c^1.8 in this regime,
confirmed by sweep from 9.7-30 kg/s), driving `p_out` well under the 75 bar
target unless `p_in` is raised to compensate. Unlike the water case, N2
is NOT wall-temperature-constrained here -- `T_wg_max` stays under ~520 K
across the whole sweep (INCO718's characterized ceiling is 1033 K), so
there is no material-limit conflict to trade against; `mass_flow_c` and
`p_in` can be tuned freely against the volumetric-flow and pressure targets
alone.

**Chosen design point**: `mass_flow_g=0.115 kg/s`, `mass_flow_c=17.5 kg/s`,
`p_in=88 bar`. Result: `Q_tot=942.5 kW` (close to the ~900 kW target),
`T_c_out=124.3 K`, `p_out=76.10 bar` (on the ~75 bar target),
`Vdot_out=30.0 L/s` (on target), `T_wg_max=475.0 K` (far under the material
ceiling), supercritical throughout (`quality` NaN).

Run: `python -m hps_combustor.validation.friday_shelltube_n2`
"""
from __future__ import annotations

import numpy as np

from ..input_data import (
    CorrelationCoefficients,
    combustorProp,
    coolantProp,
    hotgasProp,
    numericalProp,
    runProp,
    shellTubeProp,
    system_requirements,
)
from ..main_steady import run_steady
from ..physics.liquid_flow.coolprop_state_cache import coolprop_fluid_string
from ..physics.liquid_flow.regime import real_fluid_state_ph
from ..result_package import package_steady_run

# Design point found by search (see module docstring).
MASS_FLOW_G = 0.115  # kg/s, hot gas (diesel/O2) -- free variable, tuned for target power
MASS_FLOW_C = 17.5   # kg/s, N2 coolant -- free variable, tuned for outlet Vdot target
T_IN_K = 100.0
P_IN_PA = 88.0e5
VDOT_TARGET_LS = 30.0  # target outlet volumetric flow rate, L/s (delivered-gas basis)


def build_inputs():
    return {
        "coolant": coolantProp(
            coolant="Nitrogen", coolant_model="equilibrium_liquid",
            mass_flow_c=MASS_FLOW_C, T_in=T_IN_K, p_in=P_IN_PA,
        ),
        "hotgas": hotgasProp(mass_flow_g=MASS_FLOW_G),
        "combustor": combustorProp(HX_config="shellntube", flow_config="co"),
        "shelltube": shellTubeProp(),
        "numerical": numericalProp(chemistry_model="finite_rate"),
        "system": system_requirements(),
        "correlations": CorrelationCoefficients(),
        "run": runProp(run_name="friday_shellntube_n2_30Ls"),
    }


def run_case():
    inputs = build_inputs()
    solver, summary = run_steady(inputs)

    sl = solver.shell_liquid
    Twg = np.asarray(solver.tube["T_wg"])

    cool_cp = coolprop_fluid_string("Nitrogen", solver._liquid_backend)
    outlet_state = real_fluid_state_ph(cool_cp, sl["p_out"], sl["h_out"])
    Vdot_out_Ls = MASS_FLOW_C / outlet_state.rho_kg_m3 * 1000.0

    summary = dict(summary)
    summary.update({
        "mass_flow_g_kg_s": MASS_FLOW_G,
        "mass_flow_c_kg_s": MASS_FLOW_C,
        "p_in_bar": P_IN_PA / 1e5,
        "p_out_bar": float(sl["p_out"] / 1e5),
        "quality_min": float(sl["quality"].min()),
        "quality_max": float(sl["quality"].max()),
        "supercritical_throughout": bool(np.all(np.isnan(sl["quality"]))),
        "rho_out_kg_m3": float(outlet_state.rho_kg_m3),
        "Vdot_out_Ls": float(Vdot_out_Ls),
        "Vdot_target_Ls": VDOT_TARGET_LS,
        "T_wg_max_vs_inco718_data_ceiling_K": (float(Twg.max()), 1033.15),
    })
    return solver, summary


def main():
    solver, summary = run_case()
    package = package_steady_run(solver, build_inputs(), summary)

    print()
    print("=" * 60)
    print("FRIDAY DESIGN POINT — shell-and-tube, LN2/supercritical N2")
    print("=" * 60)
    print(f"  mass_flow_g = {MASS_FLOW_G:.3f} kg/s   mass_flow_c = {MASS_FLOW_C:.2f} kg/s")
    print(f"  p_in = {P_IN_PA/1e5:.1f} bar   T_in = {T_IN_K:.1f} K")
    print(f"  Q_tot        = {summary['Q_tot_kW']:.1f} kW   (target ~900 kW)")
    print(f"  T_c_out      = {summary['T_c_out_K']:.1f} K")
    print(f"  p_out        = {summary['p_out_bar']:.2f} bar  (target ~75 bar)")
    print(f"  Vdot_out     = {summary['Vdot_out_Ls']:.2f} L/s  (target {summary['Vdot_target_Ls']:.0f} L/s, "
          f"delivered-gas basis)")
    print(f"  T_wg_max     = {summary['T_wg_max_vs_inco718_data_ceiling_K'][0]:.1f} K "
          f"(INCO718 characterized data ceiling: {summary['T_wg_max_vs_inco718_data_ceiling_K'][1]:.0f} K)")
    print(f"  collapse_margin = {summary['collapse_margin']:.4f}  (want << 1)")
    print(f"  supercritical throughout = {summary['supercritical_throughout']} "
          f"(quality is NaN when true -- expected, no two-phase dome crossed)")
    print(f"  n_sweeps     = {summary['n_sweeps']}")
    print()
    print(f"  Saved: {package['folder']}")
    if package["archive"]:
        print(f"  Archive: {package['archive']}")
    return solver, summary


if __name__ == "__main__":
    main()

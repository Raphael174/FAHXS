"""Shell-and-tube (Inconel 718) steady design point - LN2/supercritical N2.

Target: from a 100 K inlet, deliver of order 0.9 MW at ~75 bar exit pressure
on the CURRENT shell-and-straight-tube Inconel 718 design.

**Coolant mass flow is now a FIXED INPUT, not a free variable** (given
2026-08-20): `mass_flow_c = 1.916 kg/s`, set by independent system needs
outside this model. Only `mass_flow_g` and `p_in` are tuned here. This
supersedes the earlier 17.5 kg/s point, which had been found by treating the
coolant flow as free in order to chase a 30 L/s outlet volumetric target.

N2's critical pressure is ~34 bar, so any operating point able to hold ~75
bar at exit is supercritical throughout and `quality` reports NaN (confirmed,
not a gap). This reads as "cryogenic-N2-supplied, operating supercritically"
rather than "boiling LN2": there is no two-phase-outlet constraint to satisfy
the way there is for water.

**Chosen design point**: `mass_flow_g=0.125 kg/s`, `mass_flow_c=1.916 kg/s`
(fixed), `p_in=101.5 bar`. Result: `Q_tot=918.4 kW` (on the ~0.9 MW target),
`T_c_out=400.7 K`, `p_out=73.86 bar` (on the ~75 bar target),
`Vdot_out=31.7 L/s`, `T_wg_max=1032.2 K` -- 0.95 K inside INCO718's
characterized-data ceiling of 1033.15 K, and far inside the ~1500 K
allowance confirmed 2026-08-19. Supercritical throughout.

That wall-temperature margin is thin. Backing `mass_flow_g` off to 0.120
kg/s at `p_in=101.0 bar` trades 5% of the duty for 13 K of wall margin:
Q_tot=867 kW, p_out=75.65 bar, T_wg_max=1019.6 K. Use that point instead if
the characterized-data boundary matters more than hitting 0.9 MW exactly.

**The character of this case changed with the lower coolant flow.** At 17.5
kg/s the N2 stayed cryogenic end to end (T_c_out = 124 K, liquid-like dense
fluid). At 1.916 kg/s the same duty is carried by roughly a ninth of the
mass, so the coolant leaves at **402 K** -- warm, gas-like, and no longer
cryogenic. Anything downstream that assumed a cold N2 outlet needs
rechecking against this point. It also pushes T_wg_max from 475 K up to
1031 K, because less coolant mass flow means a hotter wall for the same duty.

**p_in is 101.5 bar, not ~88 bar**, for two reasons that compound: the
shell-side pressure-drop model was corrected (see below), and the coolant
leaves hot and therefore much less dense, so it accelerates through the
bundle and the drop through the back half of the exchanger is large
(27.6 bar total here, of which 6.96 bar -- 23% -- is the momentum
/acceleration term that had previously been omitted altogether).

Sensitivity to `mass_flow_g` at the fixed coolant flow, `p_in` retuned each
time to hold the exit near 75 bar:

    mass_flow_g   Q_tot     T_c_out   p_out    T_wg_max
      0.070 kg/s   429 kW    189 K     74.8 b    858 K
      0.090 kg/s   588 kW    249 K     71.2 b    929 K
      0.115 kg/s   819 kW    354 K     75.0 b   1005 K
      0.125 kg/s   922 kW    402 K     74.6 b   1031 K   <- chosen
      0.130 kg/s   969 kW    425 K     74.6 b   1042 K

Duty and wall temperature rise together; 0.125 kg/s is the point that meets
the 0.9 MW target while still sitting inside the characterized material data.

**Shell-side pressure drop model (corrected 2026-08-20)**: the pressure march
now uses Bell-Delaware rather than the closure's own friction gradient. The
Gungor-Winterton/MSH and supercritical-registry gradients are straight-TUBE
correlations - they model axial flow along one L_tube-long channel with wall
skin friction, whereas the real shell-side path crosses the bundle
N_baffles+1 times through N_tcc rows each: about 7.5x the path length on this
geometry, with form drag rather than skin friction. The straight-tube
gradient under-predicts by roughly 25x here, which is far too large to treat
as an acceptable extrapolation. Note also that Bell-Delaware's own baffle-
leakage ratio r_lm is ~6.5 on this geometry against a fitted range of about
r_lm <= 1, so its leakage corrections are themselves extrapolated; the solver
reports this at runtime.

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
MASS_FLOW_G = 0.125  # kg/s, hot gas (diesel/O2) -- free variable, tuned for target power
MASS_FLOW_C = 1.916  # kg/s, N2 coolant -- FIXED by independent system needs (2026-08-20)
T_IN_K = 100.0
P_IN_PA = 101.5e5
VDOT_TARGET_LS = 30.0  # reference only: coolant flow is now fixed, so Vdot is an outcome


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
        "run": runProp(run_name="friday_shellntube_n2_0p9MW"),
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

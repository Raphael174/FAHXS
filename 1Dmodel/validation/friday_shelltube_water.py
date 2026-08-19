"""Shell-and-tube (Inconel 718) steady design-point search — Water coolant.

Target (given 2026-08-19, updated 2026-08-19, due 2026-08-21): from 300 K
inlet, find the operating point on the CURRENT shell-and-straight-tube
Inconel 718 design (`shellTubeProp()` defaults -- 235 tubes, grooved,
INCO718 -- confirmed this IS the Inconel design; no separate preset file
exists) that delivers up to 700 K outlet / 2.53 MW, at ~75 bar exit
pressure. **Updated requirement**: the outlet must actually be superheated
STEAM (single-phase vapor, to use as a pressurant), not a two-phase
liquid/vapor mixture -- and the coolant flow rate is constrained by a
target ~30 L/s volumetric flow rate MEASURED AT THE OUTLET (delivered-gas
basis, confirmed via AskUserQuestion: pressurant systems are sized by
delivered volumetric flow, not feed-liquid flow).

Found by direct 3-variable search this session (`mass_flow_g`,
`mass_flow_c`, `p_in`; geometry/material fixed): at the ORIGINAL water point
(`mass_flow_c=2.0 kg/s`) the coolant only reaches quality=0.058 (barely
started boiling, no steam) at `T_c_out=565 K`. Driving `mass_flow_c` down
far enough to fully vaporize AND superheat the coolant (quality > 1, single
phase) directly trades against the tube gas-wall temperature: less coolant
mass flow means a hotter wall for the same duty.

The first pass of this search capped `T_wg_max` at INCO718's
characterized-data ceiling (1033.15 K, `R02_INCO718["x"]` in
`mechanical/material_specs/`), which forced a compromise point short of
both the 2.53 MW and 30 L/s targets. **Updated 2026-08-19**: confirmed
acceptable for `T_wg` to run as high as ~1500 K (well past the
characterized-data range -- material behavior there is extrapolated, not
validated, but explicitly accepted for this study) -- re-opening the
search over `mass_flow_g` up to ~0.42 kg/s and `mass_flow_c` up to ~0.95
kg/s.

**Chosen design point**: `mass_flow_g=0.39 kg/s`, `mass_flow_c=0.86 kg/s`,
`p_in=81 bar`. Result: `Q_tot=2563 kW` (101% of the 2530 kW target --
matched, not just approached), `T_c_out=656.5 K` (still short of 700 K --
see below), `p_out=76.30 bar` (on the ~75 bar target), `quality_out=1.23`
(genuinely superheated single-phase steam, confirmed via
`real_fluid_state_ph` on the converged outlet (p,h) -- quality > 1 is this
codebase's convention for post-dome superheat, not clamped to [0,1]),
`Vdot_out=29.96 L/s` (essentially exact match to the 30 L/s target),
`T_wg_max=1070.9 K` (well inside the ~1500 K allowance, though still above
the material's characterized-data ceiling -- flagged, not hidden).
`T_c_out=656.5 K` remains short of 700 K: pushing further needs even less
coolant mass flow, which now trades against the *duty and volumetric-flow*
match instead of the wall ceiling (confirmed by sweep: matching Q and Vdot
simultaneously pins `mass_flow_c` near 0.86 kg/s at this `mass_flow_g`) --
a different, still-real, trade, not a stale constraint.

Run: `python -m hps_combustor.validation.friday_shelltube_water`
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
MASS_FLOW_G = 0.39   # kg/s, hot gas (diesel/O2) -- free variable, tuned for target power
MASS_FLOW_C = 0.86   # kg/s, water coolant -- free variable, tuned for outlet Vdot target
T_IN_K = 300.0
P_IN_PA = 81.0e5
INCO718_DATA_CEILING_K = 1033.15  # R02_INCO718["x"] characterized upper bound
T_WG_ALLOWANCE_K = 1500.0  # confirmed 2026-08-19: acceptable even past the data ceiling
VDOT_TARGET_LS = 30.0  # target outlet volumetric flow rate, L/s (delivered-gas basis)


def build_inputs():
    return {
        "coolant": coolantProp(
            coolant="Water", coolant_model="equilibrium_liquid",
            mass_flow_c=MASS_FLOW_C, T_in=T_IN_K, p_in=P_IN_PA,
        ),
        "hotgas": hotgasProp(mass_flow_g=MASS_FLOW_G),
        "combustor": combustorProp(HX_config="shellntube", flow_config="co"),
        "shelltube": shellTubeProp(),
        "numerical": numericalProp(chemistry_model="finite_rate"),
        "system": system_requirements(),
        "correlations": CorrelationCoefficients(),
        "run": runProp(run_name="friday_shellntube_water_steam"),
    }


def run_case():
    inputs = build_inputs()
    solver, summary = run_steady(inputs)

    sl = solver.shell_liquid
    Twg = np.asarray(solver.tube["T_wg"])
    chf = np.asarray(solver.tube["chf_margin"], dtype=float)

    cool_cp = coolprop_fluid_string("Water", solver._liquid_backend)
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
        "quality_out": float(outlet_state.quality),
        "outlet_is_superheated_steam": bool(outlet_state.quality > 1.0),
        "rho_out_kg_m3": float(outlet_state.rho_kg_m3),
        "Vdot_out_Ls": float(Vdot_out_Ls),
        "Vdot_target_Ls": VDOT_TARGET_LS,
        "chf_margin_finite_nodes": int((~np.isnan(chf)).sum()),
        "chf_margin_min_where_finite": float(np.nanmin(chf)) if np.any(~np.isnan(chf)) else None,
        "T_wg_max_vs_inco718_data_ceiling_K": (float(Twg.max()), INCO718_DATA_CEILING_K),
    })
    return solver, summary


def main():
    solver, summary = run_case()
    package = package_steady_run(solver, build_inputs(), summary)

    print()
    print("=" * 60)
    print("FRIDAY DESIGN POINT — shell-and-tube, Water (steam outlet)")
    print("=" * 60)
    print(f"  mass_flow_g = {MASS_FLOW_G:.3f} kg/s   mass_flow_c = {MASS_FLOW_C:.2f} kg/s")
    print(f"  p_in = {P_IN_PA/1e5:.1f} bar   T_in = {T_IN_K:.1f} K")
    print(f"  Q_tot        = {summary['Q_tot_kW']:.1f} kW   (target 2530 kW)")
    print(f"  T_c_out      = {summary['T_c_out_K']:.1f} K   (target 700 K -- see docstring: this point")
    print(f"                 prioritizes matching Q_tot and Vdot_out simultaneously instead)")
    print(f"  p_out        = {summary['p_out_bar']:.2f} bar  (target ~75 bar)")
    print(f"  quality_out  = {summary['quality_out']:.3f}  (>1 = superheated single-phase steam: "
          f"{summary['outlet_is_superheated_steam']})")
    print(f"  Vdot_out     = {summary['Vdot_out_Ls']:.2f} L/s  (target {summary['Vdot_target_Ls']:.0f} L/s, "
          f"delivered-gas basis)")
    print(f"  T_wg_max     = {summary['T_wg_max_vs_inco718_data_ceiling_K'][0]:.1f} K "
          f"(INCO718 characterized data ceiling: {summary['T_wg_max_vs_inco718_data_ceiling_K'][1]:.2f} K; "
          f"user-confirmed allowance: {T_WG_ALLOWANCE_K:.0f} K)")
    print(f"  collapse_margin = {summary['collapse_margin']:.4f}  (want << 1)")
    print(f"  quality range (full profile) = [{summary['quality_min']:.3f}, {summary['quality_max']:.3f}]")
    print(f"  CHF margin computed at {summary['chf_margin_finite_nodes']}/{len(solver.tube['chf_margin'])} nodes "
          f"(only where 0<=quality<=1; NaN elsewhere is expected once past full vaporization)")
    print(f"  n_sweeps     = {summary['n_sweeps']}")
    print()
    print(f"  Saved: {package['folder']}")
    if package["archive"]:
        print(f"  Archive: {package['archive']}")
    return solver, summary


if __name__ == "__main__":
    main()

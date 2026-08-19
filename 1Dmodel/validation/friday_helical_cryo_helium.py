"""Shell-and-helical-tube transient — cold gas-mode (cryo) Helium.

Target (given 2026-08-19, due 2026-08-21): a transient result for the
helical-coil design with cryogenic-supply Helium coolant. Confirmed (via
AskUserQuestion) as "cold gas-mode Helium" — same solver path already fixed
and regression-tested this session for the CFL instability
(`tests/test_transient_core_helical.py`), just at a cold inlet temperature,
not the near-critical/supercritical real-fluid path.

**T_in/p_in are a placeholder pending your confirmation**: you confirmed
"80-100 K range" is acceptable without a stricter number, so this uses
`T_in=80 K`, `p_in=80 bar` (matching the project's usual Helium pressure
convention) — already verified working standalone before this script was
written. Below `T_in=60 K` the path currently hits a hardcoded,
non-configurable guard constant in `core/coolant.py::enforce_internal_energy_bounds`
(`T_floor=60 K`, unrelated to Helium's real CoolProp range of 2.18 K) — not
attempted here since 80-100 K was confirmed sufficient; flag if a colder
point turns out to be needed later, the fix is small and understood but not
built.

**Transient scenario is a placeholder too**: no specific schedule (startup
ramp, bang-bang, etc.) was given, so this uses a simple linear mass-flow
ramp from near-zero to nominal over the first half of the run, then holds
steady — enough to show real transient wall/coolant response without
depending on a specific valve-schedule assumption. Replace
`schedule_mass_flow_c`/`schedule_mass_flow_g` below if a specific duty cycle
is actually required.

CFL-safe subcycling (this session's fix) means wall-clock cost scales with
how far `max_step` exceeds the coolant's cell residence time — reported
after the run, not promised in advance.

Run: `python -m hps_combustor.validation.friday_helical_cryo_helium`
"""
from __future__ import annotations

import time

from ..input_data import coolantProp
from ..main_transient import build_inputs as _base_build_inputs
from ..main_transient import run_transient
from ..result_package import package_transient_run

# Design point — see module docstring for the "placeholder, confirm with
# user" caveats on both the thermal point and the schedule.
T_IN_K = 80.0
P_IN_PA = 80.0e5
T_END_S = 2.0
MAX_STEP_S = 0.25
N_SAVE = 9

# Linear ramp 0 -> t_end/2, then hold at nominal for the rest of the run.
NOMINAL_MDOT_C = 0.06   # kg/s, matches this project's usual helical He flow scale
NOMINAL_MDOT_G = 0.10   # kg/s, matches shellTubeProp-adjacent hot-gas default scale
SCHEDULE_MDOT_C = ((0.0, 1e-3), (T_END_S / 2, NOMINAL_MDOT_C), (T_END_S, NOMINAL_MDOT_C))
SCHEDULE_MDOT_G = ((0.0, 5e-3), (T_END_S / 2, NOMINAL_MDOT_G), (T_END_S, NOMINAL_MDOT_G))


def build_inputs():
    inputs = _base_build_inputs()
    inputs["coolant"] = coolantProp(
        coolant="Helium", coolant_model="single_phase_coolprop",
        T_in=T_IN_K, p_in=P_IN_PA,
    )
    inputs["combustor"].HX_config = "shellnHelicalTube"
    inputs["combustor"].flow_config = "counter"
    inputs["transient"].fluid_model = "transient_coolant"
    inputs["transient"].t_end = T_END_S
    inputs["transient"].max_step = MAX_STEP_S
    inputs["transient"].n_save = N_SAVE
    inputs["transient"].schedule_mass_flow_c = SCHEDULE_MDOT_C
    inputs["transient"].schedule_mass_flow_g = SCHEDULE_MDOT_G
    inputs["run"].run_name = "friday_helical_cryo_helium_80K"
    return inputs


def run_case():
    inputs = build_inputs()
    t0 = time.perf_counter()
    solver, summary = run_transient(inputs)
    elapsed = time.perf_counter() - t0
    summary = dict(summary)
    summary["wall_clock_s"] = elapsed
    summary["T_in_K"] = T_IN_K
    summary["p_in_bar"] = P_IN_PA / 1e5
    return solver, summary, inputs


def main():
    solver, summary, inputs = run_case()
    package = package_transient_run(solver, inputs, summary)

    print()
    print("=" * 60)
    print("FRIDAY DESIGN POINT — shell-and-helical-tube, cryo He (transient)")
    print("=" * 60)
    print(f"  T_in = {T_IN_K:.1f} K   p_in = {P_IN_PA/1e5:.1f} bar")
    print(f"  t_end = {T_END_S:.1f} s   max_step = {MAX_STEP_S:.2f} s")
    print(f"  wall-clock = {summary['wall_clock_s']:.1f} s")
    print(f"  T_c_out_final    = {summary.get('T_c_out_final'):.1f} K")
    print(f"  T_wall_max_final = {summary.get('T_wall_max_final'):.1f} K")
    print(f"  energy_residual_J_final = {summary.get('energy_residual_J_final'):.3e} "
          f"(machine-precision expected — this session's CFL fix)")
    print()
    print(f"  Saved: {package['folder']}")
    if package["archive"]:
        print(f"  Archive: {package['archive']}")
    return solver, summary


if __name__ == "__main__":
    main()

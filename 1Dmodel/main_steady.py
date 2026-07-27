"""User entry point for steady HX simulations.

Choose the HX configuration in input_data.py:
  combustorProp.HX_config = "shellnHelicalTube" or "shellntube"

Run:
  python -m hps_combustor.main_steady
"""

from __future__ import annotations

import numpy as np

from .input_data import (
    CorrelationCoefficients,
    combustorProp,
    coolantProp,
    hotgasProp,
    numericalProp,
    runProp,
    shellTubeProp,
    system_requirements,
)
from .main_solve import (
    main_solver,
    solve_counterflow_liquid_reference,
    solve_counterflow_physical_reference,
)
from .main_solve_shellntube import shellntube_solver
from .result_package import package_steady_run


def main():
    inputs = build_inputs()
    solver, summary = run_steady(inputs)
    package = package_steady_run(solver, inputs, summary)

    print()
    print("Saved steady run:")
    print(f"  folder:  {package['folder']}")
    if package["archive"]:
        print(f"  archive: {package['archive']}")
    return solver


def build_inputs():
    return {
        "coolant": coolantProp(),
        "hotgas": hotgasProp(),
        "combustor": combustorProp(),
        "shelltube": shellTubeProp(),
        "numerical": numericalProp(),
        "system": system_requirements(),
        "correlations": CorrelationCoefficients(),
        "run": runProp(),
    }


def run_steady(inputs):
    hx_config = inputs["combustor"].HX_config
    if hx_config == "shellnHelicalTube":
        if inputs["combustor"].flow_config == "counter":
            # Always shoot on the physical cold-end inlet (T_in/p_in) for the
            # steady entry point, regardless of
            # numericalProp.counterflow_physical_steady_reference: the
            # alternative plain march starts from a guessed T_out/p_out at the
            # hot end and never actually enforces the user's T_in/p_in. For
            # coolant_model="equilibrium_liquid" that guess can't even
            # represent a two-phase state (T,P aren't independent inside the
            # dome), so it silently produces wrong ducts or a hard crash. This
            # exact failure mode recurred repeatedly with the guess-based
            # path, so the entry point no longer offers it for either fluid.
            # (numericalProp.counterflow_physical_steady_reference remains a
            # separate opt-in knob for the transient solver's settle-check
            # probe in main_solve_transient.py - unaffected by this.)
            reference_fn = (
                solve_counterflow_liquid_reference
                if inputs["coolant"].coolant_model == "equilibrium_liquid"
                else solve_counterflow_physical_reference
            )
            solver = reference_fn(
                coolantProp=inputs["coolant"],
                hotgasProp=inputs["hotgas"],
                combustorProp=inputs["combustor"],
                numericalProp=inputs["numerical"],
                system_requirements=inputs["system"],
                corrCoeffs=inputs["correlations"],
            )
        else:
            solver = main_solver(
                coolantProp=inputs["coolant"],
                hotgasProp=inputs["hotgas"],
                combustorProp=inputs["combustor"],
                numericalProp=inputs["numerical"],
                system_requirements=inputs["system"],
                corrCoeffs=inputs["correlations"],
            )
            solver.solver()
        summary = solver.compute_performance()
        solver.print_summary()
        return solver, summary

    if hx_config == "shellntube":
        solver = shellntube_solver(
            coolantProp=inputs["coolant"],
            hotgasProp=inputs["hotgas"],
            shellTubeProp=inputs["shelltube"],
            numericalProp=inputs["numerical"],
            system_requirements=inputs["system"],
            corrCoeffs=inputs["correlations"],
            N_axial=inputs["run"].shelltube_steady_nodes,
            flow_config=inputs["combustor"].flow_config,
        )
        solver.solve(verbose=True)
        stress = solver.compute_stress()
        solver.print_summary()
        return solver, _shelltube_summary(solver, stress)

    raise ValueError(
        "Unsupported combustorProp.HX_config for steady run: "
        f"{hx_config!r}. Use 'shellnHelicalTube' or 'shellntube'."
    )


def _shelltube_summary(solver, stress):
    tube = solver.tube
    return {
        "hx_config": "shellntube",
        "Q_tot_kW": float(solver.Q_tot / 1e3),
        "T_g_in_K": float(tube["T_g"][0]),
        "T_g_out_K": float(tube["T_g_out"]),
        "T_c_in_K": float(solver.coolantProp.T_in),
        "T_c_out_K": float(solver.T_c_out),
        "T_wg_max_K": float(np.max(tube["T_wg"])),
        "T_wc_max_K": float(np.max(tube["T_wc"])),
        "n_sweeps": int(solver.n_sweeps),
        **{key: float(value) for key, value in stress.items()},
    }


if __name__ == "__main__":
    main()

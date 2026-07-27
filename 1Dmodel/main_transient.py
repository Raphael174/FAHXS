"""User entry point for transient HX simulations.

Choose the HX configuration in input_data.py:
  combustorProp.HX_config = "shellnHelicalTube" or "shellntube"

Optional time schedules are set with runProp.schedule_file. CSV works without
extra dependencies; XLSX works when pandas/openpyxl are installed.

Run:
  python -m hps_combustor.main_transient
"""

from __future__ import annotations

from .input_data import (
    CorrelationCoefficients,
    combustorProp,
    coolantProp,
    hotgasProp,
    numericalProp,
    runProp,
    shellTubeProp,
    system_requirements,
    transientProp,
)
from .main_solve_shellntube_transient import shellntube_transient_solver
from .main_solve_transient import transient_solver
from .result_package import get_transient_time_series, package_transient_run
from .schedule_inputs import apply_schedule_file


def main():
    inputs = build_inputs()
    solver, summary = run_transient(inputs)
    package = package_transient_run(solver, inputs, summary)

    print()
    print("Saved transient run:")
    print(f"  folder:  {package['folder']}")
    if package["archive"]:
        print(f"  archive: {package['archive']}")
    return solver


def build_inputs():
    coolant = coolantProp()
    hotgas = hotgasProp()
    transient = transientProp()
    run = runProp()
    apply_schedule_file(transient, hotgas, run.schedule_file, coolant=coolant)

    return {
        "coolant": coolant,
        "hotgas": hotgas,
        "combustor": combustorProp(),
        "shelltube": shellTubeProp(),
        "numerical": numericalProp(),
        "transient": transient,
        "system": system_requirements(),
        "correlations": CorrelationCoefficients(),
        "run": run,
    }


def run_transient(inputs):
    hx_config = inputs["combustor"].HX_config
    if hx_config == "shellnHelicalTube":
        solver = transient_solver(
            coolantProp=inputs["coolant"],
            hotgasProp=inputs["hotgas"],
            combustorProp=inputs["combustor"],
            numericalProp=inputs["numerical"],
            system_requirements=inputs["system"],
            transientProp=inputs["transient"],
            corrCoeffs=inputs["correlations"],
        )
        if getattr(inputs["transient"], "fluid_model", "quasi_steady") == "transient_coolant":
            solver.solve_transient_core(verbose=True)
        else:
            solver.solve_transient(verbose=True)
        return solver, _summary_from_time_series(get_transient_time_series(solver))

    if hx_config == "shellntube":
        solver = shellntube_transient_solver(
            coolantProp=inputs["coolant"],
            hotgasProp=inputs["hotgas"],
            shellTubeProp=inputs["shelltube"],
            numericalProp=inputs["numerical"],
            system_requirements=inputs["system"],
            transientProp=inputs["transient"],
            corrCoeffs=inputs["correlations"],
            N_axial=inputs["run"].shelltube_transient_nodes,
            flow_config=inputs["combustor"].flow_config,
        )
        if getattr(inputs["transient"], "fluid_model", "quasi_steady") == "transient_coolant":
            solver.solve_transient_core(verbose=True)
        else:
            solver.solve_transient(verbose=True)
        return solver, _summary_from_time_series(get_transient_time_series(solver))

    raise ValueError(
        "Unsupported combustorProp.HX_config for transient run: "
        f"{hx_config!r}. Use 'shellnHelicalTube' or 'shellntube'."
    )


def _summary_from_time_series(time_series):
    scalars = time_series["scalars"]
    summary = {}
    for key, values in scalars.items():
        if len(values):
            summary[f"{key}_final"] = float(values[-1])
    return summary


if __name__ == "__main__":
    main()

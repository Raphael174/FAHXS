"""
@ author : Raphaël Aubry

Single-objective optimisation: minimise HX+combustor mass subject to wall-temperature,
pressure-drop, and Mach-number constraints.

Run with: python -m optimization.quick_optimize  (from project root, after pip install -e .)
"""

from scipy.optimize import minimize
import numpy as np

from hps_combustor.input_data import numericalProp, system_requirements
from hps_combustor.main_solve import main_solver
from .optimizer_desVar import design_variable_space
from .optimizer_dataMap import build_input_dataclasses


# Active variables and their bounds (continuous only, in insertion order)
_active = {name: props for name, props in design_variable_space.items()
           if props.get("activated", 0) == 1 and props["type"] == "continuous"}
bounds = [props["range"] for props in _active.values()]
x0     = [np.mean(b) for b in bounds]


def evaluate_design(x):
    coolant, hotgas, combustor, numerical, sysreqs = build_input_dataclasses(x, design_variable_space)
    numerical.radiation_ON = True
    numerical.debug_verbose = False

    try:
        solver = main_solver(
            coolantProp=coolant,
            hotgasProp=hotgas,
            combustorProp=combustor,
            numericalProp=numerical,
            system_requirements=sysreqs,
        )
        solver.solver()
        perf = solver.compute_performance()
    except Exception as e:
        print("Solver failed:", e)
        return 1e6, 5000.0, 1e6, 10.0

    return (
        perf["mass_HX_kg"] + perf["mass_combustor_kg"],
        perf["T_wg_max"],
        abs(perf["dp_c_bar"]) * 1e5,
        perf["Mach_g_max"],
    )


def objective(x):
    mass, *_ = evaluate_design(x)
    return mass

def T_constraint(x):
    _, T_wg, *_ = evaluate_design(x)
    return 950 - T_wg   # K, T_wg < 950 K

def dP_constraint(x):
    _, _, dP, _ = evaluate_design(x)
    return 5e5 - dP     # Pa, dP < 5 bar

def Mach_constraint(x):
    _, _, _, Mach = evaluate_design(x)
    return 0.5 - Mach   # Mach_g < 0.5


constraints = [
    {'type': 'ineq', 'fun': T_constraint},
    {'type': 'ineq', 'fun': dP_constraint},
    {'type': 'ineq', 'fun': Mach_constraint},
]

if __name__ == "__main__":
    result = minimize(
        fun=objective,
        x0=x0,
        bounds=bounds,
        constraints=constraints,
        method='trust-constr',
        options={'disp': True, 'maxiter': 200},
    )
    print("\nOptimal variables:")
    for name, val in zip(_active.keys(), result.x):
        print(f"  {name}: {val:.4g}")
    mass, T_wg, dP, Mach = evaluate_design(result.x)
    print(f"\nMass:       {mass:.3f} kg")
    print(f"T_wg_max:   {T_wg:.1f} K")
    print(f"dP_c:       {dP/1e5:.3f} bar")
    print(f"Mach_g_max: {Mach:.4f}")

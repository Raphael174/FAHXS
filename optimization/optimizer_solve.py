""" 
@ author : Raphaël Aubry
"""

#%%
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.soo.nonconvex.de import DE
from pymoo.optimize import minimize
from pymoo.termination import get_termination
from pymoo.core.callback import Callback
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.core.evaluator import ParallelEvaluator
import numpy as np
import multiprocessing
from .optimizer_desVar import design_variable_space
from hps_combustor.main_solve import main_solver
from .optimizer_dataMap import build_input_dataclasses

from hps_combustor.input_data import numericalProp

#%%

def evaluate_func(x):
    coolant, hotgas, combustor, toaster = build_input_dataclasses(x, design_variable_space)

    try:
        solver = main_solver(
            coolantProp=coolant,
            hotgasProp=hotgas,
            combustorProp=combustor,
            ToasterProp=toaster,
            numericalProp=numericalProp
        )
        solver.solver()
        solver.HX_sizing_brief(plotON=False, printON=False)

        return (
            solver.mass_tot,
            solver.T_fin_max,
            solver.Mach_g_max,
        )

    except Exception as e:
        print("Solver failed:", e)
        return 1e6, 1e4, 1e6, 10.0, 873  # fallback values

# ---- Pymoo Problem Class ----
class HeatExchangerProblem(ElementwiseProblem):
    def __init__(self, design_variable_space, evaluate_func):
        self.active_vars = {
            name: props
            for name, props in design_variable_space.items()
            if props.get("activated", 0) == 1
        }
        self.var_names = list(self.active_vars.keys())
        xl = [props["range"][0] for props in self.active_vars.values()]
        xu = [props["range"][1] for props in self.active_vars.values()]

        self.evaluate_func = evaluate_func

        super().__init__(
            n_var=len(self.var_names),
            n_obj=1,
            n_constr=2,
            xl=np.array(xl),
            xu=np.array(xu)
        )

    def _evaluate(self, x, out, *args, **kwargs):
        input_vars = {name: val for name, val in zip(self.var_names, x)}

        try:
            mass, T_wg, dP = self.evaluate_func(input_vars)
        except Exception as e:
            print("Error in evaluation:", e)
            mass, T_wg, dP = 1e6, 1e4, 1e6

        out["F"] = mass
        out["G"] = [800 - T_wg, 5e5 - dP]

# ---- Create Problem Instance ----
problem = HeatExchangerProblem(design_variable_space, evaluate_func)

# ---- Parallel Evaluation Setup ----
n_threads = multiprocessing.cpu_count() - 1
runner = ParallelEvaluator(n_procs=n_threads)

# ---- Differential Evolution Algorithm ----
algorithm = DE(
    pop_size=40,
    sampling=FloatRandomSampling(),
    variant="DE/rand/1/bin",
    CR=0.9,
    F=0.8,
    dither="vector",
    callback=None,
    eliminate_duplicates=True
)

# ---- Solve the Optimization Problem ----
res = minimize(
    problem,
    algorithm,
    termination=("n_gen", 50),
    seed=1,
    save_history=True,
    verbose=True,
    evaluator=runner
)

# ---- Display Results ----
print("\nOptimal Inputs:")
for name, val in zip(problem.var_names, res.X):
    print(f"{name}: {val:.4g}")

print("\nOutputs:")
mass, T_wg, dP = evaluate_func({name: val for name, val in zip(problem.var_names, res.X)})
print(f"Mass: {mass:.4f} kg")
print(f"Wall Temp: {T_wg:.2f} K")
print(f"Pressure Loss: {dP:.2f} Pa")

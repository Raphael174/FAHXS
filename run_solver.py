"""Compatibility wrapper for the steady runner.

Prefer:
    python -m hps_combustor.main_steady
"""

from hps_combustor.main_steady import build_inputs, run_steady
from hps_combustor.model_data_process.data_plotting import HXDashboard

solver, summary = run_steady(build_inputs())
db = HXDashboard(solver.data_master, coolant_name=solver.coolantProp.coolant)
db.all()       # thermal, coolant, combustion, mechanical, radiation, boiling (if applicable)


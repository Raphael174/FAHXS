"""Calibrate how Q_tot scales with hot-gas mass flow, for the shell-and-tube
geometry with helium coolant at a fixed representative coolant flow -
needed to estimate the hot-gas flow rate required for the 30 L/s steam
design-feasibility study (see STEAM_DESIGN_STUDY.md)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_case import run_shellntube

results = []
for mdot_hotgas in [0.10, 0.30, 0.60, 1.0]:
    solver = run_shellntube("helium", 0.10, mdot_hotgas=mdot_hotgas)
    r = dict(mdot_hotgas=mdot_hotgas, Q_tot_kW=float(solver.Q_tot / 1e3))
    results.append(r)
    print("POINT", json.dumps(r))

out = Path(__file__).parent / "hotgas_scaling_calibration.json"
out.write_text(json.dumps(results, indent=2))
print("SAVED", out)

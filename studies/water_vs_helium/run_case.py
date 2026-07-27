"""Water vs Helium comparison study - single-case runner.

Usage:
    python run_case.py --geometry helical --fluid water --mdot 0.10
    python run_case.py --geometry shellntube --fluid helium --mdot 0.05

Fixed for the whole study (see the study report for rationale):
  - T_in = 303.15 K, p_in = 80 bar, co-flow, finite_rate chemistry, hotgasProp() defaults.
  - Helical: combustorProp() defaults except HX_config/flow_config, with
    numericalProp.L_HX_max calibrated so the coil arc length is exactly 6.000 m
    (see L_HX_MAX_6M below - default geometry gives ~7.69 m, not 6 m).
  - Shell-and-tube: shellTubeProp() defaults (matches
    docs/context/shell_and_tube_architecture_target.png exactly), N_axial=200
    (production default), Phase 3 (p,h) shell-side liquid coupling.

Saves raw data (.npz) + a summary (.json) + a 4-panel plot (.png), all zipped
into a single labeled archive under raw_data/.
"""
import argparse
import json
import zipfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from hps_combustor.input_data import (
    CorrelationCoefficients, combustorProp, coolantProp, hotgasProp,
    numericalProp, shellTubeProp, system_requirements,
)
from hps_combustor.main_solve import main_solver
from hps_combustor.main_solve_shellntube import shellntube_solver

LUT_PATH = "docs/reference/external/2006LUTdata.txt"
T_IN = 303.15
P_IN = 80e5
# Calibrated via bisection (see study report) so the helical coil's arc
# length (L_ch_max, NOT the axial L_HX_max itself) is exactly 6.000 m.
L_HX_MAX_6M = 0.47922164350748064

STUDY_DIR = Path(__file__).parent
RAW_DIR = STUDY_DIR / "raw_data"
PLOT_DIR = STUDY_DIR / "plots"


def _coolant(fluid, mdot):
    if fluid == "helium":
        return coolantProp(coolant="Helium", coolant_model="single_phase_coolprop",
                            mass_flow_c=mdot, T_in=T_IN, p_in=P_IN, T_out=650, p_out=13e5)
    return coolantProp(coolant="Water", coolant_model="equilibrium_liquid",
                        mass_flow_c=mdot, T_in=T_IN, p_in=P_IN,
                        liquid_chf_lut_path=LUT_PATH)


def run_helical(fluid, mdot, mdot_hotgas=None):
    hg = hotgasProp(mass_flow_g=mdot_hotgas) if mdot_hotgas is not None else hotgasProp()
    solver = main_solver(
        coolantProp=_coolant(fluid, mdot),
        hotgasProp=hg,
        combustorProp=combustorProp(HX_config="shellnHelicalTube", flow_config="co"),
        numericalProp=numericalProp(L_HX_max=L_HX_MAX_6M, chemistry_model="finite_rate"),
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
    )
    solver.solver()
    solver.compute_performance()
    return solver


def run_shellntube(fluid, mdot, mdot_hotgas=None):
    hg = hotgasProp(mass_flow_g=mdot_hotgas) if mdot_hotgas is not None else hotgasProp()
    solver = shellntube_solver(
        coolantProp=_coolant(fluid, mdot),
        hotgasProp=hg,
        shellTubeProp=shellTubeProp(),
        numericalProp=numericalProp(chemistry_model="finite_rate"),
        system_requirements=system_requirements(),
        corrCoeffs=CorrelationCoefficients(),
        N_axial=200,
        flow_config="co",
    )
    solver.solve(verbose=True)
    return solver


def _arr(d, k):
    v = d.get(k) if hasattr(d, "get") else None
    return np.asarray(v, dtype=float) if v is not None else None


def make_plot(label, x, T_wg, T_wc, T_c, p_c, p_g, Res_g, Res_c, Res_w, out_path):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    fig.suptitle(label, fontweight="bold")

    ax = axes[0, 0]
    ax.plot(x, T_wg - 273.15, color="darkorange", label=r"$T_{wg}$")
    ax.plot(x, T_wc - 273.15, color="mediumblue", label=r"$T_{wc}$")
    ax.set_ylabel("Wall temperature [°C]"); ax.set_xlabel("x [m]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[0, 1]
    ax.plot(x, T_c - 273.15, color="cornflowerblue", label=r"$T_c$")
    ax.set_ylabel("Coolant temperature [°C]"); ax.set_xlabel("x [m]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1, 0]
    lines = []
    if p_c is not None:
        l1, = ax.plot(x, p_c / 1e5, color="steelblue", label=r"$p_c$")
        lines.append(l1)
    ax.set_ylabel("Coolant pressure [bar]"); ax.set_xlabel("x [m]")
    if p_g is not None:
        ax2 = ax.twinx()
        l2, = ax2.plot(x, p_g / 1e5, color="salmon", linestyle="--", label=r"$p_g$")
        ax2.set_ylabel("Hot gas pressure [bar]", color="salmon")
        lines.append(l2)
    ax.legend(handles=lines, fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1, 1]
    if Res_g is not None:
        ax.plot(x, Res_g, color="red", label=r"$R_{gas}$")
    if Res_c is not None:
        ax.plot(x, Res_c, color="cornflowerblue", label=r"$R_{coolant}$")
    if Res_w is not None and np.any(np.isfinite(Res_w)):
        ax.plot(x, Res_w, color="gray", label=r"$R_{wall}$")
    ax.set_yscale("log")
    ax.set_ylabel("Thermal resistance [K/W]"); ax.set_xlabel("x [m]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_case(geometry, fluid, mdot, solver):
    label = f"{geometry}_{fluid}_{int(round(mdot * 1000))}gs_co"

    if geometry == "helical":
        d = solver.data_master
        x = np.asarray(d["L_HX"], dtype=float)
        T_wg, T_wc, T_c = _arr(d, "T_wg"), _arr(d, "T_wc"), _arr(d, "T_c")
        p_c, p_g = _arr(d, "p_c"), _arr(d, "p_g")
        Res_g, Res_c, Res_w = _arr(d, "Res_g"), _arr(d, "Res_c"), _arr(d, "Res_w")
        npz_data = {k: np.asarray(v) for k, v in d.items()}
        summary = dict(
            geometry=geometry, fluid=fluid, mdot_kg_s=mdot,
            Q_tot_kW=float(solver.Q_tot), n_nodes=int(len(x)),
            T_wg_max_K=float(np.max(T_wg)), T_wc_max_K=float(np.max(T_wc)),
            T_c_in_K=float(T_c[0]), T_c_out_K=float(T_c[-1]),
            dp_c_bar=float(abs(p_c[-1] - p_c[0]) / 1e5),
        )
        if solver._liquid_mode:
            q, cm = _arr(d, "quality_c"), _arr(d, "chf_margin_c")
            finite_cm = cm[np.isfinite(cm)]
            summary.update(
                quality_min=float(np.min(q)), quality_max=float(np.max(q)),
                min_chf_margin=(float(np.min(finite_cm)) if finite_cm.size else None),
            )
    else:
        t = solver.tube
        x = np.arange(solver.N) * solver.dx
        T_wg, T_wc = _arr(t, "T_wg"), _arr(t, "T_wc")
        T_c = np.asarray(solver.T_shell, dtype=float)
        if solver.shell_liquid is not None and "p" in solver.shell_liquid:
            p_c = np.asarray(solver.shell_liquid["p"], dtype=float)  # real lagged march
        else:
            p_c = np.full_like(x, solver.coolantProp.p_in)  # gas mode - no shell dp/dx eq yet
        p_g = _arr(t, "p_g")
        h_g, h_c = _arr(t, "h_g"), _arr(t, "h_c")
        P_tube_i = np.pi * solver.D_tube_i
        P_tube_o = np.pi * solver.stp.D_tube_outer
        Res_g = 1.0 / np.maximum(h_g * P_tube_i * solver.dx, 1e-30)
        Res_c = 1.0 / np.maximum(h_c * P_tube_o * solver.dx, 1e-30)
        Res_w = np.full_like(x, np.nan)  # not tracked per-node for shell-and-tube
        npz_data = {k: np.asarray(v) for k, v in t.items()}
        npz_data["T_shell"] = T_c
        summary = dict(
            geometry=geometry, fluid=fluid, mdot_kg_s=mdot,
            Q_tot_kW=float(solver.Q_tot / 1e3), n_nodes=int(len(x)), n_sweeps=int(solver.n_sweeps),
            T_wg_max_K=float(np.max(T_wg)), T_wc_max_K=float(np.max(T_wc)),
            T_c_in_K=float(solver.coolantProp.T_in), T_c_out_K=float(solver.T_c_out),
        )
        if solver.shell_liquid is not None:
            q = np.asarray(solver.shell_liquid["quality"], dtype=float)
            cm = _arr(t, "chf_margin")
            finite_cm = cm[np.isfinite(cm)]
            summary.update(
                quality_min=float(np.min(q)), quality_max=float(np.max(q)),
                min_chf_margin=(float(np.min(finite_cm)) if finite_cm.size else None),
            )

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    plot_path = PLOT_DIR / f"{label}.png"
    make_plot(label, x, T_wg, T_wc, T_c, p_c, p_g, Res_g, Res_c, Res_w, plot_path)

    npz_path = RAW_DIR / f"{label}.npz"
    np.savez_compressed(npz_path, x=x, **npz_data)

    summary_path = RAW_DIR / f"{label}_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))

    zip_path = RAW_DIR / f"{label}.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(npz_path, npz_path.name)
        zf.write(summary_path, summary_path.name)
        zf.write(plot_path, plot_path.name)
    npz_path.unlink()

    print(f"SAVED {zip_path}")
    print("SUMMARY_JSON", json.dumps(summary))
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", choices=["helical", "shellntube"], required=True)
    ap.add_argument("--fluid", choices=["helium", "water"], required=True)
    ap.add_argument("--mdot", type=float, required=True)
    args = ap.parse_args()

    if args.geometry == "helical":
        solver = run_helical(args.fluid, args.mdot)
    else:
        solver = run_shellntube(args.fluid, args.mdot)
    save_case(args.geometry, args.fluid, args.mdot, solver)


if __name__ == "__main__":
    main()

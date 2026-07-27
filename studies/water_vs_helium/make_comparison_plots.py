"""Build cross-run comparison plots from the 12 case summaries."""
import json
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STUDY_DIR = Path(__file__).parent
RAW_DIR = STUDY_DIR / "raw_data"
PLOT_DIR = STUDY_DIR / "plots"

COLORS = {
    ("helical", "helium"): "firebrick",
    ("helical", "water"): "steelblue",
    ("shellntube", "helium"): "darkorange",
    ("shellntube", "water"): "mediumseagreen",
}
MARKERS = {"helical": "o", "shellntube": "s"}


def load_all():
    data = {}
    for f in sorted(glob.glob(str(RAW_DIR / "*_summary.json"))):
        d = json.load(open(f))
        key = (d["geometry"], d["fluid"])
        data.setdefault(key, []).append(d)
    for key in data:
        data[key].sort(key=lambda d: d["mdot_kg_s"])
    return data


def main():
    data = load_all()

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    fig.suptitle("Water vs Helium - cross-geometry comparison (co-flow, finite-rate chemistry)",
                 fontweight="bold")

    # Q_tot vs mdot_c
    ax = axes[0, 0]
    for (geom, fluid), runs in data.items():
        x = [r["mdot_kg_s"] * 1000 for r in runs]
        y = [r["Q_tot_kW"] for r in runs]
        ax.plot(x, y, marker=MARKERS[geom], color=COLORS[(geom, fluid)],
                label=f"{geom} / {fluid}")
    ax.set_xlabel("Coolant mass flow [g/s]"); ax.set_ylabel("Q_tot [kW]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("Absorbed duty")

    # T_wg_max vs mdot_c
    ax = axes[0, 1]
    for (geom, fluid), runs in data.items():
        x = [r["mdot_kg_s"] * 1000 for r in runs]
        y = [r["T_wg_max_K"] - 273.15 for r in runs]
        ax.plot(x, y, marker=MARKERS[geom], color=COLORS[(geom, fluid)],
                label=f"{geom} / {fluid}")
    ax.set_xlabel("Coolant mass flow [g/s]"); ax.set_ylabel("Max hot-side wall temperature [°C]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("Peak wall temperature")

    # CHF margin vs mdot_c (water only)
    ax = axes[1, 0]
    for (geom, fluid), runs in data.items():
        if fluid != "water":
            continue
        x = [r["mdot_kg_s"] * 1000 for r in runs]
        y = [r.get("min_chf_margin") for r in runs]
        ax.plot(x, y, marker=MARKERS[geom], color=COLORS[(geom, fluid)], label=f"{geom} / {fluid}")
    ax.axhline(1.0, color="red", linestyle="--", linewidth=1, label="CHF limit")
    ax.set_yscale("log")
    ax.set_xlabel("Coolant mass flow [g/s]"); ax.set_ylabel("Minimum CHF margin [-]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("CHF margin (water only)")

    # Exit quality vs mdot_c (water only)
    ax = axes[1, 1]
    for (geom, fluid), runs in data.items():
        if fluid != "water":
            continue
        x = [r["mdot_kg_s"] * 1000 for r in runs]
        y = [r.get("quality_max") for r in runs]
        ax.plot(x, y, marker=MARKERS[geom], color=COLORS[(geom, fluid)], label=f"{geom} / {fluid}")
    ax.axhline(1.0, color="red", linestyle=":", linewidth=1, label="x = 1 (dryout)")
    ax.set_xlabel("Coolant mass flow [g/s]"); ax.set_ylabel("Max quality reached [-]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("Peak exit quality (water only)")

    out = PLOT_DIR / "comparison_summary.png"
    fig.savefig(out, dpi=140)
    print("SAVED", out)


if __name__ == "__main__":
    main()

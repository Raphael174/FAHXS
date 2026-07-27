"""
run_grid_study.py
─────────────────────────────────────────────────────────────────────────────
2D grid sweep: N_arc_steps_per_turn (50–120, step 5)
             × Nusselt_correction   (0.01–0.26, step 0.05)

All other settings from the experimental test case (same as run_correlation_study.py):
  Co-flow, shellnHelicalTube, ahmed_toroid + Gnielinski
  Chamber 136 mm, coil 7/2.4 mm, O/F=2.04, 88 g/s diesel, 98 g/s He @ 30C 70 bar

Targets:
  T_c_out ≈ 420 °C  (693 K)
  Q_He    ≈ 200 kW

─────────────────────────────────────────────────────────────────────────────
Runtime estimate (equilibrium chemistry, single thread):
  90 runs × ~2400 steps × ~4 ms/step  ≈  14 min
  Set  FAST_MODE = True  to use frozen chemistry (no equilibrate in loop):
  90 runs × ~2400 steps × ~0.4 ms/step ≈  1.5 min
  Frozen mode changes absolute Q by ≲5 % but preserves trends perfectly.
─────────────────────────────────────────────────────────────────────────────
"""

# ── Package bootstrap ────────────────────────────────────────────────────────
if __name__ == "__main__" and __package__ is None:
    import sys, os, importlib.util, runpy
    _pkg_dir = os.path.dirname(os.path.abspath(__file__))
    _parent  = os.path.dirname(_pkg_dir)
    _alias   = "_hps"
    if _alias not in sys.modules:
        _spec = importlib.util.spec_from_file_location(
            _alias, os.path.join(_pkg_dir, "__init__.py"),
            submodule_search_locations=[_pkg_dir],
        )
        _pkg = importlib.util.module_from_spec(_spec)
        _pkg.__path__ = [_pkg_dir]
        _pkg.__package__ = _alias
        sys.modules[_alias] = _pkg
        _spec.loader.exec_module(_pkg)
    if _parent not in sys.path:
        sys.path.insert(0, _parent)
    runpy.run_module(f"{_alias}.run_grid_study", run_name="__main__", alter_sys=True)
    raise SystemExit(0)

from .main_solve import main_solver
from .input_data import (
    coolantProp, hotgasProp, combustorProp,
    numericalProp, system_requirements,
)
import numpy as np
import matplotlib.pyplot as plt
import time

# ── Toggle ───────────────────────────────────────────────────────────────────
FAST_MODE = False   # True = frozen chemistry (~1.5 min);  False = equilibrium (~14 min)

# ── Targets ──────────────────────────────────────────────────────────────────
T_TARGET_C = 420.0   # °C  He exit temperature
Q_TARGET   = 200.0   # kW  He absorbed power

# ── Sweep axes ───────────────────────────────────────────────────────────────
N_ARC_VALUES  = list(range(50, 125, 10))          # [50, 55, …, 120]  (15 values)
NU_CORR_VALUES = np.round(
    np.arange(0.1, 0.1, 0.05), 4).tolist()     # [0.01, 0.06, …, 0.26]  (6 values)

print(f"N_arc_steps values ({len(N_ARC_VALUES)}): {N_ARC_VALUES}")
print(f"Nusselt_correction ({len(NU_CORR_VALUES)}): {[f'{v:.2f}' for v in NU_CORR_VALUES]}")
print(f"Total runs: {len(N_ARC_VALUES) * len(NU_CORR_VALUES)}")
print(f"Chemistry: {'FROZEN (fast)' if FAST_MODE else 'EQUILIBRIUM (accurate)'}")
print()

# ── Fixed base case ──────────────────────────────────────────────────────────
def make_case(n_arc: int, nu_corr: float):
    cool = coolantProp(
        mass_flow_c = 0.098,
        T_in        = 303.15,   # 30 °C
        T_out       = 700.0,    # placeholder (co-flow)
        p_in        = 70e5,
        p_out       = 68e5,     # placeholder
    )
    hot = hotgasProp(
        mixing_ratio = 2.04,
        mass_flow_g  = 0.088,
    )
    comb = combustorProp(
        flow_config         = "co",
        inner_diameter      = 136e-3,
        mixing_length       = 50e-3,
        Dh_coil             = 7e-3,
        thickness_coil_wall = 2.4e-3,
        gap_shell2coil      = 22.5e-3,
        coil_gap            = 4e-3,
        material_HX         = "ST316L",
        material_CC         = "ST316L",
        Nusselt_shell       = "ahmed_toroid",
        Nusselt_coil        = "Gnielinski",
        Nusselt_correction  = nu_corr,
    )
    num = numericalProp(
        L_HX_max             = 0.517,
        N_arc_steps_per_turn = n_arc,
        radiation_ON         = True,
        chemistry_model      = "frozen" if FAST_MODE else "equilibrium",
        # suppress verbose per-node checks to reduce output
        check_energy_balance     = False,
        check_temperature_ordering = False,
        check_mach_limits        = False,
        check_stress_limits      = False,
        check_Re_regime          = False,
        check_Z_deviation        = False,
        debug_verbose            = False,
    )
    return cool, hot, comb, num, system_requirements()


# ── Results storage ──────────────────────────────────────────────────────────
NR = len(N_ARC_VALUES)
NC = len(NU_CORR_VALUES)

grid_T_c_out  = np.full((NR, NC), np.nan)   # He exit temperature [°C]
grid_Q_He     = np.full((NR, NC), np.nan)   # He absorbed power [kW]
grid_Q_tot    = np.full((NR, NC), np.nan)   # total heat transfer [kW]
grid_T_wg_max = np.full((NR, NC), np.nan)   # max hot-wall temperature [°C]
grid_ok       = np.zeros((NR, NC), dtype=bool)

t0 = time.time()
run_n = 0
total_runs = NR * NC

for i, n_arc in enumerate(N_ARC_VALUES):
    for j, nu_corr in enumerate(NU_CORR_VALUES):
        run_n += 1
        elapsed = time.time() - t0
        eta = (elapsed / run_n) * (total_runs - run_n) if run_n > 1 else 0
        print(f"  [{run_n:3d}/{total_runs}]  N_arc={n_arc:3d}  Nu_corr={nu_corr:.2f}"
              f"   elapsed {elapsed:5.0f}s   ETA {eta:5.0f}s", end="  ", flush=True)
        try:
            cool, hot, comb, num, sys_req = make_case(n_arc, nu_corr)
            s = main_solver(
                coolantProp=cool, hotgasProp=hot,
                combustorProp=comb, numericalProp=num,
                system_requirements=sys_req,
            )
            s.solver()
            m = s.compute_performance()

            grid_T_c_out [i, j] = float(s.data_master['T_c'][-1]) - 273.15
            grid_Q_He    [i, j] = m['Q_He_kW']
            grid_Q_tot   [i, j] = m['Q_tot_kW']
            grid_T_wg_max[i, j] = m['T_wg_max'] - 273.15
            grid_ok      [i, j] = True

            print(f"T_c_out={grid_T_c_out[i,j]:6.1f} °C   Q_He={grid_Q_He[i,j]:6.1f} kW")
        except Exception as exc:
            print(f"FAILED: {exc}")

total_time = time.time() - t0
print(f"\nTotal elapsed: {total_time:.1f} s  ({total_time/60:.1f} min)")


# ── 2D Tables ────────────────────────────────────────────────────────────────
def print_grid(title, grid, fmt, unit, targets=None):
    """targets: list of (row_i, col_j) of closest-to-target cells to mark."""
    W = 8 * NC + 36
    print(f"\n{'='*W}")
    print(f"  {title}  [{unit}]")
    print('='*W)
    hdr = f"  {'N_arc':>6s}  |" + "".join(f" {v:>7.2f}" for v in NU_CORR_VALUES)
    print(hdr)
    print(f"  {'':>6s}  | " + "  Nu_correction  →")
    print('-'*W)
    for i, n_arc in enumerate(N_ARC_VALUES):
        row = f"  {n_arc:>6d}  |"
        for j in range(NC):
            if np.isnan(grid[i, j]):
                row += f"  {'FAIL':>7s}"
            else:
                val = grid[i, j]
                cell = fmt.format(val)
                if targets and (i, j) in targets:
                    cell = f"[{cell}]"
                else:
                    cell = f" {cell} "
                row += f"{cell:>9s}"
        print(row)
    print('='*W)

# Find closest-to-target cells
def find_closest(grid, target):
    diff = np.abs(grid - target)
    diff[np.isnan(grid)] = np.inf
    idx = np.unravel_index(np.argmin(diff), grid.shape)
    return idx  # (row, col)

best_T = find_closest(grid_T_c_out, T_TARGET_C)
best_Q = find_closest(grid_Q_He,    Q_TARGET)

print_grid(
    f"He EXIT TEMPERATURE  (target ≈ {T_TARGET_C:.0f} °C)",
    grid_T_c_out, "{:6.1f}", "°C",
    targets=[best_T],
)
print_grid(
    f"He ABSORBED POWER    (target ≈ {Q_TARGET:.0f} kW)",
    grid_Q_He, "{:6.1f}", "kW",
    targets=[best_Q],
)
print_grid("HOT WALL MAX TEMPERATURE", grid_T_wg_max, "{:6.1f}", "°C")

# ── Best match report ─────────────────────────────────────────────────────────
print()
print("  Best match for T_c_out target:")
ni, nj = best_T
if not np.isnan(grid_T_c_out[ni, nj]):
    print(f"    N_arc={N_ARC_VALUES[ni]}, Nusselt_correction={NU_CORR_VALUES[nj]:.2f}")
    print(f"    T_c_out = {grid_T_c_out[ni, nj]:.1f} °C  (target {T_TARGET_C:.0f} °C)")
    print(f"    Q_He    = {grid_Q_He[ni, nj]:.1f} kW")

print()
print("  Best match for Q_He target:")
qi, qj = best_Q
if not np.isnan(grid_Q_He[qi, qj]):
    print(f"    N_arc={N_ARC_VALUES[qi]}, Nusselt_correction={NU_CORR_VALUES[qj]:.2f}")
    print(f"    Q_He    = {grid_Q_He[qi, qj]:.1f} kW  (target {Q_TARGET:.0f} kW)")
    print(f"    T_c_out = {grid_T_c_out[qi, qj]:.1f} °C")


# ── Plots ─────────────────────────────────────────────────────────────────────
nu_arr = np.array(NU_CORR_VALUES)
n_arr  = np.array(N_ARC_VALUES)
chem_tag = "frozen" if FAST_MODE else "equilibrium"

fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)
fig.suptitle(
    f"Grid study: N_arc_steps/turn × Nusselt_correction\n"
    f"ahmed_toroid + Gnielinski, co-flow 136mm, 7/2.4mm coil, diesel O/F=2.04, "
    f"98 g/s He @ 30°C 70 bar  [{chem_tag}]",
    fontsize=10, fontweight='bold',
)

# ── Heatmap helper ─────────────────────────────────────────────────────────────
def heatmap(ax, data, title, unit, target, cmap, fmt="{:.0f}"):
    im = ax.imshow(data, aspect='auto', origin='lower', cmap=cmap,
                   extent=[nu_arr[0]-0.025, nu_arr[-1]+0.025,
                           n_arr[0]-2.5,    n_arr[-1]+2.5])
    ax.set_xlabel("Nusselt_correction", fontsize=10)
    ax.set_ylabel("N_arc_steps / turn", fontsize=10)
    ax.set_xticks(nu_arr);  ax.set_xticklabels([f"{v:.2f}" for v in nu_arr], fontsize=7)
    ax.set_yticks(n_arr);   ax.set_yticklabels(n_arr, fontsize=7)
    ax.set_title(title, fontsize=10)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label(unit, fontsize=9)
    # Target contour
    try:
        cs = ax.contour(nu_arr, n_arr, data, levels=[target],
                        colors='red', linewidths=2)
        ax.clabel(cs, fmt=f"target={target:.0f}", fontsize=8)
    except Exception:
        pass
    # Cell values
    for i in range(NR):
        for j in range(NC):
            if not np.isnan(data[i, j]):
                ax.text(nu_arr[j], n_arr[i], fmt.format(data[i, j]),
                        ha='center', va='center', fontsize=6,
                        color='white' if abs(data[i, j] - target) < 0.15 * target else 'black')


heatmap(axes[0, 0], grid_T_c_out,  f"He exit temperature [°C]",        "°C", T_TARGET_C, "RdYlGn")
heatmap(axes[0, 1], grid_Q_He,     f"He absorbed power [kW]",           "kW", Q_TARGET,   "RdYlGn")
heatmap(axes[1, 0], grid_T_wg_max, f"Max hot-wall temperature [°C]",    "°C", 800.0,      "YlOrRd", fmt="{:.0f}")

# ── Line plot: T_c_out vs N_arc for each Nu_corr ─────────────────────────────
ax = axes[1, 1]
cmap_lines = plt.get_cmap('viridis')
for j, nu_corr in enumerate(NU_CORR_VALUES):
    color = cmap_lines(j / max(NC - 1, 1))
    mask  = grid_ok[:, j]
    if mask.any():
        ax.plot(np.array(N_ARC_VALUES)[mask], grid_T_c_out[mask, j],
                'o-', color=color, linewidth=1.5, markersize=4,
                label=f"Nu_c={nu_corr:.2f}")
ax.axhline(T_TARGET_C, color='red', linestyle='--', linewidth=1.5,
           label=f"Target {T_TARGET_C:.0f} °C")
ax.set_xlabel("N_arc_steps / turn", fontsize=10)
ax.set_ylabel("T_c_out [°C]", fontsize=10)
ax.set_title("He exit temperature vs N_arc (convergence check)", fontsize=10)
ax.legend(fontsize=7, ncol=2)
ax.grid(linewidth=0.5, color='#DDDDDD')
ax.set_xlim(N_ARC_VALUES[0] - 3, N_ARC_VALUES[-1] + 3)

plt.show()

"""
run_diameter_pressure_study.py
-----------------------------------------------------------------------------
Helium hydraulic-diameter -> pressure-drop sweep, under a PRESCRIBED linear
coolant temperature profile (90 K -> 650 K).

Goal
----
Deliver a coil HX that heats helium 90 K -> ~650 K while keeping the coolant
pressure drop below 10 bar.  The hot-gas Nusselt number is unvalidated, so we
sidestep it entirely: prescribe the helium temperature rise and ask the purely
mechanical question — *which coil hydraulic diameters keep dp under budget?*

As helium heats, its density collapses and it accelerates, so small diameters
choke on friction + acceleration losses.  This script sweeps Dh and reports the
admissible band.  With --friction-error N, it also overlays +/-N% friction
uncertainty bands to show how the admissible limit shifts.

Standalone: does NOT modify the core solver.  `python main_solve.py` unchanged.
Run with:
  python run_diameter_pressure_study.py --mdot 75
  python run_diameter_pressure_study.py --mdot 150 --friction-error 20
-----------------------------------------------------------------------------
"""

# -- Package bootstrap --------------------------------------------------------
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
    runpy.run_module(f"{_alias}.run_diameter_pressure_study", run_name="__main__", alter_sys=True)
    raise SystemExit(0)

import argparse
import numpy as np
import matplotlib.pyplot as plt

from .mechanical.geometry.helix_geometry import HelixGeometryRadiusCST
from .physics.pressure_drop_prescribedT import march_pressure_prescribedT
from .input_data import combustorProp, numericalProp, CorrelationCoefficients

# -- CLI overrides -------------------------------------------------------------
def _parse_args():
    p = argparse.ArgumentParser(description="He diameter -> pressure-drop sweep")
    p.add_argument("--mdot",           type=float, default=96.0,  help="mass flow [g/s]")
    p.add_argument("--T-in",           type=float, default=30+273,  help="He inlet temperature [K]")
    p.add_argument("--T-out",          type=float, default=420+273, help="He outlet temperature [K]")
    p.add_argument("--p-in",           type=float, default=70,  help="He inlet pressure [bar]")
    p.add_argument("--dp-limit",       type=float, default=10.0,  help="pressure-drop budget [bar]")
    p.add_argument("--N-ch",           type=int,   default=1,     help="number of parallel coil channels")
    p.add_argument("--Dh-min",         type=float, default=5,   help="min hydraulic diameter [mm]")
    p.add_argument("--Dh-max",         type=float, default=14,  help="max hydraulic diameter [mm]")
    p.add_argument("--Dh-step",        type=float, default=0.5,   help="diameter step [mm]")
    p.add_argument("--friction-error", type=float, default=20.0,
                   help="+/-%% friction-correlation prefactor uncertainty (e.g. 20 for +/-20%%)")
    p.add_argument("--rough-lo",       type=float, default=1.5,
                   help="optimistic inner roughness [um] (ideal drawn tube)")
    p.add_argument("--rough-nom",      type=float, default=15.0,
                   help="nominal inner roughness [um] (cold-drawn / worked steel)")
    p.add_argument("--rough-hi",       type=float, default=50.0,
                   help="pessimistic inner roughness [um] (welded/machined/bent steel)")
    p.add_argument("--no-plot",        action="store_true",       help="skip matplotlib figure")
    return p.parse_args()

_args = _parse_args()

# -- Study configuration (CLI-overridable) ------------------------------------
T_IN     = _args.T_in
T_OUT    = _args.T_out
P_IN     = _args.p_in     * 1e5   # bar -> Pa
MDOT     = _args.mdot     * 1e-3  # g/s -> kg/s
N_CH     = _args.N_ch
DP_LIMIT = _args.dp_limit * 1e5   # bar -> Pa
FERR     = _args.friction_error / 100.0   # fractional (0.20 for +/-20 %)
EPS_LO   = _args.rough_lo  * 1e-6   # um -> m
EPS_NOM  = _args.rough_nom * 1e-6
EPS_HI   = _args.rough_hi  * 1e-6

DH_RANGE = np.arange(_args.Dh_min * 1e-3,
                      _args.Dh_max * 1e-3 + 1e-9,
                      _args.Dh_step * 1e-3)

comb = combustorProp()
num  = numericalProp()
cc   = CorrelationCoefficients()


def build_geometry(Dh):
    """Return (D_coil, Rc, L_pipe) for a given coil inner diameter, or None if infeasible."""
    wall       = comb.thickness_coil_wall
    coil_pitch = Dh + 2 * wall + comb.coil_gap
    D_coil     = comb.inner_diameter - 2 * comb.gap_shell2coil - Dh - 2 * wall
    if D_coil <= 0:
        return None
    Rc     = D_coil / 2 * (1 + (coil_pitch / (np.pi * D_coil)) ** 2)
    L_coil = (num.L_HX_max - comb.mixing_length) - 2 * comb.length_2_coil - (Dh + 2 * wall)
    if L_coil <= 0:
        return None
    _, _, L_pipe = HelixGeometryRadiusCST(coil_pitch=coil_pitch, D_coil=D_coil, L_coil=L_coil)
    return D_coil, Rc, L_pipe


def run_sweep(friction_error_factor=1.0, roughness_height=EPS_LO):
    """Run one full Dh sweep at a given friction-prefactor factor and inner
    roughness.  Returns list of result dicts."""
    rows = []
    for Dh in DH_RANGE:
        geom = build_geometry(Dh)
        if geom is None:
            continue
        D_coil, Rc, L_pipe = geom
        r = march_pressure_prescribedT(
            Dh=Dh, mdot=MDOT, N_ch=N_CH, T_in=T_IN, T_out=T_OUT, p_in=P_IN,
            D_coil=D_coil, Rc=Rc, L_pipe=L_pipe,
            roughness=comb.channel_roughness,
            friction_selector=comb.friction_coil,
            corrCoeffs=cc,
            friction_error_factor=friction_error_factor,
            roughness_height=roughness_height,
        )
        r["Dh"] = Dh
        rows.append(r)
    return rows


def admissible_Dh(rows):
    """Return smallest Dh [m] that passes the dp budget, or None."""
    passing = [r for r in rows if r["valid"] and r["dp_total"] <= DP_LIMIT]
    if not passing:
        return None
    return min(passing, key=lambda r: r["Dh"])


def print_sweep(rows, label="nominal"):
    """Print the detailed table for one sweep."""
    header = (f"{'Dh':>6} {'turns':>6} {'L_pipe':>7} {'Re_in':>9} {'Re_out':>9} "
              f"{'U_out':>8} {'Mach':>6} {'dp_fric':>8} {'dp_acc':>8} {'dp_tot':>8}  ")
    units  = (f"{'[mm]':>6} {'[-]':>6} {'[m]':>7} {'[-]':>9} {'[-]':>9} "
              f"{'[m/s]':>8} {'[-]':>6} {'[bar]':>8} {'[bar]':>8} {'[bar]':>8}")
    print(f"\n-- {label} --")
    print(header); print(units); print("-" * 104)
    for r in rows:
        Dh = r["Dh"]
        if not r["valid"]:
            print(f"{Dh*1e3:6.1f} {r['n_turns']:6.1f} {r['L_pipe']:7.3f} "
                  f"{r['Re_in']:9.0f} {'--':>9} {'--':>8} {'--':>6} "
                  f"{'--':>8} {'--':>8} {'>'+f'{P_IN/1e5:.0f}':>7}  CHOKED "
                  f"(reached {r['frac_reached']*100:.0f}% of coil)")
            continue
        flag = "PASS" if r["dp_total"] <= DP_LIMIT else "fail"
        mw   = " (!Mach)" if r["Mach_out"] > 0.3 else ""
        print(f"{Dh*1e3:6.1f} {r['n_turns']:6.1f} {r['L_pipe']:7.3f} "
              f"{r['Re_in']:9.0f} {r['Re_out']:9.0f} {r['U_out']:8.1f} "
              f"{r['Mach_out']:6.3f} {r['dp_friction']/1e5:8.3f} {r['dp_accel']/1e5:8.3f} "
              f"{r['dp_total']/1e5:8.3f}  {flag}{mw}")
    print("-" * 104)
    best = admissible_Dh(rows)
    if best:
        marker = " (!Mach)" if best["Mach_out"] > 0.3 else "  (subsonic)"
        print(f" Admissible: Dh >= {best['Dh']*1e3:.1f} mm  "
              f"(dp={best['dp_total']/1e5:.2f} bar, Mach_out={best['Mach_out']:.3f}{marker})")
    else:
        print(f" No diameter in range meets the {DP_LIMIT/1e5:.0f} bar budget.")


def dp_curve(rows):
    """Return (Dh_mm, dp_bar) arrays over the valid points of a sweep."""
    v = [r for r in rows if r["valid"]]
    return (np.array([r["Dh"] * 1e3 for r in v]),
            np.array([r["dp_total"] / 1e5 for r in v]))


def main():
    print("=" * 104)
    print(f" Helium dp sweep  |  T: {T_IN:.0f}->{T_OUT:.0f} K  |  "
          f"mdot={MDOT*1e3:.0f} g/s  |  p_in={P_IN/1e5:.0f} bar  |  N_ch={N_CH}")
    print(f" Friction: {comb.friction_coil} (Darcy)  |  budget dp < {DP_LIMIT/1e5:.0f} bar")
    print(f" Roughness eps [um]:  lo={EPS_LO*1e6:.1f} (drawn)  nom={EPS_NOM*1e6:.1f} (worked steel)"
          f"  hi={EPS_HI*1e6:.1f} (welded/bent)   |   correlation prefactor +/-{FERR*100:.0f}%")
    print("=" * 104)

    # ---- scenarios: (label, eps, friction_factor) -----------------------------
    # roughness-only axis at nominal correlation
    rows_nom  = run_sweep(1.0,        EPS_NOM)        # central estimate
    rows_smth = run_sweep(1.0,        EPS_LO)         # smooth tube
    rows_rgh  = run_sweep(1.0,        EPS_HI)         # rough tube
    # global envelope corners (roughness AND correlation combined)
    rows_best = run_sweep(1.0 - FERR, EPS_LO)         # optimistic: smooth + f low
    rows_wrst = run_sweep(1.0 + FERR, EPS_HI)         # pessimistic: rough  + f high

    # full per-Dh table only for the nominal scenario (keeps output readable)
    print_sweep(rows_nom, label=f"nominal  (eps={EPS_NOM*1e6:.0f} um, f x1.00)")

    # ---- Dh_min summary across scenarios --------------------------------------
    print("\n" + "=" * 78)
    print(" Minimum admissible Dh  (smallest Dh with dp <= budget)")
    print("=" * 78)
    print(f" {'scenario':<40} {'eps[um]':>8} {'f x':>6} {'Dh_min[mm]':>11} {'dp[bar]':>8}")
    print("-" * 78)
    scen = [
        ("optimistic  (smooth + low friction)",  EPS_LO,  1.0 - FERR, rows_best),
        ("smooth tube",                          EPS_LO,  1.0,        rows_smth),
        ("NOMINAL",                              EPS_NOM, 1.0,        rows_nom),
        ("rough tube",                           EPS_HI,  1.0,        rows_rgh),
        ("pessimistic (rough + high friction)",  EPS_HI,  1.0 + FERR, rows_wrst),
    ]
    for label, eps, ff, rows in scen:
        best = admissible_Dh(rows)
        if best:
            print(f" {label:<40} {eps*1e6:>8.1f} {ff:>6.2f} "
                  f"{best['Dh']*1e3:>11.1f} {best['dp_total']/1e5:>8.2f}")
        else:
            print(f" {label:<40} {eps*1e6:>8.1f} {ff:>6.2f} {'none':>11} {'--':>8}")
    print("-" * 78)
    bn, bw = admissible_Dh(rows_best), admissible_Dh(rows_wrst)
    if bn and bw:
        print(f" GLOBAL ENVELOPE on Dh_min:  {bn['Dh']*1e3:.1f} mm (best) ... "
              f"{bw['Dh']*1e3:.1f} mm (worst)  -> design to >= {bw['Dh']*1e3:.1f} mm")
    print("=" * 78 + "\n")

    # -- plot -----------------------------------------------------------------
    if _args.no_plot:
        return

    Dh_nom, dp_nom   = dp_curve(rows_nom)
    Dh_b,   dp_best  = dp_curve(rows_best)
    Dh_w,   dp_wrst  = dp_curve(rows_wrst)
    Dh_s,   dp_smth  = dp_curve(rows_smth)
    Dh_r,   dp_rgh   = dp_curve(rows_rgh)
    if Dh_nom.size == 0:
        return
    mach_nom = np.array([r["Mach_out"] for r in rows_nom if r["valid"]])

    fig, ax = plt.subplots(figsize=(9.2, 5.8))

    # global envelope band (best -> worst), aligned on the worst-case grid
    bmap = {round(d, 4): v for d, v in zip(Dh_b, dp_best)}
    common = [d for d in Dh_w if round(d, 4) in bmap]
    if common:
        xc = np.array(common)
        lo = np.array([bmap[round(d, 4)] for d in common])
        hi = np.array([dp_wrst[list(Dh_w).index(d)] for d in common])
        ax.fill_between(xc, lo, hi, color="C3", alpha=0.15,
                        label=f"global envelope (eps {EPS_LO*1e6:.0f}-{EPS_HI*1e6:.0f}um, "
                              f"f +/-{FERR*100:.0f}%)")

    ax.plot(Dh_s, dp_smth, "--", color="C0", lw=1.1, alpha=0.8,
            label=f"smooth tube (eps={EPS_LO*1e6:.0f}um)")
    ax.plot(Dh_nom, dp_nom, "o-", color="C3", lw=2,
            label=f"nominal (eps={EPS_NOM*1e6:.0f}um)")
    ax.plot(Dh_r, dp_rgh, "--", color="C1", lw=1.3, alpha=0.9,
            label=f"rough tube (eps={EPS_HI*1e6:.0f}um)")
    ax.plot(Dh_w, dp_wrst, ":", color="darkred", lw=1.3,
            label=f"worst case (rough + f+{FERR*100:.0f}%)")
    ax.axhline(DP_LIMIT / 1e5, color="k", ls=":", lw=1.5, label=f"budget {DP_LIMIT/1e5:.0f} bar")

    # admissible region governed by the worst case
    if common:
        below = hi <= DP_LIMIT / 1e5
        if below.any():
            ax.fill_between(xc, 0, DP_LIMIT / 1e5 * 6, where=below,
                            color="green", alpha=0.07, label="admissible (worst-case)")

    ax.set_xlabel("coil inner hydraulic diameter  Dh  [mm]")
    ax.set_ylabel("helium pressure drop  [bar]")
    ax.set_title(f"He {T_IN:.0f}->{T_OUT:.0f} K, {MDOT*1e3:.0f} g/s, p_in={P_IN/1e5:.0f} bar, "
                 f"single coil  |  Reprise Essai GOx/Diesel 22 mai 26")
    ax.set_ylim(0, DP_LIMIT / 1e5 * 4)
    ax.set_xlim(Dh_nom.min(), Dh_nom.max())
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)

    ax2 = ax.twinx()
    ax2.plot(Dh_nom, mach_nom, ".-", color="C2", alpha=0.5, lw=1)
    ax2.set_ylabel("outlet Mach [-] (nominal)", color="C2")
    ax2.tick_params(axis="y", labelcolor="C2")
    ax2.set_ylim(0, max(0.35, float(np.nanmax(mach_nom)) * 1.1))
    ax2.axhline(0.3, color="C2", ls=":", lw=0.8, alpha=0.6)

    fig.tight_layout()
    out_png = f"diameter_pressure_study_mdot{MDOT*1e3:.0f}gs_rough_fricEr{FERR*100:.0f}_EssaiGOxDiesel22mai26.png"
    fig.savefig(out_png, dpi=130)
    print(f" Figure saved to {out_png}")
    plt.show()


if __name__ == "__main__":
    main()

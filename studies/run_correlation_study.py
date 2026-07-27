"""
run_correlation_study.py
────────────────────────────────────────────────────────────────────────────
Sweep all shell × coil Nusselt combinations for a 200 kW experimental case.
Prints a comparison table sorted by error, a resistance breakdown, and
proposes CorrelationCoefficients / Nusselt_correction tuning if needed.

Geometry / conditions:
  Co-flow, shellnHelicalTube
  Chamber ID       = 136 mm
  Coil ID / wall   = 7 mm / 2.4 mm  (OD = 11.8 mm)
  gap_shell2coil   = 22.5 mm  (outer coil face to inner shell wall)
  coil_gap         = 4 mm     (axial gap between turns)
  mixing_length    = 50 mm
  L_HX_max         = 517 mm  (= 445 mm HX + 50 mm mixing + 2x5 mm end clearance + 11.8 mm init offset)
  Fuel             = diesel-C16H34, O/F = 2.04, mdot_g = 88 g/s, p0 = 5 bar
  Helium           = mdot_c = 98 g/s, T_in = 30 C, p_in = 70 bar  (co-flow)
  Materials        = ST316L (HX + CC)
  N_arc_steps/turn = 50, radiation ON (Ehlme-2025 WSGGM)

Experimental reference: Q_He = 200 kW
"""

# ── Package bootstrap (folder '1Dmodel' is not a valid Python identifier) ───
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
    runpy.run_module(f"{_alias}.run_correlation_study", run_name="__main__", alter_sys=True)
    raise SystemExit(0)

from .main_solve import main_solver
from .input_data import (
    coolantProp, hotgasProp, combustorProp,
    numericalProp, system_requirements, CorrelationCoefficients,
)
import numpy as np
import matplotlib.pyplot as plt

# ────────────────────────────────────────────────────────────────────────────
Q_EXP_KW = 200.0   # kW — experimental reference

SHELL_CORRS = [
    "ahmed_toroid",
    "salimpour2008",
    "churchill_bernstein_tightcoil",
    "churchill_bernstein",
]
COIL_CORRS = ["mori1967", "Gnielinski"]

SHELL_LABELS = {
    "ahmed_toroid":                  "Ahmed-1997 toroid",
    "salimpour2008":                 "Salimpour-2008",
    "churchill_bernstein_tightcoil": "C-B tight-coil",
    "churchill_bernstein":           "Churchill-Bernstein",
}
COIL_LABELS = {
    "mori1967":    "Mori-1967",
    "Gnielinski":  "Gnielinski",
}
# ────────────────────────────────────────────────────────────────────────────


def make_case(nu_shell: str, nu_coil: str):
    cool = coolantProp(
        mass_flow_c = 0.098,
        T_in        = 303.15,   # 30 °C
        T_out       = 700.0,    # placeholder — not used for co-flow init
        p_in        = 70e5,
        p_out       = 68e5,     # placeholder
    )
    hot = hotgasProp(
        mixing_ratio = 2.04,
        mass_flow_g  = 0.088,   # 88 g/s
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
        Nusselt_shell       = nu_shell,
        Nusselt_coil        = nu_coil,
    )
    num = numericalProp(
        L_HX_max             = 0.517,   # 445 mm coil + 50 mm mixing + offsets (see header)
        N_arc_steps_per_turn = 50,
        radiation_ON         = True,
        chemistry_model      = "equilibrium",
    )
    return cool, hot, comb, num, system_requirements()


# ── Sweep ────────────────────────────────────────────────────────────────────
results  = []
solvers  = {}

for nu_sh in SHELL_CORRS:
    for nu_co in COIL_CORRS:
        sh_lbl = SHELL_LABELS[nu_sh]
        co_lbl = COIL_LABELS[nu_co]
        print(f"\n{'='*70}")
        print(f"  Shell: {sh_lbl}   x   Coil: {co_lbl}")
        print('='*70)
        try:
            cool, hot, comb, num, sys_req = make_case(nu_sh, nu_co)
            s = main_solver(coolantProp=cool, hotgasProp=hot,
                            combustorProp=comb, numericalProp=num,
                            system_requirements=sys_req)
            s.solver()
            m = s.compute_performance()
            d = s.data_master

            err_pct = (m['Q_He_kW'] - Q_EXP_KW) / Q_EXP_KW * 100.0

            # Mean thermal resistances for resistance-breakdown analysis
            Rg = float(np.mean(d['Res_g']))
            Rc = float(np.mean(d['Res_c']))
            Rw = float(np.mean(d['Res_w']))
            Rt = Rg + Rc + Rw

            results.append(dict(
                shell   = nu_sh, coil    = nu_co,
                sh_lbl  = sh_lbl, co_lbl = co_lbl,
                Q_He    = m['Q_He_kW'],
                Q_tot   = m['Q_tot_kW'],
                T_wg_max = m['T_wg_max'] - 273.15,
                T_g_out  = float(d['T_g'][-1]) - 273.15,
                T_c_out  = float(d['T_c'][-1]) - 273.15,
                eta      = m['eta_HX'],
                dp_c     = m['dp_c_bar'],
                Nu_g_mean = float(np.mean(d['Nu_g'])),
                Nu_c_mean = float(np.mean(d['Nu_c'])),
                h_g_mean  = float(np.mean(d['h_g_conv'])),
                h_c_mean  = float(np.mean(d['h_c'])),
                frac_g = Rg / Rt * 100,
                frac_c = Rc / Rt * 100,
                frac_w = Rw / Rt * 100,
                err = err_pct, ok = True,
            ))
            solvers[(nu_sh, nu_co)] = s
            print(f"  --> Q_He = {m['Q_He_kW']:.1f} kW  ({err_pct:+.1f}% vs {Q_EXP_KW:.0f} kW exp.)")

        except Exception as exc:
            import traceback
            traceback.print_exc()
            results.append(dict(
                shell=nu_sh, coil=nu_co,
                sh_lbl=sh_lbl, co_lbl=co_lbl,
                ok=False, err_msg=str(exc),
            ))


# ── Summary table ─────────────────────────────────────────────────────────────
ok_r  = sorted([r for r in results if r['ok']], key=lambda r: abs(r['err']))
W = 122
print(f"\n\n{'='*W}")
print(f"  CORRELATION SWEEP   |   Target: Q_He = {Q_EXP_KW:.0f} kW  (experiment)")
print('='*W)
print(
    f"  {'Shell':28s}  {'Coil':14s}  "
    f"{'Q_He':>7s}  {'Q_tot':>7s}  {'Err%':>7s}  "
    f"{'T_wg':>7s}  {'T_g_out':>8s}  {'T_c_out':>7s}  "
    f"{'eta':>5s}  {'dP_c':>6s}  "
    f"{'Nu_g':>6s}  {'Nu_c':>6s}  "
    f"{'Rg%':>5s}  {'Rc%':>5s}  {'Rw%':>5s}"
)
print('-'*W)
for r in ok_r:
    flag = "  <--" if abs(r['err']) <= 10 else ""
    print(
        f"  {r['sh_lbl']:28s}  {r['co_lbl']:14s}  "
        f"{r['Q_He']:>7.1f}  {r['Q_tot']:>7.1f}  {r['err']:>+7.1f}%  "
        f"{r['T_wg_max']:>7.1f}  {r['T_g_out']:>8.1f}  {r['T_c_out']:>7.1f}  "
        f"{r['eta']:>5.3f}  {r['dp_c']:>6.3f}  "
        f"{r['Nu_g_mean']:>6.1f}  {r['Nu_c_mean']:>6.1f}  "
        f"{r['frac_g']:>5.1f}  {r['frac_c']:>5.1f}  {r['frac_w']:>5.1f}"
        f"{flag}"
    )
for r in [r for r in results if not r['ok']]:
    print(f"  {r['sh_lbl']:28s}  {r['co_lbl']:14s}  FAILED: {r['err_msg'][:50]}")
print('='*W)
print(
    f"  {'Units':28s}  {'':14s}  "
    f"{'[kW]':>7s}  {'[kW]':>7s}  {'':>7s}  "
    f"{'[C]':>7s}  {'[C]':>8s}  {'[C]':>7s}  "
    f"{'':>5s}  {'[bar]':>6s}  "
    f"{'[-]':>6s}  {'[-]':>6s}  "
    f"{'[%]':>5s}  {'[%]':>5s}  {'[%]':>5s}"
)
print('='*W)
print("  Sorted by |error| ascending.  Rg/Rc/Rw = fraction of total thermal resistance.")
print("  <-- = within +/-10 % of experiment.")


# ── Tuning proposals ──────────────────────────────────────────────────────────
if ok_r:
    best = ok_r[0]
    print(f"\n{'─'*W}")
    print("  TUNING ANALYSIS")
    print(f"{'─'*W}")
    print(f"  Best match : {best['sh_lbl']}  x  {best['co_lbl']}")
    print(f"  Q_He       = {best['Q_He']:.1f} kW  ({best['err']:+.1f} % vs {Q_EXP_KW:.0f} kW)")
    print()
    print(f"  Resistance breakdown (mean over HX length):")
    print(f"    Hot-gas side  Rg = {best['frac_g']:5.1f} %   h_g_conv = {best['h_g_mean']:6.0f} W/m2/K")
    print(f"    Wall          Rw = {best['frac_w']:5.1f} %")
    print(f"    Coolant side  Rc = {best['frac_c']:5.1f} %   h_c      = {best['h_c_mean']:6.0f} W/m2/K")
    print()

    if abs(best['err']) <= 5:
        print("  Best correlation is within +/-5 % — no tuning required.")
    else:
        ratio = Q_EXP_KW / best['Q_He']  # > 1 if model under-predicts
        direction = "under-predicts" if ratio > 1 else "over-predicts"
        print(f"  Model {direction} by {abs(best['err']):.1f} %.  Required scaling: x{ratio:.3f}")
        print()

        # Identify dominant resistance and suggest the right knob
        dom_name, dom_pct = max(
            [("gas-side", best['frac_g']),
             ("wall",     best['frac_w']),
             ("coolant",  best['frac_c'])],
            key=lambda x: x[1],
        )
        print(f"  Dominant thermal resistance: {dom_name} ({dom_pct:.1f} %)")
        print()

        if dom_name == "gas-side":
            print("  --> Gas-side limited: tune Nu_shell.")
            print()
            # Option A: Nusselt_correction (works for all correlations)
            print(f"  Option A  combustorProp.Nusselt_correction = {ratio:.4f}")
            print(f"            (universal multiplier on top of any Nu_shell correlation)")
            print()

            if best['shell'] == "salimpour2008":
                a_new = 0.317 * ratio
                print(f"  Option B  CorrelationCoefficients(salimpour_a={a_new:.4f})")
                print(f"            (default salimpour_a = 0.317; scales Nu_shell ~ salimpour_a)")
            elif best['shell'] == "ahmed_toroid":
                print(f"  Option B  ahmed_toroid has no dedicated corrCoeffs prefactor.")
                print(f"            Use Nusselt_correction (Option A).")
            elif "churchill_bernstein" in best['shell']:
                print(f"  Option B  Churchill-Bernstein has no dedicated corrCoeffs prefactor.")
                print(f"            Use Nusselt_correction (Option A).")

        elif dom_name == "coolant":
            print("  --> Coolant-side limited: tune Nu_coil.")
            print()
            if best['coil'] == "mori1967":
                # Mori Nu ~ 1/mori_a_lo, so to increase Nu by ratio: decrease a_lo by ratio
                a_new = 26.2 / ratio
                print(f"  Option A  CorrelationCoefficients(mori_a_lo={a_new:.2f})")
                print(f"            (default mori_a_lo = 26.2; Nu_coil inversely proportional)")
            print(f"  Option B  numericalProp.artificial_error_Nu_cold = {ratio - 1.0:+.4f}")
            print(f"            (additive relative error on Nu_c; 0 = nominal)")

        elif dom_name == "wall":
            print("  --> Wall-resistance limited: Nu tuning has limited effect.")
            print("  Check wall thickness (2.4 mm) and ST316L conductivity (~15-24 W/m/K).")
            print("  At 2.4 mm wall on 7 mm ID, R_wall is significant — this is expected.")

        print()

        # Combined suggestion for the best case
        if abs(best['err']) > 30:
            print(f"  NOTE: Error is {abs(best['err']):.0f} % — very large mismatch.")
            print("  Before tuning coefficients, verify:")
            print("    - gap_shell2coil = 22.5 mm  (outer coil face to shell inner wall)")
            print("    - Dh_coil = 7 mm  (inner diameter), thickness = 2.4 mm")
            print("    - mass_flow_g = 0.088 kg/s,  O/F = 2.04")
            print("    - p_combustor and T_g_init (affects hot-gas properties)")
            print()

        print("  Ready-to-run tuned configuration:")
        print()
        corr_line = f"Nusselt_correction={ratio:.4f},"
        if best['shell'] == "salimpour2008":
            a_new = 0.317 * ratio
            corr_extra = f"\n        # CorrelationCoefficients(salimpour_a={a_new:.4f}) is equivalent"
        else:
            corr_extra = ""

        print(f"    solver = main_solver(")
        print(f"        coolantProp = coolantProp(mass_flow_c=0.098, T_in=303.15, p_in=70e5),")
        print(f"        hotgasProp  = hotgasProp(mixing_ratio=2.04, mass_flow_g=0.088),")
        print(f"        combustorProp = combustorProp(")
        print(f"            flow_config='co', inner_diameter=136e-3, mixing_length=50e-3,")
        print(f"            Dh_coil=7e-3, thickness_coil_wall=2.4e-3,")
        print(f"            gap_shell2coil=22.5e-3, coil_gap=4e-3,")
        print(f"            material_HX='ST316L', material_CC='ST316L',")
        print(f"            Nusselt_shell='{best['shell']}', Nusselt_coil='{best['coil']}',")
        print(f"            {corr_line}{corr_extra}")
        print(f"        ),")
        print(f"        numericalProp = numericalProp(L_HX_max=0.495, N_arc_steps_per_turn=50),")
        print(f"        system_requirements = system_requirements(),")
        print(f"    )")
    print()


# ── Bar chart ─────────────────────────────────────────────────────────────────
if ok_r:
    q_vals   = np.array([r['Q_He']  for r in ok_r])
    err_vals = np.array([r['err']   for r in ok_r])
    frac_g   = np.array([r['frac_g'] for r in ok_r])
    frac_c   = np.array([r['frac_c'] for r in ok_r])
    frac_w   = np.array([r['frac_w'] for r in ok_r])
    labels   = [f"{r['sh_lbl']}\nx {r['co_lbl']}" for r in ok_r]
    bar_colors = [
        'tab:green'  if abs(e) <= 10 else
        'tab:orange' if abs(e) <= 25 else 'tab:steelblue'
        for e in err_vals
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    fig.suptitle(
        f"Correlation sweep  |  Target Q_He = {Q_EXP_KW:.0f} kW  (experiment)\n"
        f"Co-flow · 136 mm ID · 7/2.4 mm coil · O/F=2.04 · 88 g/s diesel · 98 g/s He @ 30°C 70 bar",
        fontsize=10, fontweight='bold',
    )

    # ---- (left) Q_He bar chart ----
    ax = axes[0]
    bars = ax.bar(range(len(labels)), q_vals, color=bar_colors, edgecolor='k', linewidth=0.7)
    ax.axhline(Q_EXP_KW, color='red', linewidth=2, linestyle='--',
               label=f"Experiment {Q_EXP_KW:.0f} kW")
    ax.axhspan(Q_EXP_KW * 0.90, Q_EXP_KW * 1.10, color='red', alpha=0.10, label="±10 %")
    for b, v, e in zip(bars, q_vals, err_vals):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5,
                f"{v:.0f}\n({e:+.0f}%)", ha='center', va='bottom', fontsize=7.5,
                linespacing=1.3)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Q_He [kW]", fontsize=11)
    ax.set_title("He absorbed power (sorted by |error|)", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(axis='y', linewidth=0.5, color='#DDDDDD')
    ax.set_ylim(0, max(q_vals) * 1.30)

    # ---- (middle) error horizontal bar ----
    ax = axes[1]
    err_colors = [
        'tab:green'  if abs(e) <= 10 else
        'tab:orange' if abs(e) <= 25 else 'tab:red'
        for e in err_vals
    ]
    ax.barh(range(len(labels)), err_vals, color=err_colors, edgecolor='k', linewidth=0.7)
    ax.axvline(0, color='k', linewidth=0.8)
    ax.axvspan(-10, 10, color='green', alpha=0.10, label="±10 %")
    for i, (e, v) in enumerate(zip(err_vals, q_vals)):
        off = 0.6 if e >= 0 else -0.6
        ax.text(e + off, i, f"{e:+.1f}%  ({v:.0f} kW)",
                va='center', ha='left' if e >= 0 else 'right', fontsize=7.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Error vs experiment [%]", fontsize=11)
    ax.set_title("Relative error", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(axis='x', linewidth=0.5, color='#DDDDDD')

    # ---- (right) stacked resistance breakdown ----
    ax = axes[2]
    x = np.arange(len(labels))
    w = 0.55
    b1 = ax.bar(x, frac_g, w, label="R_g (hot gas)", color='tomato',      edgecolor='k', linewidth=0.5)
    b2 = ax.bar(x, frac_w, w, bottom=frac_g,          label="R_w (wall)",  color='lightgray', edgecolor='k', linewidth=0.5)
    b3 = ax.bar(x, frac_c, w, bottom=frac_g + frac_w, label="R_c (He)",   color='cornflowerblue', edgecolor='k', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("% of total thermal resistance", fontsize=10)
    ax.set_title("Resistance breakdown", fontsize=9)
    ax.legend(fontsize=8)
    ax.set_ylim(0, 105)
    ax.grid(axis='y', linewidth=0.5, color='#DDDDDD')
    for xi, (g, w_, c) in enumerate(zip(frac_g, frac_w, frac_c)):
        ax.text(xi, g / 2,       f"{g:.0f}%", ha='center', va='center', fontsize=7, color='white', fontweight='bold')
        ax.text(xi, g + w_ / 2,  f"{w_:.0f}%", ha='center', va='center', fontsize=7, color='k')
        ax.text(xi, g + w_ + c / 2, f"{c:.0f}%", ha='center', va='center', fontsize=7, color='white', fontweight='bold')

    plt.show()

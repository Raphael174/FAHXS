"""
Bell-Delaware shell-side method for a segmentally-baffled shell-and-tube HX.

@ author : Raphaël Aubry  (shell-and-tube extension, WP1)

Computes the shell-side heat-transfer coefficient and pressure drop for the new
baffled shell-and-tube config (hot combustion gas inside straight tubes, coolant
in shell-side zig-zag cross-flow around segmental baffles — the EchTherm-style
geometry). See DESIGN_PLAN_shellntube_transient.md section 2.2.

    h_shell = h_ideal * Jc * Jl * Jb * Js * Jr

Formulation follows the standard, self-consistent Bell-Delaware presentation in
  R.W. Serth, "Process Heat Transfer: Principles and Applications" (2007), ch. 6,
which itself derives from Bell (1963) and the Delaware final report. The ideal
tube-bank j- and f-factor correlations here are the ones the J/R correction
factors were calibrated against — do NOT substitute a different bank correlation
(e.g. Zukauskas) without re-deriving the corrections.

Unit-testable in isolation: `bell_delaware_shell()` takes a geometry dict (from
`shelltube_geometry.compute_bell_delaware_geometry`) plus fluid Re/Pr/props and
returns h and the pressure-drop breakdown. A worked example in the module test
reproduces textbook J-factor magnitudes.
"""
import numpy as np

# ---------------------------------------------------------------------------
# Ideal tube-bank j and f correlation coefficients (Bell-Delaware / Serth Table 6.1).
# Keyed by layout; each entry is a list of (Re_max, a1, a2, a3, a4) rows with the
# applicable Reynolds-number ceiling. a (the j exponent modifier) uses a3,a4.
# ---------------------------------------------------------------------------
_J_COEFFS = {
    "triangular30": [
        (1e1,  1.4,   -0.667, 1.450, 0.519),
        (1e2,  1.36,  -0.657, 1.450, 0.519),
        (1e3,  0.593, -0.477, 1.450, 0.519),
        (1e4,  0.321, -0.388, 1.450, 0.519),
        (1e6,  0.321, -0.388, 1.450, 0.519),
    ],
    "square90": [
        (1e1,  0.970, -0.667, 1.187, 0.370),
        (1e2,  0.900, -0.631, 1.187, 0.370),
        (1e3,  0.408, -0.460, 1.187, 0.370),
        (1e4,  0.107, -0.266, 1.187, 0.370),
        (1e6,  0.370, -0.395, 1.187, 0.370),
    ],
    "rotated45": [
        (1e1,  1.550, -0.667, 1.930, 0.500),
        (1e2,  0.498, -0.656, 1.930, 0.500),
        (1e3,  0.730, -0.500, 1.930, 0.500),
        (1e4,  0.370, -0.396, 1.930, 0.500),
        (1e6,  0.370, -0.396, 1.930, 0.500),
    ],
}
_F_COEFFS = {
    "triangular30": [
        (1e1,  48.0,  -1.0,   7.00,  0.500),
        (1e2,  45.1,  -0.973, 7.00,  0.500),
        (1e3,  4.570, -0.476, 7.00,  0.500),
        (1e4,  0.486, -0.152, 7.00,  0.500),
        (1e6,  0.372, -0.123, 7.00,  0.500),
    ],
    "square90": [
        (1e1,  35.0,  -1.0,   6.30,  0.378),
        (1e2,  32.1,  -0.963, 6.30,  0.378),
        (1e3,  6.090, -0.602, 6.30,  0.378),
        (1e4,  0.0815, 0.022, 6.30,  0.378),
        (1e6,  0.391, -0.148, 6.30,  0.378),
    ],
    "rotated45": [
        (1e1,  32.0,  -1.0,   6.59,  0.520),
        (1e2,  26.2,  -0.913, 6.59,  0.520),
        (1e3,  3.500, -0.476, 6.59,  0.520),
        (1e4,  0.333, -0.136, 6.59,  0.520),
        (1e6,  0.303, -0.126, 6.59,  0.520),
    ],
}


def _pick_row(table, Re):
    for (Re_ceiling, c1, c2, c3, c4) in table:
        if Re <= Re_ceiling:
            return c1, c2, c3, c4
    return table[-1][1:]


def ideal_bank_j(layout, Re, pitch_ratio):
    """Colburn j-factor for an ideal tube bank (Bell-Delaware)."""
    a1, a2, a3, a4 = _pick_row(_J_COEFFS[layout], Re)
    a = a3 / (1.0 + 0.14 * Re ** a4)
    return a1 * (1.33 / pitch_ratio) ** a * Re ** a2


def ideal_bank_f(layout, Re, pitch_ratio):
    """Fanning-type friction factor for an ideal tube bank (Bell-Delaware)."""
    b1, b2, b3, b4 = _pick_row(_F_COEFFS[layout], Re)
    b = b3 / (1.0 + 0.14 * Re ** b4)
    return b1 * (1.33 / pitch_ratio) ** b * Re ** b2


# ---------------------------------------------------------------------------
# Correction factors
# ---------------------------------------------------------------------------
def J_c(Fc):
    """Baffle-configuration (window) correction. Fc = fraction of tubes in cross-flow.
    Jc ~ 1 for ~25% cut; ranges ~0.53 (large cut) to ~1.15 (small cut)."""
    return 0.55 + 0.72 * Fc


def J_l(r_lm, r_s):
    """Baffle-leakage correction. r_lm=(S_sb+S_tb)/S_m, r_s=S_sb/(S_sb+S_tb)."""
    p = 0.44 * (1.0 - r_s)
    return p + (1.0 - p) * np.exp(-2.2 * r_lm)


def J_b(Fsbp, r_ss, laminar=False):
    """Bundle-bypass correction. Fsbp = bypass area / S_m, r_ss = N_ss / N_tcc.
    C = 1.35 (laminar, Re<100) else 1.25."""
    if r_ss >= 0.5:
        return 1.0
    C = 1.35 if laminar else 1.25
    return np.exp(-C * Fsbp * (1.0 - (2.0 * r_ss) ** (1.0 / 3.0)))


def J_s(Nb, Lsi_ratio, Lso_ratio, n1):
    """Unequal end baffle-spacing correction.
    Lsi_ratio = Ls_inlet/Ls_central, Lso_ratio = Ls_outlet/Ls_central.
    n1 = 0.6 (turbulent) typically."""
    if Nb <= 1:
        return 1.0
    num = (Nb - 1) + Lsi_ratio ** (1 - n1) + Lso_ratio ** (1 - n1)
    den = (Nb - 1) + Lsi_ratio + Lso_ratio
    return num / den


def J_r(Re, N_c_total):
    """Laminar adverse-temperature-gradient correction. Only < Re~100; 1.0 above."""
    if Re >= 100:
        return 1.0
    Jr_lam = (10.0 / max(N_c_total, 1.0)) ** 0.18
    if Re <= 20:
        return Jr_lam
    # linear interpolate between Re=20 (Jr_lam) and Re=100 (1.0)
    return Jr_lam + (1.0 - Jr_lam) * (Re - 20.0) / 80.0


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------
def bell_delaware_shell(geom, Re_s, Pr_s, k_s, cp_s, mu_s, mdot_s,
                        mu_ratio=1.0, corrCoeffs=None):
    """
    Shell-side h [W/m2K] and pressure-drop breakdown via Bell-Delaware.

    Parameters
    ----------
    geom : dict from shelltube_geometry.compute_bell_delaware_geometry(), providing
        S_m, S_w, S_tb, S_sb, S_b, D_otl, N_tcc, N_tcw, N_baffles, Fc, Fsbp,
        r_ss, layout, pitch_ratio, D_tube_outer, Lsi_ratio, Lso_ratio, D_w_eq.
    Re_s : shell Reynolds = D_o * (mdot_s/S_m) / mu_s   (max-velocity basis).
    Pr_s, k_s, cp_s, mu_s : shell-fluid properties.
    mdot_s : shell mass flow [kg/s].
    mu_ratio : (mu_bulk/mu_wall) for the ^0.14 property correction (1.0 if unknown).
    corrCoeffs : optional; uses .bell_Jl_factor, .bell_Jb_factor, .zukauskas_C_factor
        (reused as an overall ideal-h prefactor knob) if present.

    Returns dict: h_shell, Jc, Jl, Jb, Js, Jr, h_ideal, dp_shell, and dp components.
    """
    layout = geom["layout"]; pr = geom["pitch_ratio"]; Do = geom["D_tube_outer"]
    Sm = geom["S_m"]
    G_s = mdot_s / Sm                      # max mass flux [kg/m2s]

    # --- ideal-bank heat transfer ---
    j = ideal_bank_j(layout, Re_s, pr)
    h_ideal = j * cp_s * G_s * Pr_s ** (-2.0 / 3.0) * mu_ratio ** 0.14
    if corrCoeffs is not None:
        h_ideal *= getattr(corrCoeffs, "zukauskas_C_factor", 1.0)

    # --- correction factors ---
    laminar = Re_s < 100
    Jc = J_c(geom["Fc"])
    r_lm = (geom["S_sb"] + geom["S_tb"]) / Sm
    r_s = geom["S_sb"] / max(geom["S_sb"] + geom["S_tb"], 1e-30)
    Jl = J_l(r_lm, r_s)
    Jb = J_b(geom["Fsbp"], geom["r_ss"], laminar=laminar)
    Js = J_s(geom["N_baffles"], geom["Lsi_ratio"], geom["Lso_ratio"], n1=0.6)
    N_c_total = (geom["N_tcc"] + geom["N_tcw"]) * (geom["N_baffles"] + 1)
    Jr = J_r(Re_s, N_c_total)
    if corrCoeffs is not None:
        Jl *= getattr(corrCoeffs, "bell_Jl_factor", 1.0)
        Jb *= getattr(corrCoeffs, "bell_Jb_factor", 1.0)

    h_shell = h_ideal * Jc * Jl * Jb * Js * Jr

    # --- pressure drop (three zones) ---
    f = ideal_bank_f(layout, Re_s, pr)
    rho_s = mu_s * Re_s / (Do * G_s / mu_s * mu_s)  # placeholder; caller may pass rho
    # dp per ideal cross-flow row group: dp_bk = 2 f (G_s^2/rho) N_tcc (mu ratio)^-0.14
    # rho must be supplied via geom['rho_s'] for accuracy:
    rho = geom.get("rho_s", None)
    if rho is None:
        rho = 1.0  # caller should set geom['rho_s']; avoids div-by-zero
    dp_ideal_cross = 2.0 * f * (G_s ** 2 / rho) * geom["N_tcc"] * mu_ratio ** (-0.14)

    # leakage/bypass corrections on dp (Rl, Rb) — standard Bell-Delaware forms
    Rl = np.exp(-1.33 * (1.0 + r_s) * r_lm ** (0.8 - 0.15 * (1.0 + r_s)))
    Cb = 4.5 if laminar else 3.7
    Rb = np.exp(-Cb * geom["Fsbp"] * (1.0 - (2.0 * geom["r_ss"]) ** (1.0 / 3.0))
                if geom["r_ss"] < 0.5 else 0.0)
    Rs = (geom["Lsi_ratio"]) ** (-1.0) + (geom["Lso_ratio"]) ** (-1.0)  # end-zone factor (approx)

    Nb = geom["N_baffles"]
    # window dp per baffle
    Sw = geom["S_w"]
    G_w = mdot_s / np.sqrt(Sm * Sw)
    dp_window_ideal = (2.0 + 0.6 * geom["N_tcw"]) * G_w ** 2 / (2.0 * rho)

    dp_cross = dp_ideal_cross * (Nb - 1) * Rb * Rl
    dp_window = dp_window_ideal * Nb * Rl
    dp_ends = dp_ideal_cross * (1.0 + geom["N_tcw"] / max(geom["N_tcc"], 1e-9)) * Rb * Rs
    dp_shell = dp_cross + dp_window + dp_ends

    return dict(h_shell=float(h_shell), h_ideal=float(h_ideal),
                Jc=float(Jc), Jl=float(Jl), Jb=float(Jb), Js=float(Js), Jr=float(Jr),
                j_ideal=float(j), f_ideal=float(f),
                dp_shell=float(dp_shell), dp_cross=float(dp_cross),
                dp_window=float(dp_window), dp_ends=float(dp_ends),
                Rl=float(Rl), Rb=float(Rb))


if __name__ == "__main__":
    # Smoke test: J-factor magnitudes for a representative geometry (25% cut, p/D=1.25).
    # Reference sanity: Jc~1.0 at Fc~0.65, Jl~0.7-0.8, Jb~0.9, Js~1, Jr=1 (turbulent).
    geom = dict(layout="triangular30", pitch_ratio=1.25, D_tube_outer=19.05e-3,
                S_m=0.0344, S_w=0.0129, S_tb=0.00358, S_sb=0.00198, S_b=0.00518,
                D_otl=0.483, N_tcc=9, N_tcw=4, N_baffles=14, Fc=0.65, Fsbp=0.15,
                r_ss=0.0, Lsi_ratio=1.0, Lso_ratio=1.0, rho_s=1000.0)
    r = bell_delaware_shell(geom, Re_s=1.0e4, Pr_s=6.0, k_s=0.6, cp_s=4180.0,
                            mu_s=1.0e-3, mdot_s=30.0)
    print("Bell-Delaware smoke test (25% cut, p/D=1.25, Re=1e4):")
    for k in ("h_ideal", "Jc", "Jl", "Jb", "Js", "Jr", "h_shell", "dp_shell"):
        print(f"  {k:9s} = {r[k]:.4g}")
    assert 0.9 < r["Jc"] < 1.15, r["Jc"]
    assert 0.6 < r["Jl"] < 0.95, r["Jl"]
    assert 0.8 < r["Jb"] <= 1.0, r["Jb"]
    assert r["Jr"] == 1.0
    print("  J-factor ranges OK.")

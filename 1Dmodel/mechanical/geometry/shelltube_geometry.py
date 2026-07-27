"""
Bell-Delaware geometry for the baffled shell-and-tube config.

@ author : Raphaël Aubry  (shell-and-tube extension, WP1)

Computes all the flow areas and tube-count factors the Bell-Delaware shell-side
method (physics/bell_delaware.py) needs, from the physical inputs on the EchTherm
GEOMETRY screen (shell ID, tube OD/pitch/layout, baffle cut/spacing/clearances).
Formulas follow R.W. Serth, "Process Heat Transfer" (2007), ch. 6.

Single entry point: compute_bell_delaware_geometry(...) -> dict consumed directly
by bell_delaware_shell(geom, ...).
"""
import numpy as np


def _pitch_parallel(layout, Pt):
    """Tube pitch parallel to cross-flow (row-to-row spacing in the flow direction)."""
    if layout == "triangular30":
        return Pt * np.sqrt(3.0) / 2.0
    if layout == "rotated45":
        return Pt / np.sqrt(2.0)
    return Pt  # square90


def _pitch_normal(layout, Pt):
    """Tube pitch normal to cross-flow."""
    if layout == "rotated45":
        return Pt * np.sqrt(2.0)
    return Pt  # triangular30, square90


def estimate_tube_count(D_shell_inner, D_tube_outer, pitch_ratio, layout,
                        bundle_shell_clearance=0.0):
    """Rough tube count from bundle circle packing (Serth's CTP/CL area method)."""
    Pt = pitch_ratio * D_tube_outer
    CL = 0.87 if layout in ("triangular30", "rotated45") else 1.0  # layout constant
    CTP = 0.93  # one tube pass
    D_otl = D_shell_inner - bundle_shell_clearance
    A_bundle = np.pi / 4.0 * D_otl ** 2
    Nt = CTP * A_bundle / (CL * Pt ** 2)
    return max(int(round(Nt)), 1)


def compute_bell_delaware_geometry(
    D_shell_inner, D_tube_outer, pitch_ratio, layout,
    N_tubes, N_baffles, baffle_cut, L_tube,
    clearance_tube_baffle=0.8e-3,
    clearance_baffle_shell=1.6e-3,
    clearance_bundle_shell=0.0,
    N_sealing_strip_pairs=0,
    baffle_spacing=None, L_inlet_spacing=None, L_outlet_spacing=None,
):
    """
    Returns the Bell-Delaware geometry dict.

    Parameters (SI)
    ----------
    D_shell_inner : shell inside diameter [m]
    D_tube_outer  : tube outside diameter [m]
    pitch_ratio   : Pt / D_o
    layout        : "triangular30" | "square90" | "rotated45"
    N_tubes       : number of tubes
    N_baffles     : number of baffles
    baffle_cut    : fraction of D_shell (e.g. 0.20 for 20% cut)
    L_tube        : tube length [m]
    clearance_*   : diametral clearances [m]
    N_sealing_strip_pairs : sealing strip pairs
    baffle_spacing : central baffle spacing [m] (default = L_tube/(N_baffles+1))
    L_inlet_spacing, L_outlet_spacing : end baffle spacings [m] (default = central)
    """
    Ds = D_shell_inner
    Do = D_tube_outer
    Pt = pitch_ratio * Do
    Bc = baffle_cut
    Pp = _pitch_parallel(layout, Pt)
    Pn = _pitch_normal(layout, Pt)

    # central and end baffle spacings. EchTherm exposes these separately from
    # front/rear end-zone lengths, so keep them independent when provided.
    B = baffle_spacing if baffle_spacing else L_tube / (N_baffles + 1)
    Lsi = L_inlet_spacing if L_inlet_spacing else B
    Lso = L_outlet_spacing if L_outlet_spacing else B

    # bundle outer-tube-limit diameter and centreline tube-limit diameter
    D_otl = Ds - clearance_bundle_shell
    Dctl = D_otl - Do

    # --- crossflow area at the shell centreline (S_m) ---
    S_m = B * ((Ds - D_otl) + (D_otl - Do) / Pn * (Pt - Do))

    # --- fraction of tubes in cross-flow (Fc) and window (Fw) ---
    # angle argument, clipped to valid arccos domain
    arg = (Ds - 2.0 * Bc * Ds) / Dctl
    arg = float(np.clip(arg, -1.0, 1.0))
    theta_ctl = 2.0 * np.arccos(arg)
    Fc = (1.0 / np.pi) * (np.pi + 2.0 * arg * np.sin(np.arccos(arg)) - 2.0 * np.arccos(arg))
    Fc = float(np.clip(Fc, 0.0, 1.0))
    Fw = (1.0 - Fc) / 2.0

    # --- window area (gross minus tube-occupied) ---
    theta_ds = 2.0 * np.arccos(1.0 - 2.0 * Bc)
    S_wg = (Ds ** 2 / 4.0) * (0.5 * (theta_ds - np.sin(theta_ds)))  # gross window
    S_wt = N_tubes * Fw * (np.pi / 4.0) * Do ** 2                    # occupied by tubes
    S_w = max(S_wg - S_wt, 1e-9)

    # --- leakage areas ---
    Ltb = clearance_tube_baffle
    Lsb = clearance_baffle_shell
    S_tb = (np.pi / 4.0) * ((Do + Ltb) ** 2 - Do ** 2) * N_tubes * (1.0 - Fw)
    S_sb = Ds * Lsb / 2.0 * (np.pi - 0.5 * theta_ds) / np.pi * np.pi  # = Ds*Lsb/2*(pi - 0.5 theta_ds)
    S_sb = Ds * Lsb / 2.0 * (np.pi - 0.5 * theta_ds)

    # --- bundle bypass area ---
    S_b = B * (Ds - D_otl)
    Fsbp = S_b / S_m

    # --- tube rows crossed ---
    N_tcc = int(round(Ds * (1.0 - 2.0 * Bc) / Pp))              # in one cross-flow section
    N_tcw = int(round(0.8 / Pp * (Bc * Ds - (Ds - Dctl) / 2.0)))  # in one window
    N_tcw = max(N_tcw, 0)
    r_ss = N_sealing_strip_pairs / max(N_tcc, 1)

    # window equivalent hydraulic diameter (for window-zone dp, optional)
    D_w_eq = 4.0 * S_w / (np.pi * Do * N_tubes * Fw + np.pi * Ds * theta_ds / (2 * np.pi))

    return dict(
        layout=layout, pitch_ratio=pitch_ratio, D_tube_outer=Do,
        S_m=float(S_m), S_w=float(S_w), S_tb=float(S_tb), S_sb=float(S_sb), S_b=float(S_b),
        D_otl=float(D_otl), N_tcc=int(N_tcc), N_tcw=int(N_tcw), N_baffles=int(N_baffles),
        Fc=float(Fc), Fw=float(Fw), Fsbp=float(Fsbp), r_ss=float(r_ss),
        Lsi_ratio=float(Lsi / B), Lso_ratio=float(Lso / B),
        B_central=float(B), L_inlet_spacing=float(Lsi), L_outlet_spacing=float(Lso),
        D_w_eq=float(D_w_eq),
    )


if __name__ == "__main__":
    # EchTherm-style reference case (from the screenshot geometry)
    Nt = estimate_tube_count(D_shell_inner=110e-3, D_tube_outer=5e-3,
                             pitch_ratio=1.3, layout="triangular30")
    print(f"estimated tube count (Ø110 shell, Ø5 tube, p/D=1.3, tri): {Nt}  (EchTherm shows 235)")
    geom = compute_bell_delaware_geometry(
        D_shell_inner=110e-3, D_tube_outer=5e-3, pitch_ratio=1.3, layout="triangular30",
        N_tubes=235, N_baffles=15, baffle_cut=0.20, L_tube=235e-3,
        clearance_tube_baffle=1e-3, clearance_baffle_shell=1e-3, clearance_bundle_shell=0.0,
        baffle_spacing=12e-3, L_inlet_spacing=8e-3, L_outlet_spacing=8e-3)
    print("Bell-Delaware geometry for the EchTherm case:")
    for k, v in geom.items():
        print(f"  {k:12s} = {v}")
    # sanity checks
    assert geom["S_m"] > 0 and geom["S_w"] > 0
    assert 0 < geom["Fc"] < 1
    assert geom["S_tb"] > 0 and geom["S_sb"] > 0
    assert geom["N_tcc"] > 0
    print("geometry sanity OK.")

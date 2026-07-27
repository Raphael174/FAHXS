"""
Digitized curve scaffolding for Yang (2005) Ch.2 — Jet Injector loss coefficients
(Figs. 22–25). This file defines:

• Axis specs (min, max, scale) for each figure exactly as plotted.
• Empty numpy arrays you can paste (x,y) digitized points into.
• Interpolator helpers (linear or log-x) that enforce axis ranges.
• Public functions that return ξ-values from the digitized curves.

Intended use: digitize the curves with your preferred tool (e.g. WebPlotDigitizer),
then paste the (x,y) arrays below. Keep x strictly increasing in the appropriate
space (linear for Figs. 23–25; logarithmic-x for Fig. 22).

All units SI/deg as labeled on the figures. No external dependencies beyond numpy.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, Optional
import numpy as np

# ------------------------------
# Axis specifications
# ------------------------------



""" 
Fig 22 : 
{
"x": [108.0751,122.7183,134.6602,162.6355,201.941,239.4472,286.2645,344.0926,429.5709,550.9923,662.2979,817.923,1054.8077,1405.1916,1871.9654,2349.6736,3013.8261,4058.6437,5436.175,7440.5811,10294.8671,13420.9138,17214.4335,24471.3418,31388.346,39612.1474,50261.8618,64120.8084,75423.3702,88239.4206,96742.339],
"y": [0.8786,0.8368,0.8013,0.7442,0.691,0.6484,0.6136,0.5753,0.5333,0.489,0.4626,0.4278,0.3931,0.3607,0.3235,0.296,0.2684,0.2396,0.2109,0.1833,0.1557,0.1341,0.1186,0.0934,0.079,0.0634,0.0478,0.0346,0.0275,0.0203,0.0167]
}

Fig 23 
lin/di=0.6
{
"x": [0.2823,0.8709,1.6068,2.6369,3.6671,4.8445,6.169,7.935,9.701,12.0557,14.2632,17.3537,21.4744,26.1838,31.0404,37.2214,42.078,45.3157,48.259,51.791,55.3231,60.1796,64.1532,69.1569,73.8662,79.3115,84.4623,89.6132,93.4396,97.8546,103.0055,108.0092,112.1299,117.2808,123.0203,129.4957,134.7938,139.5032,146.5672,165.8462,180.7102],
"y": [0.4935,0.4634,0.4356,0.4093,0.3836,0.358,0.3317,0.3002,0.2738,0.2431,0.216,0.1889,0.1633,0.1391,0.123,0.1106,0.1062,0.1054,0.1069,0.1084,0.1113,0.1171,0.123,0.134,0.1464,0.1625,0.1786,0.1962,0.2094,0.224,0.2401,0.2562,0.2709,0.287,0.3046,0.3273,0.3448,0.3602,0.3822,0.4444,0.4913]
}
lin/di=0.15
{
"x": [0.5766,3.8143,5.8746,7.935,10.7312,12.6444,15.8821,18.6783,21.3273,24.7121,27.9498,31.7762,36.0441,40.1648,44.4326,48.7005,52.2326,55.9117,60.0325,63.7117,68.8625,72.5417,76.9568,81.3718,84.9038,87.8472,91.9679,96.6773,101.0923,106.3904,110.8054,117.5751,124.1977,139.5032,150.9823,180.563],
"y": [0.4949,0.4407,0.4137,0.3873,0.3573,0.3346,0.3053,0.276,0.2526,0.2277,0.2065,0.1867,0.1699,0.1559,0.1486,0.145,0.1435,0.1435,0.1464,0.1508,0.1589,0.1662,0.1786,0.1896,0.1984,0.2072,0.2196,0.2306,0.2445,0.2606,0.2767,0.2987,0.3199,0.3697,0.4049,0.4905]
}
lin/di=0.1
{
"x": [0.4294,4.9916,7.6407,10.7312,14.5576,18.5311,22.2103,25.5951,29.2743,34.1309,38.1044,42.078,45.61,48.8477,52.2326,56.2061,60.0325,63.2701,67.0965,71.3644,75.4851,78.1341,81.8133,86.0812,90.9377,96.5301,101.8282,109.3337,115.2204,120.3713,127.2882,134.7938,142.7408,149.8049,180.4158],
"y": [0.4949,0.4444,0.4159,0.388,0.3551,0.3243,0.2936,0.2724,0.2519,0.2292,0.2174,0.2028,0.1926,0.1874,0.1808,0.1757,0.1742,0.1757,0.1801,0.1882,0.1955,0.2043,0.2138,0.2248,0.2365,0.2548,0.2702,0.2914,0.3097,0.3243,0.3441,0.3661,0.3902,0.4093,0.4913]
}

Fig 24
{
"x": [0.0183,0.0551,0.092,0.1329,0.178,0.2251,0.2926,0.3766,0.4564,0.5322,0.5957,0.651,0.737,0.8189,0.9008,0.9868,1.0748,1.169,1.257,1.3451,1.4249,1.5171,1.601,1.6706,1.7423,1.7935],
"y": [0.4929,0.4674,0.4399,0.4113,0.3827,0.349,0.3133,0.2725,0.2429,0.2164,0.199,0.1857,0.1653,0.147,0.1306,0.1174,0.1051,0.0939,0.0847,0.0755,0.0694,0.0633,0.0572,0.0531,0.0469,0.0449]
}
"""


"""
Simple, readable digitized curves for Yang (2005) — Jet Injector loss coefficients
(Figs. 22–25). Uses only the Python standard library. Paste x,y arrays and get
interpolated values with linear (or log-x) interpolation.

Provided by user (high priority data):
- Fig. 22: ξ_{1→c}(Re) — full dataset
- Fig. 23: ξ_in(β) for l_in/d_i = 0.60, 0.15, 0.10 — full datasets
- Fig. 24: ξ_in(τ/d_i) — full dataset
Fig. 25 currently contains placeholder points (until user supplies data).
"""

from typing import List, Tuple
import math

# ------------------------------
# Helpers (simple & readable)
# ------------------------------

def clamp(v: float, vmin: float, vmax: float) -> float:
    if v < vmin:
        return vmin
    if v > vmax:
        return vmax
    return v


def interp1_linear(x: float, xs: List[float], ys: List[float]) -> float:
    """Linear 1D interpolation. Assumes xs strictly increasing and same length as ys.
    Clamps x to [xs[0], xs[-1]]."""
    n = len(xs)
    if n == 0 or n != len(ys):
        raise ValueError("interp1_linear: invalid input lengths")
    if n == 1:
        return ys[0]
    # clamp
    x = clamp(x, xs[0], xs[-1])
    # find interval
    lo, hi = 0, n - 1
    # small linear search is fine for modest n; keeps code clear
    for i in range(1, n):
        if x <= xs[i]:
            lo, hi = i - 1, i
            break
    x0, x1 = xs[lo], xs[hi]
    y0, y1 = ys[lo], ys[hi]
    if x1 == x0:
        return y0
    t = (x - x0) / (x1 - x0)
    return y0 * (1.0 - t) + y1 * t


def interp1_logx(x: float, xs: List[float], ys: List[float]) -> float:
    """Interpolation with logarithmic x-axis (base-10).
    Clamps x to [xs[0], xs[-1]]. All xs must be > 0."""
    if any(xx <= 0.0 for xx in xs) or x <= 0.0:
        raise ValueError("interp1_logx: x values must be positive")
    lxs = [math.log10(xx) for xx in xs]
    lx = math.log10(clamp(x, xs[0], xs[-1]))
    # reuse linear interpolator on log10(x)
    return interp1_linear(lx, lxs, ys)

# ------------------------------
# Axis ranges (for reference & optional validation)
# ------------------------------
RE_MIN, RE_MAX = 1.0e2, 1.0e5            # Fig. 22
BETA_MIN, BETA_MAX = 0.0, 180.0           # Fig. 23
TAU_DI_MIN, TAU_DI_MAX = 0.0, 1.8         # Fig. 24 (book ≈0–1.7)
ALPHA_MIN, ALPHA_MAX = 30.0, 85.0         # Fig. 25

# ------------------------------
# Data (user-supplied)
# ------------------------------
# Fig. 22 — ξ_{1→c}(Re)
RE_fig22 = [108.0751,122.7183,134.6602,162.6355,201.941,239.4472,286.2645,344.0926,429.5709,550.9923,662.2979,817.923,1054.8077,1405.1916,1871.9654,2349.6736,3013.8261,4058.6437,5436.175,7440.5811,10294.8671,13420.9138,17214.4335,24471.3418,31388.346,39612.1474,50261.8618,64120.8084,75423.3702,88239.4206,96742.339]
XI1C_fig22 = [0.8786,0.8368,0.8013,0.7442,0.691,0.6484,0.6136,0.5753,0.5333,0.489,0.4626,0.4278,0.3931,0.3607,0.3235,0.296,0.2684,0.2396,0.2109,0.1833,0.1557,0.1341,0.1186,0.0934,0.079,0.0634,0.0478,0.0346,0.0275,0.0203,0.0167]

# Fig. 23 — ξ_in(β) at fixed l_in/d_i
#   l_in/d_i = 0.60
BETA_fig23_060 = [0.2823,0.8709,1.6068,2.6369,3.6671,4.8445,6.169,7.935,9.701,12.0557,14.2632,17.3537,21.4744,26.1838,31.0404,37.2214,42.078,45.3157,48.259,51.791,55.3231,60.1796,64.1532,69.1569,73.8662,79.3115,84.4623,89.6132,93.4396,97.8546,103.0055,108.0092,112.1299,117.2808,123.0203,129.4957,134.7938,139.5032,146.5672,165.8462,180.7102]
XIIN_fig23_060 = [0.4935,0.4634,0.4356,0.4093,0.3836,0.358,0.3317,0.3002,0.2738,0.2431,0.216,0.1889,0.1633,0.1391,0.123,0.1106,0.1062,0.1054,0.1069,0.1084,0.1113,0.1171,0.123,0.134,0.1464,0.1625,0.1786,0.1962,0.2094,0.224,0.2401,0.2562,0.2709,0.287,0.3046,0.3273,0.3448,0.3602,0.3822,0.4444,0.4913]
#   l_in/d_i = 0.15
BETA_fig23_015 = [0.5766,3.8143,5.8746,7.935,10.7312,12.6444,15.8821,18.6783,21.3273,24.7121,27.9498,31.7762,36.0441,40.1648,44.4326,48.7005,52.2326,55.9117,60.0325,63.7117,68.8625,72.5417,76.9568,81.3718,84.9038,87.8472,91.9679,96.6773,101.0923,106.3904,110.8054,117.5751,124.1977,139.5032,150.9823,180.563]
XIIN_fig23_015 = [0.4949,0.4407,0.4137,0.3873,0.3573,0.3346,0.3053,0.276,0.2526,0.2277,0.2065,0.1867,0.1699,0.1559,0.1486,0.145,0.1435,0.1435,0.1464,0.1508,0.1589,0.1662,0.1786,0.1896,0.1984,0.2072,0.2196,0.2306,0.2445,0.2606,0.2767,0.2987,0.3199,0.3697,0.4049,0.4905]
#   l_in/d_i = 0.10
BETA_fig23_010 = [0.4294,4.9916,7.6407,10.7312,14.5576,18.5311,22.2103,25.5951,29.2743,34.1309,38.1044,42.078,45.61,48.8477,52.2326,56.2061,60.0325,63.2701,67.0965,71.3644,75.4851,78.1341,81.8133,86.0812,90.9377,96.5301,101.8282,109.3337,115.2204,120.3713,127.2882,134.7938,142.7408,149.8049,180.4158]
XIIN_fig23_010 = [0.4949,0.4444,0.4159,0.388,0.3551,0.3243,0.2936,0.2724,0.2519,0.2292,0.2174,0.2028,0.1926,0.1874,0.1808,0.1757,0.1742,0.1757,0.1801,0.1882,0.1955,0.2043,0.2138,0.2248,0.2365,0.2548,0.2702,0.2914,0.3097,0.3243,0.3441,0.3661,0.3902,0.4093,0.4913]

# Fig. 24 — ξ_in(τ/d_i)
TAU_over_DI_fig24 = [0.0183,0.0551,0.092,0.1329,0.178,0.2251,0.2926,0.3766,0.4564,0.5322,0.5957,0.651,0.737,0.8189,0.9008,0.9868,1.0748,1.169,1.257,1.3451,1.4249,1.5171,1.601,1.6706,1.7423,1.7935]
XIIN_fig24 = [0.4929,0.4674,0.4399,0.4113,0.3827,0.349,0.3133,0.2725,0.2429,0.2164,0.199,0.1857,0.1653,0.147,0.1306,0.1174,0.1051,0.0939,0.0847,0.0755,0.0694,0.0633,0.0572,0.0531,0.0469,0.0449]

# Fig. 25 — placeholder until provided (kept simple & monotonic)
ALPHA_deg_fig25 = [30, 40, 50, 60, 70, 80, 85]
XIIN_fig25 = [0.88, 0.80, 0.72, 0.66, 0.60, 0.54, 0.50]

# ------------------------------
# Evaluators
# ------------------------------

def xi1c_from_fig22(Re: float) -> float:
    """Return ξ_{1→c}(Re) using Fig. 22 with log-x interpolation."""
    # Guard domain
    Re = clamp(Re, RE_MIN, RE_MAX)
    return interp1_logx(Re, RE_fig22, XI1C_fig22)


def _interp_beta_on_curve(beta_deg: float, beta_curve: List[float], xi_curve: List[float]) -> float:
    # Clamp β to data range of this specific curve
    b = clamp(beta_deg, beta_curve[0], beta_curve[-1])
    return interp1_linear(b, beta_curve, xi_curve)


def xiin_from_fig23(beta_deg: float, lin_over_di: float) -> float:
    """Return ξ_in(β, l_in/d_i) via bilinear interpolation across the three curves.
    Clamps β to [min,max] of the data and clamps l_in/d_i to [0.10, 0.60]."""
    # Interpolate along each curve at requested β
    y_010 = _interp_beta_on_curve(beta_deg, BETA_fig23_010, XIIN_fig23_010)
    y_015 = _interp_beta_on_curve(beta_deg, BETA_fig23_015, XIIN_fig23_015)
    y_060 = _interp_beta_on_curve(beta_deg, BETA_fig23_060, XIIN_fig23_060)

    # Interpolate across l_in/d_i
    keys = [0.10, 0.15, 0.60]
    ys = [y_010, y_015, y_060]

    L = clamp(lin_over_di, keys[0], keys[-1])

    # find bracket
    if L <= keys[1]:
        # between 0.10 and 0.15
        k0, k1 = keys[0], keys[1]
        y0, y1 = ys[0], ys[1]
    else:
        # between 0.15 and 0.60
        k0, k1 = keys[1], keys[2]
        y0, y1 = ys[1], ys[2]

    t = (L - k0) / (k1 - k0)
    return y0 * (1.0 - t) + y1 * t


def xiin_from_fig24(tau_over_di: float) -> float:
    """Return ξ_in(τ/d_i) from Fig. 24 (linear interpolation)."""
    x = clamp(tau_over_di, TAU_DI_MIN, TAU_DI_MAX)
    return interp1_linear(x, TAU_over_DI_fig24, XIIN_fig24)


def xiin_from_fig25(alpha_deg):
    """Return ξ_in(α) from Fig. 25 (placeholder until real data)
    linear negative function
    """
    # y = mx + b
    m = (0.50-0.90)/(85-30)
    return m*(alpha_deg-30) + 0.9

    return interp1_linear(a, ALPHA_deg_fig25, XIIN_fig25)

# ------------------------------
# Callables to plug into the sizing code
# ------------------------------

def xi_1c_callable_from_fig22():
    return lambda Re: xi1c_from_fig22(float(Re))


def xi_in_callable_conical(beta_deg: float, lin_over_di: float):
    return lambda geom: xiin_from_fig23(beta_deg, lin_over_di)


def xi_in_callable_rounding(tau_over_di: float):
    return lambda geom: xiin_from_fig24(tau_over_di)


def xi_in_callable_tilt(alpha_deg: float):
    return lambda geom: xiin_from_fig25(alpha_deg)


if __name__ == "__main__":
    # quick prints to sanity-check
    print("Fig22 ξ1→c at Re=1e4:", xi1c_from_fig22(1.0e4))
    print("Fig23 ξin at β=40°, l_in/d_i=0.3:", xiin_from_fig23(40.0, 0.3))
    print("Fig24 ξin at τ/d_i=0.5:", xiin_from_fig24(0.5))
    print("Fig25 ξin at α=60° (placeholder):", xiin_from_fig25(60))

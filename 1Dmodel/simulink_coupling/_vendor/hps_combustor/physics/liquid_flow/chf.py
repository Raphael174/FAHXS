"""Critical heat flux (CHF) / dryout lookup for the liquid coolant path.

Independent of ``correlations.py``: this module only interpolates a supplied
CHF lookup table and applies a diameter correction, with no fluid-property
calls of its own.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np

GROENEVELD_2006_PRESSURES_MPA = np.array(
    [0.10, 0.30, 0.50, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 21.0],
    dtype=float,
)
GROENEVELD_2006_MASS_FLUXES = np.array(
    [0, 50, 100, 300, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000],
    dtype=float,
)
GROENEVELD_2006_QUALITIES = np.array(
    [-0.50, -0.40, -0.30, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00],
    dtype=float,
)


def chf_regime(quality: float) -> str:
    """Classify the CHF mechanism from local equilibrium quality.

    DNB (departure from nucleate boiling, low/negative quality) is a sudden
    collapse of the liquid film into vapor blanketing at the wall; dryout
    (high quality, annular-flow liquid-film depletion) is a gradual thinning
    of an already-established liquid film. The 2006 LUT does not need this
    distinction to interpolate a CHF value (quality is just a lookup axis),
    but the two mechanisms have different practical implications (DNB is a
    sudden, often more severe transition), so this is a diagnostic label, not
    a modeling branch. The x=0.1 cutoff is a fixed, simplified boundary for
    labeling purposes -- the real DNB/dryout transition quality is fluid-,
    pressure-, and mass-flux-dependent and is not a sharp line.
    """
    return "DNB" if quality < 0.1 else "dryout"


def local_chf_diameter_correction(q_chf_8mm_W_m2: float, diameter_m: float) -> float:
    """Groeneveld 2006 local CHF diameter correction from an 8 mm tube basis."""
    if q_chf_8mm_W_m2 <= 0.0 or diameter_m <= 0.0:
        raise ValueError("CHF and diameter must be positive")
    return q_chf_8mm_W_m2 * (diameter_m / 0.008) ** -0.5


def interpolate_chf_table(
    *,
    p_MPa: float,
    mass_flux_kg_m2_s: float,
    quality: float,
    pressures_MPa: np.ndarray,
    mass_fluxes_kg_m2_s: np.ndarray,
    qualities: np.ndarray,
    chf_kW_m2: np.ndarray,
) -> float:
    """Trilinear interpolation for a Groeneveld-style CHF lookup table.

    The full 2006 table is not encoded in the codebase yet. This helper makes
    the data dependency explicit and is validated with a small synthetic table.
    """
    axes = (
        np.asarray(pressures_MPa, dtype=float),
        np.asarray(mass_fluxes_kg_m2_s, dtype=float),
        np.asarray(qualities, dtype=float),
    )
    values = np.asarray(chf_kW_m2, dtype=float)
    if values.shape != tuple(len(axis) for axis in axes):
        raise ValueError("CHF table shape must match pressure, mass-flux, and quality axes")
    point = (p_MPa, mass_flux_kg_m2_s, quality)
    idx = []
    weights = []
    for axis, coordinate in zip(axes, point):
        if coordinate < axis[0] or coordinate > axis[-1]:
            raise ValueError("CHF lookup point is outside the supplied table")
        upper = int(np.searchsorted(axis, coordinate, side="right"))
        upper = min(max(upper, 1), len(axis) - 1)
        lower = upper - 1
        span = axis[upper] - axis[lower]
        w = 0.0 if span == 0.0 else (coordinate - axis[lower]) / span
        idx.append((lower, upper))
        weights.append(w)
    out = 0.0
    for i in (0, 1):
        wi = weights[0] if i else 1.0 - weights[0]
        for j in (0, 1):
            wj = weights[1] if j else 1.0 - weights[1]
            for k in (0, 1):
                wk = weights[2] if k else 1.0 - weights[2]
                out += wi * wj * wk * values[idx[0][i], idx[1][j], idx[2][k]]
    return out * 1000.0


def load_groeneveld_2006_lut(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load the open-source 2006 Groeneveld CHF LUT text table (cached).

    The table layout follows the MATLAB implementation from
    ``greenwoodms06/2006_Groeneveld_CriticalHeatFlux_LUT``:
    consecutive pressure blocks, each with mass-flux rows and quality columns.

    Results are cached per resolved path (``functools.lru_cache``) since a
    coupled march calls this once per node — re-parsing the text file on
    every node would otherwise be a per-node file-I/O cost (e.g. thousands of
    re-reads across a shell-and-tube sweep). Callers must treat the returned
    arrays as read-only; they are shared across calls.

    Returns
    -------
    pressures_MPa, mass_fluxes_kg_m2_s, qualities, chf_kW_m2
        ``chf_kW_m2`` has shape ``(pressure, mass_flux, quality)``.
    """
    return _load_groeneveld_2006_lut_cached(str(Path(path)))


@lru_cache(maxsize=8)
def _load_groeneveld_2006_lut_cached(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    raw = np.loadtxt(Path(path), dtype=float)
    n_p = len(GROENEVELD_2006_PRESSURES_MPA)
    n_g = len(GROENEVELD_2006_MASS_FLUXES)
    n_x = len(GROENEVELD_2006_QUALITIES)
    if raw.shape != (n_p * n_g, n_x):
        raise ValueError(
            f"expected Groeneveld LUT shape {(n_p * n_g, n_x)}, got {raw.shape}"
        )
    table = np.empty((n_p, n_g, n_x), dtype=float)
    for i_p in range(n_p):
        table[i_p, :, :] = raw[i_p * n_g : (i_p + 1) * n_g, :]
    return (
        GROENEVELD_2006_PRESSURES_MPA.copy(),
        GROENEVELD_2006_MASS_FLUXES.copy(),
        GROENEVELD_2006_QUALITIES.copy(),
        table,
    )


def groeneveld_2006_chf(
    *,
    p_Pa: float,
    mass_flux_kg_m2_s: float,
    quality: float,
    diameter_m: float,
    lut_path: str | Path,
) -> float:
    """Interpolate the 2006 Groeneveld LUT and apply the tube diameter correction."""
    p_axis, g_axis, x_axis, table = load_groeneveld_2006_lut(lut_path)
    q_8mm = interpolate_chf_table(
        p_MPa=p_Pa / 1.0e6,
        mass_flux_kg_m2_s=mass_flux_kg_m2_s,
        quality=quality,
        pressures_MPa=p_axis,
        mass_fluxes_kg_m2_s=g_axis,
        qualities=x_axis,
        chf_kW_m2=table,
    )
    return local_chf_diameter_correction(q_8mm, diameter_m)

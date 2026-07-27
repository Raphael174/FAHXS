"""
@ author : Raphaël Aubry
"""

# --- radiation_helpers.py ---
import numpy as np
from scipy import constants

sigma = constants.Stefan_Boltzmann  # 5.670374419e-8

def qrad_net_mbl(Tg, Ts, eps_g_emit, eps_g_abs, eps_s):
    """
    Net gas->surface radiative flux [W/m^2] using a compact two-gray-surface form.
    eps_g_emit : gas emissivity at Tg   (includes MBL via Le)
    eps_g_abs  : gas absorptivity at Ts (≈ emissivity at Ts with same Le)
    eps_s      : tube surface total hemispherical emissivity (0–1)
    """
    # effective gas "grayness" for multiple-reflection correction
    eps_g_eff = 0.5 * (np.clip(eps_g_emit, 1e-6, 0.999999) +
                       np.clip(eps_g_abs,  1e-6, 0.999999))
    eps_s = np.clip(eps_s, 1e-6, 0.999999)
    denom = (1.0/eps_g_eff) + (1.0/eps_s) - 1.0
    return sigma * (Tg**4 - Ts**4) / denom

def hrad_from_q(Tg, Ts, qpp):
    """Equivalent radiative h on the hot side [W/m^2-K] from a flux and ΔT."""
    dT = Tg - Ts
    if abs(dT) < 1e-6:  # avoid divide-by-zero; consistent secant linearization
        # small-ΔT fallback: local slope ~ derivative at mean temperature
        Tm = 0.5*(Tg + Ts)
        return max(0.0, 4.0 * sigma * (Tm**3))
    return max(0.0, qpp / dT)

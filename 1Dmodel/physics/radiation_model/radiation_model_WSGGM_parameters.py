"""
@ author : Raphaël Aubry 
"""

"""
#!
#! TO BE RUN ONCE TO EXTRACT PARAMETERS FOR H2O/CO2 RADIATION
#!

Ehlmé et al. (2025) — Updated WSGGM parameters up to 5000 K

This script tabulates the published coefficients for the 4-gray-gas WSGGM covering
(H2O, CO2) mixtures and the single-species limits, and writes them to JSON files.

It also provides small helper functions to evaluate κ_j(MR) and a_j(T, MR).

Notation (matches the paper’s eqs. 3–5):
  a_j(T) = sum_{i=1..4} c_{i,j}(MR) * (T/T_ref)^{i-1}
  κ_j(MR) = poly_4(MR; K_{j,1..5})
  MR = Y_H2O / Y_CO2 (molar ratio by mole fraction)
  T_ref = 1200 K

Validity (metadata below): T=500–5000 K; mixtures MR=0.4–4.0; pressure-pathlength p·L = 0.01–60 atm·m.

Outputs (in the chosen output directory):
  - ehlme2025_mixture.json
  - ehlme2025_pure_H2O.json
  - ehlme2025_pure_CO2.json
"""


import json
from typing import Dict, List, Tuple

# -------------------------------
# Canonical data structure
# -------------------------------
EHLME2025_WSGGM: Dict = {
    "metadata": {
        "model": "Ehlmé-2025 WSGGM (4+1 gray gases)",
        "T_ref_K": 1200.0,
        "weight_poly_degree": 4,  # in temperature (powers 0..3)
        "ratio_poly_degree": 4,    # in MR (powers 0..4)
        "MR_definition": "MR = Y_H2O / Y_CO2 (molar ratio)",
        "validity": {
            "temperature_K": [500.0, 5000.0],
            "mixture_MR": [0.4, 4.0],
            "pressure_pathlength_atm_m": [0.01, 60.0],
            "pressure_note": "Coefficients were fit at ~1 atm; validity stated as p·L in atm·m.",
        },
        "notes": [
            "Weights include only the 4 absorbing gray gases. Clear-gas weight = 1 - sum(a_j).",
            "For single-species (H2O in clear gas or CO2 in clear gas), MR polynomial is not used.",
        ],
    },
    "mixture": {
        # κ_j(MR) = K1 + K2*MR + K3*MR^2 + K4*MR^3 + K5*MR^4
        "K": {
            "1": [0.0304, 0.0015, -0.0014, 4.00e-04, 0.0],
            "2": [12.590, -1.5369, 1.2160, -0.3610, 0.0365],
            "3": [0.1673, 0.0057, 0.0016, -6.00e-04, 1.00e-04],
            "4": [0.9753, 0.1190, -0.0171, -5.00e-04, 2.00e-04],
        },
        # c_{i,j}(MR) = C1 + C2*MR + C3*MR^2 + C4*MR^3 + C5*MR^4
        # Stored as C[str(i)][str(j)] -> [C1..C5]
        "C": {
            "1": {  # i = 1
                "1": [-0.0445, -0.0060, 0.0199, -0.0068, 7.00e-04],
                "2": [0.3040, 0.3128, -0.1528, 0.0353, -0.0031],
                "3": [0.2204, -0.1699, 0.0850, -0.0203, 0.0018],
                "4": [0.3645, -0.1980, 0.0799, -0.0162, 0.0013],
            },
            "2": {  # i = 2
                "1": [0.4789, -0.1145, 0.0376, -0.0066, 5.00e-04],
                "2": [-0.2090, -0.2939, 0.1428, -0.0329, 0.0029],
                "3": [0.0341, 0.2025, -0.0970, 0.0228, -0.0020],
                "4": [-0.1936, 0.4754, -0.2242, 0.0513, -0.0045],
            },
            "3": {  # i = 3
                "1": [-0.1672, 0.0296, -0.0069, 8.00e-04, 0.0],
                "2": [0.0563, 0.0855, -0.0418, 0.0096, -8.00e-04],
                "3": [-0.0393, -0.0256, 0.0100, -0.0022, 2.00e-04],
                "4": [0.0177, -0.1941, 0.0936, -0.0216, 0.0019],
            },
            "4": {  # i = 4
                "1": [0.0150, 7.00e-04, -0.0016, 5.00e-04, -1.00e-04],
                "2": [-0.0056, -0.0079, 0.0039, -9.00e-04, 1.00e-04],
                "3": [0.0045, -0.0017, 0.0014, -3.00e-04, 0.0],
                "4": [0.0022, 0.0222, -0.0108, 0.0025, -2.00e-04],
            },
        },
    },
    "pure": {
        # For single-species in a clear gas: a_j(T) = sum_i c_{i,j} (T/T_ref)^{i-1}; κ_j = constant (no MR dependence)
        "H2O": {
            "kappa": [0.0468, 12.065, 0.2612, 1.3838],
            "c": {
                "1": [0.0958, 0.4709, 0.0098, 0.1402],
                "2": [0.2559, -0.4102, 0.3383, 0.2131],
                "3": [-0.1009, 0.1193, -0.1267, -0.1310],
                "4": [0.0105, -0.0115, 0.0123, 0.0176],
            },
        },
        "CO2": {
            "kappa": [0.0309, 68.594, 0.2221, 2.2780],
            "c": {
                "1": [0.0645, 0.1220, 0.1421, 0.0596],
                "2": [0.1870, -0.0476, -0.0656, 0.0646],
                "3": [-0.0880, 0.0014, 0.0176, -0.0403],
                "4": [0.0102, 8.00e-04, -0.0021, 0.0054],
            },
        },
    },
}

# -------------------------------
# Helper evaluation utilities
# -------------------------------

def _poly_eval(x: float, coeffs: List[float]) -> float:
    """Evaluate sum_{n=0..N-1} coeffs[n] * x^n (coeffs[n] = constant for n=0)."""
    p = 0.0
    xn = 1.0
    for cn in coeffs:
        p += cn * xn
        xn *= x
    return p


def kappa_mixture(MR: float) -> List[float]:
    """Return [κ1, κ2, κ3, κ4] for a mixture at MR = Y_H2O/Y_CO2.
       Clamps MR to metadata validity range for evaluation robustness.
    """
    lo, hi = EHLME2025_WSGGM["metadata"]["validity"]["mixture_MR"]
    MRc = max(lo, min(hi, MR))
    K = EHLME2025_WSGGM["mixture"]["K"]
    return [
        _poly_eval(MRc, K["1"]),
        _poly_eval(MRc, K["2"]),
        _poly_eval(MRc, K["3"]),
        _poly_eval(MRc, K["4"]),
    ]


def c_ij_mixture(MR: float, i: int, j: int) -> float:
    """Return c_{i,j}(MR) for mixture.
    i in {1..4}, j in {1..4}.
    """
    lo, hi = EHLME2025_WSGGM["metadata"]["validity"]["mixture_MR"]
    MRc = max(lo, min(hi, MR))
    C = EHLME2025_WSGGM["mixture"]["C"][str(i)][str(j)]
    return _poly_eval(MRc, C)


def a_j_mixture(T: float, MR: float, j: int) -> float:
    """Weight a_j(T, MR) for mixture.
    j in {1..4}.
    """
    Tref = EHLME2025_WSGGM["metadata"]["T_ref_K"]
    theta = T / Tref
    # sum_{i=1..4} c_{i,j}(MR) * theta^{i-1}
    return (
        c_ij_mixture(MR, 1, j) * (theta ** 0)
        + c_ij_mixture(MR, 2, j) * (theta ** 1)
        + c_ij_mixture(MR, 3, j) * (theta ** 2)
        + c_ij_mixture(MR, 4, j) * (theta ** 3)
    )


def a_and_kappa_mixture(T: float, MR: float) -> Tuple[List[float], List[float], float]:
    """Return (a_j list, kappa_j list, a_clear) for a mixture at (T, MR).
    a_clear = 1 - sum(a_j). No further clipping is applied.
    """
    a = [a_j_mixture(T, MR, j) for j in range(1, 5)]
    k = kappa_mixture(MR)
    a_clear = 1.0 - sum(a)
    return a, k, a_clear


def a_and_kappa_pure(T: float, species: str) -> Tuple[List[float], List[float], float]:
    """Return (a_j list, kappa_j list, a_clear) for a single species in clear gas.
    species = 'H2O' or 'CO2'.
    """
    Tref = EHLME2025_WSGGM["metadata"]["T_ref_K"]
    theta = T / Tref
    data = EHLME2025_WSGGM["pure"][species]
    k = list(data["kappa"])  # copy
    c = data["c"]
    a = [
        c["1"][0] * (theta ** 0) + c["2"][0] * (theta ** 1) + c["3"][0] * (theta ** 2) + c["4"][0] * (theta ** 3),
        c["1"][1] * (theta ** 0) + c["2"][1] * (theta ** 1) + c["3"][1] * (theta ** 2) + c["4"][1] * (theta ** 3),
        c["1"][2] * (theta ** 0) + c["2"][2] * (theta ** 1) + c["3"][2] * (theta ** 2) + c["4"][2] * (theta ** 3),
        c["1"][3] * (theta ** 0) + c["2"][3] * (theta ** 1) + c["3"][3] * (theta ** 2) + c["4"][3] * (theta ** 3),
    ]
    a_clear = 1.0 - sum(a)
    return a, k, a_clear


# -------------------------------
# JSON writers
# -------------------------------

def write_json(out_dir: str = ".") -> Tuple[str, str, str]:
    """Write three JSON files; return their paths."""
    # mixture JSON
    mix = {
        "metadata": EHLME2025_WSGGM["metadata"],
        "mixture": EHLME2025_WSGGM["mixture"],
    }
    pure_h2o = {
        "metadata": EHLME2025_WSGGM["metadata"],
        "pure": {"H2O": EHLME2025_WSGGM["pure"]["H2O"]},
    }
    pure_co2 = {
        "metadata": EHLME2025_WSGGM["metadata"],
        "pure": {"CO2": EHLME2025_WSGGM["pure"]["CO2"]},
    }

    path_mix = f"{out_dir.rstrip('/')}/ehlme2025_mixture.json"
    path_h2o = f"{out_dir.rstrip('/')}/ehlme2025_pure_H2O.json"
    path_co2 = f"{out_dir.rstrip('/')}/ehlme2025_pure_CO2.json"

    with open(path_mix, "w", encoding="utf-8") as f:
        json.dump(mix, f, indent=2)
    with open(path_h2o, "w", encoding="utf-8") as f:
        json.dump(pure_h2o, f, indent=2)
    with open(path_co2, "w", encoding="utf-8") as f:
        json.dump(pure_co2, f, indent=2)

    return path_mix, path_h2o, path_co2


# -------------------------------
# Small self-test / example usage
# -------------------------------
if __name__ == "__main__":
    paths = write_json(".")
    print("Wrote:")
    for p in paths:
        print("  ", p)

    # Example evaluations
    T_example = 2500.0  # K
    MR_example = 1.0    # Y_H2O / Y_CO2

    a_mix, k_mix, a0_mix = a_and_kappa_mixture(T_example, MR_example)
    print(f"\nMixture example (T={T_example} K, MR={MR_example}):")
    print("a_j   =", [round(x, 6) for x in a_mix], "  sum(a_j)=", round(sum(a_mix), 6), "  a_clear=", round(a0_mix, 6))
    print("kappa =", [round(x, 6) for x in k_mix])

    a_h2o, k_h2o, a0_h2o = a_and_kappa_pure(T_example, "H2O")
    print(f"\nPure H2O in clear gas (T={T_example} K):")
    print("a_j   =", [round(x, 6) for x in a_h2o], "  sum(a_j)=", round(sum(a_h2o), 6), "  a_clear=", round(a0_h2o, 6))
    print("kappa =", [round(x, 6) for x in k_h2o])

    a_co2, k_co2, a0_co2 = a_and_kappa_pure(T_example, "CO2")
    print(f"\nPure CO2 in clear gas (T={T_example} K):")
    print("a_j   =", [round(x, 6) for x in a_co2], "  sum(a_j)=", round(sum(a_co2), 6), "  a_clear=", round(a0_co2, 6))
    print("kappa =", [round(x, 6) for x in k_co2])



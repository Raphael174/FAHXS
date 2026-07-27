""" 
@ author : Raphaël Aubry
"""

# -----------------------------------------------------------------------------
# Ehlmé-2025 WSGGM support (polynomial-in-MR coefficients, up to 5000 K)
# -----------------------------------------------------------------------------

import json
from dataclasses import dataclass
import numpy as np 

@dataclass
class Ehlme2025Coeffs:
    """Container for Ehlmé et al. (2025) coefficient polynomials.

    JSON schema expected (from the companion coeff builder):
      {
        "metadata": {"T_ref_K": 1200.0, ...},
        "mixture": {
           "K": {"1": [K11..K15], ..., "4": [...]},
           "C": {"1": {"1": [C111..C115], ..., "4": [...]},
                   "2": {...}, "3": {...}, "4": {...}}
        },
        "pure": {"H2O": {"kappa": [..4..], "c": {"1":[..4..],..,"4":[..4..]}},
                  "CO2": {"kappa": [..4..], "c": {...}} }
      }
    All K/C coefficients are dimensionless polynomials in MR and θ=T/Tref.
    """
    Tref: float
    mix_K: dict
    mix_C: dict
    pure: dict

    @classmethod
    def from_json(cls, path: str) -> "Ehlme2025Coeffs":
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        md = data.get("metadata", {})
        mixture = data.get("mixture", {})
        pure = data.get("pure", {})
        return cls(
            Tref=float(md.get("T_ref_K", 1200.0)),
            mix_K=mixture.get("K", {}),
            mix_C=mixture.get("C", {}),
            pure=pure,
        )

    # ---- polynomial evaluators ----
    @staticmethod
    def _poly(x: float, coeffs: list[float]) -> float:
        p = 0.0; xn = 1.0
        for c in coeffs:
            p += c * xn
            xn *= x
        return p

    def kappa_mix(self, MR: float) -> np.ndarray:
        # returns (4,) array for gray gases j=1..4
        K = self.mix_K
        return np.array([
            self._poly(MR, K["1"]),
            self._poly(MR, K["2"]),
            self._poly(MR, K["3"]),
            self._poly(MR, K["4"]),
        ], float)

    def a_mix(self, T: float, MR: float) -> np.ndarray:
        # a_j(T,MR) = sum_{i=1..4} C_{i,j}(MR) * (T/Tref)^{i-1}
        C = self.mix_C
        th = float(T)/float(self.Tref)
        thp = np.array([1.0, th, th*th, th*th*th])
        a = np.zeros(4, float)
        for j in range(1,5):
            # gather C_i_j(MR) for i=1..4
            cvec = np.array([
                self._poly(MR, C["1"][str(j)]),
                self._poly(MR, C["2"][str(j)]),
                self._poly(MR, C["3"][str(j)]),
                self._poly(MR, C["4"][str(j)]),
            ])
            a[j-1] = float(np.dot(cvec, thp))
        return a

    def kappa_pure(self, species: str) -> np.ndarray:
        return np.array(self.pure[species]["kappa"], float)

    def a_pure(self, T: float, species: str) -> np.ndarray:
        # a_j(T) = sum_{i=1..4} c_i_j * (T/Tref)^{i-1}
        th = float(T)/float(self.Tref)
        thp = np.array([1.0, th, th*th, th*th*th])
        c = np.array(self.pure[species]["c"])  # shape (4,)
        # stored as dict-of-lists; rebuild 4x4
        cmat = np.vstack([
            np.array(self.pure[species]["c"]["1"], float),
            np.array(self.pure[species]["c"]["2"], float),
            np.array(self.pure[species]["c"]["3"], float),
            np.array(self.pure[species]["c"]["4"], float),
        ])  # (4,4) with columns j=1..4
        return (thp @ cmat).astype(float)


class WSGGM_Ehlme:
    """Evaluator implementing Ehlmé-2025 equations.

    Mixture case uses MR = yH2O / yCO2 (by mole fraction). If either species
    is ~0, falls back to the corresponding single-species set.
    κ_j multiplies p_a = (p_H2O + p_CO2); a_j depends on (T, MR) or (T) for pure.
    """
    def __init__(self, coeffs: Ehlme2025Coeffs, *, eps_clip=(1e-9, 0.999999), MR_bounds=(0.4, 4.0)):
        self.c = coeffs
        self.eps_min, self.eps_max = eps_clip
        self.MR_lo, self.MR_hi = MR_bounds

    def _mr(self, pw: float, pc: float) -> float:
        pw = max(0.0, float(pw)); pc = max(0.0, float(pc))
        if pc <= 0.0:
            return float(np.inf) if pw > 0.0 else 1.0
        return pw/pc

    # def emissivity(self, T: float, pw_Pa: float, pc_Pa: float, Le_m: float) -> float:
    #     pw = float(pw_Pa); pc = float(pc_Pa)
    #     pa = pw + pc
    #     if pa <= 0.0 or Le_m <= 0.0:
    #         return 0.0
    #     if pw <= 0.0 and pc <= 0.0:
    #         return 0.0
    #     if pw <= 0.0 and pc > 0.0:
    #         a = self.c.a_pure(T, "CO2"); k = self.c.kappa_pure("CO2")
    #     elif pc <= 0.0 and pw > 0.0:
    #         a = self.c.a_pure(T, "H2O"); k = self.c.kappa_pure("H2O")
    #     else:
    #         MR = self._mr(pw, pc)
    #         MRc = max(self.MR_lo, min(self.MR_hi, MR))
    #         a = self.c.a_mix(T, MRc)
    #         k = self.c.kappa_mix(MRc)
    #     tau = np.clip(k * pa * float(Le_m), 0.0, 700.0)
    #     eps = float(np.sum(a * (1.0 - np.exp(-tau))))
        
    def emissivity(self, T: float, pw_Pa: float, pc_Pa: float, Le_m: float) -> float:
        # partial pressures in Pa -> convert to atm for use with kappa in 1/(atm·m)
        pw = float(pw_Pa); pc = float(pc_Pa)
        pa_atm = (pw + pc) / 101325.0
        if pa_atm <= 0.0 or Le_m <= 0.0:
            return 0.0

        # choose coefficient set
        if pw <= 0.0 and pc > 0.0:
            a = self.c.a_pure(T, "CO2"); k = self.c.kappa_pure("CO2")
        elif pc <= 0.0 and pw > 0.0:
            a = self.c.a_pure(T, "H2O"); k = self.c.kappa_pure("H2O")
        else:
            MR  = self._mr(pw, pc)           # pw/pc is fine -> equals yH2O/yCO2
            MRc = max(self.MR_lo, min(self.MR_hi, MR))
            a   = self.c.a_mix(T, MRc)
            k   = self.c.kappa_mix(MRc)

        # optical thickness; clip to avoid exp underflow (exp(-50) ~ 2e-22)
        tau = np.clip(k * pa_atm * float(Le_m), 0.0, 50.0)
        eps = float(np.sum(a * (1.0 - np.exp(-tau))))
        return float(np.clip(eps, self.eps_min, self.eps_max))
        

    def absorptivity(self, Ts: float, pw_Pa: float, pc_Pa: float, Le_m: float) -> float:
        return self.emissivity(Ts, pw_Pa, pc_Pa, Le_m)


class RadiativeBackendEhlme:
    """Adapter matching your solver API for the Ehlmé model."""
    def __init__(self, model: WSGGM_Ehlme):
        self.model = model

    def __call__(self, *, T_eval: float, p: float, yH2O: float, yCO2: float, Le: float, **_):
        pw = float(p) * float(yH2O)
        pc = float(p) * float(yCO2)
        return self.model.emissivity(T_eval, pw, pc, Le)


def make_ehlme_backend(json_path_or_dir: str) -> RadiativeBackendEhlme:
    """Load Ehlmé-2025 coefficients and return a backend.

    Accepts either:
      - path to a *combined* JSON (mixture+pure), or
      - a *directory* containing the three files emitted by the coeff builder:
        ehlme2025_mixture.json, ehlme2025_pure_H2O.json, ehlme2025_pure_CO2.json
    """
    import os
    path = json_path_or_dir
    if os.path.isdir(path):
        # Merge the three files into one dict shape expected by Ehlme2025Coeffs
        with open(os.path.join(path, "ehlme2025_mixture.json"), "r", encoding="utf-8") as f:
            mix = json.load(f)
        with open(os.path.join(path, "ehlme2025_pure_H2O.json"), "r", encoding="utf-8") as f:
            h2o = json.load(f)
        with open(os.path.join(path, "ehlme2025_pure_CO2.json"), "r", encoding="utf-8") as f:
            co2 = json.load(f)
        merged = {
            "metadata": mix.get("metadata", {}),
            "mixture": mix.get("mixture", {}),
            "pure": {"H2O": h2o.get("pure", {}).get("H2O", {}),
                     "CO2": co2.get("pure", {}).get("CO2", {})},
        }
        tmp_json = os.path.join(path, "_ehlme2025_merged.json")
        with open(tmp_json, "w", encoding="utf-8") as f:
            json.dump(merged, f, indent=2)
        coeffs = Ehlme2025Coeffs.from_json(tmp_json)
    else:
        coeffs = Ehlme2025Coeffs.from_json(path)
    model = WSGGM_Ehlme(coeffs)
    return RadiativeBackendEhlme(model)


# ---- Example wiring (disabled by default) ----
if __name__ == "__main__":
    from pathlib import Path
    from radiation_model.radiation_equations import qrad_net_mbl, hrad_from_q


    # try:
    # point this to your emitted JSON (mixture + pure combined file)
    # For convenience, you can also embed pure into the same file as in the builder.
    backend = make_ehlme_backend(Path(__file__).parent/"ehlme2025_mixture.json")  # if only mixture is present
    # If you built a single JSON that also includes 'pure', the loader will use it automatically.
    p = 5e5; yH2O=0.12; yCO2=0.09; Le=0.08; Tg=1900.0; Ts=950.0; eps_s=0.8
    eps_emit = backend(T_eval=Tg, p=p, yH2O=yH2O, yCO2=yCO2, Le=Le)
    eps_abs  = backend(T_eval=Ts, p=p, yH2O=yH2O, yCO2=yCO2, Le=Le)
    qpp      = qrad_net_mbl(Tg, Ts, eps_emit, eps_abs, eps_s)
    hrad     = hrad_from_q(Tg, Ts, qpp)
    print(f"Ehlmé backend: eps_emit={eps_emit:.4f}, eps_abs={eps_abs:.4f}, q''={qpp:.1f} W/m^2, h_rad={hrad:.1f} W/m^2-K")
    # except FileNotFoundError:
    #     print("?")
    #     pass

# Latent heats of vaporization (approximate), molar basis at ~298 K unless noted
# Sources are NIST Chemistry WebBook pages (see citations below).

#%%
import numpy as np 

# Units: kJ/mol
DHVAP_KJMOL = {
    # Alcohol
    "ETHANOL": 42.3,   # Ethanol ΔHvap°(298 K) = 42.34 ± 0.08 kJ/mol (Green 1960, compiled by NIST)

    # Paraffins
    "IC8H18": 35.0,    # iso-octane (2,2,4-trimethylpentane). NIST gives Majer & Svoboda (1985) correlation:
                       # ΔHvap = A*exp(-β*Tr)*(1-Tr)^β with A=50.28 kJ/mol, β=0.2668, Tc=543.9 K.
                       # Evaluated near 298 K → ~35 kJ/mol.

    "NC7H16": 36.0,    # n-heptane. NIST Majer & Svoboda (1985) correlation (Tc=540.2 K) → ~36–37 kJ/mol at 298 K.

    # Olefins
    "C6H12-1": 31.5,   # 1-hexene. NIST page lists ΔHvap° and correlation; value near 298 K ≈ 31–32 kJ/mol.
    "C5H10-1": 29.7,   # 1-pentene. NIST correlation (Majer & Svoboda, Tc=464.7 K) gives ≈ 29–30 kJ/mol at ~298 K.

    # Light olefin (gas at 298 K; value is at saturation, close to normal bp)
    "IC4H8": 21.3,     # isobutene (2-methylpropene). NIST lists ΔHvap ≈ 21.3 kJ/mol near 261 K (at 1 atm).
}

# Molecular weights for conversion to kJ/kg (Cantera has these, but include here for clarity), g/mol
MW = {
    "ETHANOL": 46.06844,
    "IC8H18": 114.231,
    "NC7H16": 100.205,
    "C6H12-1": 84.1595,   # 1-hexene
    "C5H10-1": 70.133,    # 1-pentene
    "IC4H8": 56.1063,     # isobutene
}

def dHvap_kJ_per_kg(species: str) -> float:
    """Convert ΔHvap from kJ/mol to kJ/kg using molecular weight."""
    return DHVAP_KJMOL[species] / (MW[species] / 1000.0)

# Example: mixture latent heat by mass-fraction weighting
def mixture_dHvap_kJ_per_kg(mass_fracs: dict) -> float:
    """
    mass_fracs: dict {species: Y_i} that sums to 1 over the species present here.
    Returns Σ Y_i * (ΔHvap_i in kJ/kg).
    """
    return sum(mass_fracs[s] * dHvap_kJ_per_kg(s) for s in mass_fracs if s in DHVAP_KJMOL)

#%%

Y_fuel_E5 = {'C2H5OH': np.float64(0.022521123344674744), 'C6H12-1': np.float64(0.05760037978051016),
            'C5H10-1': np.float64(0.020571564207325047),'IC4H8': np.float64(0.010971500910573358),
            'IC8H18': np.float64(0.6121661440853741), 'NC7H16': np.float64(0.2761692876715427)}

Y_fuel_E10 = {'C2H5OH': np.float64(0.04648289377205732), 'C6H12-1': np.float64(0.059442690615179886),
            'C5H10-1': np.float64(0.021229532362564234), 'IC4H8': np.float64(0.01132241726003426),
            'IC8H18': np.float64(0.5936889008300161), 'NC7H16': np.float64(0.2678335651601482)}

#%%

print(f"dH SP95 E5 {mixture_dHvap_kJ_per_kg(Y_fuel_E5)}")
print(f"dH SP95 E10 {mixture_dHvap_kJ_per_kg(Y_fuel_E10)}")

#%%
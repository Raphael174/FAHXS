""" 
Use this script just to build the initial fuel mass composition to then 
save as pre-computed for code runs. 

Here producing an approimate Gasoline (SP95) E5/E10 fuel 
"""

import cantera as ct

def build_e5_e10_like_fuel_X(
    gas: ct.Solution,
    X_ethanol: float = 0.10,  # overall mole fraction of ethanol (e.g., 0.05 for E5, 0.10 for E10)
    olefin_split: dict[str, float] | None = None,  # relative split within the olefin bucket (any positive weights)
    paraffin_split: dict[str, float] | None = None # relative split within the paraffin bucket (any positive weights)
):
    """
    Build a gasoline-like surrogate composition (mole fractions) using three buckets:
      Ethanol, Olefins, Paraffins.

    - X_ethanol: total mole fraction of ethanol (C2H5OH).
    - olefin_split: dict of {species: weight} for the olefin bucket (weights need not sum to 1).
        Defaults: {"C6H12-1": 7, "C5H10-1": 3, "IC4H8": 2}  # totals 12 parts -> 12% if using defaults below
    - paraffin_split: dict of {species: weight} for the paraffin bucket.
        Defaults (scaled to fill remainder): {"IC8H18": 35, "NC7H16": 18}

    Returns:
        X_dict: normalized mole-fraction dict
        X_string: Cantera composition string (mole basis)
    """
    # ---- defaults (chosen to mirror the earlier example buckets) ----
    if olefin_split is None:
        olefin_split = {"C6H12-1": 7.0, "C5H10-1": 3.0, "IC4H8": 2.0}  # relative parts
        X_olefins_total_default = 0.12  # overall mole fraction for the olefin bucket
    else:
        X_olefins_total_default = 0.12  # you can tweak this line or expose as an argument if you like

    if paraffin_split is None:
        paraffin_split = {"IC8H18": 35.0, "NC7H16": 18.0}  # relative parts

    # ---- sanity checks ----
    if not (0.0 <= X_ethanol < 1.0):
        raise ValueError("X_ethanol must be in [0,1).")

    # Compute olefin bucket normalized split
    ole_sum = sum(max(0.0, w) for w in olefin_split.values())
    if ole_sum <= 0.0:
        raise ValueError("olefin_split must contain positive weights.")
    olefin_norm = {s: max(0.0, w) / ole_sum for s, w in olefin_split.items()}

    # Total mole fraction assigned to the olefin bucket
    X_olefins = X_olefins_total_default

    # Remainder goes to paraffins
    X_paraffins = 1.0 - X_ethanol - X_olefins
    if X_paraffins <= 1e-12:
        raise ValueError("Paraffin bucket got non-positive remainder. Reduce X_ethanol and/or olefin total.")

    # Normalize paraffin split
    para_sum = sum(max(0.0, w) for w in paraffin_split.values())
    if para_sum <= 0.0:
        raise ValueError("paraffin_split must contain positive weights.")
    para_norm = {s: max(0.0, w) / para_sum for s, w in paraffin_split.items()}

    # ---- build raw mole fractions ----
    X = {}

    # Ethanol
    ethanol_sp = "C2H5OH"
    X[ethanol_sp] = X_ethanol

    # Olefins
    for sp, f in olefin_norm.items():
        X[sp] = X.get(sp, 0.0) + X_olefins * f

    # Paraffins
    for sp, f in para_norm.items():
        X[sp] = X.get(sp, 0.0) + X_paraffins * f

    # ---- validate species exist in mechanism & normalize ----
    missing = [s for s in X if s not in gas.species_names]
    if missing:
        raise KeyError(f"Species not found in mechanism: {missing}")

    # Clean small negatives, normalize
    for s,v in list(X.items()):
        if v < 0.0:
            X[s] = 0.0
    S = sum(X.values())
    if S <= 0:
        raise ValueError("Resulting composition is empty.")
    X = {s: v/S for s, v in X.items()}

    # Make a Cantera composition string (mole basis)
    X_string = ", ".join(f"{s}:{X[s]:.8g}" for s in X)

    return X, X_string

#%%

def mole_to_mass_fractions(gas: ct.Solution, X_dict: dict[str, float]):
    # filter out zeros and normalize X
    X = {s: x for s, x in X_dict.items() if x > 0.0}
    ssum = sum(X.values())
    if ssum <= 0:
        raise ValueError("Empty or zero mole-fraction set.")
    X = {s: x/ssum for s, x in X.items()}

    # get MWs from the mechanism
    MW = {}
    missing = []
    for s in X:
        if s in gas.species_names:
            MW[s] = gas.molecular_weights[gas.species_index(s)]  # g/mol
        else:
            missing.append(s)
    if missing:
        raise KeyError(f"Species not in mechanism: {missing}")

    # denominator: sum_j X_j * W_j
    denom = sum(X[s] * MW[s] for s in X)

    # mass fractions
    Y = {s: (X[s] * MW[s]) / denom for s in X}

    # tiny numerical cleanup & renormalize
    Ysum = sum(Y.values())
    Y = {s: max(0.0, y) for s, y in Y.items()}
    Y = {s: y / Ysum for s, y in Y.items()}  # ensure sum to 1.0

    # Make a Cantera composition string (mass-fraction form)
    Y_string = ", ".join(f"{s}:{Y[s]:.8g}" for s in Y)
    return Y, Y_string

# --- example usage ---
# gas = ct.Solution("gasoline_surrogate.yaml")
# X_in = {"C2H5OH":0.10, "C6H5CH3":0.18, "XYLENE":0.065, "C6H6":0.005,
#         "C6H12-1":0.07, "C5H10-1":0.03, "IC4H8":0.02, "IC8H18":0.35, "NC7H16":0.18}
# Y_dict, Y_str = mole_to_mass_fractions(gas, X_in)
# print(Y_str)
# gas.TPY = 300, ct.one_atm, Y_str

#%%
def set_mix_at_OF_mass(gas, Y_fuel, OF, T_ign=800, P=ct.one_atm, o2_name="O2"):
    """Mix fuel (mass fractions) with pure O2 at target O/F (mass)."""
    s = sum(Y_fuel.values())
    if s <= 0: raise ValueError("Y_fuel is empty.")
    Yf = {k:v/s for k,v in Y_fuel.items()}
    Y_mix = {k: (1.0/(1.0+OF))*Yf[k] for k in Yf}
    Y_mix[o2_name] = OF/(1.0+OF)
    gas.TPY = T_ign, P, Y_mix
    gas.equilibrate("HP")
    return gas

#%%
# ----------------
# Example usage:
gas = ct.Solution("C:/Users/raubry/Desktop/Combustor-HX/combustion_chemistry/gasoline_surrogate/llnl_gasoline_Detailed.yaml")
X_dict, X_str = build_e5_e10_like_fuel_X(
    gas,
    X_ethanol=0.1,  # E10-like
    olefin_split={"C6H12-1":7, "C5H10-1":3, "IC4H8":2},  # defaults
    paraffin_split={"IC8H18":35, "NC7H16":18}            # defaults
)

Y_dict, Y_string = mole_to_mass_fractions(gas, X_dict)

#%%
gas = set_mix_at_OF_mass(gas, Y_dict, OF=2.8)

#%%
"""
SP95 E10
C2H5OH:0.1, C6H12-1:0.07, C5H10-1:0.03, IC4H8:0.02, IC8H18:0.51509434, NC7H16:0.26490566

SP95 E5 
C2H5OH:0.05, C6H12-1:0.07, C5H10-1:0.03, IC4H8:0.02, IC8H18:0.54811321, NC7H16:0.28188679
"""



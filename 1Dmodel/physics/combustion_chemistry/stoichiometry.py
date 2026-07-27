""" 
@ author : Raphaël Aubry
"""

""" 
Script to estimate Kerosene combustion chemistry using cantera
"""
#%%
import cantera as ctr 
from pathlib import Path
import numpy as np 
import matplotlib.pyplot as plt

ctr.suppress_thermo_warnings()

#%%

base_path = Path(__file__).parent

JETA_chem_path = base_path / "A2highT.yaml"
phase = ctr.Solution(JETA_chem_path)


def compute_OF_st(fuel_species):
    """Mass-based stoich O/F for pure O2 from species elemental makeup."""
    iF = phase.species_index(fuel_species)  # will error if not in mech
    a = phase.n_atoms(fuel_species, "C")
    b = phase.n_atoms(fuel_species, "H")

    nu_O2 = (a + b/4)                        # mol O2 per mol fuel
    MW_O2 = phase.molecular_weights[phase.species_index("O2")]
    MW_F  = phase.molecular_weights[iF]
    return (nu_O2 * MW_O2) / MW_F

#%%

fuel = "POSF10325"

print(f"OF_st for {fuel}/O2 mixture = {compute_OF_st(fuel)}")


#%%
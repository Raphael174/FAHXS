""" 
Misc functions

"""
import numpy as np

def print_mass_of_N_most_abundant_species(GAS, N):
    """ 
    Upon user request, print the N most abundant species with name and mass fraction
    """

    unsorted_names = np.array(GAS.species_names)
    unsorted_mass_fractions = GAS.Y

    # python sorts array directly in ascending order
    mass_fractions_indices_ascending = np.argsort(unsorted_mass_fractions)
    # Take the first element as the most abundant one
    # This is the full list of indices for decreasing Y

    mass_fractions_indices_DESCENDING = mass_fractions_indices_ascending[::-1]
    # This is the list of indices for the top N most abundant
    mass_fractions_N_indices_DESCENDING = mass_fractions_indices_DESCENDING[:N]

    species_names_sorted = unsorted_names[mass_fractions_N_indices_DESCENDING]
    species_massFractions_sorted = unsorted_mass_fractions[mass_fractions_N_indices_DESCENDING]

    for i in range(N):
        print(
            f'    Y_{species_names_sorted[i]} : {np.round(species_massFractions_sorted[i], 3)}')

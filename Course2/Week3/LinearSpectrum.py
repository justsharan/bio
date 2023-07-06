from week3 import INTEGER_MASSES
from display import print_list

def linear_spectrum(peptide: str) -> list[int]:
    prefix_mass = []
    for i in range(0, len(peptide)+1):
        # Sum up the masses of each amino acid in peptide
        prefix_mass.append(sum(INTEGER_MASSES[aa] for aa in peptide[:i]))
    linear_spectrum = [0]
    for i in range(len(prefix_mass)):
        for j in range(i+1, len(prefix_mass)):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(linear_spectrum)

peptide_seq = 'MVDFELFRVMIYKRKGYLKFVYSILDHFIKLSTPNH'
print_list(linear_spectrum(peptide_seq))

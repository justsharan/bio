from week3 import INTEGER_MASSES
from display import print_list

def cyclic_spectrum(peptide: str) -> list[int]:
    prefix_mass = []
    for i in range(0, len(peptide)+1):
        prefix_mass.append(sum(INTEGER_MASSES[aa] for aa in peptide[:i]))
    peptide_mass = prefix_mass[-1]
    cyclic_spectrum = [0]
    for i in range(len(prefix_mass)):
        for j in range(i+1, len(prefix_mass)):
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if (i > 0) and (j < len(peptide)):
                cyclic_spectrum.append(peptide_mass - cyclic_spectrum[-1])
    return sorted(cyclic_spectrum)

peptide_seq = 'ELGPWPKRTMRCH'
print_list(cyclic_spectrum(peptide_seq))

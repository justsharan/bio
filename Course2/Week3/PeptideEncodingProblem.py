from helpers import window, reverse_comp, DNA_to_peptide
from display import print_list

def peptide_encoding_problem(dna: str, peptide: str) -> list[str]:
    fragments = []
    for frag in window(dna, len(peptide)*3):
        if (peptide == DNA_to_peptide(frag)) or (peptide == DNA_to_peptide(reverse_comp(frag))):
            fragments.append(frag)
    return fragments

with open('dataset_96_7.txt', 'r') as file:
    dna = file.readline().strip()
    peptide = file.readline().strip()
    print_list(peptide_encoding_problem(dna, peptide))

with open('Bacillus_brevis.txt', 'r') as file:
    genome = file.read().strip().replace('\n', '')
    tyrocidine_b1 = 'VKLFPWFNQY'
    for fragment in window(genome, len(tyrocidine_b1)*3):
        if DNA_to_peptide(fragment) == tyrocidine_b1:
            print(fragment)
        elif DNA_to_peptide(reverse_comp(fragment)) == tyrocidine_b1:
            print(fragment)

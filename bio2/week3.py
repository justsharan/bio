from typing import Generator

codon_table = {}
with open('bio2/RNA_codon_table.txt', 'r') as file:
    for line in file:
        items = line.strip().split(' ')
        codon_table[items[0]] = items[1] if len(items) == 2 else ''

integer_mass_table = {}
with open('bio2/integer_mass_table.txt', 'r') as file:
    for line in file:
        aa, mass = line.strip().split(' ')
        integer_mass_table[aa] = int(mass)

def chunks(seq: str, n: int = 3) -> Generator[str, None, None]:
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

def sliding_window(seq: str, n: int = 3) -> Generator[str, None, None]:
    for i in range(len(seq)):
        curr = seq[i:i+n]
        if len(curr) != n:
            break
        else:
            yield seq[i:i+n]

def complement(seq: str) -> str:
    pairs = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' }
    return ''.join(pairs[nuc] for nuc in seq)

def reverse_complement(seq: str) -> str:
    return complement(seq)[::-1]

def transcribe(seq: str) -> str:
    return seq.replace('T', 'U')

def translate(seq: str) -> str:
    return ''.join(codon_table[codon] for codon in chunks(seq, n=3))

def DNA_to_peptide(seq: str) -> str:
    return translate(transcribe(seq))

# with open('dataset_96_4.txt', 'r') as file:
#     seq = file.read().strip()
#     print(translate(seq))

# with open('dataset_96_7.txt', 'r') as file:
#     input1 = file.readline().strip()
#     input2 = file.readline().strip()
#     for fragment in sliding_window(input1, len(input2)*3):
#         if input2 in [DNA_to_peptide(fragment), DNA_to_peptide(reverse_complement(fragment))]:
#             print(fragment)

with open('Bacillus_brevis.txt', 'r') as file:
    genome = file.read().strip().replace('\n', '')
    tyrocidine_b1 = 'VKLFPWFNQY'
    for fragment in sliding_window(genome, len(tyrocidine_b1)*3):
        if DNA_to_peptide(fragment) == tyrocidine_b1:
            print(fragment)
        elif DNA_to_peptide(reverse_complement(fragment)) == tyrocidine_b1:
            print(fragment)

def linear_spectrum(peptide: str) -> list[int]:
    prefix_mass = []
    for i in range(0, len(peptide)+1):
        prefix_mass.append(sum(integer_mass_table[aa] for aa in peptide[:i]))
    linear_spectrum = [0]
    for i in range(len(prefix_mass)):
        for j in range(i+1, len(prefix_mass)):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(linear_spectrum)

# peptide_seq = 'MVDFELFRVMIYKRKGYLKFVYSILDHFIKLSTPNH'
# spectrum = map(str, linear_spectrum(peptide_seq))
# print(' '.join(spectrum))

def cyclic_spectrum(peptide: str) -> list[int]:
    prefix_mass = []
    for i in range(0, len(peptide)+1):
        prefix_mass.append(sum(integer_mass_table[aa] for aa in peptide[:i]))
    peptide_mass = prefix_mass[-1]
    cyclic_spectrum = [0]
    for i in range(len(prefix_mass)):
        for j in range(i+1, len(prefix_mass)):
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if (i > 0) and (j < len(peptide)):
                cyclic_spectrum.append(peptide_mass - cyclic_spectrum[-1])
    return sorted(cyclic_spectrum)

# peptide_seq = 'ELGPWPKRTMRCH'
# spectrum = map(str, cyclic_spectrum(peptide_seq))
# print(' '.join(spectrum))
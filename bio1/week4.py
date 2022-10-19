import itertools
import numpy as np
import random
from tools import window
from week3 import profile, profitable_kmer, score

def motifs(profile_matrix: np.ndarray, dna: list[str]) -> np.ndarray:
  res = []
  for strand in dna:
    res.append(profitable_kmer(strand, profile_matrix.shape[1], profile_matrix))
  return np.array([list(i) for i in res])

def randomized_motif_search(dna: list[str], k: int, t: int):
  # make a motif matrix by picking a random k-mer from each strand
  # rest of the function tries to improve Score(Motif) by altering it
  _motifs = np.array([list(random.choice(list(itertools.islice(window(strand, k), t)))) for strand in dna])
  best_motifs = np.copy(_motifs)
  while True:
    _motifs = motifs(profile(_motifs, pseudocount=1), dna)
    if score(_motifs) < score(best_motifs):
      best_motifs = _motifs
    else:
      return best_motifs

(k, t) = (15, 20)
Dna = [
  'CTTAATTACGATCATTGCATCCTTCACAAATGCCCACGTGGTGATTTACCTTTATAGGTCCCATCATGGCAGGGAATGGTCCGACTCCAGGGGTTATTATTTCCACCGTCGAGTAGTTAACTCGGCCCATCGCCGTAGAAGCTTTCATGTTGCAGGATTATTTGTGTAAAGGCTGACAGGGACAGCCTTAATTACGATCAT',
  'TGCATCCTTCACAAATGCCCACGTGGTGATTTACCTTTATAGGTCCCATCATGGCAGGGAATGGTCCGACTCCAGGGGTTATTATTTCCACCGTCGAGTAGTTAACTCGGCCCATCGCCGTAGAAGCTTTCATGTTGCAGGATTATTTGTGTAAAGGCTGACAGGGACAGACAAGAAGACCTCAACCTTAATTACGATCAT',
  'GAAACTGGGGTAGGGCCAGAGCTAGCCGAAAATCGTCTTCTTACCCACTGCCGGGGCTGGTGTGAAAATAGCGGATCGGACGCTGGGAGTAAACTGCGGATGGCACACTAACTAATTGATAATGCCAACGAACTCCGTCACCCGCTCCTAAGATGGCGACACAATTACGGATAGACCTCAAATATACTACGGCGGGTTATG',
  'AGAATCCGCGCACGGAATGCATTAATGCTCGTGGGGTTGGGCTAGGGTCAATAAGCGTGTTACACCGGTACCTCAAGCTCGCCGCCATCCAGCCATAGGTGTGTTGACCTGAGTCTGGTGGCGCGAAGCGCTTCGGCAGAGGTTACTCCCTGGCCCTACGTGAGAAGGCGTCATACAATCGGATGGTGTTCAGGACCACAA',
  'GAACCGGTACACCCTAGTCTTGAATTCGGCAGGGACCGCCGCAGTAAGGATTTCTAGACCTCAACCCGTTCGGCCTTCCCCCTATAACTAGTGACCGGGCTAAACATTTCGACTCTATCTACAGTCTAGTTCCGTTATGATGATTACTGGCCAGCACACAATTCGTCTCGTAGCTGGGATCAGGGTTTTGCAAAATACTGG',
  'CAGGGCCAATGGCACATATCCCGAGCATTGAATTTATATGTTTGAACGACTACTAAAATCAAGGACAGGGTTATCATTTCTTGAGGAGATCTTCTTCGCTCGCCGTATCTGAGGACCAGCAGCAGGGGTCGATTAGAGTAGGCTTAGATAGTGCAAAGATAGACTATGTTCACTGACGAACACCTCCCCCTCAATAAACGT',
  'CGCTCATCGTGAGGACTAATCGGCGCGAAGCGATAGACTCAAGACTGGCTTCACGACTGCATTAGTTAATCTCCTCGTTCTACTAAGTGCGGGAGCTCAAATTCGCGACTGCACCTTGAGTCTCACCTAGACCTCTTGTGTTGAGCGCGCTGGTTGAACATGCTTCGTTGTACCACACTGTAGGGCACATGAAATGCGAGA',
  'CATCGCCGCAAGCCGTCGTTCTATCCCACCTAGTTACTTGTCGCGAAGTGACTAGAATGTAGTAGTGGGTGCATCGACACCTAGGTATCAATATCGGGGTTGTTTTAGGAGCCTCGCGGTGATCAAGTATCGCTCCGCTAGCAATACGCCGCATTCGAACCACATACAGTGATCATCAGCCCGGCCGATATGTTCCAAGAT',
  'TCGAGGGACGCGCTAATTGTAACATCATAGGTACGGAGTGTATACACCTAGACGAGAAGTCTGTCACGAGATTTGAATACTAGAAGCCGTCAGTGGATCCGCGCCGAACAAAAGTAAGTCATAAGCGCATTTAATCCGTAACAGTGTCGAATAGAAATTTATACAGAAATTGTACTTAGGCCTCATTGAACTCTGACAATT',
  'TACAGAAAGACCTCAACGGCTATATTATCGGCCTCGATGGTGCTCACGTGACGCAATCATTAGTCTTCACCGATGTAATATTTTATCACCACGGCGATCCAAGAGGTCCTAAGGATCCGAGCCGGCTTCCGGTCGCTAGGAGGTGGGTAGAACCCCGGTCTGTATCAATGATCAAAGTAAAGACTTGGTATTAGTAATGTG',
  'ATGTATCATGAGATGCATGACTGCGGTGGGGTAGAAATCCGTATATGATCAACTACTTAACGGTACTGAGAAAATCTGGCTATCGCATTGTCTGTTAGGAGGGTAGCCCAACGAAGCCGTCACTCCCCCTATGGGGATCTACGAACGGAACATACTAAACACCTGAGCCTCAAGTTTCCTACCGGTGTGCCCCGATCTAGT',
  'TCTCGGAAACAGCTCTTACCGAAACGTCGCCAAGACTCTGCGGATAGAACAACTTTCAAAAGTAACATGGCGGTTTCGTAATTGGATAATGGAAAGACACGGGGACCTCAAGAAAAACAAGGGGCATTGCTGCGTGGCTTATCATGCCACTTGATAGGACTTATTTAGTTGATAAACACCAATTAATCTGTCTGCTAAGCT',
  'CTTGACTCATTGGGTGCCAGCTGAGGAATCACTTACGTACGTAACTATAGTGTTCTCCCACACCGACAGGCACGCACTCACTGCGTTAAAAAAGCATGCAGCTTAAGCGACCAGTCCTGGCGCTAACTCCTTGTATACTCGGGTGATTGGGGCACCCCGCTAAGTTGAATCTAGAAAAACACCTAGACCGTTACACTTATC',
  'GCATGAATTTATACTATATGAAATACAGGCGAACTACCTAGGACAATTAGAATTTCCCCGATACAAAGAGAGCTCTCACACCTATTGCTCAAGATCCCGCATGTATCTCCAGGGACAGTGCAAGGTTGCTGGAATTTAAGCTGTCTACGAGACCCGCAAGCGCTCTTACTATGAATAGCAGTTAGAACATCCGTATATATT',
  'GGTGCGATCGTGTCGATGGTTTAACTGGAGCTCTTCTTCGGAGATGATTTAAACACCAGAATTGCGAGAATCGCCAACTCATATACACTGTAGGAGGTAACACCTAGGGTCACAACCCCATATGTCACAGTCCGGGTCAAACCCCTAAATATTATGAGATGTGCACGTCCTAGACCTCAAGGTTGAGAACACCGCATGTTC',
  'ATTCCAAACACCCGGAAAGGACAACCGGCTAATGCTTAAAAGTTGAGCATCACCCTGTTTACAAACGCACCCACTGAGAGCAGGGCGGATCCACACCGTAACCTCAAATATAGAAGCTGCCTCGTACTTACGGTTTATTAGCAGGTTTGCCGTCTAGTCTTGAACCACCAGCGGAGCAAGGCAAGCTCGGGTCTGCCAGAG',
  'GAACAAAAACATCGGATATGTAGACACTAAAACTCTAAGTTGTATTACAACGTCGCCTCGTTCGGAGGTCCCAATAGAATCAAGCAACGAATATGGCTGACATTGTTTAGCTGGGCTGTCCCTCTATTTGAGGCAGAGTACCAGTATGTATTTGTTCACACAACGACCTCAAGGATTTTCGCGCTGGATAACTCTTATAGC',
  'ATTCGTCACGTACTAATTATGCCCCAGCTGATGTTACACCGGACACATGTCTAGTCTTTCAGTCGCCTTAACGGCCGATACTCTCTCTAACCCATTGACGACCGAGCGTGCCGTTTAAGACACCTAGACCTGCCTATTAACTATATGCTGAACGCGCGTACTATCGGAAGAAATTAGTAACGGGCTAAGGGGATGACTACT',
  'GATCCAAGAAAGTGGTGCATCCTCATTGGCGAGAGGGATCGCATTCTGGGTTTTGCTCGGGGGTTTGAACTCTAATGTGACCAAGGTGGCCGGGCTCACGGAAGGCCGGATGGCTTCAACACCTAGAAAGCAACATGTGCGACCTCAACGTGGTGCACCATAACGTGAAAAAATCCCGTACACGTCCAAGGGAAGCGACAA',
  'TTAAGCAACCTAGACCTCACCGTTTCCACCTGGTCTCCTATCGGCGCGACGACCCCTGCCTGGAACCTTCCGAGCATCCGGCTGCGCGACATTCATAGGGGTATATCGGGCGCGAGCGTTCAAGGCCCTGCTCTATATATTACACGCTCGTTGACATTTTGAGGGTCAGATAGTTAAGCCCTCAGCCCACAAACTGACGCT'
]

result = np.empty(shape=(t, k), dtype=str)
for i in range(1000):
  new_res = randomized_motif_search(Dna, k, t)
  if score(new_res) < score(result) or score(result) == 0:
    result = new_res
  print(i) # for keeping track of progress while the loop is running

# print result
print('\n'.join([''.join(i) for i in result]))
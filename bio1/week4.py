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
  # Generate a random k-mer for each strand in dna
  _motifs = np.array([list(random.choice(list(itertools.islice(window(strand, k), t)))) for strand in dna])
  best_motifs = np.copy(_motifs)
  while True:
    _motifs = motifs(profile(_motifs, pseudocount=1), dna)
    if score(_motifs) < score(best_motifs):
      best_motifs = np.copy(_motifs)
    else:
      return best_motifs

import itertools
import numpy as np
import random
from tools import window
from week3 import profile, profile_prob, profitable_kmer, score

def motifs(profile_matrix: np.ndarray, dna: list[str]) -> np.ndarray:
  res = []
  for strand in dna:
    res.append(profitable_kmer(strand, profile_matrix.shape[1], profile_matrix))
  return np.array([list(i) for i in res])

def randomized_motif_search(dna: list[str], k: int, t: int):
  # Generate a random k-mer for each strand in dna
  chosen_motifs = np.array([
    list(random.choice(list(itertools.islice(window(strand, k), t))))
    for strand in dna
  ])
  best_motifs = np.copy(chosen_motifs)
  while True:
    chosen_motifs = motifs(profile(chosen_motifs, pseudocount=1), dna)
    if score(chosen_motifs) < score(best_motifs):
      best_motifs = np.copy(chosen_motifs)
    else:
      return best_motifs

def profile_random_kmer(profile_matrix: np.ndarray, text: str, k: int) -> dict[str, float]:
  return { kmer: profile_prob(profile_matrix, kmer) for kmer in window(text, k) }

def gibbs_sampler(dna: list[str], k: int, t: int, n: int):
  # Generate a random k-mer for each strand in dna
  chosen_motifs = np.array([
    list(random.choice(list(itertools.islice(window(strand, k), t))))
    for strand in dna
  ])
  best_motifs = np.copy(chosen_motifs)
  for _ in range(n):
    i = random.randrange(t)
    pfile = profile(np.delete(chosen_motifs, i, axis=0), pseudocount=1)
    distribution = profile_random_kmer(pfile, dna[i], k)
    chosen_motifs[i,:] = list(random.choices(list(distribution.keys()), weights=list(distribution.values()), k=1)[0])
    if score(chosen_motifs) < score(best_motifs):
      best_motifs = np.copy(chosen_motifs)
  return best_motifs

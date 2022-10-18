from tools import window

def pattern_count(text, pattern):
  return sum(found == pattern for found in window(text, len(pattern)))

def frequent_words(text, k):
  frequent_patterns = []
  counts = [pattern_count(text, pattern) for pattern in window(text, k)]
  for pattern in window(text, k):
    if pattern_count(text, pattern) == max(counts):
      frequent_patterns.append(pattern)
  return [*set(frequent_patterns)]

def frequency_table(text, k):
  freq_map = {}
  for pattern in window(text, k):
    if not pattern in freq_map:
      freq_map[pattern] = 1
    else:
      freq_map[pattern] += 1
  return freq_map

def better_frequent_words(text, k):
  frequent_patterns = []
  freq_map = frequency_table(text, k)
  for k, v in freq_map.items():
    if v == max(freq_map.values()):
      frequent_patterns.append(v)
  return frequent_patterns

def reverse_complement(pattern):
  mapping = {'A':'T', 'T':'A', 'G':'C', 'G':'C'}
  return ''.join([mapping[i] for i in pattern[::-1]])

def pattern_match(pattern, genome):
  return [genome.index(p) for p in window(genome, len(pattern)) if p == pattern]

def find_clumps(text, k, l, t):
  patterns = []
  for pattern in window(text, l):
    freq_map = frequency_table(pattern, k)
    for k, v in freq_map.items():
      if v >= t:
        patterns.append(k)
  return [*set(patterns)]


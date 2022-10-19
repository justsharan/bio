from tools import window

def skew(pattern):
  count = 0
  out = []
  for nuc in pattern:
    if nuc == 'C':
      count -= 1
    elif nuc == 'G':
      count += 1
    out.append(count)
  return out

def min_skew(pattern):
  skew_list = skew(pattern)
  return [i for i in range(len(pattern)) if skew_list[i] == min(skew_list)]

def hamming_distance(pattern1: str, pattern2: str):
  count = 0
  for a, b in zip(pattern1, pattern2):
    count = count if a == b else count + 1
  return count

def approx_pattern_match(pattern, text, d):
  for i, curr_pattern in enumerate(window(text, len(pattern))):
    if hamming_distance(pattern, curr_pattern) <= d:
      yield i

def approx_pattern_count(text, pattern, d):
  count = 0
  for curr_pattern in window(text, len(pattern)):
    if hamming_distance(pattern, curr_pattern) <= d:
      count += 1
  return count

def immediate_neighbors(pattern, d):
  if d == 0:
    return [pattern]
  if len(pattern == 1):
    return list('ACGT')
  neighborhood = []
  suffix_neighbors = immediate_neighbors(pattern[1:], d)
  for text in suffix_neighbors:
    if hamming_distance(pattern[1:], text) < d:
      for nuc in 'ACGT':
        neighborhood.append(nuc+text)
    else:
      neighborhood.append(pattern[0]+text)
  return neighborhood

def freq_words_mismatches(text, k, d):
  freq_map = {}
  for i in range(len(text)-k):
    pattern = text[i:i+k]
    for neighbor in immediate_neighbors(pattern, d):
      freq_map[neighbor] = freq_map[neighbor]+1 if neighbor in freq_map else 1
  return [k for k, v in freq_map.items() if v == max(freq_map.values())]


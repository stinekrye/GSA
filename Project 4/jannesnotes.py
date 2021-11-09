import numpy as np

def c_table(x, sa):
    counts, buckets = {}, {}

    for i in range(0, len(sa)):
        char = x[sa[i]]
        if char not in counts:
            buckets[char] = i
            counts[char] = 0
        counts[char] += 1

    return counts, buckets

def o_table(x, sa):
    alphabet = sorted(set(x))
    alphabet_size = len(alphabet)
    length = len(sa) + 1

    O = np.zeros((length, alphabet_size))
    
    for i in range(1, length):
        for a in range(0, alphabet_size):
            if x[sa[i - 1] - 1] == alphabet[a]:
                O[i][a] = O[i - 1][a] + 1
            else: 
                O[i][a] = O[i - 1][a]

    return O

def fm_search(x, p, sa):
    
    return



sa = [11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]
x = 'mississippi$'
y = 'GATGCGAGAGATG$'
p = 'iss'

c_table(x, sa)
o_table(x, sa)
fm_search(x, p, sa)
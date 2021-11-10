import numpy as np
import gen_lcp, parsers
import sys, argparse

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
    n = len(sa)

    O = np.zeros((n + 1, alphabet_size))
    
    for i in range(1, n + 1):
        for a in range(0, alphabet_size):
            if x[sa[i - 1] - 1] == alphabet[a]:
                O[i][a] = O[i - 1][a] + 1
            else: 
                O[i][a] = O[i - 1][a]

    return O

def fm_search(x, p, sa):
    alphabet = {a:i for i, a in enumerate(sorted(set(x)))} # To find the correct column in the O table
    O = o_table(x, sa)
    _, c_buckets = c_table(x, sa)
    matches = None
    L, R = 0, len(sa)
    
    for a in reversed(p):
        if L == R: break
        L = c_buckets[a] + O[int(L)][alphabet[a]]
        R = c_buckets[a] + O[int(R)][alphabet[a]]
    
    if L != R:
        matches = sa[int(L):int(R)]
    
    return matches

# Wrapper function
def search_fm(sa, fastq):

    if len(sa) < 0 or len(fastq) < 0:
        return "Problems with either fastq file or the SA and LCP"

    flag, mapq, pnext, tlen = 0,0,0,0
    rnext = "*"

    for x in sa.items():
        rname = x[0]
        y = x[1][0] + "$"
        sa = x[1][1]

        for p in fastq.items():
            qname = p[0]
            substring = p[1][0]
            cigar = str(len(substring)) + "M"
            qual = p[1][1]

            matches = fm_search(y, substring, sa)

            if matches is not None:
                for match in matches:
                    pos = int(match) + 1
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)

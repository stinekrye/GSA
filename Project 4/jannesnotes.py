import numpy as np
import gen_sa, parsers
import sys, argparse

def c_table(x, sa):
    buckets = {}

    for i in range(0, len(sa)):
        char = x[sa[i]]
        if char not in buckets:
            buckets[char] = i

    return buckets

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

def fm_search(O, C, p, sa, alpha):
    matches = None
    L, R = 0, len(sa)
    
    for a in reversed(p):
        if L == R: break
        L = C[a] + O[int(L)][alpha[a]]
        R = C[a] + O[int(R)][alpha[a]]
    
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

            alpha = {a:i for i, a in enumerate(sorted(set(y)))}
            O = o_table(y, sa)
            C = c_table(y, sa)

            matches = fm_search(O, C, substring, sa, alpha)

            if matches is not None:
                for match in matches:
                    pos = int(match) + 1
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)


sa = {
    "one":["mississippi", [11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]], 
    "two":["mississippimississippi", [22, 21, 10, 18, 7, 15, 4, 12, 1, 11, 0, 20, 9, 19, 8, 17, 6, 14, 3, 16, 5, 13, 2]]}

fastq_dict = {
    "iss":["iss", "~~~"], 
    "mis":["mis", "~~~"], 
    "ssi":["ssi", "~~~"], 
    "ssippi":["ssippi", "~~~~~~"]}

# x = "mississippi$"
# p = "iss"
# sa = [11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]
# alpha = {a:i for i, a in enumerate(sorted(set(x)))}
# O = o_table(x, sa)
# C = c_table(x, sa)
# fm_search(O, C, p, sa, alpha)
search_fm(sa, fastq_dict)

import numpy as np
import gen_sa, parsers
import sys, argparse, json

def c_table(sa_dict, fastaname):
    f = open(f"{fastaname}.c-table", "w")
    buckets = None

    for key, value in sa_dict.items():
        buckets = {}
        sa = value[1]
        seq = value[0] + "$"
        
        for i in range(0, len(sa)):
            char = seq[sa[i]]
            if char not in buckets:
                buckets[char] = i
        
        f.write(">" + str(key) + "\n")
        f.write(str(buckets) + "\n")
    f.close()

    return "Done"

def o_table(sa_dict, fastaname):
    f = open(f"{fastaname}.o-table", "w")

    for key, value in sa_dict.items():
        sa = value[1]
        seq = value[0] + "$"

        alphabet = sorted(set(seq))
        alphabet_size = len(alphabet)
        n = len(sa)

        O = np.zeros((n + 1, alphabet_size))
        
        for i in range(1, n + 1):
            for a in range(0, alphabet_size):
                if seq[sa[i - 1] - 1] == alphabet[a]:
                    O[i][a] = O[i - 1][a] + 1
                else: 
                    O[i][a] = O[i - 1][a]
        
        f.write(">" + str(key) + "\n")
        f.write(str(O.tolist()) + "\n")

    f.close()

    return "Done"

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

def search_fm(sa, fastq, o_dict, c_dict):

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
            O = o_dict[rname]
            C = c_dict[rname]

            matches = fm_search(O, C, substring, sa, alpha)

            if matches is not None:
                for match in matches:
                    pos = int(match) + 1
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)
    
    return 


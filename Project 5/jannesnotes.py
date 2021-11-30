import numpy as np
import gen_sa, parsers
import sys, argparse, itertools

def c_table(sa):
    buckets = None

    for key, value in sa.items():
        buckets = {}
        sa = value[1]
        seq = value[0] + "$"
        
        for i in range(0, len(sa)):
            char = seq[sa[i]]
            if char not in buckets:
                buckets[char] = i

    return buckets

def o_table(sa):

    for key, value in sa.items():
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

    return O

def d_table(RO, C, sa, p, alpha):
    min_edits = 0
    m = len(p)
    L, R = 0, len(sa)
    d = np.zeros(m, dtype=int)

    for i in range(m):
        a = p[i]
        L = C[a] + RO[int(L)][alpha[a]]
        R = C[a] + RO[int(R)][alpha[a]]
        if L >= R: 
            min_edits += 1
            L = 0
            R = len(sa)
        d[i] = min_edits

    return d

def bw_approx(O, C, p, d, sa, alpha, max_edits):
    matches = []
    L, R = 0, len(sa)
    i = len(p) - 1
    cigar = np.zeros(len(p)+max_edits, dtype=str)
    
    # Match
    c_index = 0
    for a in list(alpha.keys())[1:]:
        new_L = C[a] + O[int(L)][alpha[a]]
        new_R = C[a] + O[int(R)][alpha[a]]
        
        # Match or mismatch
        if a == p[i]: 
            edit_cost = 0
        else:
            edit_cost = 1
        
        if max_edits - edit_cost < 0: continue
        if new_L >= new_R: continue

        cigar[c_index] = 'M'
        rec_approx(
            d, sa, O, C, cigar, new_L, new_R, 
            i - 1, max_edits - edit_cost, c_index + 1)
        

    # Insertion
    cigar[c_index] = 'I'
    return rec_approx(
        d, sa, O, C, cigar, L, R, 
        i - 1, max_edits - 1, c_index + 1)
    
    return matches, cigar

def rec_approx(d, sa, O, C, cigar, L, R, i, edits_left, c_index):
    lower_limit = d[i]
    if edits_left < lower_limit:
        return None
    if i < 0: # Means we have a match
        matches = sa[int(L):int(R)]
        print(matches, cigar[:c_index])
        return
    
    
    for a in list(alpha.keys())[1:]:
        new_L = C[a] + O[int(L)][alpha[a]]
        new_R = C[a] + O[int(R)][alpha[a]]

        if a == p[i]:
            edit_cost = 0
        else:
            edit_cost = 1
        if edits_left - edit_cost < 0: continue
        if new_L >= new_R: continue

        cigar[c_index] = 'M'
        rec_approx(
            d, sa, O, C, cigar, new_L, new_R, 
            i - 1, edits_left - edit_cost, c_index + 1)
    
    # Insertion
    cigar[c_index] = 'I'
    rec_approx(
        d, sa, O, C, cigar, L, R, 
        i - 1, edits_left - 1, c_index + 1)

    #  Deletion
    cigar[c_index] = 'D'
    for a in list(alpha.keys())[1:]:
        new_L = C[a] + O[int(L)][alpha[a]]
        new_R = C[a] + O[int(R)][alpha[a]]
        if new_L >= new_R: continue

        rec_approx(
            d, sa, O, C, cigar, new_L, new_R, 
            i, edits_left - 1, c_index + 1)

    return

x = "mississippi$"
sa_dict = {"one":["mississippi", [11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]]}
rsa_dict = {"one":["mississippi", [0, 10, 1, 7, 4, 11, 3, 2, 9, 6, 8, 5]]}

p = "is"
sa = [11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]
alpha = {a:i for i, a in enumerate(sorted(set(x)))}
O = o_table(sa_dict)
C = c_table(sa_dict)
RO = o_table(rsa_dict)
d = d_table(RO, C, sa, p, alpha)

bw_approx(O, C, p, d, sa, alpha, 1)


def search_bw(sa, fastq, o_dict, c_dict, ro_dict):

    if len(sa) < 0 or len(fastq) < 0:
        return "Problems with either fastq file or the SA and LCP"

    flag, mapq, pnext, tlen = 0,0,0,0
    rnext = "*"

    for x in sa.items():
        rname = x[0]
        y = x[1][0] + "$"
        sa = x[1][1]

        O = o_dict[rname]
        RO = ro_dict[rname]
        C = c_dict[rname]

        for p in fastq.items():
            qname = p[0]
            substring = p[1][0]
            cigar = str(len(substring)) + "M"
            qual = p[1][1]

            alpha = {a:i for i, a in enumerate(sorted(set(y)))}
            
            d = d_table(RO)
            matches, cigars = bw_approx(O, C, d, substring, sa, alpha)

            if matches is not None:
                for match, cigar in itertools.zip_longest(matches, cigars):
                    pos = int(match) + 1
                    cigar = ''.join([f'{value}{key}'for key, value in cigars.items()])
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)
    
    return 

# # Creating first parser
# parser1 = argparse.ArgumentParser(description='SA computation using SAIS')
# parser1.add_argument('-p', help='Create SA from fastafile')
# args1 = parser1.parse_known_args()

# if args1[0].p:
#     fastafile = parsers.read_fasta_file(args1[0].p)
#     fastaname = f"{args1[0].p}"

#     # Make suffix array and read file
#     gen_sa.gen_sa(fastafile, fastaname) # Use Stine's new function here instead
#     sa = parsers.read_SA(f"{fastaname}.sa")
#     rsa = parsers.read_SA(f"{fastaname}.rsa") # Make sure that it also creates the sa for the reversed string

#     # Make C, O, and RO table
#     c_table(sa, fastaname)
#     o_table(sa, fastaname)
#     o_table(rsa, fastaname)


# else:
#     # Creating second parser if -p is not given
#     parser2 = argparse.ArgumentParser(description='Pattern matching using suffix tree')
#     parser2.add_argument('fastafile', help="Input fasta file")
#     parser2.add_argument('fastqfile', help="Input fastq file")
#     args2 = parser2.parse_args()
#     fastq = parsers.read_fastq_file(args2.fastqfile)
#     sa = parsers.read_SA(f"{args2.fastafile}.sa")
    
#     o_dict = parsers.read_o(f"{args2.fastafile}.o-table")
#     ro_dict = parsers.read_o(f"{args2.fastafile}.ro-table")
#     c_dict = parsers.read_c(f"{args2.fastafile}.c-table")
    
#     search_bw(sa, fastq, o_dict, c_dict, ro_dict)
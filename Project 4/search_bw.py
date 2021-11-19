
import numpy as np
import gen_sa, parsers
import sys, argparse

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

#### RUNNING THE SCRIPT
# If the -p option is given with a fastafile (e.g. "python search_st2.py -p test.fa"), 
# the script will construct a suffix tree, create the SA and LCP from it, and output the 
# two arrays to a textfile named "test.fa.sa-lcp".
# If the script is not given the -p option, it takes two input files: a fastq file and a 
# fasta file name. From the fasta file name it finds a file that has the fasta name as 
# prefix and a suffix of .sa-lcp. From this it searches for the pattern found in the fastq 
# file.

# Creating first parser
parser1 = argparse.ArgumentParser(description='SA computation from suffix tree')
parser1.add_argument('-p', help='Create SA from fastafile')
args1 = parser1.parse_known_args()

if args1[0].p:
    fastafile = parsers.read_fasta_file(args1[0].p)
    fastaname = f"{args1[0].p}"

    # Make suffix array and read file
    gen_sa.gen_sa(fastafile, fastaname)
    sa = parsers.read_SA(f"{fastaname}.sa")

    # Make C and O table
    c_table(sa, fastaname)
    o_table(sa, fastaname)


else:
    # Creating second parser if -p is not given
    parser2 = argparse.ArgumentParser(description='Pattern matching using suffix tree')
    parser2.add_argument('fastafile', help="Input fasta file")
    parser2.add_argument('fastqfile', help="Input fastq file")
    args2 = parser2.parse_args()
    fastq = parsers.read_fastq_file(args2.fastqfile)
    sa = parsers.read_SA(f"{args2.fastafile}.sa")
    
    o_dict = parsers.read_o(f"{args2.fastafile}.o-table")
    c_dict = parsers.read_c(f"{args2.fastafile}.c-table")
    
    search_fm(sa, fastq, o_dict, c_dict)


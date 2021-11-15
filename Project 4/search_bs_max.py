import sys
import argparse
import gen_lcp, parsers

def binary_search(x, p, sa):
    low, high = 0, len(sa)
    while low < high:
        mid = low + (high - low) // 2
        if x[sa[mid]:(sa[mid]+len(p))] == p:
            lower_bound = lower_bound_search(x, p, sa, low, mid)
            upper_bound = upper_bound_search(x, p, sa, high, mid)
            occ = sa[int(lower_bound):int(upper_bound)]
            return occ
        elif x[sa[mid]:(sa[mid]+len(p))] > p:
            high = mid - 1
        else:
            low = mid + 1
    mid = low + (high - low) // 2
    lower_bound = lower_bound_search(x, p, sa, low, mid)
    upper_bound = upper_bound_search(x, p, sa, high, mid)
    occ = sa[int(lower_bound):int(upper_bound)]
    return occ

def lower_bound_search(x, p, sa, low, mid):
    while low < mid:
        middle = (low + mid) // 2
        if x[sa[middle]:(sa[middle] + len(p))] >= p:
            mid = middle
        elif x[sa[middle]:(sa[middle] + len(p))] < p:
            low = middle + 1
    return low

def upper_bound_search(x, p, sa, high, mid):
    while mid < high:
        middle = (mid + high) // 2
        if x[sa[middle]:(sa[middle] + len(p))] > p:
            high = middle
        elif x[sa[middle]:(sa[middle] + len(p))] <= p:
            mid = middle + 1
    if high == len(p):
        return high
    else:
        return high + 1


def search_bs(sa, fastq):

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

            matches = binary_search(y, substring, sa)

            if matches is not None:
                for match in matches:
                    pos = int(match) + 1
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)

#### RUNNING THE SCRIPT
# If the -p option is given with a fastafile (e.g. "python search_bs.py -p test.fa"),
# the script will construct a suffix tree, create the SA and LCP from it, and output the
# two arrays to a textfile named "test.fa.sa-lcp".
# If the script is not given the -p option, it takes two input files: a fasta file and a
# fastq file name. From the fasta file name it finds a file that has the fasta name as
# prefix and a suffix of .sa-lcp. From this it searches for the pattern found in the fastq
# file.

# Creating first parser
parser1 = argparse.ArgumentParser(description='SA computation from suffix tree')
parser1.add_argument('-p', help='Create SA from fastafile')
args1 = parser1.parse_known_args()

if args1[0].p:
    fastafile = parsers.read_fasta_file(args1[0].p)
    fastaname = f"{args1[0].p}"
    gen_lcp.gen_lcp(fastafile, fastaname)

else:
    # Creating second parser if -p is not given
    parser2 = argparse.ArgumentParser(description='Pattern matching using suffix tree')
    parser2.add_argument('fastafile', help="Input fasta file")
    parser2.add_argument('fastqfile', help="Input fastq file")
    args2 = parser2.parse_args()
    sa = parsers.read_SA_LCP(f"{args2.fastafile}.sa-lcp")
    fastq = parsers.read_fastq_file(args2.fastqfile)
    search_bs(sa, fastq)
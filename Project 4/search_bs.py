import gen_lcp, parsers
import sys, argparse

def binary1(p,x,sa,k, l, u): # returns a match in the SA
    low = l
    high = u
    mid = low + ((high-low)//2)

    while high-low != 1:
        if p[k] == x[sa[mid]+k:sa[mid]+k+1]:
            return mid, high, low

        elif p[k] > x[sa[mid]+k:sa[mid]+k+1]:
            low = mid
            mid = low + (high-low)//2

        else:
            high = mid
            mid = low + ((high-low)//2)
def binary2(p,x,sa,k, mid, l, u, upper):
    if upper == False: # If we are searching for lower bound
        low = l
        high = mid
        mid = low + ((high-low) // 2)
        while high - low != 1:
            if p[k] == x[sa[mid] + k:sa[mid] + k + 1]:
                high = mid
                mid = low + ((high-low)//2)
            else:  # p[k] > x[sa[mid] + k:sa[mid] + k + 1] #Try to draw this. We do not need to consider the case where p[k] < x[....]
                low = mid
                mid = low + ((high-low)//2)
            # else:
            #     high = mid
            #     mid = low + ((high-low)//2)


    else:              # If we are searching for upper bound
        low = mid
        high = u
        mid = low + ((high-low) // 2)
        while high - low != 1:
            if p[k] == x[sa[mid] + k:sa[mid] + k + 1]:
                low = mid
                mid = low + ((high-low)//2)
            # elif p[k] > x[sa[mid] + k:sa[mid] + k + 1]:
            #     low = mid
            #     mid = low + ((high-low)//2)
            else: # p[k] < x[.......]
                high = mid
                mid = low + ((high-low)//2)
    return mid
def binary3(p,x,sa):
    l = 0
    u = len(sa)
    k = 0
    while k < len(p):
        first_match, high, low = binary1(p,x,sa,k, l, u)
        if first_match:
            l = binary2(p, x, sa, k, first_match, low, high, upper=False)        # Points one past the interval
            u = binary2(p, x, sa, k, first_match, low, high, upper=True) + 1     # Points one past the interval
            if k == len(p)-1:
                return sa[l+1:u]
            else:
                k += 1
        else:
            return "No match"

    else:
        return "No match"
def search_bs(sa, fastq, test = False):

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

            matches = binary3(substring, y, sa)
            # matches = sorted(matches)
            if test == False:
                if matches is not None:
                    for match in matches:
                        pos = int(match) + 1
                        print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)
            if test == True:
                if matches is not None:
                    for match in matches:
                        pos = int(match) + 1
                        yield f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}"

def print_test(iter):
    iter = sorted(iter)
    for i in iter:
        print(i, file = sys.stdout)

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

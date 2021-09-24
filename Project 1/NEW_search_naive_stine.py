from read_fasta import read_fasta_file
from read_fastq import read_fastq_file
import argparse, sys

# Initialize parser
parser = argparse.ArgumentParser(description='A naive approach for exact pattern matching')

# Add arguments
parser.add_argument(
    'fastafile',
    help="Input fasta file"
)
parser.add_argument(
    'fastqfile',
    help="Input fastq file"
)

args = parser.parse_args()

# Read in files
fasta = read_fasta_file(args.fastafile)
fastq = read_fastq_file(args.fastqfile)



# Search function
def naive_matching(x, p):
    i, j = 0, 0
    while j <= (len(x) - len(p)):  # break when we reach the end of the sequence
        if p[i] == x[j]:
            while i < len(p) and p[i] == x[j]:
                j += 1
                i += 1

            if i == len(p):  # if a match is found
                yield j-i,x[j-i:j+1]
                j = j-i+1  # We want to look at the next letter in j, from where we found the previous match
                i = 0

            else:
                j = j-i+1  # We want to look at the next letter in j, from where we found the previous match
                i = 0

        else:
            j += 1



# Wrapper function
def search_naive(fasta,fastq):

    if len(fasta) < 0 or len(fastq) < 0:
        return "Problems with either fasta or fastq file"

    flag, mapq, pnext, tlen = 0,0,0,0
    rnext = "*"

    for p in fastq.items():
        for x in fasta.items():
            qname = p[0]
            rname = x[0]
            substring = p[1][0]
            cigar = str(len(substring)) + "M"
            qual = p[1][1]
            matches = naive_matching(x[1], substring)

            for match in matches:
                pos = match[0]+1
                print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}\t", file = sys.stdout)

search_naive(fasta,fastq)
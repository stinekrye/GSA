from read_fasta import read_fasta_file
from read_fastq import read_fastq_file
import re, argparse


# Initialize parser
parser = argparse.ArgumentParser(description='A naive approach for exact pattern matching')

# Add arguments
fasta = parser.add_argument(
    'fastafile',
    help="Input fasta file"
)
fastq = parser.add_argument(
    'fastqfile',
    help="Input fastq file"
)

args = parser.parse_args()

# Read in files
fasta = read_fasta_file(args.fastafile)
fastq = read_fastq_file(args.fastqfile)

# Search function
def naive_matching_ez(x, p):

    # Definitions for the SAM format
    qname = p[0]
    flag = 0
    rname = x[0]
    pos = 0
    mapq = 0
    substring = p[1][0]
    cigar = str(len(substring)) + 'M'
    rnext = '*'
    pnext = 0
    tlen = 0
    qual = p[1][1]

    # Other definitions
    rseq = x[1]

    # Find exact match
    for match in re.finditer(substring, rseq):
        print(qname, flag, rname, match.start() + 1, mapq, cigar, rnext, pnext, tlen, substring, qual)

for item in fastq.items():
    for seqs in fasta.items():
        naive_matching_ez(seqs, item)
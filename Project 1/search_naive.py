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

# # Read in files
fasta = read_fasta_file(args.fastafile)
fastq = read_fastq_file(args.fastqfile)

# Search function
def naive_matching(x, p):

    if len(x) > 0 and len(p) > 0:
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
        j, i = 0, 0
        currentstr = ''

        # Find exact match
        while j < len(rseq):
            if substring[i] != rseq[j]:

                # Reset i, currentstr, and pos, but increment j if we haven't found any matches yet
                if len(currentstr) == 0:
                    j += 1
                    i = 0
                    currentstr = ''
                    pos = 0
                # If we have found some matches, don't increment j
                else: 
                    i = 0
                    currentstr = ''
                    pos = 0

            elif substring[i] == rseq[j]:

                # If it's the first match, update the pos variable
                if pos == 0:
                    pos = j + 1

                # Update currentstr
                currentstr = currentstr + rseq[j]

                if currentstr == substring:

                    # If we've matched the full substring, print the SAM information
                    print(
                        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, substring, qual,
                        file = sys.stdout
                    )

                    # And reset i, currentstr, and pos
                    j = pos
                    i = 0
                    pos = 0
                    currentstr = ''
                
                else:
                    # If we haven't yet found the full substring, increment i and j and continue matching
                    i += 1
                    j += 1
    else:
        print('One or both input files are empty')

for item in fastq.items():
    for seq in fasta.items():
        naive_matching(seq, item)
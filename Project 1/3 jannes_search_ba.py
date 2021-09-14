from read_fasta import read_fasta_file
from read_fastq import read_fastq_file
import argparse
 
# Initialize parser
parser = argparse.ArgumentParser(description='Using border arrays for exact pattern matching')

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

# Border array function
def borderarray(x):
    return

# Search function
def search_ba(x):
    return

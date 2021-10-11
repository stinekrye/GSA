from parsers.read_fasta import read_fasta_file
from parsers.read_fastq import read_fastq_file
import argparse, sys

# Definitions for argparse
parser = argparse.ArgumentParser(description='Pattern matching using suffix tree construction')

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

# Function to remap sequences from bases to numbers
def remap(x):
    m = {a:i+1 for i,a in enumerate(sorted(set(x)))}
    n = [m[a] for a in x]
    return n


class Node():

    def __init__(self, range_start = None, range_end = None, string_label = None, children = [None]*5, parent = None):
        self.range_start = range_start
        self.range_end = range_end
        self.string_label = string_label
        self.children = children # Is defined later
        self.parent = parent

    def __repr__(self):
        rep = f"Node({self.range_start},{self.range_end})"
        return rep

class SuffixTree():

    def __init__(self, seq = None):
        self.seq = seq
        self.root = Node()
        self.root.parent = self.root
    
    # Fix this
    def find_outgoing_edge(self, v, x, i):
        return

    def split_edge(self, v, i, k):
        # Make U
        u = Node(range_start = v.range_start, range_end = i + k, parent = v.parent)
        
        # Update p(v)
        v.parent.children = u

        # Update v
        v.parent = u
        v.range_start = i + k

        # Add v as child of u                                   
        u.children = v  

        return u

    def insert_child(self, v, i, x, xend):
        new = Node(range_start = x + i, range_end = xend, string_label = i, parent = v)
        return new

    # Stine's function
    def naive_insert(self, c, x, i, n):
        child = self.find_outgoing_edge(c, x, i)
        if child:  # If we have an outgoing edge
            while x + i < child.range_end and self.seq[child.range_start + x] == self.seq[i+x]:  # If we have have not reached the end of the edge yet and we can extend the match
                x += 1
            if x + i == child.range_end:  # if we reach the end of the edge, we have to find the new outgoing edge to search from
                self.naive_insert(child, x, i, n)
            else:
                u = self.split_edge(child, x, i, n)
                self.insert_child(u, x, i, n)


        else:  # If we have not got an outgoing edge
            self.insert_child(c, x, i, n)
    
    # Function to insert the first suffix manually, and then the rest using naive_insert 
    def naive_st(self, str):
        # Add sentinel

        # Remap sequence, define x and xend
        self.seq = remap(str)
        x = 0
        xend = len(str)

        # Inserting first suffix
        first_node = Node(range_start = x, range_end = xend, parent = self.root, string_label = 0) # Not sure about label

        # Inserting the rest
        for i in range(1, xend):
            self.naive_insert(self.root, x, i, xend)

        return


x = "ATGGACTTACGCGGTTACCAATCAAGTA"
st = SuffixTree()
st.naive_st(x)
print(st)

def st_matching(x, p):

    """
    A function which searches for patterns in a string using a suffix tree construction

    Args:
        x: A string of length n
        p: A string of length m
    
    Yields:
        The start position of the match in x
    """

    st = SuffixTree()
    j = 0

    while j < len(p):
        if "p[j] matches whereever we are in the st":
            j += 1
        # Something more
        if j == len(p):  # if a match is found
            yield j - "something" # To get the position in x where p starts

# Stine's match function 
def st_search(fasta, fastq):

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
            matches = st_matching(x[1], substring)

            for match in matches:
                pos = match[0]+1
                print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)
    
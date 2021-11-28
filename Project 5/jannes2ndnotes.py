import gen_lcp, parsers
import sys, argparse, itertools

def remap(x):
    m = {a:i for i,a in enumerate(sorted(set(x)))}
    n = [m[a] for a in x]
    return n, m

def return_labels(node):
    stack = [node]
    # stack.append(node)

    while len(stack) > 0:
        if stack[-1].string_label != None:
            yield stack[-1].string_label
            stack.pop(-1)
            # p.parent = None
            # p.children = None

        else:
            parent = stack.pop(-1)
            for child in parent.children:
                if child != None:
                    stack.append(child)
            # parent.parent = None
            # parent.children = None

class Node():
    def __init__(self, range_start = None, range_end = None, string_label = None, children = [], parent = None):
        self.range_start = range_start
        self.range_end = range_end
        self.string_label = string_label
        self.children = children
        self.parent = parent

    def __repr__(self):
        rep = f"Node({self.range_start},{self.range_end})"
        return rep

class SuffixTree2():

    def __init__(self, seq = None):
        self.seq = seq
        self.root = Node(range_start = 0, range_end = 0)
        self.n = 0
        self.alpha = None

    def append_child(self, v, sa, lcp, length):
        new_leaf = Node(
            range_start=sa+lcp,
            range_end=length,
            string_label=sa,
            parent = v,
            children=[]
        )
        v.children.append(new_leaf)

        return new_leaf

    def split_edge(self, v, length_up):                       
        # New node
        u = Node(
            range_start = v.range_start, 
            range_end = v.range_end - length_up, 
            parent = v.parent, 
            children=[],
            string_label=None
            )

        # Adding u to the children list of v's parent and removing v from it
        v.parent.children.append(u)
        v.parent.children.remove(v)

        # Update v  
        v.parent = u                                                                    
        v.range_start = v.range_end - length_up

        # Add v as child of u                                                                                         
        u.children.append(v)

        return u
    
    def lcp_insert(self, i, sa, lcp, v, length):
        length_up = length - sa[i-1] - lcp[i]
        v_edge_length = v.range_end - v.range_start

        while length_up >= v_edge_length and length_up != 0:
            length_up -= v_edge_length
            v = v.parent
            v_edge_length = v.range_end - v.range_start
        
        if length_up == 0:
            new_leaf = self.append_child(v, sa[i], lcp[i], length)
        else: 
            u = self.split_edge(v, length_up)
            new_leaf = self.append_child(u, sa[i], lcp[i], length)

        return new_leaf

    def lcp_suffix_tree(self, sa, lcp):
        seq_length = len(sa)

        # Initialize with SA[0]
        v = Node(
            range_start=sa[0], 
            range_end=seq_length, 
            string_label=sa[0], 
            children=[], 
            parent=self.root)

        # Make it child of root
        self.root.children.append(v)

        # Insert the rest
        for i in range(1, seq_length):
            v = self.lcp_insert(i, sa, lcp, v, seq_length)

    def find_edge_search(self, c, char):
        for child in c.children:
            if self.seq[child.range_start] == char:
                return child

    def recurse_children():
        return

    def approx_search(self, node, x, k, n, cigar, edits):
        
        child = self.find_edge_search(node, n[x])
        if child:  # If we have an outgoing edge
            k = child.range_start
            while k < child.range_end and x < len(n) and self.seq[k] == n[x]:  # self.root.range_start has to be the beginning of the match. In other words: We have to begin at each node instead of the root
                x += 1
                k += 1

                if 'M' not in cigar:
                    cigar['M'] = 1
                else:
                    cigar['M'] += 1

            if edits < 0: # We've used too many edits
                return None
            
            elif x == len(n):  # We've found a match
                res = return_labels(child)
                return res, cigar  # Return both the interval and the cigar
            
            elif k >= child.range_end and x < len(n):  # if we reach the end of the edge, and we have not found a match, we have to search down a new node
                return self.approx_search(child, x, k, n, cigar, edits) # We need a recurse_children here
            
            elif self.seq[k] != n[x]:
                if 'D' not in cigar:
                    cigar['D'] = 1
                else:
                    cigar['D'] += 1

                return self.approx_search(child, x, k + 1, n, cigar, edits - 1)
            elif self.seq[k] != n[x]:
                if 'I' not in cigar:
                    cigar['I'] = 1
                else:
                    cigar['I'] += 1

                return self.approx_search(child, x + 1, k, n, cigar, edits - 1)
            else:
                return None
        else:
            return None

# Wrapper function
def search_suffix(sa_lcp, fastq, no_edits):

    if len(sa_lcp) < 0 or len(fastq) < 0:
        return "Problems with either fastq file or the SA and LCP"

    flag, mapq, pnext, tlen = 0,0,0,0
    rnext = "*"

    for x in sa_lcp.items():
        rname = x[0]
        tree = None
        tree = SuffixTree2()
        tree.root.children = []
        y = x[1][0] + "0"
        tree.seq, tree.alpha = remap(y)
        sa = x[1][1]
        lcp = x[1][2]
        tree.lcp_suffix_tree(sa, lcp)


        for p in fastq.items():
            qname = p[0]
            substring = p[1][0]
            qual = p[1][1]

            seq = [tree.alpha[a] for a in substring]
            cigar = {}

            matches, cigars = tree.approx_search(tree.root, 0, 0, seq, cigar, no_edits)


            if matches is not None:
                for match, cigar in itertools.zip_longest(matches, cigars):
                    pos = int(match) + 1
                    cigar = ''.join([f'{value}{key}'for key, value in cigars.items()])
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)
sa_lcp = {"seq1": ["mississippi", [11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2], [0, 0, 1, 1, 4, 0, 0, 1, 0, 2, 1, 3]]}
fastq = {"isi": ["isi", "~~~"]}
search_suffix(sa_lcp, fastq, 1)

#### RUNNING THE SCRIPT
# If the -p option is given with a fastafile (e.g. "python search_st2.py -p test.fasta"), 
# the script will construct a suffix tree, create the SA and LCP from it, and output the 
# two arrays to a textfile named "SA_LCP_name". Where the name is the sequence name from 
# the fastafile.
# If the script is not given the -p option, it takes two input files: a fastq file and a 
# text file containing the SA and LCP. From this it constructs a suffix tree and searches
# for the pattern found in the fastq file.

# # Creating first parser
# parser1 = argparse.ArgumentParser(description='SA and LCP computation from suffix tree')
# parser1.add_argument('-p', help='Create SA and LCP from fastafile')
# args1 = parser1.parse_known_args()

# if args1[0].p:
#     fastafile = parsers.read_fasta_file(args1[0].p)
#     fastaname = f"{args1[0].p}"
#     gen_lcp.gen_lcp(fastafile, fastaname)

# else:
#     # Creating second parser if -p is not given
#     parser2 = argparse.ArgumentParser(description='Pattern matching using suffix tree')
#     parser2.add_argument('fastafile', help="Input fasta file")
#     parser2.add_argument('fastqfile', help="Input fastq file")
#     args2 = parser2.parse_args()
#     sa_lcp = parsers.read_SA_LCP(f"{args2.fastafile}.sa-lcp")
#     fastq = parsers.read_fastq_file(args2.fastqfile)
#     search_suffix(sa_lcp, fastq)

from parsers.read_fasta import read_fasta_file
from parsers.read_fastq import read_fastq_file
import sys, argparse

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

class SuffixTree():

    def __init__(self, seq = None):
        self.seq = seq
        self.root = Node(range_start = 0, range_end = 0)
        self.n = 0
        # self.alpha = None

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
            self.seq[child.range_start] == char
            return child

    def search_rec(self,c,x,k, n):

        child = self.find_edge_search(c, n[x])
        if child:  # If we have an outgoing edge
            k = child.range_start
            while k < child.range_end and x < len(n) and self.seq[k] == n[x]:  # self.root.range_start has to be the beginning of the match. In other words: We have to begin at each node instead of the root
                x += 1
                k += 1
            if k >= child.range_end and x < len(n):  # if we reach the end of the edge, and we have not found a match, we have to search down a new node
                return self.search_rec(child, x, k, n)
            elif x == len(n):  # if we reach the end of an edge and we find a match (also works on edges
                res = return_labels(child)
                return res  # Be carefull with this recursion
            else:
                return None
        else:
            return None


# Wrapper function
def search_suffix(fasta,fastq):

    if len(fasta) < 0 or len(fastq) < 0:
        return "Problems with either fasta or fastq file"

    flag, mapq, pnext, tlen = 0,0,0,0
    rnext = "*"

    for x in fasta.items():
        rname = x[0]
        tree = None
        tree = SuffixTree()
        tree.root.children = [None, None, None, None, None]
        tree.insert(x[1])


        for p in fastq.items():
            qname = p[0]
            substring = p[1][0]
            cigar = str(len(substring)) + "M"
            qual = p[1][1]

            seq = [tree.alpha[a] for a in substring]

            matches = tree.search_rec(tree.root, 0, 0, seq)

            if matches is not None:
                for match in matches:
                    pos = int(match) + 1
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)

# fasta = read_fasta_file("fasta_test.fasta")
# fastq = read_fastq_file("fastq_test.fastq")

# sa = [8, 7, 4, 0, 5, 2, 1, 6, 3]
# lcp = [0, 0, 1, 0, 1, 3, 0, 0, 2]
sa = [11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]
lcp = [0, 0, 1, 1, 4, 0, 0, 1, 0, 2, 1, 3]
st = SuffixTree()
st.lcp_suffix_tree(sa, lcp)
print(st.find_edge_search(Node(10, 11), 8))

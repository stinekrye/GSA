from parsers.read_fasta import read_fasta_file
from parsers.read_fastq import read_fastq_file
import argparse, sys

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
    def __init__(self, range_start = None, range_end = None, string_label = None, children = [None,None,None,None,None], parent = None):
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
        self.root = Node(range_start = 0, range_end = 0)
        self.n = 0
        self.alpha = None
       # self.root.parent = self.root

    def find_edge(self,c,x,i):
        char = self.seq[x+i]
        return c.children[char]

    def insert_child(self,c,x,i):
        char = self.seq[x + i]
        new = Node(range_start = i + x, range_end = len(self.seq), string_label=i, parent = c, children=[None, None, None, None, None])
        c.children[char] = new
        self.n += 1
        return


    def split_edge(self,c,x):
        c_char = self.seq[c.range_start]                                                # first character in c
        u = Node(range_start = c.range_start, range_end = c.range_start + x, parent = c.parent, children=[None, None, None, None, None])    # create u node
        c.parent.children[c_char] = u                                                   # update the child of P(c) to be u. The first character in u is c_char
        c.parent = u                                                                    # set the parent of c to be the new node u
        c.range_start = u.range_start + x                                                           # update the range of c to start where u ends
        c_char = self.seq[c.range_start]                                                # find the new beginning character of c
        u.children[c_char] = c                                                          # set c to be the child of u at the right place in the list using the first character of c
        self.n += 1
        return u

    def naive_insert(self, c, x, k, i):
        child = self.find_edge(c, x, i)
        if child:  # If we have an outgoing edge
            while k + child.range_start < child.range_end and self.seq[child.range_start + k] == self.seq[i+x]:  # If we have have not reached the end of the edge yet and we can extend the match
                x += 1
                k += 1
            if k + child.range_start == child.range_end:  # if we reach the end of the edge, we have to find the new outgoing edge to search from
                k = 0
                self.naive_insert(child, x,k, i) # Be carefull with this recursion
            else:
                u = self.split_edge(child, k)
                self.insert_child(u, x, i)

        else:  # If we have not got an outgoing edge
            self.insert_child(c, x, i)


    def insert(self, y):
        y = y + "0"
        self.seq,self.alpha = remap(y)

        x = 0
        i = 0
        k = 0

        # insert first sequence manually
        char = self.seq[x + i]
        new = Node(range_start = i + x, range_end = len(self.seq), string_label = i, parent = self.root, children = [None, None, None, None, None])
        self.root.children[char] = new
        i += 1


        while i < len(y):                                                                    # Change this to iterate through the full sequence when done
            self.naive_insert(self.root,x,k,i)                                                  # input = suffix tree, root to search from, x (incrementor), i = string label. n+1 = length of sequence
            i += 1


    def find_edge_search(self, c, char):
        return c.children[char]

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





fasta = read_fasta_file("fasta_test.fasta")
fastq = read_fastq_file("fastq_test.fastq")

search_suffix(fasta,fastq)



# #
# res = SuffixTree()
# # res.insert("MISSISSIPPI")
# res.insert("GTAGTA")
# print(res)
# n = res.search("GTA")
# for i in n:
#     print(i)
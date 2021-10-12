# Where to put this function?
import graphviz

def remap(x):
    m = {a:i for i,a in enumerate(sorted(set(x)))}
    n = [m[a] for a in x]
    return n, m

def return_labels(node):
    if node.string_label != None:
        a = node.string_label
        yield a
    else:
        for child in node.children:
            if child != None:
                if child.string_label != None:
                    yield child.string_label
                else:
                    return_labels(child)

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
 #       self.root.parent = self.root

    def find_edge(self,c,x,i):
        char = self.seq[x+i]
        return c.children[char]

    def insert_child(self,c,x,i,n):
        char = self.seq[x + i]
        new = Node(range_start=i + x, range_end=n, string_label=i, parent = c, children=[None, None, None, None, None])
        c.children[char] = new
        self.n += 1
        return


    def split_edge(self,c,x,i,n):
        c_char = self.seq[c.range_start]                                                # first character in c
        u = Node(range_start = c.range_start, range_end = c.range_start + x, parent = c.parent, children=[None, None, None, None, None])    # create u node
        c.parent.children[c_char] = u                                                   # update the child of P(c) to be u. The first character in u is c_char
        c.parent = u                                                                    # set the parent of c to be the new node u
        c.range_start = u.range_start + x                                                           # update the range of c to start where u ends
        c_char = self.seq[c.range_start]                                                # find the new beginning character of c
        u.children[c_char] = c                                                          # set c to be the child of u at the right place in the list using the first character of c
        self.n += 1
        return u

    def naive_insert(self, c, x, i, n):
        child = self.find_edge(c, x, i)
        if child:  # If we have an outgoing edge
            while x + i < child.range_end and self.seq[child.range_start + x] == self.seq[i+x]:  # If we have have not reached the end of the edge yet and we can extend the match
                x += 1
            if x + i == child.range_end:  # if we reach the end of the edge, we have to find the new outgoing edge to search from
                self.naive_insert(child, x, i, n) # Be carefull with this recursion
            else:
                u = self.split_edge(child, x, i, n)
                self.insert_child(u, x, i, n)


        else:  # If we have not got an outgoing edge
            self.insert_child(c, x, i, n)


    def insert(self, y):
        y = y + "0"                                                                     # Can be done smarter
        self.seq,self.alpha = remap(y)
        n = len(self.seq)

        x = 0
        i = 0

        # insert first sequence manually
        char = self.seq[x + i]
        new = Node(range_start = i + x, range_end = n, string_label = i, parent = self.root, children = [None, None, None, None, None])
        self.root.children[char] = new
        i += 1


        while i < len(y):                                                                    # Change this to iterate through the full sequence when done
            self.naive_insert(self.root,x,i,n)                                                  # input = suffix tree, root to search from, x (incrementor), i = string label. n+1 = length of sequence
            i += 1


    def find_edge_search(self, c, char):
        return c.children[char]

    def search_rec(self,c,x,n):
        child = self.find_edge_search(c, n[x])
        if child:  # If we have an outgoing edge
            while x < child.range_end and x < len(n) and self.seq[child.range_start + x] == n[x]:  # self.root.range_start has to be the beginning of the match. In other words: We have to begin at each node instead of the root
                x += 1
            if x == len(n):  # if we reach the end of an edge and we find a match (also works on edges
                res = return_labels(child)
                return res  # Be carefull with this recursion
            elif x + child.range_start == child.range_end and x < len(n):  # if we reach the end of the edge, and we have not found a match, we have to search down a new node
                return self.search_rec(child,x = 0, n = n[x:])
            else:
                return "Match not found in branch"
        else:
            return "Match not found in node"

    def search(self,y):
        n = [self.alpha[a] for a in y]
        c = self.root
        x = 0
        res = self.search_rec(c,x,n)
        return res




x = "GTCGTA"


res = SuffixTree()
res.insert(x)
# print(res.n)
n = res.search("GTA")
for i in n:
    print(i)
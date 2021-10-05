# Where to put this function?
def remap(x):
    m = {a:i for i,a in enumerate(sorted(set(x)))}
    n = [m[a] for a in x]
    return n


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
        self.root = Node()
 #       self.root.parent = self.root

    def find_edge(self,c,x,i):
        char = self.seq[x+i]
        return c.children[char]

    def insert_child(self,c,x,i,n):
        char = self.seq[x + i]
        new = Node(range_start=i + x, range_end=n, string_label=i, parent = c, children=[None, None, None, None, None])
        c.children[char] = new


    def split_edge(self,c,x,i,n):
        c_char = self.seq[c.range_start]                                                # first character in c
        u = Node(range_start = c.range_start, range_end = i + x, parent = c.parent)    # create u node
        c.parent.children[c_char] = u                                                   # update the child of P(c) to be u. The first character in u is c_char
        c.parent = u                                                                    # set the parent of c to be the new node u
        c.range_start = i + x                                                           # update the range of c to start where u ends
        c_char = self.seq[c.range_start]                                                # find the new beginning character of c
        u.children[c_char] = c                                                          # set c to be the child of u at the right place in the list using the first character of c

    def naive_insert(self, c, x, i, n):
        child = self.find_edge(c, x, i)
        if child:  # If we have an outgoing edge
            while x + i < child.range_end and self.seq[child.range_start + x] == self.seq[i+k]:  # If we have have not reached the end of the edge yet and we can extend the match
                x += 1
            if x + i == child.range_end:  # if we reach the end of the edge, we have to find the new outgoing edge to search from
                self.naive_insert(child, x, i, n)
            else:
                u = self.split_edge(child, x, i, n)
                self.insert_child(u, x, i, n)


        else:  # If we have not got an outgoing edge
            self.insert_child(c, x, i, n)


    def insert(self, y):
        y = y + "0"                                                                     # Can be done smarter
        self.seq = remap(y)
        n = len(self.seq)

        x = 0
        i = 0

        # insert first sequence manually
        char = self.seq[x + i]
        new = Node(range_start = i + x, range_end = n, string_label = i, parent = self.root, children = [None, None, None, None, None])
        self.root.children[char] = new
        i += 1

        c = self.root
        while i < 3:                                                                    # Change this to iterate through the full sequence when done
            self.naive_insert(c,x,i,n)                                                  # input = suffix tree, root to search from, x (incrementor), i = string label. n+1 = length of sequence
            i += 1


x = "GTCGCA"


res = SuffixTree()
res.insert(x)
print(res)
# Where to put this function?
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

    def insert(self, x):
        self.seq = remap(x)
        self.seq = self.seq + "0"
        n = len(self.seq)

        x = 0
        i = 0
        c = self.root

        while i < 2:                # Change this to iterate through the full sequence when done
            naive_insert(self,c,x,i,n) # input = suffix tree, root to search from, x (incrementor), i = string label. n+1 = length of sequence
            i + = 1

    def naive_insert(self,c,x,i,n):

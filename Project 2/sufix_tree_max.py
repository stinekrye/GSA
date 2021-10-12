def remap(x): # Converts the alphabet to numbers, so for instance ACGT would be [1234]
    m = {a:i for i,a in enumerate(sorted(set(x)))}
    n = [m[a] for a in x]
    return n


class Node():

    def __init__(self, range_start = None, range_end = None, string_label = None, parent = None, children = [None, None, None, None, None]):
        self.range_start = range_start
        self.range_end = range_end
        self.string_label = string_label
        self.parent = parent
        self.children = children # Is defined later


    def __repr__(self):
        rep = f"Node({self.range_start},{self.range_end})"
        return rep


class SuffixTree():

    def __init__(self, seq = None):
        self.seq = seq
        self.root = Node()
        # self.root.parent = self.root


    def insert_suffix(self, string):
        self.seq = remap(string) + "0" #add sentinel

        start = 0
        end = len(self.seq)
        i = 0

        #First suffix
        first = Node(range_start = start, range_end = end, string_label = i, parent = self.root, children = [None, None, None, None, None])
        self.root.children = first

        #Insert the rest
        for i in range(1, end):
            self.naive_insert()


    def naive_insert(self, st, v, start, end): #suffix tree, node to search from (root), start, end
        w = self.find_outgoing_edge(v, start)
        if w:

        else: #If there is no outgoing edge that matches, we insert child here
            self.insert_child(v, start, end)


    def find_outgoing_edge(self, v, start):
        w = v.children
        while w:
            if w.range_start == start:
                break
            else:
                w = w. #????????


    def insert_child(self, v, start, end):
        child = Node(range_start = start, range_end = end, string_label = 1, parent = v) # LAbels??
        return child


    def split_edge(self, st, v, w):
        u = Node(range_start = v.range_start, range_end = v.range_end, string_label = v.label + 1, parent = v.parent) #New node

        # Update w and v
        u.parent = v
        u.children = w
        w.range_start =
        w. parent = u

        return u

# Things to figure out:
# - Labels
# - Indexes, need for almost every function
# - Find outgoing edge
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

        start = 0 #x for stine
        end = len(self.seq) #n for stine
        i = 0

        #First suffix
        char = self.seq[start + i] #First character in the sequence
        first = Node(range_start = start + i, range_end = end, string_label = i, parent = self.root, children = [None, None, None, None, None])
        self.root.children[char] = first #Assigns the first Node as a Child to the root for the corresponding character

        #Insert the rest
        for i in range(1, end):
            self.naive_insert(self.root, v, start, end)


    def naive_insert(self, st, v, start, end, i): #suffix tree, node to search from (root), start, end
        w = self.find_outgoing_edge(st, v, start, i)
        if w: #If there is an outgoing edge
            while start + i < w.range_end and self.seq[w.range_start + start] == self.seq[start + i] #While scanning through the out edge
                start += 1
                if start + i == w.range_end: #If we reach the end of the edge, we must continue from there
                    self.naive_insert(w, start, end, i )


        else: #If there is no outgoing edge that matches, we insert child here
            self.insert_child(v, start, end, i)


    def find_outgoing_edge(self, v, start, i):
        char = self.seq[start + i]
        return v.children[char]


    def insert_child(self, v, start, end, i):
        char = self.seq[start + i]
        child = Node(range_start = start + i, range_end = end, string_label = i, parent = v, children = [None, None, None, None, None])
        v.children[char] = child
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
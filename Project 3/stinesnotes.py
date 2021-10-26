# Try to build the SA and LCP array from a suffix tree
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


def SA(node): #Rename/build to deal with LCP also
    for child in node.children:
        if child != None:
            if child.string_label != None:
                print(child.string_label)
            else:
                SA(child)





# Construct a tree for testing

x = "GTGTC"
tree = SuffixTree()
tree.insert(x)
nodes = SA(tree.root)
print(nodes)

test = [5,4,2,0,3,1]
for i in test:
    print(x[i:])

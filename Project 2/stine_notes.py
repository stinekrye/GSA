# Where to put this function?
def remap(x):
    m = {a:i+1 for i,a in enumerate(sorted(set(x)))}
    n = [m[a] for a in x]
    return n, len(m)


class Node():
    def __init__(self, range_start = None, range_end = None, string_label = None, children = None, parent = None):
        self.range_start = range_start
        self.range_end = range_end
        self.string_label = string_label
        self.children = children # Is defined later
        self.parent = parent

    def __repr__(self):
        rep = f"Node({self.range_start},{self.range_end})"
        return rep


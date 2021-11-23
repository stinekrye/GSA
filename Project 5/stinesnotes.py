import numpy as np
import pandas as pd

## Helper functions

# Remaps sequence
def remap(x):
    m = {a:i for i,a in enumerate(sorted(set(x)))}
    n = [m[a] for a in x]
    return n, m

# Produce accumulated counts
def counts(seq):
    values, counts = np.unique(seq, return_counts = True)
    acc_sum = np.cumsum(counts)
    acc_sum = np.concatenate((np.array([0]),acc_sum[:-1]))
    return acc_sum



# SAIS algorithm

### Step 1: Construct t table
def t_table(seq):
    t = np.empty([len(seq)], dtype = str)
    t[len(seq)-1] = "S"
    for i in reversed(range(len(seq)-1)):
        if x[i] < x[i+1]:
            t[i] = "S"
        elif x[i] > x[i+1]:
            t[i] = "L"
        elif x[i] == x[i+1]:
            t[i] = t[i+1]
    return t

### Step 2: Construct LMS array
def LMS_array(t_table):
    LMS = np.zeros([len(t_table)], dtype=int)
    for i in range(len(t_table)):
        if t_table[i] == "S" and t_table[i-1] == "L":
            LMS[i] = i
    return LMS

### Step 3: Create u
def make_u(LMS):
    u = np.nonzero(LMS)[0]
    return u


### Step 3: sort the LMS strings

def induced_sorting(seq, LMS):
    u = make_u(LMS)
    nlms = np.count_nonzero(LMS)
    bins = np.zeros([len(LMS)-1+nlms], dtype=int)
    u_plus = np.concatenate((u,acc_sum[:-1]))
    test = seq[u][:-1]
    print(test)
    acc_sum = counts(seq) # This does not take the LMS strings into account



    # Place references in bins
    return len(bins)



def construct_sa_lcp(seq):
    x = remap(seq)[0]
    x = np.array(list(x), dtype=str)
    tTable = t_table(x)
    LMS = LMS_array(tTable)
    u = make_u(LMS) # The individual sequences can be found by iterating through u and printing x[u[i:i+2]]

    ### No prints until here ####
    test = induced_sorting(x,LMS)

    return "test"

### Goal: make bin size with the LMS strings as well







#####
x = "mississippi0"
res = construct_sa_lcp(x)
print(res)


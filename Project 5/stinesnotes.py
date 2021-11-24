import numpy as np
import pandas as pd

## Helper functions

# Remaps sequence
def remap(x):
    m = {a:i for i,a in enumerate(sorted(set(x)))}
    n = [m[a] for a in x]
    return n, m

# Produce accumulated counts for the bins. Could be better!
def counts(seq):
    values, counts = np.unique(seq, return_counts = True)
    acc_sum = np.cumsum(counts)
    acc_sum = np.concatenate((np.array([0]),(acc_sum)))
    return acc_sum[:-1], len(values)



# SAIS algorithm

### Step 1: Construct t table
def t_table(seq):
    t = np.empty([len(seq)], dtype = str)
    t[len(seq)-1] = "S"
    for i in reversed(range(len(seq)-1)):
        if seq[i] < seq[i+1]:
            t[i] = "S"
        elif seq[i] > seq[i+1]:
            t[i] = "L"
        elif seq[i] == seq[i+1]:
            t[i] = t[i+1]
    return t

### Step 2: Construct LMS array
def LMS_array(t_table):
    LMS = np.zeros([len(t_table)], dtype=int)
    for i in range(len(t_table)):
        if t_table[i] == "S" and t_table[i-1] == "L":
            LMS[i] = i
    return LMS



### Step 3: Helper functions
#### Create alpha: Helper function
def lms_eq(lms1,lms2,seq, LMS): # Will break of it gets two identical LMS sequences which both are the last LMS index
    non_zero_LMS = np.nonzero(LMS)[0]

    ## Preprocess so we can find length of each LMS string
    if lms1 != non_zero_LMS[-1]:
        lms1_end = non_zero_LMS[np.where(non_zero_LMS == lms1)[0]+1]
    else:
        lms1_end = non_zero_LMS[-1]
    if lms2 != non_zero_LMS[-1]:
        lms2_end = non_zero_LMS[np.where(non_zero_LMS == lms2)[0]+1]
    else:
        lms2_end = non_zero_LMS[-1]

    ## Compare the length of the lms strings
    if lms1_end-lms1 != lms2_end-lms2:
        return False

    ## Compare the sequences of the LMS strings
    else:
        i = 0
        while i < len(lms1_end-lms1):
            if seq[lms1 + i] != seq[lms2 + i]:
                return False
            i += 1
    return True


def create_alpha(bins,LMS,seq):
    alpha = np.zeros([len(bins)], dtype = int)
    a = 1
    first = True
    for i in bins:
        if first:
            if i in np.nonzero(LMS)[0]:
                first = False
                alpha[i] = a
                prev = i


        elif i in np.nonzero(LMS)[0]:
            if not lms_eq(i,prev,seq, LMS):# function to compare two strings
                a += 1
            alpha[i] = a
            prev = i

    return alpha



### Step 3: sort the LMS strings

def induced_sorting(seq, LMS, t_table):
    nlms = np.count_nonzero(LMS)
    bins = np.zeros([len(LMS)], dtype=int)
    bin_start, counter_len = counts(seq) # The sequence is remapped so the index = the remapped letter

    ### Insert the LMS strings
    counter = np.zeros([counter_len], dtype = int)
    for i in np.nonzero(LMS)[0]: # Insert LMS strings in the bottom of their bin
        location = bin_start[seq[int(i)]+1]-1-counter[seq[int(i)]] # gives the location in the bin # We want the LMS strings to be found in the bottom of their bin
        bins[location] = int(i)
        counter[seq[int(i)]] += 1

    ### Insert the runs of L
    counter = np.zeros([counter_len], dtype = int)
    for i in bins:
        if i != 0 and t_table[i-1] == "L":
            location = bin_start[seq[int(i-1)]]+counter[seq[int(i-1)]]
            bins[location] = int(i-1)
            counter[seq[int(i-1)]] += 1

    ### Insert the runs of S
    counter = np.zeros([counter_len], dtype = int)
    for i in reversed(bins):
        if i != 0 and t_table[i-1] == "S":
            location = bin_start[seq[int(i-1)]+1]-1-counter[seq[int(i-1)]]
            bins[location] = int(i-1)
            counter[seq[int(i-1)]] += 1



    ### Create alpha
    a = create_alpha(bins,LMS,seq)
    print(a)
    print(np.nonzero(a))
    return "test"







def construct_sa_lcp(seq):
    x = remap(seq)[0]
    x = np.array(list(x), dtype=int)
    tTable = t_table(x)
    LMS = LMS_array(tTable)

    ### No prints until here ####
    test = induced_sorting(x,LMS, tTable)

    return test

### Goal: make bin size with the LMS strings as well







#####
x = "mississippi0"
res = construct_sa_lcp(x)
print(res)


import numpy as np
## Helper functions
# Optimize lms_eq and create alpha

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
    t = np.zeros([len(seq)], dtype = str)
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
    return np.nonzero(LMS)[0]

### Step 3: Helper functions
#### Create alpha: Helper function
def lms_eq(lms1,lms2,seq, LMS): # Will break of it gets two identical LMS sequences which both are the last LMS index
    ## Preprocess so we can find length of each LMS string
    if lms1 != LMS[-1]:
        lms1_end = LMS[np.where(LMS == lms1)[0]+1]
    else:
        lms1_end = LMS[-1]
    if lms2 != LMS[-1]:
        lms2_end = LMS[np.where(LMS == lms2)[0]+1]
    else:
        lms2_end = LMS[-1]

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
            if i in LMS:
                first = False
                alpha[i] = a
                prev = i


        elif i in LMS:
            if not lms_eq(i,prev,seq, LMS):# function to compare two strings
                a += 1
            alpha[i] = a
            prev = i

    return alpha[np.nonzero(alpha)]-1

### Step 3: sort the LMS strings

def induced_sorting(seq, LMS, ttable, map):
    bins = np.zeros([len(seq)], dtype=int)
    bin_start, counter_len = counts(seq) # The sequence is remapped so the index = the remapped letter

    ### Insert the LMS strings
    counter = np.zeros([counter_len], dtype = int)
    for i in reversed(LMS): # Insert LMS strings in the bottom of their bin. Has to be reverse because it starts at the bottom
        location = bin_start[seq[int(i)]+1]-1-counter[seq[int(i)]] # gives the location in the bin # We want the LMS strings to be found in the bottom of their bin
        bins[location] = int(i)
        counter[seq[int(i)]] += 1

    ### Insert the runs of L
    counter = np.zeros([counter_len], dtype = int)
    for i in bins:
        if i != 0 and ttable[i-1] == "L":
            location = bin_start[seq[int(i-1)]]+counter[seq[int(i-1)]]
            bins[location] = int(i-1)
            counter[seq[int(i-1)]] += 1

    ### Insert the runs of S
    counter = np.zeros([counter_len], dtype = int)
    for i in reversed(bins):
        if i != 0 and ttable[i-1] == "S":
            location = bin_start[seq[int(i-1)]+1]-1-counter[seq[int(i-1)]]
            bins[location] = int(i-1)
            counter[seq[int(i-1)]] += 1


    if map == True:
        ### Create alpha
        a = create_alpha(bins,LMS,seq)
        a.size
        if len(np.unique(a)) == len(a):
            return bins
        else:
            a_tTable = t_table(a)
            a_LMS = LMS_array(a_tTable)
            # if len(np.unique(a_LMS)) == len(a_LMS):
            #     return bins
            # else:
            return induced_sorting(a,a_LMS, a_tTable, map = True) # return or not

    return bins

def LMS_order(map_LMS, LMS):
    nLMS = LMS
    sa_LMS = np.zeros([len(nLMS)], dtype= int)
    for i in range(len(map_LMS)):
        sa_LMS[i] = nLMS[map_LMS[i]]
    return sa_LMS

def construct_sa(seq):
    x = remap(seq)[0]
    x = np.array(list(x), dtype=int)
    tTable = t_table(x)
    LMS = LMS_array(tTable)

    ### No prints until here ####
    map_LMS = induced_sorting(x, LMS, tTable, map=True)
    # if len(map_LMS) == len(x):
    #     return map_LMS
    #
    # else:
    sa_LMS = LMS_order(map_LMS, LMS)
    sorted = induced_sorting(x,sa_LMS, tTable, map = False)
    return sorted

def gen_sa(fastadict, fastaname, reverse = False):
    if reverse == True:
        f = open(f"{fastaname}.rev.sa", "w")
    else:
        f = open(f"{fastaname}.sa", "w")

    # # Generate tree, SA and LCP arrays
    for key, value in fastadict.items():
        if reverse == True:
            value = value[::-1]
        seq = value+"0" # Find a better way to do this
        SA = construct_sa(seq)

    # Write file
        f.write(">" + str(key) + "\t" + str(value) + "\n")
        for i in range(len(SA)):
            f.write(f"{SA[i]} \n")

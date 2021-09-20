from test_scripts.search_ba_stine import ba_search
import random
random.seed(2)

#%%

def dna(length=int(), letters="CGTA"):
    return ["seq",''.join(random.choices(letters, k=length))]

def fastq(dna): # Third row with + is omitted, because it is excluded by fastq reader
    res = []
    first = '@Seq'
    second = "{}".format(dna)
    fourth = "{}".format('~' * len(dna))
    res.append(first)
    res.append([second,fourth])
    res.append(fourth)
    return res


#%%
# Test 1
x = dna(20)
seq1 = x[1][:5]
p = fastq(seq1)
real_res1 = "@Seq 0 seq 1 0 5M * 0 0 {} {}".format(seq1, "~" * len(seq1))


test_res = ba_search(x, p)
print(test_res)




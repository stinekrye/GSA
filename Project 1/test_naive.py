from search_naive_stine import naive_matching
import random

#%%

def dna(length=int(), letters="CGTA"):
    return ["seq",''.join(random.choices(letters, k=length))]

def fastq(dna): # Third row with + is omitted, because it is excluded by fastq reader
    res = []
    first = '@Seq'
    second = "{}".format(dna[1])
    fourth = "{}".format('~' * len(dna[1]))
    res.append(first)
    res.append([second,fourth])
    res.append(fourth)
    return res


#%%

fasta = dna(20)
idx = random.randrange(5,len(fasta[1]),5)
seq1 = fasta[:5]
seq2 = fasta[-5:]
seq3 = fasta[idx:idx+6]

p = fastq(seq1)
# print(p[1][1])

res = naive_matching(fasta,p)
print(res)
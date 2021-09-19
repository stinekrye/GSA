from search_naive_stine import naive_matching
import unittest
import random
random.seed(2)

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
# Test 1
fasta = dna(20)
idx = random.randrange(5,len(fasta[1]),5)
seq1 = fasta[:5]
p = fastq(seq1)
real_res1 = "@Seq 0 seq 1 0 20M * 0 0 AACCATTGTTTCGGTAATGG ~~~~~~~~~~~~~~~~~~~~"

seq2 = fasta[-5:]
seq3 = fasta[idx:idx+6]




# print(p[1][1])

test = []
res = naive_matching(fasta,p)
for value in res:
    test.append(value)
test = "\n".join(test)
print(test)


class TestNaive(unittest.TestCase):
    def TestBeginning(self,real_res1):
        iterator = naive_matching(fasta,p)
        test_res = []
        for value in iterator:
            test_res.append(value)
        test_res = "\n".join(test_res)
        self.assertEqual(test_res, real_res1)

if __name__ == '__main__':
    unittest.main()

from unittest import TestCase
from search_naive_stine import naive_search
import random
random.seed(367)
# Fails when has to deal with two or more matches. Problems with the way I have defined the result. This will be easier, when I can use this alg as "expected result"

#%%  Helper functions
def dna(length=int(), letters="CGTA"):
    return ["seq", ''.join(random.choices(letters, k=length))]


def fastq(dna):  # Third row with + is omitted, because it is excluded by fastq reader
    res = []
    first = '@Seq'
    second = "{}".format(dna)
    fourth = "{}".format('~' * len(dna))
    res.append(first)
    res.append([second, fourth])
    res.append(fourth)
    return res

#%% Test 1 data BEGINNING                              # make this into a table to test 5 times each case
x1 = dna(20)
seq1 = x1[1][:5]
p1 = fastq(seq1)
real_res1 = "@Seq 0 seq {} 0 {}M * 0 0 {} {}".format(1, len(seq1), seq1, "~" * len(seq1))

#%% Test 2 data END
x2 = dna(20)
seq2 = x2[1][-5:]
p2 = fastq(seq2)
real_res2 = "@Seq 0 seq {} 0 {}M * 0 0 {} {}".format(len(x2[1])-5+1, len(seq2), seq2, "~" * len(seq2))

#%% Test 3 data RANDOM
x3 = dna(20)
idx = random.randrange(5,len(x3[1]),5)
seq3 = x3[1][idx:idx+6]
p3 = fastq(seq3)
real_res3 = "@Seq 0 seq {} 0 {}M * 0 0 {} {}".format(idx+1, len(seq3), seq3, "~" * len(seq3))


#%% Test 4 data fasta is empty
x4 = ""
p4 = fastq(seq3)

#%% Test 5 data fasta is empty
x5 = dna(20)
p5 = ""

#%% Test 6 multiple hits
x6 = ["seq", "AGCTAGCTCCCCCCAGCT"]
seq6 = "AGCT"
p6 = fastq(seq6)
hit = [0,4,14]
real_res6 = f"@Seq 0 seq {hit[0]+1} 0 {len(seq6)}M * 0 0 {seq6} {'~' * len(seq6)}\n" \
            f"@Seq 0 seq {hit[1]+1} 0 {len(seq6)}M * 0 0 {seq6} {'~' * len(seq6)}\n" \
            f"@Seq 0 seq {hit[2]+1} 0 {len(seq6)}M * 0 0 {seq6} {'~' * len(seq6)}"

#%% Test 7 multiple overlapping hits
x7 = ["seq", "AGAGAGACCCCCCAGAG"]
seq7 = "AGAG"
p7 = fastq(seq7)
hit = [0,2,13]
real_res7 = f"@Seq 0 seq {hit[0]+1} 0 {len(seq7)}M * 0 0 {seq7} {'~' * len(seq7)}\n" \
            f"@Seq 0 seq {hit[1]+1} 0 {len(seq7)}M * 0 0 {seq7} {'~' * len(seq7)}\n" \
            f"@Seq 0 seq {hit[2]+1} 0 {len(seq7)}M * 0 0 {seq7} {'~' * len(seq7)}"


class TestRandom(TestCase):
    def testBeginning(self):
        test_res1 = naive_search(x1, p1)
        self.assertEqual(test_res1, real_res1)
    def testEnd(self):
        test_res2 = naive_search(x2, p2)
        self.assertEqual(test_res2, real_res2)
    def testRandom(self):
        test_res3 = naive_search(x3, p3)
        self.assertEqual(test_res3, real_res3)
    def testEmptyFasta(self):
        with self.assertRaises(IndexError):
            naive_search(x4, p4)
    def testEmptyFastq(self):
        with self.assertRaises(IndexError):
            naive_search(x5, p5)
    def testMultiple(self):
        test_res6 = naive_search(x6,p6)
        self.assertEqual(test_res6,real_res6)
    def testOverlap(self):
        test_res7 = naive_search(x7, p7)
        self.assertEqual(test_res7, real_res7)
from unittest import TestCase
from search_naive import search_naive
from stinesnotes import search_bs
from stinesnotes import print_test
import parsers, random, sys

random.seed(367)



#%% Test 1 data compare output of testfile                           # make this into a table to test 5 times each case
fasta = parsers.read_fasta_file(f"fasta_test.fa")
fastq = parsers.read_fastq_file(f"fastq_test.fq")
sa = parsers.read_SA_LCP(f"fasta_test.fa.sa-lcp")




class TestRandom(TestCase):
    def testMissisipi(self):
        search_naive(fasta,fastq)
        right_answer = sys.stdout.getvalue().strip()
        iter = search_bs(sa,fastq, test = True)
        print_test(iter)
        t = sys.stdout.getvalue().strip()
        self.assertEqual(t, right_answer)
    # def testOverlapping(self):
    #     search_naive(fasta1,fastq1)
    #     t = sys.stdout.getvalue().strip()
    #     search_bs(fasta1,fastq1)
    #     right_answer = sys.stdout.getvalue().strip()
    #     self.assertEqual(t, right_answer)

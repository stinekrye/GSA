from unittest import TestCase
from search_naive import search_naive
from stinesnotestest import search_bs
from stinesnotestest import print_test
import parsers, random, sys

random.seed(367)



#%% Test 1 data compare output of testfile                           # make this into a table to test 5 times each case
fasta = parsers.read_fasta_file(f"fasta_test.fa")
fastq = parsers.read_fastq_file(f"fastq_test.fq")
sa = parsers.read_SA(f"fasta_test.fa.sa")

fasta1 = parsers.read_fasta_file(f"fasta_test1.fasta")
fastq1 = parsers.read_fastq_file(f"fastq_test1.fastq")
sa1 = parsers.read_SA(f"fasta_test1.fasta.sa")



class TestRandom(TestCase):
    def testMissisipi(self):
        search_naive(fasta,fastq)
        right_answer = sys.stdout.getvalue().strip()
        iter = search_bs(sa,fastq, test = True)
        print_test(iter)
        t = sys.stdout.getvalue().strip()
        self.assertEqual(t, right_answer)
    def testOverlapping(self):
        search_naive(fasta1,fastq1)
        right_answer = sys.stdout.getvalue().strip()
        iter = search_bs(sa1,fastq1, test = True)
        print_test(iter)
        t = sys.stdout.getvalue().strip()
        self.assertEqual(t, right_answer)

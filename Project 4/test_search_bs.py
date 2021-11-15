from unittest import TestCase
from search_bs import search_bs
from ../Project_1/
import parsers
import random, sys
random.seed(367)

#%% Test 1 data compare output of testfile                           # make this into a table to test 5 times each case
fasta = parsers.read_fasta_file(f"fasta_test.fa")
fastq = parsers.read_fastq_file(f"fastq_test.fq")



class TestRandom(TestCase):
    def testMissisipi(self):
        naive_matching(fasta,fastq)
        test = sys.stdout.getvalue().strip()
        search_ba_wrapper(fasta,fastq)
        right_answer = sys.stdout.getvalue().strip()
        self.assertEqual(test, right_answer)
    def testOverlapping(self):
        naive_matching(fasta1,fastq1)
        test = sys.stdout.getvalue().strip()
        search_ba_wrapper(fasta1,fastq1)
        right_answer = sys.stdout.getvalue().strip()
        self.assertEqual(test, right_answer)

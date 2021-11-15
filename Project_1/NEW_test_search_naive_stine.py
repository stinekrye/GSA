from unittest import TestCase
from search_ba import search_ba_wrapper
from read_fasta import read_fasta_file
from read_fastq import read_fastq_file
from NEW_search_naive_stine import search_naive
import random, sys
random.seed(367)

#%% Test 1 data compare output of testfile                           # make this into a table to test 5 times each case
fasta = read_fasta_file(f"resubmission/fasta_test.fasta")
fastq = read_fastq_file(f"resubmission/fastq_test.fastq")
fasta1 = read_fasta_file(f"fasta_test1.fasta")
fastq1 = read_fastq_file(f"fastq_test1.fastq")



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

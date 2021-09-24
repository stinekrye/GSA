from unittest import TestCase
from NEW_search_naive_stine import naive_matching
from search_ba import search_ba_wrapper
from read_fasta import read_fasta_file
from read_fastq import read_fastq_file
import random, sys
random.seed(367)

#%% Test 1 data compare output of testfile                           # make this into a table to test 5 times each case
fasta = read_fasta_file(f"../fasta_test.fasta")
fastq = read_fastq_file(f"../fastq_test.fastq")



class TestRandom(TestCase):
    def testBeginning(self):
        naive_matching(fasta,fastq)
        test = sys.stdout.getvalue().strip()
        search_ba_wrapper(fasta,fastq)
        right_answer = sys.stdout.getvalue().strip()
        self.assertEqual(test, right_answer)

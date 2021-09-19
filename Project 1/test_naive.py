from search_naive import naive_matching
from read_fasta import read_fasta_file
from read_fastq import read_fastq_file

# x er fasta file, p er fastq file

x = read_fasta_file("stinetest.fasta")
p = read_fastq_file("fastqtest.fastq")

res = naive_matching()
print(p["Seq1"][0])
print(x)
# naive_matching(x,p)
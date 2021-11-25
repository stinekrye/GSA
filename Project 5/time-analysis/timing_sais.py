#%% md

## Timing our algorithms

#%%

import time, random
import pandas as pd
import seaborn as sns
import matplotlib as plt
import sys

#%%

import parsers, stinesnotes

#%%

def dna(length=int(), letters="CGTA"):
    return''.join(random.choices(letters, k=length))

def create_one_fasta(n):
    name = f"test_files/n_{n}.fasta"
    f = open(name, "w")

    f.write('>Seq' + str(1) + '\n' + dna(n) + '\n')

    f.close()

    return name

def create_many_fasta(start, stop, step):
    name = f"test_files/n_{stop}.fasta"
    f = open(name, "w")
    counter = 1

    for i in range(start, stop, step):
        f.write('>Seq' + str(counter) + '\n' + dna(i) + '\n')
        counter += 1
    f.close()

    return name

def create_one_fastq(m):
    name = f"test_files/m_{m}.fastq"
    f = open(name, "w")
    f.write(
        '@' + 'Seq' + str(1) + '\n' +
        dna(m) + '\n' +
        '+' + '\n' +
        '~' * m + '\n')
    f.close()

    return name

def create_many_fastq(start, stop, step):
    name = f"test_files/m_{stop}.fastq"
    f = open(name, "w")

    counter = 1

    for i in range(start, stop, step):
        f.write(
            '@' + 'Seq' + str(counter) + '\n' +
            dna(i) + '\n' +
            '+' + '\n' +
            '~' * i + '\n'
        )
        counter += 1

    f.close()

    return name



#%%

def time_construction(n):
    df = pd.DataFrame(range(20, n, int(n/10)), columns=['n'])
    fasta_file = create_many_fasta(20, n, int(n/10))
    fasta = parsers.read_fasta_file(fasta_file)
    fastaname = "fasta_test.fa"
    times = []

    if len(fasta) < 0:
        return "Problems with either fasta or fastq file"

    stinesnotes.gen_sa(fasta,fastaname)


    df['Time'] = times

    df['Time/Expected time'] = df['Time']/df['n']
    return df

#%%

n = 100000
m = 1000

#%%

df_stc = time_construction(n)

#%%

h = sns.lineplot(x = 'n', y = 'Time/Expected time', marker = '.',
                    data = df_stc)
h.set_title('Suffix tree construction')
h.figure.savefig('../figures/st_construction.pdf')

# %% md

## Timing our algorithms

# %%

import time, random
import pandas as pd
import seaborn as sns
import matplotlib as plt
import sys

# %%

import parsers, stinesnotes


# %%

def dna(length=int(), letters="CGTA"):
    return ''.join(random.choices(letters, k=length))


def create_one_fasta(n):
    name = f"test_files/n_{n}.fa"
    f = open(name, "w")

    f.write('>Seq' + str(1) + '\n' + dna(n) + '\n')

    f.close()

    return name


def create_many_fasta(start, stop, step):
    name = f"test_files/n_{stop}.fa"
    f = open(name, "w")
    counter = 1

    for i in range(start, stop, step):
        f.write('>Seq' + str(counter) + '\n' + dna(i) + '\n')
        counter += 1
    f.close()

    return name


def create_one_fastq(m):
    name = f"test_files/m_{m}.fq"
    f = open(name, "w")
    f.write(
        '@' + 'Seq' + str(1) + '\n' +
        dna(m) + '\n' +
        '+' + '\n' +
        '~' * m + '\n')
    f.close()

    return name


def create_many_fastq(start, stop, step):
    name = f"test_files/m_{stop}.fq"
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


# %%

def time_bssearch(sa, m):
    df = pd.DataFrame(range(20, m, int(m / 10)), columns=['m'])
    fastq_file = create_many_fastq(20, m, int(m / 10))

    fastq = parsers.read_fastq_file(fastq_file)
    sa = parsers.read_SA_LCP(sa)
    times = []

    flag, mapq, pnext, tlen = 0, 0, 0, 0
    rnext = "*"

    for x in sa.items():
        rname = x[0]
        y = x[1][0] + "$"
        sa = x[1][1]

        for p in fastq.items():
            qname = p[0]
            substring = p[1][0]
            cigar = str(len(substring)) + "M"
            qual = p[1][1]

            start = time.time()
            matches = stinesnotes.binary3(substring, y, sa)
            end = time.time()
            difference = end - start
            times.append(difference)

            if matches is not None:
                for match in matches:
                    pos = int(match) + 1
                    print(
                        f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}",
                        file=sys.stdout)

    df['Time'] = times

    df['Time/Expected time'] = df['Time'] / df['m']

    return df


# %%

n = 100000
m = 10000

# %%

df_search = time_bssearch("test_files/n_100000.fa.sa-lcp", m)

# %%

h = sns.lineplot(x='m', y='Time/Expected time', marker='.',
                 data=df_search)
h.set_title('Pattern matching using Binary search')
h.figure.savefig('figures/search_bw.pdf')

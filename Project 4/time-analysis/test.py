import time, random, sys
import pandas as pd
import parsers, jannesnotes

def dna(length=int(), letters="CGTA"):
    return''.join(random.choices(letters, k=length))

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

def time_bwsearch(sa, m, o_dict, c_dict):
    df = pd.DataFrame(range(20, m, int(m/10)), columns=['m'])
    fastq_file = create_many_fastq(20, m, int(m/10))

    fastq = parsers.read_fastq_file(fastq_file)
    sa = parsers.read_SA(sa)
    times = []

    flag, mapq, pnext, tlen = 0,0,0,0
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

            alpha = {a:i for i, a in enumerate(sorted(set(y)))}
            O = o_dict[rname]
            C = c_dict[rname]

            start = time.time()
            matches = jannesnotes.fm_search(O, C, substring, sa, alpha)
            end = time.time()
            difference = end - start
            times.append(difference)

            if matches is not None:
                for match in matches:
                    pos = int(match) + 1
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)
    
    df['Time'] = times

    df['Time/Expected time'] = df['Time']/df['m']

    return df

m = 10000
o_dict = parsers.read_o("G:\Mit drev\Genome-scale algorithms\GSA\Project 4\time-analysis\test_files\n_100000.fa.o-table")
c_dict = parsers.read_c("G:\Mit drev\Genome-scale algorithms\GSA\Project 4\time-analysis\test_files\n_100000.fa.c-table")
# c_ = {"Seq1": {'$': 0, 'A': 1, 'C': 24825, 'G': 49820, 'T': 75029}}
# o = {"Seq1": [[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0, 1.0], [0.0, 1.0, 1.0, 1.0, 1.0], [0.0, 1.0, 1.0, 1.0, 2.0], [0.0, 1.0, 1.0, 2.0, 2.0], [0.0, 1.0, 2.0, 2.0, 2.0], [0.0, 1.0, 2.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0, 3.0], [0.0, 2.0, 2.0, 3.0, 3.0], [0.0, 2.0, 3.0, 3.0, 3.0], [0.0, 2.0, 3.0, 4.0, 3.0], [0.0, 2.0, 3.0, 5.0, 3.0], [0.0, 2.0, 4.0, 5.0, 3.0], [0.0, 2.0, 4.0, 5.0, 4.0], [0.0, 3.0, 4.0, 5.0, 4.0], [0.0, 4.0, 4.0, 5.0, 4.0], [0.0, 5.0, 4.0, 5.0, 4.0], [0.0, 5.0, 5.0, 5.0, 4.0], [0.0, 5.0, 6.0, 5.0, 4.0], [0.0, 5.0, 6.0, 6.0, 4.0], [0.0, 5.0, 6.0, 7.0, 4.0], [0.0, 6.0, 6.0, 7.0, 4.0], [0.0, 7.0, 6.0, 7.0, 4.0], [0.0, 8.0, 6.0, 7.0, 4.0], [0.0, 9.0, 6.0, 7.0, 4.0], [0.0, 10.0, 6.0, 7.0, 4.0], [0.0, 10.0, 6.0, 7.0, 5.0]]}
df_search = time_bwsearch("test_files/n_100000.fa.sa", m, o_dict, c_dict)



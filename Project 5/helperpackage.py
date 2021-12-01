from os import read
import ast
import numpy as np

#### This file contains all our parsers and our functions for creating C, O and d tables.

# Parsers

def read_fasta_file(filename):
    sequences_lines = {}
    current_sequence_lines = None
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            if line.startswith('>'):
                line = line.replace(" ", "")
                sequence_name = line.lstrip('>')
                current_sequence_lines = []
                sequences_lines[sequence_name] = current_sequence_lines
            else:
                if current_sequence_lines is not None:
                    current_sequence_lines.append(line)
    sequences = {}
    for name, lines in sequences_lines.items():
        sequences[name] = ''.join(lines)
    return sequences

def read_fastq_file(file):
    seq = {}
    current_sequence = None
    with open(file) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('+') or not line:
                continue
            if line.startswith('@'):
                sequence_name = line.lstrip('@')
                current_sequence = []
                seq[sequence_name] = current_sequence
            else:
                if current_sequence is not None:
                    current_sequence.append(line)
    seqs = {}
    for name, lines in seq.items():
        seqs[name] = lines

    return seqs

def read_SA(file):
    seq = {}
    sequence, current_sa = None, None
    with open(file) as fp:
        for line in fp:
            line = line.rstrip()
            if line.startswith('SA') or not line:
                continue
            if line.startswith('>'):
                line = line.replace(" ", "")
                line = line.split()
                sequence_name = line[0].lstrip('>')
                sequence = line[1]
                current_sa = []
                seq[sequence_name] = sequence, current_sa
            else:
                line = line.split()
                if current_sa is not None:
                    current_sa.append(int(line[0]))
    seqs = {}
    for name, lines in seq.items():
        seqs[name] = lines

    return seqs

def read_c(file):
    c_table = {}
    with open(file) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('>'):
                sequence_name = line.lstrip('>')
            else:
                line = ast.literal_eval(line)
                c_table[sequence_name] = line
    seqs = {}
    for name, lines in c_table.items():
        seqs[name] = lines

    return seqs

def read_o(file):
    o_table = {}
    with open(file) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('>'):
                sequence_name = line.lstrip('>')
            else:
                line = ast.literal_eval(line)
                o_table[sequence_name] = line
    seqs = {}
    for name, lines in o_table.items():
        seqs[name] = lines

    return seqs

# Tables

def c_table(sa_dict, fastaname):
    f = open(f"{fastaname}.c-table", "w")
    buckets = None

    for key, value in sa_dict.items():
        buckets = {}
        sa = value[1]
        seq = value[0] + "$"
        
        for i in range(0, len(sa)):
            char = seq[sa[i]]
            if char not in buckets:
                buckets[char] = i
        
        f.write(">" + str(key) + "\n")
        f.write(str(buckets) + "\n")
    f.close()

    return

def o_table(sa_dict, fastaname, ro = False):
    if ro == True:
        f = open(f"{fastaname}.ro-table", "w")
    else:
        f = open(f"{fastaname}.o-table", "w")

    for key, value in sa_dict.items():
        sa = value[1]
        seq = value[0] + "$"

        alphabet = sorted(set(seq))
        alphabet_size = len(alphabet)
        n = len(sa)

        O = np.zeros((n + 1, alphabet_size))
        
        for i in range(1, n + 1):
            for a in range(0, alphabet_size):
                if seq[sa[i - 1] - 1] == alphabet[a]:
                    O[i][a] = O[i - 1][a] + 1
                else: 
                    O[i][a] = O[i - 1][a]
        
        f.write(">" + str(key) + "\n")
        f.write(str(O.tolist()) + "\n")

    f.close()

    return

def d_table(RO, C, sa, p, alpha):
    min_edits = 0
    m = len(p)
    L, R = 0, len(sa)
    d = np.zeros(m, dtype=int)

    for i in range(m):
        a = p[i]
        L = C[a] + RO[int(L)][alpha[a]]
        R = C[a] + RO[int(R)][alpha[a]]
        if L > R:
            min_edits += 1
            L = 0
            R = len(sa)
        d[i] = min_edits

    return d
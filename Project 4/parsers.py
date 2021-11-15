
from os import read
import ast


def read_fasta_file(filename):
    sequences_lines = {}
    current_sequence_lines = None
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            if line.startswith('>'):
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

def read_SA_LCP(file):
    seq = {}
    sequence, current_sa, current_lcp = None, None, None
    with open(file) as fp:
        for line in fp:
            line = line.rstrip()
            if line.startswith('SA') or line.startswith('LCP') or not line:
                continue
            if line.startswith('>'):
                line = line.split()
                sequence_name = line[0].lstrip('>')
                sequence = line[1]
                current_sa = []
                current_lcp = []
                seq[sequence_name] = sequence, current_sa, current_lcp
            else:
                line = line.split()
                if current_sa is not None:
                    current_sa.append(int(line[0]))
                if current_lcp is not None:
                    current_lcp.append(int(line[1]))
    seqs = {}
    for name, lines in seq.items():
        seqs[name] = lines

    return seqs

def read_c(file):
    c_table = {}
    current_c, sequence = None, None
    with open(file) as fp:
        for line in fp:
            line = line.rstrip()
            if line.startswith('>'):
                line = line.split()
                sequence_name = line[0].lstrip('>')
                sequence = line[1]
                current_c = []
                c_table[sequence_name] = current_c
                current_c.append(sequence)
            else:
                if current_c is not None:
                    line = ast.literal_eval(line)
                    current_c.append(line)
    seqs = {}
    for name, lines in c_table.items():
        seqs[name] = lines

    return seqs

def read_o(file):
    o_table = {}
    current, current_o, sequence = None, None, None
    with open(file) as fp:
        for line in fp:
            line = line.rstrip()
            if line.startswith('>'):
                line = line.split()
                sequence_name = line[0].lstrip('>')
                sequence = line[1]
                current = []
                o_table[sequence_name] = current
                current.append(sequence)
            else:
                if current is not None:
                    current.append(line)
    seqs = {}
    for name, lines in o_table.items():
        seqs[name] = lines

    return seqs

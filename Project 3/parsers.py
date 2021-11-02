import numpy as np

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

# def read_SA_LCP(fasta):
#     data = np.loadtxt(f"SA_LCP_{fasta}.txt", dtype = int, skiprows = 1)
#     SA = data[:,0]
#     LCP = data[:,1]
#     return SA,LCP

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
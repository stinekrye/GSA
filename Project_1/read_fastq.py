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
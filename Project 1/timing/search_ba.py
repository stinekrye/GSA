# Border array function
def borderarray(x):
    ba = [0] * len(x)
    for i in range(1, len(x)):
        b = ba[i - 1]
        while b > 0 and x[i] != x[b]:
            b = ba[b - 1]
        if x[i] == x[b]:
            ba[i] = b + 1

    return ba

# Search function
def search_ba(x, p):
    # Definitions for the SAM format
    qname = p[0]
    flag = 0
    rname = x[0]
    pos = 0
    mapq = 0
    substring = p[1][0]
    cigar = str(len(substring)) + 'M'
    rnext = '*'
    pnext = 0
    tlen = 0
    qual = p[1][1]

    # Other definitions
    rseq = x[1]
    m = len(substring)
    s = substring + '$' + rseq
    b = borderarray(s)

    for i in range(len(b)):
        if b[i] == m:
            pos = i - 2 * m + 1
            print(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, substring, qual)
            pos = 0
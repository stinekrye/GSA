import sys
import argparse
import gen_lcp, parsers

def binary_search(x, p, sa):
    low, high = 0, len(sa)
    while low < high:
        mid = low + (high - low) // 2
        if x[sa[mid]:(sa[mid]+len(p))] == p:
            lower_bound = lower_bound_search(x, p, sa, low, mid)
            upper_bound = upper_bound_search(x, p, sa, high, mid)
            occ = sa[int(lower_bound):int(upper_bound)]
            return occ
        elif x[sa[mid]:(sa[mid]+len(p))] > p:
            high = mid - 1
        else:
            low = mid + 1
    mid = low + (high - low) // 2
    lower_bound = lower_bound_search(x, p, sa, low, mid)
    upper_bound = upper_bound_search(x, p, sa, high, mid)
    occ = sa[int(lower_bound):int(upper_bound)]
    return occ

def lower_bound_search(x, p, sa, low, mid):
    while low < mid:
        middle = (low + mid) // 2
        if x[sa[middle]:(sa[middle] + len(p))] >= p:
            mid = middle
        elif x[sa[middle]:(sa[middle] + len(p))] < p:
            low = middle + 1
    return low

def upper_bound_search(x, p, sa, high, mid):
    while mid < high:
        middle = (mid + high) // 2
        if x[sa[middle]:(sa[middle] + len(p))] > p:
            high = middle
        elif x[sa[middle]:(sa[middle] + len(p))] <= p:
            mid = middle + 1
    if high == len(p):
        return high
    else:
        return high + 1


def search_bs(sa, fastq):

    if len(sa) < 0 or len(fastq) < 0:
        return "Problems with either fastq file or the SA and LCP"

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

            matches = binary_search(y, substring, sa)

            if matches is not None:
                for match in matches:
                    pos = int(match) + 1
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)
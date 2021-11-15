import gen_lcp, parsers
import sys, argparse

def binary1(p,x,sa,k, l, u): # returns a match in the SA
    low = l
    high = u
    mid = low + ((high-low)//2)

    while high-low != 1:
        if p[k] == x[sa[mid]+k:sa[mid]+k+1]:
            return mid, high, low

        elif p[k] > x[sa[mid]+k:sa[mid]+k+1]:
            low = mid
            mid = low + (high-low)//2

        else:
            high = mid
            mid = low + ((high-low)//2)

def binary2(p,x,sa,k, mid, l, u, upper):
    if upper == False: # If we are searching for lower bound
        low = l
        high = mid
        mid = low + ((high-low) // 2)
        while high - low != 1:
            if p[k] == x[sa[mid] + k:sa[mid] + k + 1]:
                high = mid
                mid = low + ((high-low)//2)
            else:  # p[k] > x[sa[mid] + k:sa[mid] + k + 1] #Try to draw this. We do not need to consider the case where p[k] < x[....]
                low = mid
                mid = low + ((high-low)//2)
            # else:
            #     high = mid
            #     mid = low + ((high-low)//2)


    else:              # If we are searching for upper bound
        low = mid
        high = u
        mid = low + ((high-low) // 2)
        while high - low != 1:
            if p[k] == x[sa[mid] + k:sa[mid] + k + 1]:
                low = mid
                mid = low + ((high-low)//2)
            # elif p[k] > x[sa[mid] + k:sa[mid] + k + 1]:
            #     low = mid
            #     mid = low + ((high-low)//2)
            else: # p[k] < x[.......]
                high = mid
                mid = low + ((high-low)//2)
    return mid


def binary3(p,x,sa):
    l = 0
    u = len(sa)
    k = 0
    while k < len(p):
        first_match, high, low = binary1(p,x,sa,k, l, u)
        if first_match:
            l = binary2(p, x, sa, k, first_match, low, high, upper=False)        # Points one past the interval
            u = binary2(p, x, sa, k, first_match, low, high, upper=True) + 1     # Points one past the interval
            if k == len(p)-1:
                return sa[l+1:u]
            else:
                k += 1
        else:
            return "No match"

    else:
        return "No match"

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

            matches = binary3(substring, y, sa)

            if matches is not None:
                for match in matches:
                    pos = int(match) + 1
                    print(f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}", file = sys.stdout)



###################################
p = "ISS"
x = "MISSISSIPPI0"
sa = [11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]
#
# p = "AT"
# x = "GAGAT"
# sa = [5,1,3,0,2,4]
# l = 0
# u = len(x)

#
# p = "AA"
# x = "AAAAA0"
# sa = [5,4,3,2,1,0]


# res = binary1(p,x,sa,k, l, u)
# print(res)
# res1 = binary2(p,x,sa,k,res,l, u, upper = True)
# res2 = binary3(p,x,sa)

fastq = parsers.read_fastq_file("fastq_test.fq")
sa = parsers.read_SA_LCP("fasta_test.fa.sa-lcp")
res3 = search_bs(sa, fastq)
print(res3)


# The lower bound is one past and the upper is the last one

#############
# Find a way to do the wrapper iteratively
# I can find the upper and lower bound of an interval, but I have to run through all k and return a match if any is found
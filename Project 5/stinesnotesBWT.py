import numpy as np
import SAIS, helperpackage, sys, argparse, os

def compress_cigar(string):
    res2 = np.empty(0, dtype = str)
    count = 1

    for i in range(len(string)-1):
        if string[i] == string[i+1]:
            count += 1
        else:
            res2 = np.append(res2, str(count))
            res2 = np.append(res2, string[i])
            count = 1
    res2 = np.append(res2, str(count))
    res2 = np.append(res2, string[i+1])
    return "".join(res2)

def print_sam(matches, cigar,qname,rname,substring,qual):

    flag, mapq, pnext, tlen = 0, 0, 0, 0
    rnext = "*"
    cigar = ''.join(cigar)
    cigar = cigar[::-1]
    cigar = compress_cigar(cigar)
    for match in matches:
        pos = int(match) + 1
        print(
        f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{substring}\t{qual}",
        file=sys.stdout)
    return


def bw_approx(O, C, p, d, sa, alpha, max_edits, qname, rname, qual):
    L, R = 0, len(sa)
    i = len(p) - 1
    cigar = np.zeros(len(p) + max_edits, dtype=str)

    # Match
    c_index = 0
    for a in list(alpha.keys())[1:]:
        new_L = C[a] + O[int(L)][alpha[a]]
        new_R = C[a] + O[int(R)][alpha[a]]

        # Match or mismatch
        if a == p[i]:
            edit_cost = 0
        else:
            edit_cost = 1

        if max_edits - edit_cost < 0: continue
        if new_L >= new_R: continue

        cigar[c_index] = 'M'
        rec_approx(
            d, sa, O, C, cigar, new_L, new_R,
            i - 1, max_edits - edit_cost, c_index + 1,
            qname, rname, p, qual, alpha)

    # Insertion
    cigar[c_index] = 'I'
    rec_approx(
        d, sa, O, C, cigar, L, R,
        i - 1, max_edits - 1, c_index + 1,
        qname, rname, p, qual, alpha)


def rec_approx(d, sa, O, C, cigar, L, R, i, edits_left, c_index, qname, rname, p, qual, alpha):
    lower_limit = d[i]
    if edits_left < lower_limit:
        return None
    if i < 0:  # Means we have a match
        matches = sa[int(L):int(R)]
        # print(matches, cigar[:c_index][::-1])
        print_sam(matches, cigar[:c_index], qname, rname, p, qual)
        return

    for a in list(alpha.keys())[1:]:
        new_L = C[a] + O[int(L)][alpha[a]]
        new_R = C[a] + O[int(R)][alpha[a]]

        if a == p[i]:
            edit_cost = 0
        else:
            edit_cost = 1
        if edits_left - edit_cost < 0: continue
        if new_L >= new_R: continue

        cigar[c_index] = 'M'
        rec_approx(
            d, sa, O, C, cigar, new_L, new_R,
            i - 1, edits_left - edit_cost, c_index + 1,
            qname, rname, p, qual, alpha)

    # Insertion
    cigar[c_index] = 'I'
    rec_approx(
        d, sa, O, C, cigar, L, R,
        i - 1, edits_left - 1, c_index + 1,
        qname, rname, p, qual, alpha)

    #  Deletion
    cigar[c_index] = 'D'
    for a in list(alpha.keys())[1:]:
        new_L = C[a] + O[int(L)][alpha[a]]
        new_R = C[a] + O[int(R)][alpha[a]]
        if new_L >= new_R: continue

        rec_approx(
            d, sa, O, C, cigar, new_L, new_R,
            i, edits_left - 1, c_index + 1,qname,rname,p,qual, alpha)

    return

def search_bw(sa_dict, rsa_dict, fastq, o_dict, c_dict, ro_dict, max_edits):
    if len(sa_dict) < 0 or len(fastq) < 0:
        return "Problems with either fastq file or the SA"

    for x in sa_dict.items():
        rname = x[0]
        y = x[1][0]
        sa = x[1][1]
        rsa = rsa_dict[rname][1]

        O = o_dict[rname]
        C = c_dict[rname]
        RO = ro_dict[rname]
        alpha = {a: i for i, a in enumerate(sorted(set(y + '0')))}

        for p in fastq.items():
            qname = p[0]
            substring = p[1][0]
            qual = p[1][1]

            d = helperpackage.d_table(RO, C, rsa, substring, alpha)

            bw_approx(O, C, substring, d, sa, alpha, max_edits, qname, rname, qual)
    return

#
#
# # Making a folder for the preprocessing files
# parent_dir = os.getcwd()
# if not os.path.exists('preprocessing'):
#     os.makedirs('preprocessing')
# p_path = os.path.join(parent_dir, 'preprocessing')
#
# # Creating first parser
# parser1 = argparse.ArgumentParser(description='SA and RSA computation using SAIS')
# parser1.add_argument('-p', help='Create SA from fastafile')
# args1 = parser1.parse_known_args()
#
# if args1[0].p:
#     fastafile = helperpackage.read_fasta_file(args1[0].p)
#     fastaname = f"{args1[0].p}"
#
#     # Make suffix array and read file
#     os.chdir(p_path)
#     SAIS.gen_sa(fastafile, fastaname)
#     SAIS.gen_sa(fastafile, fastaname, True)
#     sa = helperpackage.read_SA(f"{fastaname}.sa")
#     rsa = helperpackage.read_SA(f"{fastaname}.rev.sa")
#
#     # Make C, O, and RO table
#     helperpackage.c_table(sa, fastaname)
#     helperpackage.o_table(sa, fastaname)
#     helperpackage.o_table(rsa, fastaname, True)
#
# else:
#     # Creating second parser if -p is not given
#     parser2 = argparse.ArgumentParser(description='Approximate pattern matching using BWT with d-table')
#     parser2.add_argument('fastafile', help="Input fasta file")
#     parser2.add_argument('fastqfile', help="Input fastq file")
#     parser2.add_argument('-d', type=int, help='The maximum number of edits allowed')
#     args2 = parser2.parse_args()
#
#     # Read in files
#     os.chdir(parent_dir)
#     fastq = helperpackage.read_fastq_file(args2.fastqfile)
#     os.chdir(p_path)
#     sa = helperpackage.read_SA(f"{args2.fastafile}.sa")
#     rsa = helperpackage.read_SA(f"{args2.fastafile}.rev.sa")
#     o_dict = helperpackage.read_o(f"{args2.fastafile}.o-table")
#     ro_dict = helperpackage.read_o(f"{args2.fastafile}.ro-table")
#     c_dict = helperpackage.read_c(f"{args2.fastafile}.c-table")
#
#     search_bw(sa, rsa, fastq, o_dict, c_dict, ro_dict, args2.d)
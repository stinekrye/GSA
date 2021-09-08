from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="FASTA parser")
parser.add_argument("filename")
args = parser.parse_args()

with open(args.filename) as file:
    for record in SeqIO.parse(file, "fasta"):
        print('>'+record.name+'\n'+record.seq)
        

        
    #      sequences = {}
    # for record in SeqIO.parse(file, "fasta"):   
    #     sequences[record.id] = record.seq
    #     print(record.id + "")
    # print (sequences)
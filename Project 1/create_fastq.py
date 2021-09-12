import argparse
import random
 
# Initialize parser
parser = argparse.ArgumentParser(description='A script that generates a fastq file')

# Add arguments
parser.add_argument(
    'n',
    type=int,
    help="Length of the sequence(s)"
)
parser.add_argument(
    'number',
    type=int,
    help="The number of sequences"
)
parser.add_argument(
    'outputname',
    help="The name of the outputfile"
)

args = parser.parse_args()

def dna(length=int(), letters="CGTA"):
    return''.join(random.choices(letters, k=length))

# Creating a new file
f = open(f"{args.outputname}.fastq", "w")

for i in range(args.number): 
    f.write(
        '@' + 'Seq' + str(i + 1) + '\n' +
        dna(args.n) + '\n' +
        '+' + '\n' +
        '~' * args.n + '\n'

    )

f.close()
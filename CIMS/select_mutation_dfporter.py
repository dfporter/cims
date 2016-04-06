"""
Covert the mismatch file output by novoalign2bed.pl to
the required format. Selects for only one type of mutation.

Usage:
    python select_mutation_dfporter.py --mut AtoT -i <in.mismatch> -o <out.mismatch>

"""
import glob
import sys
import os
import re
import HTSeq
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--mut', dest='mutation',
                   default='TtoC',
                   help='''Type of mutation (default: TtoC),
options are NtoN, del, ins, sub. sub includes all substitions.''')
parser.add_argument('-i',
                    help="Input novoalign2bed.pl mismatch filename.")

parser.add_argument('-o',
                    help="Output .bed filename")

parser.add_argument('-t',
                    help="tags .bed file, for filtering by unique tags")
args = parser.parse_args()

print args

tags = set()
with open(args.t, 'r') as f:
    for line in f:
        s = line.split('\t')
        tags.add(s[3])


specific_sub = re.match('(\w)to(\w)', args.mutation)
if specific_sub is not None:
    from_nt = specific_sub.group(1)
    to_nt = specific_sub.group(2)
    print "Looking for {i} to {o} substitions.".format(
        i=from_nt, o=to_nt)
    args.mutation = "Specific sub"

if not os.path.exists(args.i):
    print "Input file does not exist."
    sys.exit()
# Example of input file format:
#V	13027027	13027028	HWI-D00256:355:H3KFTBCXX:1:1101:1861:2206#AACTGGCCC	3	-	3	.	+	TGG	1
#I	5788710	5788711	HWI-D00256:355:H3KFTBCXX:1:1101:1908:2226#CACTGGCTC	17	-	17	A	>	G	1

def write_line(cols, output_file):
    output_file.write("\t".join(cols[0:6]) + '\n')

output_file = open(args.o, 'w')
with open(args.i, 'r') as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if cols[3] not in tags: continue
        # col 9: - if del, + if ins, > if sub.
        if cols[8] == '-' and args.mutation == 'del':
            write_line(cols, output_file)
            continue
        if cols[8] == '+' and args.mutation == 'ins':
            write_line(cols, output_file)
            continue
        if cols[8] == '>' and args.mutation == "sub":
            write_line(cols, output_file)
            continue
        if cols[8] == '>' and args.mutation == "Specific sub":
            if cols[7] == from_nt and cols[9] == to_nt:
                write_line(cols, output_file)
output_file.close()


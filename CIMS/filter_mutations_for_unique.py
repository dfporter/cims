"""
Filter a mutations file by a collapsed tag file.

Usage:
python filter_mutations_for_unique.py <tags_file> <mutations_file> <out>

"""
import re
import glob
import sys
import os

tags_filename = sys.argv[1]
mutations_filename = sys.argv[2]

tags = set()
with open(tags_filename, 'r') as f:
    for line in f:
        s = line.split('\t')
        tags.add(s[3])

out_mutations_file = open(sys.argv[3], 'w')
with open(mutations_filename, 'r') as f:
    for line in f:
        s = line.split('\t')
        if s[3] in tags:
            out_mutations_file.write(line)
out_mutations_file.close()

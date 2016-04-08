
import glob
import sys
import os
import re
from fastqs import *
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser("FASTQ > CIMS tables/fasta")
    parser.add_argument('-i', '--fastq', default='fastq/',
        description='Folder of fastq files, required even if not mapping')
    parser.add_argument('-m', '--map', default=False,
        action='store_true',
        description='''Map the input fastq folder, and proceed to write
commands for read collapsing, then exit. Default: False''')
    parser.add_argument('-c', '--cims', default=False,
        action='store_true', description='''
Combine replicates, reformat, split mutation types, call CIMS and
create output tables in cims_tables/''')
    parser.add_argument('-l', '--lib', default='config.ini',
description='''Config ini file of the format:
Experiment class 1 name: FOG-3
Experiment class 1 files: fog_GGCA, fog_CGGA, fog_GGTT, fog_TGGC
Experiment class 2 name: Control
Experiment class 2 files: control_TTGT, control_CCGG, control_TTAA, control_AATA
''')
    args = parser.parse_args()
    fqs = fastqs(path=args.fastq, name='FBF', lib=args.lib)
    print fqs
    if args.map:
        fqs.call_novoalign()
        fqs.convert_novoalign_to_bed()
        fqs.combine_replicates()
        fqs.collapse_duplicate_tags_for_cims(just_print=True)
    if args.cims:
    #    fqs.combine_replicates(top_dir='novo_tags_collapse/')
        #fqs.reformat_collapsed_tags(
        #    input_dir='bed_collapsed/no_rrna',
            #input_dir='novo_tags_collapse',
        #    output_dir='collapsed_reformated')
        #fqs.split_types(in_dir='mismatch_tags/')
        #fqs.combine_replicates_mismatches(top_dir='mismatches_by_type')
        #fqs.call_cims(in_dir='collapsed_reformated/',
        #              mismatches_dir='mismatches_by_type/',
        #              output_dir='cims_out/')
        #fqs.create_tables(cims=True)
    fqs.get_reproducible_peaks(in_dir='cims_tables/', min_reps=2)



import glob
import sys
import os
import re
from fastqs import *

if __name__ == '__main__':
    fqs = fastqs(path=sys.argv[1], name='FOG-3')
    print fqs
    #fqs.call_novoalign()
    #fqs.convert_novoalign_to_bed()
    #fqs.combine_replicates()
    #fqs.collapse_duplicate_tags_for_cims(just_print=True)
    #fqs.combine_replicates(top_dir='novo_tags_collapse/')
    fqs.combine_replicates_mismatches(top_dir='mismatches_by_type')
    #fqs.reformat_collapsed_tags()
    #fqs.split_types()
    fqs.call_cims()
    fqs.create_tables(cims=True)

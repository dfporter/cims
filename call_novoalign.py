"""
python ../clip/src/call_novoalign.py

#Output files are now in ./novoaligned/.

perl ../clip/CIMS/novoalign2bed.pl --mismatch-file mm.bed novoaligned/in.novo out.bed

"""

import glob
import sys
import os
import re



def call_novoalign(input_filenames):
    os.system('mkdir novoaligned/')
    for input_filename in input_filenames:
        output_filename = './novoaligned/' + os.path.basename(
            input_filename).partition('.fastq')[0] + '.novo'
        if os.path.exists(output_filename):
            print "call_novoalign(): Output file %s already exists. Skipping..." % (
                output_filename)
            continue
        cmd = '/opt/novocraft/novoalign -d /scratch/indexes/ws235.novo'
        cmd += ' -f {infile} > {outfile}'.format(
            infile=input_filename, outfile=output_filename)
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    call_novoalign(glob.glob(sys.argv[1] + '/*.fastq'))

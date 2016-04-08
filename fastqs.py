
import glob
import sys
import os
import re
import collections
import HTSeq
import pandas
import HTSeq

def _mk(d):
    if not os.path.exists(d):
        os.system('mkdir ' + d)

class fastqs(object):

    def __init__(self, path='fastqs/', name='Unnamed', lib='none'):
        self.path = path
        self.name = name
        self.files = glob.glob(path + '/*.fastq')
        #self.cims_path = '/groups/Kimble/Aman\ Prasad/clip/CIMS/'
        self.cims_path = '/groups/Kimble/Common/fog_iCLIP/cims/CIMS/'
        if lib == 'none':
            self.lib = """
Experiment class 1 name: FOG-3
Experiment class 1 files: fog_GGCA, fog_CGGA, fog_GGTT, fog_TGGC
Experiment class 2 name: Control
Experiment class 2 files: control_TTGT, control_CCGG, control_TTAA, control_AATA
"""
        else:
            self.lib = "\n".join(open(lib).readlines())
        self.parse_lib(self.lib)
        self.gtf_raw = '/groups/Kimble/Common/fog_iCLIP/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf'
        self.gtf = '/groups/Kimble/Common/fog_iCLIP/lib/gtf_with_names_column.txt'
        for path in [
            self.cims_path, self.gtf_raw, self.gtf,
            '/opt/novocraft/novoalign', '/scratch/indexes/ws235.novo',
            '/scratch/indexes/WS235.fa']:
            if not os.path.exists(path):
                if raw_input(
                    """{a} could not be found. Run anyway?
""".format(a=path)).upper() != 'Y': sys.exit()

    def parse_lib(self, astr):
        self.experiments = collections.defaultdict(dict)
        for li in astr.split('\n'):
            print li,
            s = li.rstrip('\n').split(':')
            mname = re.search('Experiment class (\d+) name', s[0])
            mfiles = re.search('Experiment class (\d+) files', s[0])
            if mname:
                self.experiments[mname.group(1)]['Name'] = \
                    re.sub('\s', '', s[1])
            if mfiles:
                self.experiments[mfiles.group(1)]['Files'] = \
                    [re.sub('\s', '', x) for x in s[1].split(',')]
        self.exp = {}
        for exp in self.experiments.values():
            self.exp[exp['Name']] = set(exp['Files'])
        self.experiments = self.exp
        self.file_to_experiment = collections.defaultdict(str)
        for key in self.experiments:
            for fname in self.experiments[key]:
                self.file_to_experiment[fname] = key
        print self.experiments


    def __str__(self):
        return """
Name: {a}
Input dir: {b}
Input fastq files: {c}
""".format(
    a=self.name, b=self.path, c=self.files)

    def call_novoalign(self):
        _mk('novoaligned/')
        for fastq in self.files:
            output_filename = './novoaligned/' + os.path.basename(
                fastq).partition('.fastq')[0] + '.novo'
            if os.path.exists(output_filename):
                print "call_novoalign(): Output file %s already exists. Skipping..." % (
                    output_filename)
                continue
            cmd = '/opt/novocraft/novoalign -d /scratch/indexes/ws235.novo'
            cmd += ' -f {infile} > {outfile}'.format(
                infile=fastq, outfile=output_filename)
            print cmd
            os.system(cmd)

    def convert_novoalign_to_bed(self, indir='novoaligned/'):
        _mk('novo_tags')
        _mk('mismatch_tags')
        for input_filename in glob.glob(indir + '/*.novo'):
            tags_filename = 'novo_tags/' + os.path.basename(
                input_filename).partition('.novo')[0] + '.bed'
            mismatch_filename = 'mismatch_tags/' + os.path.basename(
                input_filename).partition('.novo')[0] + '.bed'
            cmd = '{a}/novoalign2bed.pl '.format(a=self.cims_path)
            cmd += '--mismatch-file {mm} {i} {o}'.format(
                mm=mismatch_filename, i=input_filename, o=tags_filename)
            print cmd
            os.system(cmd)

    def combine_replicates_mismatches(self, top_dir='mismatches_by_type/'):
        _mk('mismatches_by_type/')
        by_type = collections.defaultdict(dict)
        for filename in glob.glob(top_dir + '/*'):
            bname = os.path.basename(filename).partition('.bed')[0]
            for x in ['_ins', '_del', '_sub']:
                if re.search(x, bname):
                    exp = self.file_to_experiment[bname.partition(x)[0]]
                    by_type[x].setdefault(exp, [])
                    by_type[x][exp].append(bname.partition(x)[0])
        for _type in by_type:
            for exp in by_type[_type]:
                cmd = 'cat '
                cmd += ' '.join([
                    top_dir + '/' + x + _type + '.bed' for \
                    x in by_type[_type][exp]])
                cmd += ' > ' + top_dir + '/all_' + exp + _type + '.bed'
                print cmd
                os.system(cmd)

    def combine_replicates(self, top_dir='novo_tags/'):
        """Identify what the replicate files are and concatenate those bed
        files.
        """
        top_dir = top_dir.rstrip('/') + '/'
        # Combine replicates before collapsing.
        reps = collections.defaultdict(list)
        for filename in glob.glob(top_dir + '/*'):
            basename = os.path.basename(filename).partition('.bed')[0]
            name = self.file_to_experiment[basename]
            reps[name].append(filename)
        for exp in reps:
            output_filename = top_dir + 'all_' + exp + '.bed'
            if os.path.exists(output_filename):
                print output_filename + " exists, skipping..."
                continue
            cmd = 'cat '
            for filename in reps[exp]: cmd += ' %s' % filename
            cmd += ' > ' + output_filename
            print cmd
            os.system(cmd)

    def collapse_duplicate_tags_for_cims(self, just_print=True):
        print "export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/"
        os.system('mkdir novo_tags_collapse')
        #for input_filename in glob.glob('novo_tags/*'):
        for input_filename in glob.glob('novo_tags/*'):
            output_filename = 'novo_tags_collapse/'
            output_filename += os.path.basename(input_filename)
            if os.path.exists(output_filename): continue
            cmd = 'perl {a}/tag2collapse.pl '.format(a=self.cims_path)
            cmd += ' --random-barcode --seq-error-model fix=0.02 '
            cmd += ' -EM 50 --keep-cache {i} {o}'.format(
                i=input_filename, o=output_filename)
            nohup_cmd = 'nohup {c} > {o}.log 2>&1 </dev/null &'.format(
                c=cmd, o=os.path.basename(input_filename).partition('.bed')[0]
                )
            if just_print:
                print nohup_cmd
            else:
                print cmd
                os.system(cmd)

    def reformat_collapsed_tags(self, input_dir='novo_tags_collapse',
                                output_dir='collapsed_reformated'):
        """Writes in the format:
chrm left right name "num from [n=\d]" strand
The n=\d number appears to be the number of duplicates collapsed.
"""
        _mk(output_dir)
        for input_filename in glob.glob(input_dir + '/*'):
            print input_filename
            output_filename = output_dir + '/' + os.path.basename(input_filename)
            out_f = open(output_filename, 'w')
            with open(input_filename, 'r') as f:
                for line in f:
                    s = line.rstrip('\n').split('\t')
                    name = s[3].partition('[')[0]
                    m = re.search('\[n=(\d+)\]', s[3])
                    out_f.write("\t".join([
                        s[0], s[1], s[2], name, m.group(1), s[5]]) + '\n')

    def split_types(self, in_dir='mismatch_tags/'):
        _mk('mismatches_by_type')
        for filename in glob.glob(in_dir + '/*'):
            for this_mut in ['ins', 'del', 'sub']:
                cmd = 'python {a}/select_mutation_dfporter.py '.format(
                    a=self.cims_path)
                cmd += '--mut {mut} -t {tags} -i {i} -o {o}'.format(
                    tags='collapsed_reformated/' + os.path.basename(filename),
                    mut=this_mut, i=filename,
                    o='mismatches_by_type/{base}_{mm_type}.bed'.format(
                        base=os.path.basename(filename).partition('.bed')[0],
                        mm_type=this_mut)
                )
                print cmd
                os.system(cmd)

    def call_cims(self, in_dir='collapsed_reformated/',
                  mismatches_dir='mismatches_by_type/',
                  output_dir='cims_out/'):
        _mk('cims_out')
        for tag_filename in glob.glob(in_dir + '/*'):
            mut_filename = mismatches_dir.rstrip('/') + '/'
            mut_filename += os.path.basename(tag_filename).partition('.bed')[0] + '_del.bed'
            out_filename = output_dir + '/' + os.path.basename(mut_filename)
            cmd = 'perl '
            cmd += '{a}/CIMS.pl --keep-cache -v {tags} {muts} {out}'.format(
                a=self.cims_path, tags=tag_filename, muts=mut_filename,
                out=out_filename)
            print cmd
            os.system(cmd)

    def create_tables(self, cims=True):
        # Load some library files.
        print "create_tables() called."
        fasta_filename = '/scratch/indexes/WS235.fa'
        sequences = dict(
            (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
        gtf_ga = get_gtf(None, gtf_file=self.gtf_raw)
        gtf_df = pandas.read_csv(
            self.gtf, sep='\t')
        if cims:
            globstr = 'cims_out/*'
            cits_option = False
            fdr = 0.001
            table_dir = 'cims_tables/'
            _mk('cims_tables')
        else:
            create_cits_tables(sequences, gtf_df)
            return
        for filename in glob.glob(globstr):
            bname = os.path.basename(filename)
            # Apply filter also renames columns.
            print "create_tables(): %s" % filename
            print "Loading peaks..."
            peaks = pandas.read_csv(filename, sep='\t')
            print peaks
            #peaks = peaks.head()
            rename_columns(peaks)
            peaks = apply_filter(peaks, cits=cits_option, fdr=fdr)
            assign_to_gene(peaks, gtf_ga)
            assign_biotype(peaks, gtf_df)
            peaks = peaks[peaks['biotype']=='protein_coding']
            peaks.sort('height', ascending=0, inplace=True)
            print peaks
            add_sequences(peaks, sequences, expand=10)
            write_fasta(peaks, 'fasta/' + os.path.basename(filename))
            peaks.to_csv(table_dir + '/' + bname, sep='\t', index=False)
            print "Finished processing..."

    def get_reproducible_peaks(self, in_dir='cims_tables/', min_reps=2):
        by_rep = {}
        df = {}
        self.min_reps = min_reps
        for x in glob.glob(in_dir + '/*'):
            bname = os.path.basename(x)
            exp = None
            for _file in self.file_to_experiment:
                if re.search(_file, bname):
                    exp = self.file_to_experiment[_file]
            if exp is None: continue
            by_rep.setdefault(exp, {})
            df.setdefault(exp, {})
            df[exp][x] = pandas.read_csv(x, sep='\t')
            by_rep[exp][x] = zip(
                df[exp][x].chrm, df[exp][x].left, df[exp][x].strand)
        _mk('reproducible/')
        for exp in df:
            repro_df = self.get_reproducible_peaks_this_exp(
                df[exp], by_rep[exp])
            repro_df.to_csv('reproducible/' + exp, sep='\t', index=False)
            write_fasta(repro_df, 'fasta/' + os.path.basename(exp))

    def get_reproducible_peaks_this_exp(self, df, loci_by_rep):
        rep_counts = collections.defaultdict(int)
        for _file in loci_by_rep:
            for loci in set(loci_by_rep[_file]):
                rep_counts[loci] += 1
        for _file in df:
            df[_file]['n_reps'] = [
                rep_counts[x] for x in loci_by_rep[_file]]
            df[_file] = df[_file][df[_file]['n_reps']>=self.min_reps]
        df_exp = self.take_highest(df)
        return df_exp

    def take_highest(self, df):
        alist = []
        for _file in df:
            alist.extend(df[_file].to_dict('records'))
        highest_at_loci = {}
        for row in alist:
            loci = (row['chrm'], row['left'], row['strand'])
            if loci in highest_at_loci:
                if row['height'] > highest_at_loci[loci]['height']:
                    highest_at_loci[loci] = row
            else:
                highest_at_loci[loci] = row
        df['combined'] = pandas.DataFrame(highest_at_loci.values())
        return df['combined']

def assign_biotype(peaks, gtf_df):
    to_biotype = collections.defaultdict(str)
    to_biotype.update(dict(zip(gtf_df.gene_name, gtf_df.biotype)))
    peaks['biotype'] = [to_biotype[x] for x in peaks['gene_name'].tolist()]


def assign_to_gene(peaks, features):
    _ivs = zip(peaks.chrm, peaks.left, peaks.right, peaks.strand)
    _genes = []
    for tup in _ivs:
        gene_ids = set()
        for sub_iv, gene in features[HTSeq.GenomicInterval(*tup)].steps():
            gene_ids |= gene
        if len(gene_ids) == 1:
            _genes.append(list(gene_ids)[0])
        elif len(gene_ids) == 0:
            _genes.append('_no_feature')
        else:
            _genes.append('_ambiguous')
    peaks['gene_name'] = _genes


def get_gtf(lib, gtf_file=None):
    if gtf_file is not None:
        gtf_path = gtf_file
    else:
        gtf_path = lib['gtf_raw']
    print gtf_path
    gtf_file = HTSeq.GFF_Reader(gtf_path)
    features = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    for feature in gtf_file:
        if feature.type == 'exon':
            features[feature.iv] += feature.attr['gene_name']
    return features


def seq_from_iv(sequences, x, expand=0):
    left = max([x[1] - expand, 0])
    right = min([x[2] + expand, len(sequences[x[0]])])
    if x[3] == '+':
        return sequences[x[0]][left:right]
    if x[3] == '-':
        return rc(sequences[x[0]][left:right])


def add_sequences(peaks, sequences, expand=10):
    _ivs = zip(peaks.chrm, peaks.left, peaks.right, peaks.strand)
    peaks['seq'] = [seq_from_iv(sequences, x, expand=expand) for x in _ivs]


def rename_columns(peaks, cits=False):
    if cits:
        peaks.columns = ['chrm', 'left', 'right', 'name',
                         'score', 'strand']
        peaks['height'] = peaks['score']
    else:
        peaks.columns = ['chrm', 'left', 'right', 'name', 'score', 'strand',
                         'tagNumber_k', 'mutationFreq_m', 'FDR',
                         'count_at_least_m_k']
        peaks['height'] = peaks['tagNumber_k']


def apply_filter(peaks, cits=False, fdr=False,
                 min_k=False, min_score=False,
                 highest_per_gene=False):
    """Applies a filter to dataframe.
    """
    if fdr and not cits:
        peaks = peaks[peaks['FDR'] < fdr]
    if min_k and not cits:
        peaks = peaks[peaks['tagNumber_k']>=min_k]
    return peaks

def write_fasta(peaks, filename):
    _mk(os.path.dirname(filename))
    tups = zip(peaks.gene_name, peaks.name, peaks.seq)
    li = [">{i}_{a}\n{s}\n".format(i=x[0], a=x[1], s=x[2]) for \
          x in tups]
    li = ''.join(li)
    with open(filename, 'w') as f: f.write(li)


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s

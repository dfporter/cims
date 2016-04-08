import collections
import os
import sys
import re
import glob


def read_dimers(fname, as_freq=True):
    dimers = collections.defaultdict(int)
    keys = []
    bases = ['A', 'T', 'C', 'G']
    for base_a in bases:
        keys.extend([tuple([base_a, b]) for b in bases])
    for k in keys: dimers[k] = 0
    seqs = []
    with open(fname) as f:
        while True:
            id_line = f.readline()
            seq_line = f.readline()
            seqs.append(seq_line.rstrip('\n').upper())
            if not seq_line: break
    total_len = 0
    for seq in seqs:
        for i in range(0, len(seqs), 1):
            # ACG:
            # i=0 -> seqs[0:2] = AC
            # i=2 -> 3 >= len(seqs) -> break
            if i + 1 >= len(seq): break
            dimers[tuple(seq[i:i+2])] += 1
            total_len += 1
    if as_freq:
        for k in dimers:
            dimers[k] = float(dimers[k])/float(total_len)
    return dimers, total_len

fname = sys.argv[1]
dimers, total_len = read_dimers(fname, as_freq=True)
print dimers
print total_len
bases = ['A', 'T', 'C', 'G']
import matplotlib.pyplot as plt
import numpy as np
image = []
for a in ['A', 'T', 'C', 'G']:
    image.append(np.array([dimers[tuple([a, b])] for b in bases]))
image = np.array(image)
fig, ax = plt.subplots()
ax.imshow(image, interpolation='None')
x, y = np.meshgrid(np.arange(0,4,1), np.arange(0,4,1))
#for x_val, y_val in zip(x.flatten(), y.flatten()):
for x_val, a in enumerate(bases):
    for y_val, b in enumerate(bases):
        #c = 'x' if (x_val + y_val)%2 else 'o'
        c = "%.2f" % dimers[(bases[y_val], bases[x_val])]
        #c += str((bases[y_val], bases[x_val]))
        ax.text(x_val, y_val, c, va='center', ha='center')
ax.set_ylabel('First base')
ax.set_xlabel('Second base')
ax.set_xticks(range(0,4))
ax.set_xticklabels(bases)
ax.set_yticks(range(0,4))
ax.set_yticklabels(bases)
#ax.grid()
plt.savefig('%s_dinucleotide_frequencies.pdf' % os.path.basename(fname),
            format='pdf')

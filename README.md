CIMS analysis
----

|    Folder     | Description |
|    ------     | ----------- |
|cims/cims_tables | Results of cims analysis |
|cims/CIMS      | Zhang lab CIMS scripts with modifications|
|cims/fasta     | Fasta files from cims analysis |
|cims/novoaligned | Raw output from novoalign in .novo |
|cims/novo_tags | .novo files in novoaligned converted to .bed|
|cims/novo_tags_collapse | novo_tags after Zhang collapsing (still bed)|
|cims/collapsed_reformated | novo_tags_collapse reformated (still bed) for input into CIMS analysis |
|cims/mismatch_tags | bed-like, holds read and mutation |
|cims/mismatches_by_type | mismatch_tags split by ins/del/sub |
|cims_out | Raw output of CIMS scripts |

There is a cims github repo.

```bash
# This maps the fastq files and writes nohup commands for read collapsing.
python cims_master.py -i FOLDER_OF_FASTQ_FILES --lib fog_cims.ini --map
# Then run the read collapsing.
# Then do the cims analysis.
python cims_master.py -i FOLDER_OF_FASTQ_FILES --lib fog_cims.ini --cims
# Output is in cims_tables/

# To get dinucleotide frequencies in a fasta.
python dimers.py FASTA_FILE
# Outputs a heatmap.

# To get frequncy near a morif of interest:
pyhton pos_vs_motif.py -i cims_tables/ -m A_MOTIF
```


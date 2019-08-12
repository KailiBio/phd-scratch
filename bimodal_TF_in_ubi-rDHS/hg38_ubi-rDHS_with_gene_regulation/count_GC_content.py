#!/usr/bin/env python

# -- Kaili
# This script is for counting GC content of given fasta.

# INPUT: fasta file & output file
# OUTPUT: txt file with GC context (two column, uniqID&GCcontent)
# EXP: python count_GC_content.py \
#			/data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/TSS_overlapped_rOCRs.fa \
#			/data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/TSS_overlapped_rOCRs_GCcontent.txt

import re, os, sys
from collections import Counter

def calculate_GC_content(fasta_file, output):
    for line in open(fasta_file).readlines():
        if line.startswith(">"):
            name = line.rstrip("\n")[1:]
        else:
            count = Counter(line.rstrip("\n").upper())
            GC_content = round(float(count['C']+count['G'])/(count['A']+count['T']+count['C']+count['G']),2)
            print >> output, "\t".join([name, str(GC_content)])


if __name__ == "__main__":
    fasta_file = sys.argv[1]
    outDir = sys.argv[2]

    output = open(outDir, 'w')
    calculate_GC_content(fasta_file, output)
    output.close()

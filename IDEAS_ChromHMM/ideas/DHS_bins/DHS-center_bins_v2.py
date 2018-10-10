#!/usr/bin/env python

# -- Kaili
# This script is for getting version 2 DHS-center bins.
# INPUT: dhs.bed file (only with 3 columns: chr, start, end; sorted would be better.)
#        dhs_gap.bed file (only with 3 columns: chr, start, end; sorted would be better.)
# OUTPUT:dhs-center_bins file (3 columns: chr, start, end, uniq_name)
# EXP: python DHS-center_bins_v2.py /data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed
#      /data/zusers/fankaili/ideas/dhs_bins/mm10_rOCR_gap.bed
#      /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v1.bed

import re, os, sys

# dhs_file = "/data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed"
# dhs_gap_file = "/data/zusers/fankaili/ideas/dhs_bins/mm10_rOCR_gap.bed"
# out_file = "/data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v1.bed"

def cut_bins(file, output, character, with_small_bins=0):
    num = 0
    for line in file.readlines():
        chr, start, end = line.rstrip().split("\t")
        length = int(end) - int(start)
        if length >= 100:
            if length%200 > 100:
                n = length//200 + 1
            else:
                n = length//200
            s = int(start)
            for i in range(0,n):
                num = num+1
                if i!=(n-1):
                    print >> output, ("\t").join([chr, str(s), str(s+200), character+str(num)])
                    s = s+200
                else:
                    print >> output, ("\t").join([chr, str(s), end, character+str(num)])
        elif length < 100 and with_small_bins == 1:
            num = num+1
            print >> output, ("\t").join([chr, start, end, character+str(num)])


if __name__ == "__main__":
    dhs_file = sys.argv[1]
    dhs_gap_file = sys.argv[2]
    out_file = sys.argv[3]

    output = open(out_file, "w+")

    with open(dhs_file) as f1:
        cut_bins(f1, output, "O", with_small_bins=1)

    with open(dhs_gap_file) as f2:
        cut_bins(f2, output, "G", with_small_bins=1)

    output.close()

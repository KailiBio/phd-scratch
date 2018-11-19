#!/usr/bin/env python

# -- Kaili
# This script is for getting version 3 DHS-center bins.
# INPUT: dhs.bed file (only with 3 columns: chr, start, end; sorted would be better.)
#        dhs_gap.bed file (only with 3 columns: chr, start, end; sorted would be better.)
# OUTPUT:dhs-center_bins file (3 columns: chr, start, end, uniq_name)
# EXP: python DHS-center_bins_v2.py /data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed
#      /data/zusers/fankaili/ideas/dhs_bins/mm10_rOCR_gap.bed
#      /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.bed

## update at Nov19
# make each gap as odd number of bins for center.
# all input gaps are >100bp

import re, os, sys
import subprocess

# dhs_file = "/data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed"
dhs_gap_file = "/data/zusers/fankaili/ideas/dhs_bins/mm10_rOCR_gap_v3_2.bed"
out_file = "/data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3_gap.bed"

def cut_bins(file, output, character, with_small_bins=0):
    num = 0
    for line in file.readlines():
        chr, start, end = line.rstrip().split("\t")
        length = int(end) - int(start)
        # decide how many bins for this gap
        n0 = length/200
        if length<350 :
            n = 1
            double = False
        elif length>=350 and length<600:
            n = 3
            double = False
        elif n0%2 == 0:
            if (length - (n0-2)*(length/(n0-1))) < 350:
                n = n0-1
                double = False
            elif (length - (n0-2)*(length/(n0-1))) < 700:
                n = n0-1+2
                double = True
            elif (length - (n0)*(length/(n0+1))) > 350:
                n = n0+3
                double = True
            else:
                n = n0+1
                double = False
        else:
            if (length - (n0-1)*(length/(n0))) > 350:
                n = n0+2
                double = True
            else:
                n = n0
                double = False
        # get the length of each bins:
        if double == False:
            l1 = length/n
            l2 = length - (n-1)*l1
        else:
            l1 = length/(n-2)
            ll = length - (n-3)*l1
            l2 = ll/3
            l3 = ll - 2*l2
        # get bins
        s = int(start)
        for i in range(0,n):
            num = num + 1
            if double == False:
                if i == n/2:
                    print >> output, ("\t").join([chr, str(s), str(s+l2), character+str(num)])
                    s = s + l2
                else:
                    print >> output, ("\t").join([chr, str(s), str(s+l1), character+str(num)])
                    s = s + l1
            else:
                if i == n/2:
                    print >> output, ("\t").join([chr, str(s), str(s+l3), character+str(num)])
                    s = s + l3
                elif i==(n/2-1) or i==(n/2+1):
                    print >> output, ("\t").join([chr, str(s), str(s+l2), character+str(num)])
                    s = s + l2
                else:
                    print >> output, ("\t").join([chr, str(s), str(s+l1), character+str(num)])
                    s = s + l1


if __name__ == "__main__":

    output = open(out_file, "w+")
    with open(dhs_gap_file) as f2:
        cut_bins(f2, output, "G", with_small_bins=1)

    output.close()



# output = open(out_file, "w+")
# with open(dhs_gap_file) as f2:
#     cut_bins(f2, output, "G", with_small_bins=1)
#
# output.close()

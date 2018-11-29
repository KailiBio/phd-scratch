#!/usr/bin/env python

# -- Kaili
# This script is for cutting bins into 50 small windows for making aggregation plot.
# INPUT: bed file
#        output file for small windows in upsteam 500bp
#        output file for small windows in bins given
#        output file for small windows in downstream 500bp
# OUTPUT: 3 small windows bed files
# EXP: python cut_bins_into_50.py \
#       /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_signal_based.bed \
#       /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_50windows_bins.bed \
#       /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_50windows_up.bed \
#       /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_50windows_down.bed

import re, os, sys
import subprocess

def cut_bins_into_50(chr, start, end, prefix, output, m, n=50):
    '''
    This function is for cutting given bins into 50 windows.
        m is the prefix number: 0 for upstream, 50 for bins, 100 for downstream.
        n is the bins to cut, default is 50.
    '''
    length = (float(end)-float(start))/n
    s = start
    for i in range(1,51):
        x = s
        y = s+length
        name = prefix+"_"+str(m+i)
        print >> output, ("\t").join([chr, str(int(round(x))), str(int(round(y))), name])
        s = x+length


if __name__ == "__main__":

    bins = sys.argv[1]
    bins_out_file = sys.argv[2]
    upstream_out_file = sys.argv[3]
    downstream_out_file = sys.argv[4]

    output1 = open(upstream_out_file, "w+")
    output2 = open(bins_out_file, "w+")
    output3 = open(downstream_out_file, "w+")
    for line in open(bins).readlines():
        chr, start, end, prefix = line.rstrip().split("\t")
        start = int(start)
        end = int(end)
        #
        if start-500>=0:
            cut_bins_into_50(chr, start-500, start, prefix, output1, 0, 50)
        else:
            cut_bins_into_50(chr, 0, start, prefix, output1, 0, 50)
        cut_bins_into_50(chr, start, end, prefix, output2, 50, 50)
        cut_bins_into_50(chr, end, end+500, prefix, output3, 100, 50)

    output1.close()
    output2.close()
    output3.close()

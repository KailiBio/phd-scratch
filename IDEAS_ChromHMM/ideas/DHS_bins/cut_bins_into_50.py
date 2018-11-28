#!/usr/bin/env python

# -- Kaili
# This script is for cutting bins into 50 small windows for making aggregation plot.

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
    bins = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_signal_based.bed"
    out_file = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_50windows.bed"

    # bins = "/data/zusers/fankaili/ideas/run_ideas_p_value/mm10_tab.bed"
    # out_file = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_normal_bins_50windows.bed"

    output = open(out_file, "w+")
    for line in open(bins).readlines():
        chr, start, end, prefix = line.rstrip().split("\t")
        start = int(start)
        end = int(end)
        #
        cut_bins_into_50(chr, start-500, start, prefix, output, 0, 50)
        cut_bins_into_50(chr, start, end, prefix, output, 50, 50)
        cut_bins_into_50(chr, end, end+500, prefix, output, 100, 50)

    output.close()

#!/usr/bin/env python

# -- Kaili
# This script is for calculte averge TSS TS index for each gene.

import re, os, sys
import subprocess

def write_file(file, output_path):
	o= open(output_path, 'w')
	o.writelines(file)
	o.close()

if __name__ == "__main__":

ts_index_file = sys.argv[1]

ts_index_file = "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_TSS_exp_TSscore.txt"

with open(ts_index_file) as f1:
	for line in f1.readlines():
		tss, index = line.strip().split("\t")

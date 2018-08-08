#!/usr/bin/env python

# -- Kaili
# This is script for making zscore matrix based on zscore file list.
# INPUT: zscore file list
#        eg: BE2C	ENCSR000DQB-ENCFF001EQH.txt
#            Caco-2	ENCSR000DQM-ENCFF001ESW.txt
#            DOHH2	ENCSR052WRV-ENCFF985CGU.txt
# OUTPUT: zscore matrix

# EXP: python get_zscore.py /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_H3K4me3_zscore_filelist.txt /data/zusers/fankaili/ccre/tf/zscore_h3k4me3/ H3K4me3


import sys,os
import subprocess

inFile = sys.argv[1]
outDir = sys.argv[2]
mark = sys.argv[3]

filelist = open(inFile).readlines()

# out = []

for line in filelist:
    cellline = line.rstrip().split("\t")[0]
    filename = line.rstrip().split("\t")[1]
    command1 = """awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$2}}}' /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt /data/zusers/moorej3/Registry-of-ccREs/hg19/V4/signal-output/"""+filename+" > "+outDir+"hg19_ubi-rDHS_"+cellline+"_"+mark+"_zscore.txt\n"
    subprocess.call(command1, shell=True)

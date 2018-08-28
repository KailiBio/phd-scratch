#!/usr/bin/env python

# -- Kaili
# This script is for getting corresponding zscore file name from master_file: /data/zusers/moorej3/ENCODE-Registry/hg19/V4/Cell-Type-Specific/Master-Cell-List.txt
# INPUT: 1) file lines from master file
#			eg: ---	---	---	---	---	---	ENCSR000DPF	ENCFF001GVQ	A549	Group12
#				ENCSR346JWH	ENCFF891YXC	---	---	ENCSR714TJD	ENCFF743TET	---	---	A673	Group5
#        2) what kind of data: DNase, H3K4me3, H3K27ac, CTCF
# OUTPUT: file list

# EXP: python get_zscore_file_list_from_master.py /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_data_master_file.txt H3K4me3 /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_H3K4me3_zscore_filelist.txt

import sys,os

def write_file(file, output_path):
	o= open(output_path, 'w')
	o.writelines(file)
	o.close()

inFile = sys.argv[1]
mark = sys.argv[2]
outFile = sys.argv[3]
infile = open(inFile).readlines()

out = []
if mark=="DNase":
    for line in infile:
        id1 = line.rstrip().split("\t")[0]
        id2 = line.rstrip().split("\t")[1]
        if id1!="---":
            cellline = line.rstrip().split("\t")[8]
            outline = cellline + "\t" + id1 + "-" + id2 + ".txt\n"
            out.append(outline)
elif mark=="H3K4me3":
    for line in infile:
        id1 = line.rstrip().split("\t")[2]
        id2 = line.rstrip().split("\t")[3]
        if id1!="---":
            cellline = line.rstrip().split("\t")[8]
            outline = cellline + "\t" + id1 + "-" + id2 + ".txt\n"
            out.append(outline)
elif mark=="H3K27ac":
    for line in infile:
        id1 = line.rstrip().split("\t")[4]
        id2 = line.rstrip().split("\t")[5]
        if id1!="---":
            cellline = line.rstrip().split("\t")[8]
            outline = cellline + "\t" + id1 + "-" + id2 + ".txt\n"
            out.append(outline)
elif mark=="CTCF":
    for line in infile:
        id1 = line.rstrip().split("\t")[6]
        id2 = line.rstrip().split("\t")[7]
        if id1!="---":
            cellline = line.rstrip().split("\t")[8]
            outline = cellline + "\t" + id1 + "-" + id2 + ".txt\n"
            out.append(outline)

write_file(out, outFile)

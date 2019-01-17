#!/usr/bin/env python

# -- Kaili
# This script is for getting mouse rep1 ATAC data.

import re, os, sys
import subprocess

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

if __name__ == "__main__":
    datalist = "/data/zusers/fankaili/ideas/mm10_tissue_used_list_CTCF_peak_list.txt"
    outpath = "/data/zusers/fankaili/ideas/ENCODE_mouse_rep1_CTCF_filelist.txt"
    outfile = open(outpath, "w+")

    for line in open(datalist).readlines():
        line = line.split("\t")
        expID = line[3]
        myexp = Exp.fromJsonFile(expID)
        sample = "-".join(myexp.biosample_term_name.split(" ")) + "_" + myexp.age
        for myfile in myexp.files:
            if myfile.output_type=="signal p-value" and myfile.file_format=="bigWig" and myfile.bio_rep==[1] and myfile.assembly=="mm10" and myfile.file_status=="released":
                fileID = myfile.fileID
                print >> outfile, ("\t").join([sample, "CTCF", "rep1", expID, fileID])

    outfile.close()

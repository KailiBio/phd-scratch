#!/usr/bin/env python

# Kaili
# This script is for getting TF file list.
# INPUT: URL from ENCODE
# OUTPUT: file list: /data/zusers/fankaili/ccre/tf/encode_tf_file_list/

import re, os, sys
import subprocess

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

def write_file(file, output_path):
	o= open(output_path, 'w')
	o.writelines(file)
	o.close()

##################################
ref = "hg19"
outDir = "/data/zusers/fankaili/ccre/tf/encode_tf_file_list/"
urlFile = sys.argv[1]

infile = open(urlFile).readlines()

if __name__ == "__main__":
    for line in infile:
        cellline = line.rstrip().split("\t")[0]
        URL = line.rstrip().split("\t")[1]
        print(cellline)
        q = QueryDCC(auth=False)
        out = []
        for exp in q.getExps(URL):
            for f in exp.files:
                myexp = Exp.fromJsonFile(f.expID)
                for myfile in myexp.files:
                    if myfile.bio_rep==[1, 2] and myfile.file_format=="bigWig" and myfile.output_type=="fold change over control" and myfile.assembly=="hg19":
                        x = (f.expID, myfile.accession, exp.label)
                        out.append('\t'.join(x)+"\n")
        write_file(out, outDir+"encode_"+ref+"_"+cellline+"_tf_id_list0.txt")
        subprocess.call("sort "+outDir+"encode_"+ref+"_"+cellline+"_tf_id_list0.txt | uniq > "+outDir+"encode_"+ref+"_"+cellline+"_tf_id_list.txt", shell=True)

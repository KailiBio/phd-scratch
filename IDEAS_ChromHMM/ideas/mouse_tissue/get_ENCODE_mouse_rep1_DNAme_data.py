#!/usr/bin/env python

# -- Kaili
# This script is for getting mouse rep1 DNAme data.
### get CpG bedMethyl file from Bismark output here.
### can use bedMethyl2bigWig to transform to bigWig file for signal calculation.

import json
import re, os, sys

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

if __name__ == "__main__":
    jfile = json.load(open("/data/public_html_users/fankaili/IDEAS/allFiles.json"))
    data = jfile["input data"]
    outpath="/data/zusers/fankaili/ideas/ENCODE_mouse_rep1_DNAme_filelist.txt"
    outfile = open(outpath, "w+")

    for i in range(len(data)):
        if data[i]['assay_term_name']=='whole-genome shotgun bisulfite sequencing':
            sample = "-".join(data[i]['biosample_term_name'].split(" ")) + "_" + data[i]['age']
            expID = data[i]["accession"]
            myexp = Exp.fromJsonFile(expID)
            for myfile in myexp.files:
                if myfile.file_format=="bed" and myfile.output_type=="methylation state at CpG" and myfile.bio_rep==[1] and myfile.assembly=="mm10" and myfile.file_status=="released":
                    fileID = myfile.fileID
                    print >> outfile, ("\t").join([sample, "DNAme", "rep1", expID, fileID])

    outfile.close()

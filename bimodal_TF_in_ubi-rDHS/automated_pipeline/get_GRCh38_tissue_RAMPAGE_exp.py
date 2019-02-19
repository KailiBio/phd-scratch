#!/usr/bin/env python

# -- Kaili
# This script is for getting GRCh38 RNA-seq exp data in all 113 tissues from ENCODE tsv file.

# INPUT:
# OUTPUT:
# EXP: python get_GRCh38_tissue_RAMPAGE_exp.py /home/fankaili/genome/hg38_v28_comprehensive_TSS.bed
#			/data/zusers/fankaili/ccre/hg38_ubi-rDHS/tissue_rampage_v28/

import re, os, sys
import subprocess

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

##################################
if __name__ == "__main__":
    url = "https://www.encodeproject.org/search/?type=Experiment&assay_title=RAMPAGE&assay_slims=Transcription&lab.title=Thomas+Gingeras%2C+CSHL&limit=all&format=json"
    bedfile = sys.argv[1]
    signal_folder = sys.argv[2]

    q = QueryDCC(auth=False)
    expID = []
    for exp in q.getExps(url):
        if exp.biosample_type=="cell line":
	    sample = exp.biosample_term_name.replace(" ","_")
	else:
	    sample = ("_").join([exp.biosample_term_name.replace(" ","_"), exp.age_display.replace(" ","_")])
	expID = exp.files[0].expID
        plus = ""
        minus = ""
        for f in exp.files:
            if f.bio_rep==[1] and f.tech_rep==['1_1'] and f.file_format=="bigWig" and f.output_type=="plus strand signal of unique reads" and f.assembly=="GRCh38":
		plus = f.accession
		print(plus)
		subprocess.call("bigWigAverageOverBed /data/projects/encode/data/"+expID+"/"+plus+".bigWig "+bedfile+" "+signal_folder+plus+".tab", shell=True)
	    elif f.bio_rep==[1] and f.tech_rep==['1_1'] and f.file_format=="bigWig" and f.output_type=="minus strand signal of unique reads" and f.assembly=="GRCh38":
	        minus = f.accession
		print(minus)
		subprocess.call("bigWigAverageOverBed /data/projects/encode/data/"+f.expID+"/"+minus+".bigWig "+bedfile+" "+signal_folder+minus+".tab", shell=True)
            if plus != "" and minus != "" :
                break

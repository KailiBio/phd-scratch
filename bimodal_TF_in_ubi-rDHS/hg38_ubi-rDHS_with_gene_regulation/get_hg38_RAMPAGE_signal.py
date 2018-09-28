#!/usr/bin/env python

# -- Kaili
# This script is for getting TSS signal from GRCh38 155 matching RAMPAGE data (from ENCODE bigWig file).


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

if __name__ == "__main__":
    url = "https://www.encodeproject.org/search/?type=Experiment&assay_title=RAMPAGE&assay_slims=Transcription&lab.title=Thomas+Gingeras%2C+CSHL&limit=all&format=json"
    outDir = "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt"

    q = QueryDCC(auth=False)
    expID = []
    out = []
    for exp in q.getExps(url):
	sample = ("_").join([myexp.biosample_term_name.replace(" ","_"), myexp.age_display.replace(" ","_")])
	expID = exp.files[0].expID
        for f in exp.files:
            if myfile.bio_rep==[1] and myfile.tech_rep==['1_1'] and myfile.file_format=="bigWig" and myfile.output_type=="plus strand signal of unique reads" and myfile.assembly=="GRCh38":
		plus = myfile.accession
		print(plus)
		subprocess.call("bigWigAverageOverBed /data/projects/encode/data/"+f.expID+"/"+plus+".bigWig /data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered_sorted.bed /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/"+plus+".tab", shell=True)
	    elif myfile.bio_rep==[1] and myfile.tech_rep==['1_1'] and myfile.file_format=="bigWig" and myfile.output_type=="minus strand signal of unique reads" and myfile.assembly=="GRCh38":
		minus = myfile.accession
		print(minus)
		subprocess.call("bigWigAverageOverBed /data/projects/encode/data/"+f.expID+"/"+minus+".bigWig /data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered_sorted.bed /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/"+minus+".tab", shell=True)
	out.append(("\t").join([f.expID, plus, minus, sample])+"\n")
    write_file(out,outDir)

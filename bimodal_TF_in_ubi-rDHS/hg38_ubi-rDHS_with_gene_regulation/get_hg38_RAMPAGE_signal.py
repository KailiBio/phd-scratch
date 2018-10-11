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
				#subprocess.call("bigWigAverageOverBed /data/projects/encode/data/"+expID+"/"+plus+".bigWig /data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered_sorted.bed /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/"+plus+".tab", shell=True)
	    	elif f.bio_rep==[1] and f.tech_rep==['1_1'] and f.file_format=="bigWig" and f.output_type=="minus strand signal of unique reads" and f.assembly=="GRCh38":
				minus = f.accession
				print(minus)
				#subprocess.call("bigWigAverageOverBed /data/projects/encode/data/"+f.expID+"/"+minus+".bigWig /data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered_sorted.bed /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/"+minus+".tab", shell=True)
            if plus != "" and minus != "" :
                break
		out.append(("\t").join([expID, plus, minus, sample])+"\n")
    write_file(out,outDir)

#!/bin/env python

# -- Kaili
# This script is for getting mouse ATAC-seq data from ENCODE url.

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
	url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=ATAC-seq&assay_slims=DNA+accessibility&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&biosample_type=tissue&award.project=ENCODE&limit=all&format=json"
	outDir="/data/zusers/fankaili/ideas/mouse_atac_data_list.txt"

	expID = []
	out = []
	q = QueryDCC(auth=False)
	for exp in q.getExps(url):
		for f in exp.files:
			if not f.expID in expID:
				expID.append(f.expID)
				if exp.status=="released":
					sample = ("_").join([exp.biosample_term_name.replace(" ","_"), exp.age_display.replace(" ","_")])
					out.append(("\t").join([f.expID, sample, exp.assay_title])+"\n")

	out2 = []
	for line in out:
		line = line.rstrip().split("\t")
		myfile = Exp.fromJsonFile(line[0]).files
		for file in myfile:
			if file.file_type=="bigWig" and file.output_type=="signal p-value" and file.bio_rep==[1]:
				out2.append(("\t").join([line[1], line[2], line[0], file.accession])+"\n")

	write_file(out2,outDir)

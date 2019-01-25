#!/usr/bin/env python

# -- Kaili
# This script is for getting mouse RNA-seq data list from ENCODE tsv file.

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
	url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=polyA+RNA-seq&biosample_type=tissue&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&lab.title=Barbara+Wold%2C+Caltech&limit=all&format=json"
	outDir="/data/zusers/fankaili/ideas/mouse_RNA_data_list.txt"

	expID = []
	q = QueryDCC(auth=False)
	for exp in q.getExps(url):
		for f in exp.files:
			if not f.expID in expID:
				expID.append(f.expID)

	out = []
	for id in expID:
		myexp = Exp.fromJsonFile(id)
		sample = "-".join(myexp.biosample_term_name.split(" ")) + "_" + myexp.age
		for myfile in myexp.files:
			if myfile.output_type=="gene quantifications" and myfile.file_format=="tsv" and myfile.bio_rep==[1] and myfile.assembly=="mm10" and myfile.file_status=="released":
				out.append(("\t").join([sample, id, myfile.accession])+"\n")

	write_file(out,outDir)

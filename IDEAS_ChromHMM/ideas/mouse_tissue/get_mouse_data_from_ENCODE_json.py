#!/usr/bin/env python

# -- Kaili
# This script is for getting mouse tissue data list from ENCODE tsv file.

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
	url="https://www.encodeproject.org/search/?type=Experiment&assay_title=ChIP-seq&assay_title=DNase-seq&assay_title=WGBS&assay_title=ATAC-seq&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&biosample_type=tissue&assembly=mm10&award.project=ENCODE&limit=all&format=json"
	outDir="/data/zusers/fankaili/ideas/mouse_tissue_data_list.txt"

	expID = []
	out = []
	q = QueryDCC(auth=False)
	for exp in q.getExps(url):
		for f in exp.files:
			if not f.expID in expID:
				expID.append(f.expID)
				if exp.status=="released":
					sample = ("_").join([exp.biosample_term_name.replace(" ","_"), exp.age_display.replace(" ","_")])
					out.append(("\t").join([expID, sample, exp.assay_title, exp.label])+"\n")

	write_file(out,outDir)

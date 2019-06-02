#!/usr/bin/env python

# -- Kaili
# This script is for getting all mouse RNA-seq data list.

import re, os, sys
import subprocess

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

##################################

if __name__ == "__main__":
	url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_slims=Transcription&assay_title=polyA+RNA-seq&assay_title=total+RNA-seq&assembly=mm10&limit=all&format=json"
	outDir="/data/zusers/fankaili/ccre/mm10_rnaseq/mouse_RNA_data_list.txt"
	outfile = open(outDir, "w+")

	expID = []
	q = QueryDCC(auth=False)
	for exp in q.getExps(url):
		for f in exp.files:
			if not f.expID in expID:
				expID.append(f.expID)

	for id in expID:
		myexp = Exp.fromJsonFile(id)
		assay="-".join(myexp.assay_title.split(" "))
		sample = "-".join(myexp.biosample_summary.split(" "))
		print (id)
		for myfile in myexp.files:
			if myfile.output_type=="alignments" and myfile.file_format=="bam" and myfile.assembly=="mm10" and myfile.file_status=="released":
				print (myfile.accession)
				if (len(myfile.jsondata["analysis_step_version"]["analysis_step"]["pipelines"])!=0 and len(myfile.jsondata["analysis_step_version"]["analysis_step"]["pipelines"][0]["title"].split("("))>1):
					exp_type= ";".join(myfile.jsondata["analysis_step_version"]["analysis_step"]["pipelines"][0]["title"].split("(")[1].split(")")[0].split(", "))
				else:
					exp_type="-"
				outline = ("\t").join([id, myfile.accession, str(myfile.bio_rep[0]), str(myfile.tech_rep[0]), sample, assay, exp_type ])
				print >> outfile, outline.encode('ascii', 'ignore')

	outfile.close()

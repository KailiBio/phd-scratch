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
	expID_file="/data/zusers/fankaili/ccre/mm10_rnaseq/mouse_RNA_data_list_clean1_bam.txt"
	outDir="/data/zusers/fankaili/ccre/mm10_rnaseq/mouse_RNA_fastq_list.txt"
	outfile = open(outDir, "w+")

	for line in open(expID_file).readlines():
		id=line.split("\t")[0]
		#
		myexp = Exp.fromJsonFile(id)
		assay="-".join(myexp.assay_title.split(" "))
		sample = "-".join(myexp.biosample_summary.split(" "))
		print (id)
		for myfile in myexp.files:
			if myfile.output_type=="reads" and myfile.file_format=="fastq" and myfile.file_status=="released":
				print (myfile.accession)
				outline = ("\t").join([id, myfile.accession, str(myfile.bio_rep[0]), str(myfile.tech_rep[0]), sample, assay, "unstr_SE" ])
				print >> outfile, outline.encode('ascii', 'ignore')

	outfile.close()

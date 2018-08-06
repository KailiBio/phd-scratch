#!/usr/bin/env python

# Kaili
# This script is for getting the list of biosamples that from the given ENCODE URL.
# INPUT: URL from ENCODE. (Don't give it the URL of the matrix.)
# OUTPUT: /data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_tf_cellline_list0.txt

import re, os, sys

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

def write_file(file, output_path):
	o= open(output_path, 'w')
	o.writelines(file)
	o.close()

##################################
URL = "https://www.encodeproject.org/search/?type=Experiment&assay_title=ChIP-seq&target.investigated_as=transcription+factor&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assembly=hg19&biosample_type=cell+line&limit=all&format=json"
outDir = "/data/zusers/fankaili/ccre/tf/encode_tf_file_list/"

if __name__ == "__main__":
	q = QueryDCC(auth=False)
	out = []
	for exp in q.getExps(URL):
		for f in exp.files:
			myexp = Exp.fromJsonFile(f.expID)
			out.append(myexp.biosample_term_name+"\n")
	write_file(out, outDir+"encode_hg19_tf_cellline_list0.txt")

#!/usr/bin/env python

# -- Kaili
# This script is for getting TF peak file form ENCODE.
# INPUT:
# OUTPUT:

import re, os, sys
import subprocess

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

##################################

if __name__ == "__main__":
    url = sys.argv[1]
    outDir = sys.argv[2]
	out = open(outDir,"w")

    q = QueryDCC(auth=False)
    for exp in q.getExps(url):
		if exp.status=="released":
			expID = exp.files[0].expID
			TF = exp.label
			for f in exp.files:
				if f.file_format=="bed" and f.output_type=="optimal idr thresholded peaks" and f.assembly=="GRCh38":
					outline = "\t".join([TF, expID, f.accession])
					print >> out, outline

	out.close()


# url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=ChIP-seq&target.investigated_as=transcription+factor&assembly=GRCh38&biosample_ontology.term_name=K562&assay_title=ChIP-seq&limit=all&format=json"
# outDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/K562_TF_peak_filelist.txt"

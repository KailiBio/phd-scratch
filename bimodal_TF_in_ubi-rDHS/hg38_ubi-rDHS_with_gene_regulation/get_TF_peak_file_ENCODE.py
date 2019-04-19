#!/usr/bin/env python

# -- Kaili
# This script is for getting GRCh38 TF peak file.

import re, os, sys
import subprocess

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

##################################
if __name__ == "__main__":
    url = "https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=ChIP-seq&assembly=GRCh38&target.investigated_as=transcription+factor&limit=all&format=json"
    outDir = "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/TF_motif_peak_filelist.txt"
    #
    tf_list_file = "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/JASPAR_motif_list.txt"
    tf_list = {}
    for line in open(tf_list_file).readlines():
        tf_list[line.rstrip()] = 1

    q = QueryDCC(auth=False)
    expID = []
    out = open(outDir, "w")
    for exp in q.getExps(url):
		if tf_list.has_key(exp.label):
            mark = exp.label
            expID = exp.files[0].expID
            for f in exp.files:
                if f.bio_rep==[1] and f.tech_rep==['1_1'] and f.file_format=="bed" and  f.assembly=="GRCh38" and f.output_type!="peaks and background as input for IDR":
                    fileID = f.accession
                    peak_type = f.output_type
                    print >> out, ("\t").join([expID, fileID, mark, peak_type])
                elif f.bio_rep==[1,2] and f.file_format=="bed" and  f.assembly=="GRCh38" and f.output_type=="optimal idr thresholded peaks":
                    fileID = f.accession
                    peak_type = f.output_type
                    print >> out, ("\t").join([expID, fileID, mark, peak_type])

    out.close()

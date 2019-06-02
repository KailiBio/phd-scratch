#!/usr/bin/env python

# -- Kaili
# This script is for getting mouse rep2 ATAC data.

import re, os, sys
import subprocess

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

import requests
from requests.auth import HTTPBasicAuth


if __name__ == "__main__":
    expURL = "https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=ATAC-seq&assembly=mm10&biosample_ontology.classification=tissue&lab.title=Bing+Ren%2C+UCSD&limit=all&format=json"
    outpath = "/data/zusers/fankaili/ideas/ENCODE_mouse_rep2_ATAC_filelist.txt"
    outfile = open(outpath, "w+")

    # get all ATAC experiment ID
    expID = []
    q = QueryDCC(auth=False)
    for exp in q.getExps(expURL):
        for file in exp.files:
	    if not file.expID in expID:
                expID.append(file.expID)

    # get rep1 signal p-value bigWig file
    for exp in expID:
        myexp = Exp.fromJsonFile(exp)
        # sample = "-".join(myexp.biosample_term_name.split(" ")) + "_" + myexp.age
        #
        url = "https://www.encodeproject.org/experiments/"+exp+"/"
        headers = {'accept': 'application/json'}
        response = requests.get(url, headers=headers, auth=HTTPBasicAuth('WFK5V4G5', 'dult5uijwz4akf6o'))
        jfile = response.json()
        for i in range(len(jfile["original_files"])):
            # query from each file
            file = jfile["original_files"][i].split("/")[2]
            fileURL = "https://www.encodeproject.org/files/"+file+"/"
            file_response = requests.get(fileURL, headers=headers, auth=HTTPBasicAuth('WFK5V4G5', 'dult5uijwz4akf6o'))
            file_jfile = file_response.json()
            if file_jfile["file_type"]=="bigWig" and file_jfile["biological_replicates"]==[2] and file_jfile["output_type"]=="signal p-value":
                fileID = file_jfile["accession"]
                print >> outfile, ("\t").join(["ATAC", "rep2", exp, fileID])

    outfile.close()

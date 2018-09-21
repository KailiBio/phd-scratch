#!/usr/bin/env python

# -- Kaili
# This script is for getting GRCh38 RNA-seq exp data in all 113 tissues from ENCODE tsv file.
# INPUT:
# OUTPUT:

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
    url = "https://www.encodeproject.org/search/?type=Experiment&assay_title=total+RNA-seq&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&lab.title=Thomas+Gingeras%2C+CSHL&month_released=June%2C+2014&month_released=February%2C+2016&month_released=August%2C+2016&month_released=May%2C+2016&month_released=November%2C+2014&month_released=May%2C+2015&month_released=June%2C+2015&month_released=March%2C+2016&month_released=December%2C+2014&month_released=January%2C+2015&month_released=January%2C+2016&month_released=July%2C+2014&month_released=June%2C+2016&biosample_type=tissue&limit=all&format=json"
    outDir = "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_gene_exp_list.txt"

    q = QueryDCC(auth=False)
    expID = []
    out = []
    for exp in q.getExps(url):
        for f in exp.files:
            if not f.expID in expID:
                expID.append(f.expID)
                myexp = Exp.fromJsonFile(f.expID)
                for myfile in myexp.files:
                    if myfile.bio_rep==[1] and myfile.tech_rep==['1_1'] and myfile.file_format=="tsv" and myfile.output_type=="gene quantifications" and myfile.assembly=="GRCh38":
                        print(f.expID)
                        print(myfile.accession)
                        sample = ("_").join([exp.biosample_term_name.replace(" ","_"), exp.age_display.replace(" ","_")])
                        out.append(("\t").join([f.expID, myfile.accession ,exp.biosample_term_name.replace(" ","_"), sample])+"\n")
                        subprocess.call("""grep "ENSG" /data/projects/encode/data/"""+f.expID+"/"+myfile.accession+".tsv | cut -f 1,6 | sort -k1 > /data/zusers/fankaili/ccre/hg38_ubi-rDHS/tissue_gene_exp/"+f.expID+".txt", shell=True)

    write_file(out,outDir)

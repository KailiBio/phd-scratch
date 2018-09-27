#!/usr/bin/env python

# -- Kaili
# This script is for getting codes for making codes for getting signal file.

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
    url = "https://www.encodeproject.org/search/?type=Experiment&assay_title=ChIP-seq&assay_title=WGBS&assay_title=DNase-seq&assay_title=ATAC-seq&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&biosample_type=tissue&assembly=mm10&organ_slims=brain&organ_slims=embryo&organ_slims=heart&organ_slims=liver&organ_slims=limb&organ_slims=intestine&organ_slims=kidney&organ_slims=lung&organ_slims=stomach&limit=all&format=json"
    outDir = "/data/zusers/fankaili/ideas/mm10_tissue_used_list.txt"
    codeDir = "/data/zusers/fankaili/ideas/code/get_macs2_pvalue_code/get_macs2_pvalue_signal_code.sh"

    tissue_list = ["forebrain", "midbrain", "hindbrain", "neural tube", "heart", "embryonic facial prominence", "limb", "kidney", "lung", "stomach", "intestine", "liver"]

    q = QueryDCC(auth=False)
    expID = []
    out = []
    code = []
    for exp in q.getExps(url):
        if exp.biosample_term_name in tissue_list:
            # sample name
            sample = ("_").join([exp.biosample_term_name.replace(" ","_"), exp.age_display.replace(" ","_")])
            # mark
            if exp.assay_title=="ChIP-seq":
                mark = exp.label
            else:
                mark = exp.assay_title
            #out.append(sample+"\t"+mark+"\n")
            # experiment ID
            expID = exp.files[0].expID
            # file ID
            for f in exp.files:
                if mark=="WGBS" and f.output_type=="signal" and f.file_format=="bigWig" and f.bio_rep==[1] and f.tech_rep==['1_1'] and f.assembly=="mm10":
                    fileID = f.accession
                elif f.output_type=="signal p-value" and f.file_format=="bigWig" and f.bio_rep==[1] and f.tech_rep==['1_1'] and f.assembly=="mm10":
                    fileID = f.accession
            # code line for getting signal
            code.append("bigWigAverageOverBed /data/projects/encode/data/"+expID+"/"+fileID+".bw /data/zusers/fankaili/ideas/run_ideas_p_value/mm10_tab.bed /data/zusers/fankaili/ideas/signal/macs2_pvalue/"+sample+"_"+mark+"_macs2_pvalue.tab"+"\n")
            if mark=="WGBS":
                code.append("cut -f 6 /data/zusers/fankaili/ideas/signal/macs2_pvalue/"+sample+"_"+mark+"_macs2_pvalue.tab > /data/zusers/fankaili/ideas/signal/macs2_pvalue/"+sample+"_"+mark+"_macs2_pvalue.txt"+"\n")
            else:
                code.append("cut -f 5 /data/zusers/fankaili/ideas/signal/macs2_pvalue/"+sample+"_"+mark+"_macs2_pvalue.tab > /data/zusers/fankaili/ideas/signal/macs2_pvalue/"+sample+"_"+mark+"_macs2_pvalue.txt"+"\n")
            # input
            filePath="/data/zusers/fankaili/ideas/signal/macs2_pvalue/"+sample+"_"+mark+"_macs2_pvalue.txt"
            out.append(("\t").join([sample, mark, filePath, expID, fileID])+"\n")

    write_file(out,outDir)
    write_file(code,codeDir)

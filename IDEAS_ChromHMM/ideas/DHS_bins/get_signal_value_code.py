#!/usr/bin/env python

# -- Kaili
# This script is for getting 66 biosamples on 10 marks signal pvalue file for IDEAS from json and metadata.
# INPUT: bins bed file for calculating signal
#				output file.
# OUTPUT: code file
# EXP: python get_signal_pvalue_code.py /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_random_50windows_bins.bed
#	   /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/code_for_getting_signal_p_value_v3_bins.txt \
#	  /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_v3_bins/

import json
import os,sys

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

###########

def write_file(file, output_path):
	o= open(output_path, 'w')
	o.writelines(file)
	o.close()

################
jfile = json.load(open("/data/public_html_users/fankaili/IDEAS/allFiles.json"))
data = jfile["input data"]

bedfile = sys.argv[1]
output_file = sys.argv[2]
outpath = sys.argv[3]


code = []
for i in range(len(data)):
    if data[i]['assay_term_name']=='whole-genome shotgun bisulfite sequencing':
        if (data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube"):
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_DNAme"
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+bedfile+" "+outpath+signal_file_name+".tab"+"\n")
        #code.append("""awk '{print $6}' """+outpath+signal_file_name+".tab > "+outpath+signal_file_name+".txt"+"\n")
	code.append("""awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $4,a[$4]}}' """+outpath+signal_file_name+".tab "+bedfile+" > "+outpath+signal_file_name+".txt"+"\n")
    elif data[i]['assay_term_name']=='ATAC-seq':
        if data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube":
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_ATAC"
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+bedfile+" "+outpath+signal_file_name+".tab"+"\n")
        #code.append("""awk '{print $5}' """+outpath+signal_file_name+".tab > "+outpath+signal_file_name+".txt"+"\n")
	code.append("""awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $4,a[$4]}}' """+outpath+signal_file_name+".tab "+bedfile+" > "+outpath+signal_file_name+".txt"+"\n")
    elif data[i]['assay_term_name']=='ChIP-seq':
        if data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube":
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_" + data[i]['label']
        myexp = Exp.fromJsonFile(data[i]['accession'])
        for myfile in myexp.files:
        	if myfile.output_type=="signal p-value" and myfile.file_format=="bigWig" and myfile.bio_rep==[1, 2] :
        		id = myfile.fileID
        code.append("bigWigAverageOverBed /data/projects/encode/data/"+data[i]['accession']+"/"+id+".bigWig "+bedfile+" "+outpath+signal_file_name+".tab"+"\n")
        #code.append("""awk '{print $5}' """+outpath+signal_file_name+".tab > "+outpath+signal_file_name+".txt"+"\n")
	code.append("""awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $4,a[$4]}}' """+outpath+signal_file_name+".tab "+bedfile+" > "+outpath+signal_file_name+".txt"+"\n")

write_file(code, output_file)

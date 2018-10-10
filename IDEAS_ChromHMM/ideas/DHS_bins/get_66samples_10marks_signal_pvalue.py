#!/usr/bin/env python

# -- Kaili
# This script is for getting 66 biosamples on 10 marks signal pvalue file for IDEAS from json and metadata.
# INPUT: bins bed file for calculating signal (must be sorted by chromosome, bigWigAverageOverBed problem)
#		 output path.
# OUTPUT: code file and input file
# EXP: python get_66samples_10marks_signal_pvalue.py /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v1.bed
#	   /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/

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
output_path = sys.argv[2]


code = []
input = []
for i in range(len(data)):
    if data[i]['assay_term_name']=='whole-genome shotgun bisulfite sequencing':
        if (data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube"):
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_DNAme"
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+bedfile+" "+output_path+"signal/"+signal_file_name+".tab"+"\n")
        code.append("""awk '{print $6}' """+output_path+"signal/"+signal_file_name+".tab > "+output_path+"signal/"+signal_file_name+".txt"+"\n")
        input.append(biosample+" DNAme "+ output_path + "signal/" + signal_file_name + ".txt"+"\n")
    elif data[i]['assay_term_name']=='ATAC-seq':
        if data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube":
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_ATAC"
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+bedfile+" "+output_path+"signal/"+signal_file_name+".tab"+"\n")
        code.append("""awk '{print $5}' """+output_path+"signal/"+signal_file_name+".tab > "+output_path+"signal/"+signal_file_name+".txt"+"\n")
        input.append(biosample+" ATAC "+ output_path + "signal/" + signal_file_name + ".txt"+"\n")
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
        code.append("bigWigAverageOverBed /data/projects/encode/data/"+data[i]['accession']+"/"+id+".bigWig "+bedfile+" "+output_path+"signal/"+signal_file_name+".tab"+"\n")
        code.append("""awk '{print $5}' """+output_path+"signal/"+signal_file_name+".tab > "+output_path+"signal/"+signal_file_name+".txt"+"\n")
        input.append(biosample+" "+data[i]['label']+" "+ output_path + "signal/" + signal_file_name + ".txt"+"\n")

write_file(code, output_path+"codes_for_getting_signal_pvalue.sh")
write_file(input, output_path+"signal_pvalue.input")

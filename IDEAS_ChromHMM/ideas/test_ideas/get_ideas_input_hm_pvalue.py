#!/usr/bin/env python

# -- Kaili
# This script is for getting signal file for IDEAS from json and metadata.
# use signal p-value for histone data here.

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

jfile = json.load(open("/data/public_html_users/fankaili/IDEAS/allFiles.json"))
data = jfile["input data"]

path = "/data/zusers/fankaili/ideas/run_ideas_p_value/"

################

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
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"mm10_tab.bed "+path+"signal/"+signal_file_name+".tab"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+".tab > "+path+"signal/"+signal_file_name+".txt"+"\n")
        input.append(biosample+" DNAme "+ path + "signal/" + signal_file_name + ".txt"+"\n")
    elif data[i]['assay_term_name']=='ATAC-seq':
        if data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube":
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_ATAC"
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"mm10_tab.bed "+path+"signal/"+signal_file_name+".tab"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+".tab > "+path+"signal/"+signal_file_name+".txt"+"\n")
        input.append(biosample+" ATAC "+ path + "signal/" + signal_file_name + ".txt"+"\n")
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
        code.append("bigWigAverageOverBed /data/projects/encode/data/"+data[i]['accession']+"/"+id+".bigWig "+path+"mm10_tab.bed "+path+"signal/"+signal_file_name+".tab"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+".tab > "+path+"signal/"+signal_file_name+".txt"+"\n")
        input.append(biosample+" "+data[i]['label']+" "+ path + "signal/" + signal_file_name + ".txt"+"\n")


write_file(code, path+"codes_for_getting_signal_pvalue.sh")
write_file(input, path+"run_IDEAS_8hm_atac_dname_pvalue.input")

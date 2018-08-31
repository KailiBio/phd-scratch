#!/usr/bin/env python

# -- Kaili
# This script is for getting signal file for IDEAS from json and metadata.
# use bam file for histone data here.

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
data = jfile["models"]["chromhmm_15_state_8marks_200bp_more"]["bamfiles"]

path = "/data/zusers/fankaili/ideas/run_ideas_repeat/"

################

code = []
input = []
for biosample in data.keys():
	list = []
	for i in range(len(data[biosample])):
		hm = data[biosample][i]["mark"]
		case = data[biosample][i]["bam"]
		control = data[biosample][i]["control"]
		signal_file_name = biosample+"_"+hm
		if hm in list:
			code.append("bamToBed -i "+case+" > /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_case_2.bed"+"\n")
			code.append("intersectBed -a /data/zusers/fankaili/ideas/run_ideas_repeat/mm10_tab.bed -b /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_case_2.bed -wa -a > /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_case_2_count.bed"+"\n")
			code.append("bamToBed -i "+control+" > /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_control_2.bed"+"\n")
			code.append("intersectBed -a /data/zusers/fankaili/ideas/run_ideas_repeat/mm10_tab.bed -b /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_control_2.bed -wa -a > /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_control_2_count.bed"+"\n")
		else:
			list.append(hm)
			code.append("bamToBed -i "+case+" > /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_case_1.bed"+"\n")
			code.append("intersectBed -a /data/zusers/fankaili/ideas/run_ideas_repeat/mm10_tab.bed -b /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_case_1.bed -wa -a > /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_case_1_count.bed"+"\n")
			code.append("bamToBed -i "+control+" > /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_control_1.bed"+"\n")
			code.append("intersectBed -a /data/zusers/fankaili/ideas/run_ideas_repeat/mm10_tab.bed -b /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_control_1.bed -wa -a > /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_control_1_count.bed"+"\n")
			code.append("Rscript /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/test_ideas/negBinomial.R /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_case_1_count.bed /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_control_1_count.bed /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_1.txt"+"\n")
			input.append(biosample+" "+hm+" /data/zusers/fankaili/ideas/run_ideas_repeat/signal/"+signal_file_name+"_1.txt"+"\n")


write_file(code, path+"codes_for_getting_signal_bam.sh")
write_file(input, path+"run_IDEAS_8hm_atac_dname_repeat.input")

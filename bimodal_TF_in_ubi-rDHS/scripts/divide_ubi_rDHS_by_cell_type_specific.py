#!/usr/bin/env python

# -- Kaili
# This script is for divide ubi-rDHS into PLS/ELS/CTCF-only based on cell-type specific definition.
# INPUT
# OUTPUT

import re, os, sys
import subprocess

fileList=open("/data/zusers/fankaili/ccre/tf/cell_type_specific/hg19_ccREs_cell_type_specific_definition_file_list.txt").readlines()

ubirDHS="/data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt"
outDir="/data/zusers/fankaili/ccre/tf/cell_type_specific/"

if __name__ == "__main__":
    for line in fileList:
        cellline = line.rstrip().split("\t")[0]
        filename = line.rstrip().split("\t")[1]
        filename2 = filename.split(".gz")[0]
        command1 = "gunzip -c /data/projects/screen/Version-4/ver10/hg19/public_html/cts_formatted/"+filename+" > "+outDir+"temp_"+filename2
        command2 = """awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4] && $9=="255,0,0"){print $4}}}' """+ubirDHS+" "+outDir+"temp_"+filename2+" > "+outDir+"PLS/ubi-rDHS_"+cellline+"_PLS_list.txt"
        command3 = """awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4] && $9=="255,205,0"){print $4}}}' """+ubirDHS+" "+outDir+"temp_"+filename2+" > "+outDir+"ELS/ubi-rDHS_"+cellline+"_ELS_list.txt"
        command4 = """awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4] && $9=="0,176,240"){print $4}}}' """+ubirDHS+" "+outDir+"temp_"+filename2+" > "+outDir+"CTCF/ubi-rDHS_"+cellline+"_CTCF_list.txt"
        command5 = """awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4] && $9=="6,218,147"){print $4}}}' """+ubirDHS+" "+outDir+"temp_"+filename2+" > "+outDir+"DNase/ubi-rDHS_"+cellline+"_DNase_list.txt"
        command6 = """awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4] && $9=="225,225,225"){print $4}}}' """+ubirDHS+" "+outDir+"temp_"+filename2+" > "+outDir+"inactive/ubi-rDHS_"+cellline+"_inactive_list.txt"
        subprocess.call(command1, shell=True)
        subprocess.call(command2, shell=True)
        subprocess.call(command3, shell=True)
        subprocess.call(command4, shell=True)
        subprocess.call(command5, shell=True)
        subprocess.call(command6, shell=True)
        subprocess.call("rm "+outDir+"temp_"+filename2, shell=True)

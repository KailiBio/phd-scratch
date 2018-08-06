#!/usr/bin/env python

# -- Kaili
# This python scipt for calculating CTCF signal of ccREs in all biosamples.
# Get mean signal for replicates.
# INPUT: 1) a given cell line. (cell line list file: /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF.txt)
#        2) expID, fileID list files of given biosamples: /data/zusers/fankaili/ccre_old/tf/TF_list/
# OUTPUT: CTCF signal of ccREs in all celllines: /data/zusers/fankaili/ccre/tf/signal/

from __future__ import print_function

import re, os, sys
import subprocess

cellline=sys.argv[1]
infile = open("/data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_"+cellline+"_tf_id_list.txt").readlines()
outDir = "/data/zusers/fankaili/ccre/tf/signal/"

# save expID and fileID
fileDic = {}
for line in infile:
    if line.rstrip().split("\t")[2]=="CTCF":
        expID = line.rstrip().split("\t")[0]
        fileID = line.rstrip().split("\t")[1]
        if fileDic.has_key(expID):
            fileDic[expID].append(fileID)
        else:
            fileDic[expID] = [fileID]

for key,value in fileDic.items():
	print (key)
	if len(value)==1:
		subprocess.call("bigWigAverageOverBed /data/projects/encode/data/"+key+"/"+value[0]+".bigWig /data/zusers/fankaili/ccre/hg19-cREs.bed "+outDir+"CTCF_signal_ccREs_"+cellline+"_"+key+".txt", shell=True)
                print(value[0])
                subprocess.call("""awk '{FS=OFS="\t"}{print $1,$5}' """+outDir+"CTCF_signal_ccREs_"+cellline+"_"+key+".txt > "+outDir+"CTCF_signal_ccREs_"+cellline+"_"+key+"_final.txt", shell=True)
        else:
		subprocess.call("bigWigAverageOverBed /data/projects/encode/data/"+key+"/"+value[0]+".bigWig /data/zusers/fankaili/ccre/hg19-cREs.bed "+outDir+"CTCF_signal_ccREs_"+cellline+"_"+key+"_1.txt", shell=True)
                subprocess.call("bigWigAverageOverBed /data/projects/encode/data/"+key+"/"+value[1]+".bigWig /data/zusers/fankaili/ccre/hg19-cREs.bed "+outDir+"CTCF_signal_ccREs_"+cellline+"_"+key+"_2.txt", shell=True)
		subprocess.call("""awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{mean=(a[$1]+$5)/2; print $1,mean}}' """+outDir+"CTCF_signal_ccREs_"+cellline+"_"+key+"_1.txt "+outDir+"CTCF_signal_ccREs_"+cellline+"_"+key+"_2.txt > "+outDir+"CTCF_signal_ccREs_"+cellline+"_"+key+"_final.txt", shell=True)

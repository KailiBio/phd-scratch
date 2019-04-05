#!/usr/bin/env python

# -- Kaili
# This script is for getting merged TF peaks for each TFs (encode).

import re, os, sys
import subprocess

##################################
if __name__ == "__main__":
    peak_file = sys.argv[1]
    outDir = sys.argv[2]

    current_motif = ""
    subprocess.call("if [ -f "+outDir+" ];then rm "+outDir+";fi", shell=True)
    subprocess.call("if [ -f tmp.bed ];then rm tmp.bed;fi",shell=True)
    for line in open(peak_file).readlines():
        motif = line.rstrip().split("\t")[2]
        expID = line.rstrip().split("\t")[0]
        fileID = line.rstrip().split("\t")[1]

        if motif == current_motif:
            subprocess.call("if [ -f /data/projects/encode/data/"+expID+"/"+fileID+".bed ]; then cat /data/projects/encode/data/"+expID+"/"+fileID+".bed >> tmp.bed; else cp /data/projects/encode/data/"+expID+"/"+fileID+".bed.gz ./; gzip -d "+fileID+".bed.gz; cat "+fileID+".bed >> tmp.bed; rm "+fileID+".bed;fi", shell=True)
        else:
            subprocess.call('if [ -f tmp.bed ];then sort -k1,1 -k2,2n tmp.bed | cut -f 1-3 > tmp2.bed; bedtools merge -i tmp2.bed | awk -v current_motif="'+current_motif+'"'+""" '{OFS="\t"}{print $0,current_motif}' >> """+outDir+";fi", shell=True)
            current_motif = motif
            subprocess.call("if [ -f tmp.bed ];then rm tmp.bed;fi",shell=True)
            subprocess.call("if [ -f /data/projects/encode/data/"+expID+"/"+fileID+".bed ]; then cat /data/projects/encode/data/"+expID+"/"+fileID+".bed >> tmp.bed; else cp /data/projects/encode/data/"+expID+"/"+fileID+".bed.gz ./; gzip -d "+fileID+".bed.gz; cat "+fileID+".bed >> tmp.bed; rm "+fileID+".bed;fi", shell=True)

    subprocess.call("if [ -f /data/projects/encode/data/"+expID+"/"+fileID+".bed ]; then cat /data/projects/encode/data/"+expID+"/"+fileID+".bed >> tmp.bed; else cp /data/projects/encode/data/"+expID+"/"+fileID+".bed.gz ./; gzip -d "+fileID+".bed.gz; cat "+fileID+".bed >> tmp.bed; rm "+fileID+".bed;fi", shell=True)
    subprocess.call("rm tmp.bed tmp2.bed", shell=True)

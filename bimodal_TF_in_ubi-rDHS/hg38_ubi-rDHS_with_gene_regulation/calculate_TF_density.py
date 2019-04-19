#!/usr/bin/env python

# -- Kaili
# This script is for calculating TF density in rOCRs.

import re, os, sys
import subprocess

##################################
if __name__ == "__main__":
    motif_file = sys.argv[1]
    outDir = sys.argv[2]
    outDir2 = sys.argv[3]

    id = ""
    subprocess.call("if [ -f "+outDir+" ];then rm "+outDir+"; fi", shell=True)
    for line in open(motif_file).readlines():
        ocr = line.rstrip().split("\t")[3]
        if id==ocr:
            print >> out, line.rstrip()
        else:
            if id!="":
                out.close()
                subprocess.call('bedtools merge -i tmp.bed | awk -v id="'+id+'"'+""" 'BEGIN{FS=OFS="\\t";sum=0}{sum+=($3-$2)}END{print id,sum}' >> """+outDir, shell=True)
                subprocess.call('cut -f 5 tmp.bed | sort -u | wc -l | awk -v id="'+id+'"'+""" '{FS=OFS="\\t"}{print id,$1}' >> """+outDir2, shell=True)

            id = ocr
            out = open("tmp.bed", "w")
            print >> out, line.rstrip()

    out.close()
    subprocess.call('bedtools merge -i tmp.bed | awk -v id="'+id+'"'+""" 'BEGIN{FS=OFS="\\t";sum=0}{sum+=($3-$2)}END{print id,sum}' >> """+outDir, shell=True)
    subprocess.call('cut -f 5 tmp.bed | sort -u | wc -l | awk -v id="'+id+'"'+""" '{FS=OFS="\\t"}{print id,$1}' >> """+outDir2, shell=True)
    #
    subprocess.call("rm tmp.bed", shell=True)

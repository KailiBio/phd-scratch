#!/usr/bin/env python

# This script is for merging TSSs within 50bp for each gene.

# INPUT: tss file, output file, and working direactry
# OUTPUT: merged_TSSs bed file
# EXP: python merge_TSS_50bp.py hg38_v28_basic_TSS_filtered_uniq.bed merged_TSS.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

import re, os, sys
import subprocess

def merge_TSS(tss_file, outFile, workDir):
    # sort TSS by gene
    subprocess.call("sort -k7,7 -k1,1 -k2,2n "+tss_file+" > "+workDir+"tmp.tss_sorted.bed", shell=True)
    #
    gene = ""
    #
    for line in open(workDir+"tmp.tss_sorted.bed").readlines():
        if gene=="":
            gene = line.rstrip().split("\t")[6]
            output = open(workDir+"tmp.bed","w")
            a = line.rstrip().split("\t")
            print >> output, "\t".join([a[0], a[1], a[2], a[6], a[4], a[5], a[6]])
        elif gene == line.rstrip().split("\t")[6]:
            a = line.rstrip().split("\t")
            print >> output, "\t".join([a[0], a[1], a[2], a[6], a[4], a[5], a[6]])
        else:
            output.close()
            subprocess.call("bedtools merge -i "+workDir+"tmp.bed -s -d 50 > "+workDir+"tmp2.bed", shell=True)
            subprocess.call("""awk 'BEGIN{FS=OFS="\t";i=0}{if(NR==FNR){if(NR==1){a=$6;b=$7}}else{i+=1;print $0,b"_"i,".",a,b}}' """+workDir+"tmp.bed "+workDir+"tmp2.bed >> "+workDir+outFile, shell=True)
            gene = line.rstrip().split("\t")[6]
            output = open(workDir+"tmp.bed","w")
            a = line.rstrip().split("\t")
            print >> output, "\t".join([a[0], a[1], a[2], a[6], a[4], a[5], a[6]])
    #
    output.close()
    subprocess.call("bedtools merge -i "+workDir+"tmp.bed -s -d 50 > "+workDir+"tmp2.bed", shell=True)
    subprocess.call("""awk 'BEGIN{FS=OFS="\t";i=0}{if(NR==FNR){if(NR==1){a=$6;b=$7}}else{i+=1;print $0,b"_"i,".",a,b}}' """+workDir+"tmp.bed "+workDir+"tmp2.bed >> "+workDir+outFile, shell=True)
    #
    subprocess.call("rm "+workDir+"tmp.bed", shell=True)
    subprocess.call("rm "+workDir+"tmp2.bed", shell=True)
    subprocess.call("rm "+workDir+"tmp.tss_sorted.bed", shell=True)


##################################
if __name__ == "__main__":
    tss_file = sys.argv[1]
    outFile = sys.argv[2]
    workDir = sys.argv[3]

    # initiation
    subprocess.call("if [ -f "+workDir+outFile+" ]; then rm "+workDir+outFile+"; fi", shell=True)
    subprocess.call("touch "+workDir+outFile, shell=True)

    # merge
    merge_TSS(tss_file, outFile, workDir)

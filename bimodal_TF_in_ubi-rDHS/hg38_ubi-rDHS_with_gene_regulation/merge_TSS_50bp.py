#!/usr/bin/env python

# This script is for merging TSSs within 50bp for each gene.


import re, os, sys
import subprocess

##################################
if __name__ == "__main__":
    gene = ""
    output = open("/data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp.bed","w")
    subprocess.call("rm /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene.bed", shell=True)
    subprocess.call("touch /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene.bed", shell=True)
    subprocess.call("touch /data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp2.bed", shell=True)

    for line in open("/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered_sorted_gene.bed").readlines():
        if gene == "":
            gene = line.rstrip().split("\t")[6]
            output = open("/data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp.bed","w")
            a = line.rstrip().split("\t")
            print >> output, "\t".join([a[0], a[1], a[2], a[6], a[4], a[5], a[6]])
        elif gene == line.rstrip().split("\t")[6]:
            a = line.rstrip().split("\t")
            print >> output, "\t".join([a[0], a[1], a[2], a[6], a[4], a[5], a[6]])
        else:
            output.close()
            subprocess.call("bedtools merge -i /data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp.bed -s -d 50 > /data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp2.bed", shell=True)
            subprocess.call("""awk 'BEGIN{FS=OFS="\t";i=0}{if(NR==FNR){if(NR==1){a=$6;b=$7}}else{i+=1;print $0,b"_"i,".",a,b}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp.bed /data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp2.bed >> /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_merged_TSS_gene.bed""", shell=True)
            gene = line.rstrip().split("\t")[6]
            output = open("/data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp.bed","w")
            a = line.rstrip().split("\t")
            print >> output, "\t".join([a[0], a[1], a[2], a[6], a[4], a[5], a[6]])

    subprocess.call("rm /data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp.bed", shell=True)
    subprocess.call("rm /data/zusers/fankaili/ccre/hg38_ubi-rDHS/temp2.bed", shell=True)

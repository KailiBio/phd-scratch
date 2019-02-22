#!/bin/bash

# -- Kaili
# This script is for proving TSSs that overlap ubi-rOCRs are likely to cluster together.

# INPUT:
# OUTPUT:
# EXP: bash 3_TSS_cluster.sh
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

# tss_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed"
# outPath="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/"

tss_labeled_file=$1
outPath=$2

###################

# 1. get gene with singular TSS
cut -f 7 ${tss_labeled_file} | sort | uniq -u > ${outPath}gene_singular_TSS_list.txt

# 2. calculate nearest distance in gene with multiple TSSs
# get multiple TSSs file
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]!=1){print $0}}}' ${outPath}gene_singular_TSS_list.txt ${tss_labeled_file} | sort -k7,7 -k2,2n > ${outPath}TSSs_in_multiple_TSSs_genes.txt
# calculate distance
awk 'BEGIN{FS=OFS="\t";gene="";tss="";end="";distance=1000000}{if(gene==""){gene=$7;tss=$4;end=$3}else{if(gene==$7){if(distance>($2-end)){print gene,tss,$2-end}else{print gene,tss,distance};tss=$4;distance=$2-end;end=$3}else{print gene,tss,distance;gene=$7;tss=$4;end=$3;distance=1000000}}}END{print gene,tss,distance}' ${outPath}TSSs_in_multiple_TSSs_genes.txt > ${outPath}tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$3;b[$2]=1}else{if(b[$4]){print $0,a[$4]}}}' ${outPath}tmp.txt ${tss_labeled_file} | sort -u > ${outPath}nearest_TSS_distance_in_multiple_TSS_genes.txt

rm ${outPath}tmp.txt

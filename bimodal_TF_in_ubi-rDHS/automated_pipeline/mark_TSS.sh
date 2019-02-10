#!/bin/bash

# -- Kaili
# This script is for labeling TSSs as:
#          overlap_with_ubi-rOCRs, overlap_with_not-ubi_active-rOCRs, no_overlap
#          overlap_with_ubi-rOCRs, in_gene_overlapping_ubi-rOCRs, no_overlap

# INPUT: TSSs bed file and rOCRs & ubi-rOCRs bed file
# OUTPUT: TSSs bed file with labels
# EXP: bash mark_TSS.sh hg38_v28_basic_TSS_filtered_uniq.bed
#            /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

tss_bed=$1
rOCRs_file=$2
ubi_rOCRs_file=$3
workDir=$4

###################
cd ${workDir}
name=${tss_bed%.bed}

# 1. overlap_with_ubi-rOCRs, overlap_with_not-ubi_active-rOCRs, no_overlap
intersectBed -a ${tss_bed} -b ${rOCRs_file} -wa -u > tmp.overlapping_rOCR.txt
intersectBed -a ${tss_bed} -b ${ubi_rOCRs_file} -wa -u > tmp.overlapping_ubi-rOCR.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0,"overlap_with_not-ubi_active-rOCRs"}else{print $0,"no_overlap"}}}' \
tmp.overlapping_rOCR.txt ${tss_bed} > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $1,$2,$3,$4,$5,$6,$7,"overlap_with_ubi-rOCRs"}else{print $0}}}' \
tmp.overlapping_ubi-rOCR.txt tmp.txt > tmp.label1.txt

# 2. overlap_with_ubi-rOCRs, in_gene_overlapping_ubi-rOCRs, no_overlap
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$7=1]}else{if(a[$4]){print $0,"overlap_with_ubi-rOCRs"}else if(a[$4]!=1 && b[$7]){print $0,"in_gene_overlapping_ubi-rOCRs"}else{print $0,"no_overlap"}}}' \
tmp.overlapping_ubi-rOCR.txt tmp.label1.txt > ${name}_labeled.bed

rm tmp.txt tmp.*.txt

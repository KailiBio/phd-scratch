#!/bin/bash

# -- Kaili
# This script is for labeling merged-TSSs, counting number of overlapping site.
#       overlap_with_ubi-rOCRs, overlap_with_not-ubi_active-rOCRs, no_overlap

# INPUT: merged-TSSs bed file and rOCRs & ubi-rOCRs bed file
# OUTPUT: merged-TSSs bed file with labels, with counting
#                   col8: number of TSSs overlap rOCRs
#                   col9: number of TSSs overlap ubi-rOCRs
#                   col10: percentage of TSSs overlap ubi-rOCRs/rOCRs
# EXP: bash mark_merged_TSS.sh merged_TSS.bed  /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

merged_TSS_bed=$1
rOCRs_file=$2
ubi_rOCRs_file=$3
workDir=$4

###################
cd ${workDir}
name=${merged_TSS_bed%.bed}

# rOCRs
intersectBed -a ${merged_TSS_bed} -b ${rOCRs_file} -wa -c > tmp.rOCRs_num.bed
# ubi-rOCRs
intersectBed -a ${merged_TSS_bed} -b ${ubi_rOCRs_file} -wa -c > tmp.ubi-rOCRs_num.bed
# mark
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$8}else{if($8==0){print $0,a[$4],0,"no_overlap"}else if((a[$4]/$8)>=0.5){print $0,a[$4],(a[$4]/$8),"overlap_with_ubi-rOCRs"}else{print $0,a[$4],(a[$4]/$8),"overlap_with_not-ubi_active-rOCRs"}}}' \
tmp.ubi-rOCRs_num.bed tmp.rOCRs_num.bed > ${name}_labeled.bed

rm tmp.rOCRs_num.bed tmp.ubi-rOCRs_num.bed

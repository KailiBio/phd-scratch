#!/bin/bash

# -- Kaili
# This script is


cd /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -b hg38_v28_basic_TSS_filtered_uniq_labeled.bed -wa -u > ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed
intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -b hg38_v28_basic_TSS_filtered_uniq_labeled.bed -wa -u > rOCRs_overlap_hg38_v28_basic_TSS.bed

mkdir non-ubi_active_rOCR_overlap_TSS
while read line
do
    id=`awk '{print $1}' <<< ${line}`
    #
    intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/${id}_OCR.bed -b hg38_v28_basic_TSS_filtered_uniq_labeled.bed -wa -u > ./non-ubi_active_rOCR_overlap_TSS/${id}_OCR_overlap_TSS.bed
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_matched_DNase_RNA_RAMPAGE_list.txt

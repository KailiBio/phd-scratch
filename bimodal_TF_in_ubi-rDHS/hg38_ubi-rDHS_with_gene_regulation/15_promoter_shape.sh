#!/bin/bash

# -- Kaili
# This script is for analyzing ubi-rOCRs & promoter shapes.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/"
promoterShapeResult="/data/zusers/zhangx/projects/rampage/0_rampage_peak/entropy/promoter_type/"

cd ${workDir}

# 1. GC content: ubi-rOCRs vs. rOCRs with TSSs overlapping
# get TSS overlapped rOCRs
intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -b \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered.uniq.bed -wa |sort -u | sort -k1,1 -k2,2n > \
TSS_overlapped_rOCRs.bed
# get fasta
bedtools getfasta -fi /home/fankaili/genome/hg38.fa -bed TSS_overlapped_rOCRs.bed > TSS_overlapped_rOCRs.fa
python ${scriptDir}count_GC_content.py /data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/TSS_overlapped_rOCRs.fa \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/TSS_overlapped_rOCRs_GCcontent.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{id=$1":"$2"-"$3;if(b[id]){print $0,a[id]}}}' \
TSS_overlapped_rOCRs_GCcontent.txt TSS_overlapped_rOCRs.bed > TSS_overlapped_rOCRs_GCcontent_ID.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0,"ubi-rOCRs"}else{print $0,"remaining_rOCRs"}}}' \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed TSS_overlapped_rOCRs_GCcontent_ID.txt \
> TSS_overlapped_rOCRs_GCcontent_ID_annotated.txt


# 2. promoter shapes
cut -f 1-3,17 ${promoterShapeResult}ENCSR855LXJ_rampage_entropy_promoter_type.txt > ss.txt
intersectBed -a TSS_overlapped_rOCRs_GCcontent_ID_annotated.txt -b ss.txt -wa -wb > ss2.txt
cut -f 4,6,10 ss2.txt | sort -u > ss3.txt

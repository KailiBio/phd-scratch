#!/bin/bash

# -- Kaili
# This script is for labeling genes as:
#       1. genes whose TSSs overlap ubi-rOCRs
#       2. genes whose TSSs overlap other active rOCRs,
#       3. genes whose TSSs not overlap with rOCRs

# INPUT: marked uniq TSS file
#               protein-coding gene file from genome folder (with path).
#               work directory
# OUTPUT: labled gene list. (all & protein-coding only)
# EXP: bash mark_gene.sh hg38_v28_basic_TSS_filtered_uniq_labeled.bed /home/fankaili/genome/
#           hg38_v28_basic_gene_protein_coding.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

mark_tss_file=$1
genome_path=$2
protein_coding_gene_file=$3
workDir=$4

###################
cd ${workDir}
name=${mark_tss_file%_TSS_filtered_uniq_labeled.bed}

# mark genes
cut -f 7 ${mark_tss_file} | sort -u | awk '{FS=OFS="\t"}{print $1,"genes_whose_TSSs_not_overlape_rOCRs"}' > tmp.gene.txt
awk '{FS=OFS="\t"}{if(NR==FNR){if($8=="overlap_with_ubi-rOCRs"){a[$7]=1}else
if($8=="overlap_with_not-ubi_active-rOCRs"){b[$7]=1}}else{if(a[$1]){print $1,"genes_whose_TSSs_overlap_ubi-rOCRs"}
else if(b[$1]){print $1,"genes_whose_TSSs_overlap_other_active_rOCRs"}else{print $1,$2}}}' ${mark_tss_file} tmp.gene.txt > ${name}_gene_labled.bed

# get protein-coding gene list
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $0}}}' ${genome_path}${protein_coding_gene_file} \
${name}_gene_labled.bed > ${name}_gene_labled_protein_coding.bed

rm tmp.gene.txt

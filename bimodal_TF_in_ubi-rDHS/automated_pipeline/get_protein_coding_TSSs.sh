#!/bin/bash

# -- Kaili
# This script is for getting TSSs of protein-coding genes.

# INPUT:
# OUTPUT:
# EXP: bash get_protein_coding_TSSs.sh hg38_v28_basic_TSS_filtered_uniq_labeled.bed
#           hg38_v28_basic_TSS_filtered_merged_labeled.bed
#           /home/fankaili/genome/ hg38_v28_basic_TSS_protein_coding.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

tss_bed_labled=$1
merged_tss_bed_labled=$2
genome_path=$3
protein_coding_tss=$4
workDir=$5

###################
cd ${workDir}
name1=${tss_bed_labled%.bed}
name2=${merged_tss_bed_labled%.bed}

# 1. for TSSs
intersectBed -a ${tss_bed_labled} -b ${genome_path}${protein_coding_tss} -wa -wb | \
awk '{FS=OFS="\t"}{if($7==$16){print $1,$2,$3,$4,$5,$6,$16,$8,$9}}' | sort -u > ${name1}_protein_coding.bed

# 2. for merged-TSSs
intersectBed -a ${merged_tss_bed_labled} -b ${genome_path}${protein_coding_tss} -wa -wb | \
awk '{FS=OFS="\t"}{if($7==$18){print $1,$2,$3,$4,$5,$6,$18,$8,$9,$10,$11}}' | sort -u > ${name2}_protein_coding.bed

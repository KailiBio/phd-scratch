#!/bin/bash

# -- Kaili
# This script is for pre-processing gene & TSS files.

# INPUT:
# OUTPUT:
# EXP: bash 0_preprocess_gene_TSS_files.sh hg38_v28_basic_TSS_filtered.bed hg38_v28_basic_TSS_protein_coding.bed \
#           hg38_v28_basic_gene_protein_coding.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

original_tss_file="hg38_v28_basic_TSS_filtered.bed"
protein_coding_tss_file="hg38_v28_basic_TSS_protein_coding.bed"
protein_coding_gene_file="hg38_v28_basic_gene_protein_coding.txt"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/"

original_tss_file=$1
protein_coding_tss_file=$2
protein_coding_gene_file=$3
workDir=$4

## other inPath
autoScriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/automated_pipeline/"
genome_path="/home/fankaili/genome/"
hg38_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
hg38_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed"

cd ${workDir}

########################
# TSSs
name=${original_tss_file%.bed}

# 1. get coordinated TSS a uniq ID
bash ${autoScriptDir}get_uniq_TSS.sh ${original_tss_file} ${genome_path} ${workDir}

## 2) mark TSS
### overlap_with_ubi-rOCRs, overlap_with_not-ubi_active-rOCRs, no_overlap
### overlap_with_ubi-rOCRs, in_gene_overlapping_ubi-rOCRs, no_overlap
#
bash ${autoScriptDir}mark_TSS.sh ${name}_uniq.bed ${workDir}

## 3) get merged-TSSs
#### this steps take longer.
python  ${autoScriptDir}merge_TSS_50bp.py ${name}_uniq.bed ${name}_merged.bed \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

## 4) label merged-TSSs, count how many TSSs overlap ubi-rOCRs/rOCRs in merged-TSSs
bash  ${autoScriptDir}mark_merged_TSS.sh ${name}_merged.bed ${hg38_rOCRs} ${hg38_ubi_rOCRs} ${workDir}

## 5) get protein-coding gene TSSs
bash  ${autoScriptDir}get_protein_coding_TSSs.sh ${name}_uniq_labeled.bed ${name}_merged_labeled.bed \
${genome_path} ${protein_coding_tss_file} ${workDir}

########################
# genes
# label genes: genes whose TSSs overlap ubi-rOCRs, genes whose TSSs overlap other active rOCRs,
#                       and genes whose TSSs not overlap with rOCRs.
# get protein-coding gene list
bash ${autoScriptDir}mark_gene.sh ${name}_uniq_labeled.bed ${genome_path} ${protein_coding_gene_file} ${workDir}

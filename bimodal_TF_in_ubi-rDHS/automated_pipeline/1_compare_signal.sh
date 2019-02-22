#!/bin/bash

# -- Kaili
# This script is for getting signal between ubi-rOCRs and other active rOCRs.

# INPUT: file_list is the expID and fileID for DNase, RNA-seq, RAMPAGE with sampleName. (7 column)
#              rOCRs and ubi-rOCRs bed file.
# OUTPUT:
# EXP: bash 1_compare_signal.sh /data/zusers/moorej3/moorej.ghpcc.project/PsychENCODE/Registry/V0-hg38/signal-output/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_gene_labled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/comprehensive_annotation/hg38_v28_comprehensive_TSS_filtered_with_uniqID.bed
#            /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue/
#          /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

# DNase_signal_path="/data/zusers/moorej3/moorej.ghpcc.project/PsychENCODE/Registry/V0-hg38/signal-output/"
# non_ubi_outPath="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/"
# gene_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_gene_labled.bed"
# tss_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed"
# tss_with_id_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/comprehensive_annotation/hg38_v28_comprehensive_TSS_filtered_with_uniqID.bed"
# exp_signal_path="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/"
# tss_signal_path="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue/"
# outPath="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/"

DNase_signal_path=$1
non_ubi_outPath=$2
gene_labeled_file=$3
tss_labeled_file=$4
tss_with_id_file=$5
exp_signal_path=$6
tss_signal_path=$7
outPath=$8

autoScriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/automated_pipeline/"
hg38_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
hg38_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed"

matched_file_list="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_matched_DNase_RNA_RAMPAGE_list.txt"

###################
# 1. DNase
bash ${autoScriptDir}get_DNase_signal_two_groups.sh ${DNase_signal_path} ${hg38_rOCRs} ${hg38_ubi_rOCRs} ${non_ubi_outPath} ${outPath}
#
echo "finish DNase part."

# 2. RNA-seq
# get all gene expression (tsv file, use expID as filename)
# python /data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/get_GRCh38_RNA_exp.py
# get RNA-seq for genes that overlapping ubi-rOCRs and non-ubi active rOCRs
bash ${autoScriptDir}get_RNAseq_signal_two_groups.sh ${gene_labeled_file} ${tss_labeled_file} ${exp_signal_path} ${non_ubi_outPath} ${outPath}
#
echo "finish RNA-seq part."

# 3. RAMPAGE
# get all TSS RAMPAGE signal
# python /data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/get_hg38_RAMPAGE_signal.py
# get RAMPAGE for labeled TSSs
bash ${autoScriptDir}get_RAMPAGE_signal_two_groups.sh ${tss_labeled_file} ${tss_with_id_file} ${tss_signal_path} ${non_ubi_outPath} ${outPath}
#
echo "finish RAMPAGE part."

####
# make fig1

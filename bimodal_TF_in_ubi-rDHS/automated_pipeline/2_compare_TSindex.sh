#!/bin/bash

# -- Kaili
# This script is for calculating tissue-specificity index for genes (RNA-seq) and TSSs (RAMPAGE).

# INPUT:
# OUTPUT:
# EXP: bash 2_compare_TSindex.sh
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_gene_labled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_with_uniqID.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

# gene_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_gene_labled.bed"
# tss_with_id_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_with_uniqID.bed"
# tss_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed"
# outPath="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/"

gene_labeled_file=$1
tss_with_id_file=$2
tss_labeled_file=$3
outPath=$4

gene_exp_matrix="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_gene_exp_matrix.txt"
tss_exp_matrix="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_TSS_exp_matrix_v28.txt"

autoScriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/automated_pipeline/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"

###################
# 1. tissue-specificity index for genes
# get interested gene matrix
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' ${gene_labeled_file} ${gene_exp_matrix} > ${outPath}tmp.gene_matrix.txt
# quantile normalization
Rscript ${dailyCodeDir}do_quantile_normalization.R ${outPath}tmp.gene_matrix.txt 1 ${outPath}tmp.gene_matrix_quantile.txt
# calculate tissue-specificity index
python ${autoScriptDir}calculate_TSindex.py ${outPath}tmp.gene_matrix_quantile.txt ${outPath}gene_exp_TSindex.txt
sed -i 's/-0.100000/NA/g' ${outPath}gene_exp_TSindex.txt
# label
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if(FNR>1){print $0,a[$1]}}}' ${gene_labeled_file} gene_exp_TSindex.txt > ${outPath}gene_exp_TSindex_labeled.txt
#
rm ${outPath}tmp.gene_matrix.txt ${outPath}tmp.gene_matrix_quantile.txt
echo "finish TSindex for RNA-seq."

# 2. tissue-specificity index for TSSs
# get interested TSS matrix
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$8;b[$4]=1}else{if(b[$1]){id=$1;$1=a[id];print $0}}}' ${tss_with_id_file} ${tss_exp_matrix} | sort -u > ${outPath}tmp.tss_matrix.txt
# quantile normalization
Rscript ${dailyCodeDir}do_quantile_normalization.R ${outPath}tmp.tss_matrix.txt 1 ${outPath}tmp.tss_matrix_quantile.txt
# calculate tissue-specificity index
python ${autoScriptDir}calculate_TSindex.py ${outPath}tmp.tss_matrix_quantile.txt ${outPath}tss_exp_TSindex.txt
sed -i 's/-0.100000/NA/g' ${outPath}tss_exp_TSindex.txt
# label
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{if(FNR>1){print $0,a[$1]}}}' ${tss_labeled_file} ${outPath}tss_exp_TSindex.txt > ${outPath}tss_exp_TSindex_labeled.txt
#
rm ${outPath}tmp.tss_matrix.txt ${outPath}tmp.tss_matrix_quantile.txt
echo "finish TSindex for RAMPAGE."

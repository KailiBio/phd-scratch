#!/bin/bash

# -- Kaili
# This script is for repeating figs only for protein-coding genes.

# INPUT:
# OUTPUT:
# EXP: bash 4_repeat_for_protein_coding.sh /data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_gene_labled_protein_coding.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled_protein_coding.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue/
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_with_uniqID.bed
#           /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/protein_coding/

# non_ubi_outPath="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/"
# gene_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_gene_labled_protein_coding.bed"
# tss_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled_protein_coding.bed"
# exp_signal_path="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/"
# tss_signal_path="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue/"
# tss_with_id_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_with_uniqID.bed"
# outPath="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/protein_coding/"

non_ubi_outPath=$1
gene_labeled_file=$2
tss_labeled_file=$3
exp_signal_path=$4
tss_signal_path=$5
tss_with_id_file=$6
outPath=$7

autoScriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/automated_pipeline/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"

hg38_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
hg38_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed"
gene_exp_matrix="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_gene_exp_matrix.txt"
tss_exp_matrix="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_TSS_exp_matrix_v28.txt"

matched_file_list="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_matched_DNase_RNA_RAMPAGE_list.txt"

###################
# 1. fig1
bash ${autoScriptDir}get_RNAseq_signal_two_groups.sh ${gene_labeled_file} ${tss_labeled_file} ${exp_signal_path} ${non_ubi_outPath} ${outPath}
#
bash ${autoScriptDir}get_RAMPAGE_signal_two_groups.sh ${tss_labeled_file} ${tss_with_id_file} ${tss_signal_path} ${non_ubi_outPath} ${outPath}
#
echo "finish figure1 data."

# 2. fig2
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
######  TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$8]){print $0}}}' ${tss_labeled_file} ${tss_with_id_file} | sort -u > ${outPath}protein_coding_tss_with_uniqID.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$8;b[$4]=1}else{if(b[$1]){id=$1;$1=a[id];print $0}}}' ${outPath}protein_coding_tss_with_uniqID.txt ${tss_exp_matrix} | sort -u > ${outPath}tmp.tss_matrix.txt
# quantile normalization
Rscript ${dailyCodeDir}do_quantile_normalization.R ${outPath}tmp.tss_matrix.txt 1 ${outPath}tmp.tss_matrix_quantile.txt
# calculate tissue-specificity index
python ${autoScriptDir}calculate_TSindex.py ${outPath}tmp.tss_matrix_quantile.txt ${outPath}tss_exp_TSindex.txt
sed -i 's/-0.100000/NA/g' ${outPath}tss_exp_TSindex.txt
# label
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{if(FNR>1){print $0,a[$1]}}}' ${tss_labeled_file} ${outPath}tss_exp_TSindex.txt > ${outPath}tss_exp_TSindex_labeled.txt
#
rm ${outPath}tmp.tss_matrix.txt ${outPath}tmp.tss_matrix_quantile.txt
echo "finish figure2 data."

# 3. fig3
cut -f 7 ${tss_labeled_file} | sort | uniq -u > ${outPath}gene_singular_TSS_list.txt
# get multiple TSSs file
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]!=1){print $0}}}' ${outPath}gene_singular_TSS_list.txt ${tss_labeled_file} | sort -k7,7 -k2,2n > ${outPath}TSSs_in_multiple_TSSs_genes.txt
# calculate distance
awk 'BEGIN{FS=OFS="\t";gene="";tss="";end="";distance=1000000}{if(gene==""){gene=$7;tss=$4;end=$3}else{if(gene==$7){if(distance>($2-end)){print gene,tss,$2-end}else{print gene,tss,distance};tss=$4;distance=$2-end;end=$3}else{print gene,tss,distance;gene=$7;tss=$4;end=$3;distance=1000000}}}END{print gene,tss,distance}' ${outPath}TSSs_in_multiple_TSSs_genes.txt > ${outPath}tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$3;b[$2]=1}else{if(b[$4]){print $0,a[$4]}}}' ${outPath}tmp.txt ${tss_labeled_file} | sort -u > ${outPath}nearest_TSS_distance_in_multiple_TSS_genes.txt
#
rm ${outPath}tmp.txt
echo "finish figure3 data."

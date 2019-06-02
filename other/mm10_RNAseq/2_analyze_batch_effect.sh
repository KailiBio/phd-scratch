#!/bin/bash

# -- Kaili
# This script is for analyzing batch effect for all 73 RNA-seq data.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}

# 1. get expression matrix
## 1) raw signal
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_exp_matrix.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$1]}}' mouse_RNA_fastq_list.txt mouse_RNA_fastq_list_merged.txt > mouse_RNA_fastq_list_merged_SampleInfo.txt
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    fileIDs=`awk '{print $2}' <<< ${line}`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    #
    echo ${sample}
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$1]}}' ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results mm10_exp_matrix.txt > tmp.txt
    sed "s/TPM/${sample}_${rep}/" tmp.txt > mm10_exp_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt


# awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_exp_matrix_rep1.txt
# #
# while read line
# do
#     expID=`awk '{print $1}' <<< ${line}`
#     fileIDs=`awk '{print $2}' <<< ${line}`
#     rep=`awk '{print $3}' <<< ${line}`
#     sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
#     #
#     if [ "$rep" -eq "1" ];then
#         echo ${sample}
#     fi
#     #
#     awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$1]}}' ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results mm10_exp_matrix_rep1.txt > tmp.txt
#     sed "s/TPM/${sample}_${rep}/" tmp.txt > mm10_exp_matrix_rep1.txt
# done < mouse_RNA_fastq_list_merged_SampleInfo.txt

## 2) log-transformed, z-score signal
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_exp_zscore_matrix.txt
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    fileIDs=`awk '{print $2}' <<< ${line}`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    #
    echo ${sample}
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6;b[$1]=1}else{if(FNR>1 && b[$1]){print $1,a[$1]}}}' ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results mm10_exp_zscore_matrix.txt > tmp.txt
    python /data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/zscore-normalization.py tmp.txt 2 > tmp2.txt
    awk -v sample="${sample}" -v rep="${rep}" 'BEGIN{FS=OFS="\t";print "gene_id",sample"_"rep}{print $1,$2}' tmp2.txt > tmp3.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$1]}}' tmp3.txt mm10_exp_zscore_matrix.txt > tmp4.txt
    mv tmp4.txt mm10_exp_zscore_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt


# 2. detect batch effect using raw TPM
analyze_batch_effect.R

#!/bin/bash

# -- Kaili
# This script is for getting final TPM, read counts * spikeIn matrix.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}


# 1. get TPM matrix
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_M18_TPM_matrix.txt
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
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$1]}}' ./rna_exp_Jul25/${expID}/${fileIDs}/${fileIDs}.genes.results mm10_M18_TPM_matrix.txt > tmp.txt
    sed "s/TPM/${sample}_${rep}/" tmp.txt > mm10_M18_TPM_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt
rm tmp.*txt

# 2. get read counts matrix
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_M18_readCounts_matrix.txt
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
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$1]}}' ./rna_exp_Jul25/${expID}/${fileIDs}/${fileIDs}.genes.results mm10_M18_readCounts_matrix.txt > tmp.txt
    sed "s/expected_count/${sample}_${rep}/" tmp.txt > mm10_M18_readCounts_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt
rm tmp.*txt

# 3. get spike-in TPM matrix
echo "gene_id" > mm10_M18_spikeIn_TPM_matrix.txt
grep "ERCC-" gencode.vM18.basic.annotation_withSpinkIn.gtf | cut -f 1 | sort -u >> mm10_M18_spikeIn_TPM_matrix.txt
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
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$1]}}' ./rna_exp_Jul25/${expID}/${fileIDs}/${fileIDs}.genes.results mm10_M18_spikeIn_TPM_matrix.txt > tmp.txt
    sed "s/TPM/${sample}_${rep}/" tmp.txt > mm10_M18_spikeIn_TPM_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt
#
rm tmp*

# 4. get spike-in read_counts matrix
echo "gene_id" > mm10_M18_spikeIn_readCounts_matrix.txt
grep "ERCC-" gencode.vM18.basic.annotation_withSpinkIn.gtf | cut -f 1 | sort -u >> mm10_M18_spikeIn_readCounts_matrix.txt
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
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$1]}}' ./rna_exp_Jul25/${expID}/${fileIDs}/${fileIDs}.genes.results mm10_M18_spikeIn_readCounts_matrix.txt > tmp.txt
    sed "s/expected_count/${sample}_${rep}/" tmp.txt > mm10_M18_spikeIn_readCounts_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt
#
rm tmp*

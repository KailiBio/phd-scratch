#!/bin/bash

# -- Kaili
# This script is for doing gene expression regression in results have CTCF states.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/gene_regression/"

other_file_dir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/"

cd ${workDir}

# 1. get state proportion
count_state_proportion(){
    stateBedDir=$1
    prefix=$2
    state_num=$3
    #
    while read line
    do
        echo $line
        ### get state in each window
        intersectBed -a ${other_file_dir}mm10_protein_coding_promoter_bins.bed -b \
        ${stateBedDir}${line}_state_sorted.bed -wo | sort -u | \
        awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2],$5,$6,$7,$8,$10,$9}' \
        | sort -k6,6n -k5,5 -k8,8n > ./state_proportion/mm10_gene_${line}_${prefix}_state.txt
        ### calculate state count in 20 windows
        for i in {1..20}
        do
            awk -v n="$i" -v state_num="$state_num" 'BEGIN{FS=OFS="\t";gene="";for(i=0;i<state_num;i++){a[i]=0}}{l=$3-$2;if($6==n){if(gene==""){gene=$5;a[$12]+=$11/l}else if(gene==$5){a[$12]+=$11/l}else{printf gene"\t"n;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n";gene=$5;for(i=0;i<state_num;i++){a[i]=0};a[$12]+=$11/l}}}END{printf gene"\t"n;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n"}' ./state_proportion/mm10_gene_${line}_${prefix}_state.txt  | sort -k1,1 > ./state_proportion/mm10_gene_${line}_${prefix}_state_count_window_${i}.txt
        done
    done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
}

stateBedDir="/data/zusers/fankaili/ideas/dhs_ctcf/state_bed/"
prefix="dhs-ctcf"
state_num=47
count_state_proportion ${stateBedDir} ${prefix} ${state_num}

# 2. do regression
## get matched expression
head -1 ${other_file_dir}mm10_RNA_protein-coding_tpm_matrix.txt > mm10_RNA_protein-coding_tpm_matrix_matched.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$0}else{if(a[$1]){print b[$1]}}}' ${other_file_dir}mm10_RNA_protein-coding_tpm_matrix.txt ./state_proportion/mm10_gene_kidney_15.5_dhs-ctcf_state_count_window_7.txt >> mm10_RNA_protein-coding_tpm_matrix_matched.txt

# Rscript do_within_gene_regression.R

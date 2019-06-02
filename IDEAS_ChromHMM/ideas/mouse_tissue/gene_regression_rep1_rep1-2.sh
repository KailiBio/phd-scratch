#!/bin/bash

# -- Kaili
# This script is for doing gene expression regression between rep1 and rep1&2.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/"

cd ${workDir}

# 1. get state_bed file.
stateDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_result/"
prefix="DHS_reps."
stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_state_bed/"
mkdir ${stateBedDir}
nohup bash ${scriptDir}make_each_sample_state_bed_IDEAS.sh ${stateDir} ${prefix} ${stateBedDir} 132 > ./nohup/nohup.make_each_sample_state_bed_reps.out 2>&1&

# 2. get state proportion
# This is for calculating state proportion for IDEAS results using replicates.
count_state_proportion_reps(){
    stateBedDir=$1
    prefix=$2
    state_num=$3
    #
    while read line
    do
        echo $line
        ### get state in each window
        intersectBed -a ./gene_expression/mm10_protein_coding_promoter_bins.bed -b ${stateBedDir}${line}.0_state_sorted.bed -wo | sort -u | \
        awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2],$5,$6,$7,$8,$10,$9}' \
        | sort -k6,6n -k5,5 -k8,8n > ./${prefix}_state_proportion/mm10_gene_${line}.0_${prefix}_state.txt
        #
        intersectBed -a ./gene_expression/mm10_protein_coding_promoter_bins.bed -b ${stateBedDir}${line}.1_state_sorted.bed -wo | sort -u | \
        awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2],$5,$6,$7,$8,$10,$9}' \
        | sort -k6,6n -k5,5 -k8,8n > ./${prefix}_state_proportion/mm10_gene_${line}.1_${prefix}_state.txt
        #
        cat ./${prefix}_state_proportion/mm10_gene_${line}.0_${prefix}_state.txt ./${prefix}_state_proportion/mm10_gene_${line}.1_${prefix}_state.txt | sort -k6,6n -k5,5 -k8,8n > ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state.txt
        ### calculate state count in 20 windows
        for i in {1..20}
        do
            awk -v n="$i" -v state_num="$state_num" 'BEGIN{FS=OFS="\t";gene="";for(i=0;i<state_num;i++){a[i]=0}}{l=2*($3-$2);if($6==n){if(gene==""){gene=$5;a[$12]+=$11/l}else if(gene==$5){a[$12]+=$11/l}else{printf gene"\t"n;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n";gene=$5;for(i=0;i<state_num;i++){a[i]=0};a[$12]+=$11/l}}}END{printf gene"\t"n;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n"}' ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state.txt  | sort -k1,1 > ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state_count_window_${i}.txt
        done
    done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
}

mkdir DHS_reps_state_proportion

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_state_bed/"
prefix="DHS_reps"
state_num=44
count_state_proportion_reps ${stateBedDir} ${prefix} ${state_num}

# 3. do regression, make figures
# Rscript do_gene_regression_reps_DHSbins.R

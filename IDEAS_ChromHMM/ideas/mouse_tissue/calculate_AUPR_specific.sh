#!/bin/bash

# -- Kaili
# This script is for calculating PRAU for given IDEAS result.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"

# prefix="ctcf_8sample_impute_11sample"
# parafile="/data/zusers/fankaili/ideas/dhs_ctcf2/ctcf_8sample_impute_11sample_result/ctcf_8sample_impute_11sample.para0"
# workDir="/data/zusers/fankaili/ideas/dhs_ctcf2/"

prefix=$1
parafile=$2
workDir=$3

cd ${workDir}
mkdir ./state_ranked_PR/${prefix}/

if [ -f ./state_ranked_PR/${prefix}/${prefix}_PRAU.txt ];then rm ./state_ranked_PR/${prefix}/${prefix}_PRAU.txt;fi

# 1. rank states
Rscript ${scriptDir}rank_states_by_CTCF_signal.R ${workDir}state_ranked_PR/${prefix}/ ${parafile} ${prefix}_state_CTCF_ranked.txt

# 2. make PR curve
for file in `ls ./state_bed_${prefix}/*_state_sorted.bed`
do
    file0=${file%_state_sorted.bed}
    sample=${file0#./state_bed_${prefix}/}
    echo $sample
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5;b[$4]=1}else{if(b[$4]){print $4,a[$4]}}}' ./state_bed_${prefix}/${sample}_state_sorted.bed ./state_ranked_PR/state_ranked_${sample}_specific_CTCFpeak_pos.txt > tmp.${prefix}_pos.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$1}else{print $1,a[$2]}}' ./state_ranked_PR/${prefix}/${prefix}_state_CTCF_ranked.txt tmp.${prefix}_pos.txt > ./state_ranked_PR/${prefix}/${sample}_${prefix}_peak_PR_specific_pos.txt
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5;b[$4]=1}else{if(b[$4]){print $4,a[$4]}}}' ./state_bed_${prefix}/${sample}_state_sorted.bed ./state_ranked_PR/${sample}_CTCFpeak_neg.txt > tmp.${prefix}_neg.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$1}else{print $1,a[$2]}}' ./state_ranked_PR/${prefix}/${prefix}_state_CTCF_ranked.txt tmp.${prefix}_neg.txt > ./state_ranked_PR/${prefix}/${sample}_${prefix}_peak_PR_neg.txt
    #
    prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}state_ranked_PR/${prefix}/ ${sample}_${prefix}_peak_PR_specific_pos.txt ${sample}_${prefix}_peak_PR_neg.txt`
    tmp=`awk '{print $2}' <<< $prau`
    echo -e ${sample}"\t"$tmp >> ./state_ranked_PR/${prefix}/${prefix}_specific_PRAU.txt
done

rm tmp.${prefix}_pos.txt tmp.${prefix}_neg.txt

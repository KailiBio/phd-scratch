#!/bin/bash

# -- Kaili
# This script is for making PR curve by ranking CTCF states.

cd /data/zusers/fankaili/ideas/dhs_ctcf/

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"

prefix=$1
parafile=$2

# prefix="9impute11"
# parafile="/data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_impute_11sample_result/ctcf_9sample_impute_11sample.para0"

mkdir ./state_ranked_PR2/${prefix}/

# 1. rank states
# Rscript ${scriptDir}rank_states_by_CTCF_signal.R /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/${prefix}/ ${parafile} ${prefix}_state_CTCF_ranked.txt

# 2. make PR curve
for file in `ls ./state_bed_${prefix}/*_state_sorted.bed`
do
    file0=${file%_state_sorted.bed}
    sample=${file0#./state_bed_${prefix}/}
    echo $sample
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{print $4,a[$4]}}' ./state_bed_${prefix}/${sample}_state_sorted.bed ./state_ranked_PR2/${sample}_CTCFpeak_pos.txt > tmp.${prefix}_pos.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$1}else{print $1,a[$2]}}' ./state_ranked_PR/${prefix}/${prefix}_state_CTCF_ranked.txt tmp.${prefix}_pos.txt > ./state_ranked_PR2/${prefix}/${sample}_${prefix}_peak_PR_pos.txt
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{print $4,a[$4]}}' ./state_bed_${prefix}/${sample}_state_sorted.bed ./state_ranked_PR2/${sample}_CTCFpeak_neg.txt > tmp.${prefix}_neg.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$1}else{print $1,a[$2]}}' ./state_ranked_PR/${prefix}/${prefix}_state_CTCF_ranked.txt tmp.${prefix}_neg.txt > ./state_ranked_PR2/${prefix}/${sample}_${prefix}_peak_PR_neg.txt
    #
    Rscript ${scriptDir}make_PRcurve_peak.R /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR2/${prefix}/ ${sample}_${prefix}_peak_PR_pos.txt ${sample}_${prefix}_peak_PR_neg.txt ${sample}_${prefix}_PRcurve_peak_2
done

rm tmp.${prefix}_pos.txt tmp.${prefix}_neg.txt

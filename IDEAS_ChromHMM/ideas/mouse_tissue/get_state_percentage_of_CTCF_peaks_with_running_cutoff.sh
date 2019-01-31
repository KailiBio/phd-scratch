#!/bin/bash

# -- Kaili
# This script is for calculate the percentage of states with CTCF peaks with running cut-off.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

############

for i in {1..42}
do
    echo "state"${i}
    if [ -f ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_peak_count.txt ]
    then
        rm ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_peak_count.txt
    fi
    #
    total=`wc -l ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_region.bed | awk '{print $1}'`
    for j in {0..100}
    do
        echo ${j}
        awk -v signal="$j" '{FS=OFS="\t"}{if(NR==FNR && $5>signal){a[$1]=1}else{if(a[$4]){print $0}}}' \
        liver_14.5_day_CTCF_peak_signal.txt liver_14.5_day_CTCF_peak_sorted.bed > tmp.liver14.5_peak.bed
        num=`intersectBed -a ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_region.bed \
        -b tmp.liver14.5_peak.bed -f 1 -F 0.5 -e -wa | sort -u | wc -l | awk '{print $1}'`
        echo -e ${j}"\t"${total}"\t"${num} >> ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_peak_count.txt
    done
done

# Rscript make_percentage_of_state_with_CTCF_peak_cutoff.R

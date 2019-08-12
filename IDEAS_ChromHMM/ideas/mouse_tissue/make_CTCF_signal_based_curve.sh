#!/bin/bash

# -- Kaili
# This script is for making CTCF signal based curve.

cd /data/zusers/fankaili/ideas/dhs_ctcf/

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"

prefix=$1
mkdir ./CTCF_signal_based_curve/${prefix}/

for file in `ls ./state_bed_${prefix}/*_state_sorted.bed`
do
    file0=${file%_state_sorted.bed}
    sample=${file0#./state_bed_${prefix}/}
    echo $sample
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{if(a[$1]==12 || a[$1]==24 || a[$1]==26 || a[$1]==27 || a[$1]==30 || a[$1]==34 || a[$1]==37 || a[$1]==38 || a[$1]==39 || a[$1]==41 || a[$1]==45 || a[$1]==46){print $1,$5}}}' ./state_bed_${prefix}/${sample}_state_sorted.bed ${signalDir}${sample}_CTCF_dhs.tab > ./CTCF_signal_based_curve/${prefix}/${prefix}_${sample}_CTCFstate_pos.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{if(a[$1]!=12 && a[$1]!=24 && a[$1]!=26 && a[$1]!=27 && a[$1]!=30 && a[$1]!=34 && a[$1]!=37 && a[$1]!=38 && a[$1]!=39 && a[$1]!=41 && a[$1]!=45 && a[$1]!=46){print $1,$5}}}' ./state_bed_${prefix}/${sample}_state_sorted.bed ${signalDir}${sample}_CTCF_dhs.tab > ./CTCF_signal_based_curve/${prefix}/${prefix}_${sample}_CTCFstate_neg.txt
    #
    Rscript ${scriptDir}make_CTCF_signal_based_curve.R /data/zusers/fankaili/ideas/dhs_ctcf/CTCF_signal_based_curve/${prefix}/ ${prefix}_${sample}_CTCFstate_pos.txt ${prefix}_${sample}_CTCFstate_neg.txt ${prefix}_${sample}_curve
done

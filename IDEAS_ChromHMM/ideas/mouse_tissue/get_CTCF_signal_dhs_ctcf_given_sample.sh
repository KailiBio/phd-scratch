#!bin/bash

# -- Kaili
# This script is for getting CTCF signal of dhs-CTCF-imputed bins in given state in given sample.

# INPUT: state number
# OUTPUT: CTCF signal of all bins in that given state
# EXP: bash get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ 0 heart_0

workDir=$1
state=$2
sample=$3

#########
if [ -f  ${workDir}state_CTCF_signal/state_${state}_${sample}_CTCF_signal.txt ]; then
    rm ${workDir}state_CTCF_signal/state_${state}_${sample}_CTCF_signal.txt
fi
#
awk -v state="$state" -v sample="$sample" '{FS=OFS="\t"}{if(NR==FNR){if($5==state){a[$4]=1}}else{if(a[$1]){print state,sample,$1,$5}}}' ${workDir}state_bed/${sample}_state_sorted.bed /data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/${sample}_CTCF_dhs.tab >> ${workDir}state_CTCF_signal/state_${state}_${sample}_CTCF_signal.txt

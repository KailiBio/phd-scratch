#!bin/bash

# -- Kaili
# This script is for getting CTCF signal of all bins in given state.

# INPUT: state number
# OUTPUT: CTCF signal of all bins in that given state
# EXP: bash get_CTCF_signal_for_given_state.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ 0

workDir=$1
state=$2

#########
if [ -f  ${workDir}state_CTCF_signal/state_${state}_CTCF_signal.txt ]; then
    rm ${workDir}state_CTCF_signal/state_${state}_CTCF_signal.txt
fi
#
for i in {1..11}
do
    sample=`head -1 /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.chr1.state | awk -v n="$i" '{print $(n+4)}'`
    awk -v state="$state" -v sample="$sample" '{FS=OFS="\t"}{if(NR==FNR){if($5==state){a[$4]=1}}else{if(a[$1]){print state,sample,$1,$5}}}' \
    ${workDir}state_bed/${sample}_state_sorted.bed \
    /data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/${sample}_CTCF_normal.tab \
    >> ${workDir}state_CTCF_signal/state_${state}_CTCF_signal.txt
done

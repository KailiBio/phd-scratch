#!/bin/bash

# -- Kaili
# This script is for making IDEAS state output into bed file for each sample.

# INPUT:state file
#              state bed file folder.
#              number of samples.
# OUTPUT: state bed file for each sample. (4 column)
# EXP:bash make_each_sample_state_bed_tmp.sh /data/zusers/fankaili/ideas/dhs_ctcf/9-impute-11_only_100_result/9-impute-11_only_100.tmp.1.state /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed_9impute11-100/ 2

stateFile=$1
stateBedDir=$2
num=$3

#################
for ((i=1; i<=num; i++))
do
    # get sample name
    sample=`head -1 ${stateFile} | awk -v n="$i" '{print $(n+4)}'`
    #
    awk -v n="$i" '{FS=" ";OFS="\t"}{if(NR>1){print $2,$3,$4,$1,$(n+4)}}' ${stateFile} > ${stateBedDir}${sample}_state.bed
    sort -k1,1 -k2,2n ${stateBedDir}${sample}_state.bed > ${stateBedDir}${sample}_state_sorted.bed
done

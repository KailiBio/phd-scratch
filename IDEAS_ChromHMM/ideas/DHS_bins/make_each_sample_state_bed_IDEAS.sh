#!/bin/bash

# -- Kaili
# This script is for making IDEAS state output into bed file for each sample.

# INPUT:result dir of IDEAS
#              prefix name of IDEAS.
#              state bed file folder.
#              number of samples.
# OUTPUT: state bed file for each sample. (4 column)
# EXP:bash make_each_sample_state_bed_IDEAS.sh "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp_result/"
#           "DHS_v3_100-400bp." "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/" 66

stateDir=$1
prefix=$2
stateBedDir=$3
num=$4

#################
for ((i=1; i<=num; i++))
do
    # get sample name
    sample=`head -1 ${stateDir}${prefix}chr1.state | awk -v n="$i" '{print $(n+4)}'`
    if [ -f ${stateBedDir}${sample}_state.bed ]; then rm ${stateBedDir}${sample}_state.bed; fi
    # for each chromosome
    ## get chromosome ID
    for j in {1..21}
    do
        if [ ${j} -eq 20 ]; then c="X"; elif [ ${j} -eq 21 ]; then c="Y"; else c=${i}; fi
        #
        awk -v n="$i" '{FS=" ";OFS="\t"}{if(NR>1){print $2,$3,$4,$1,$(n+4)}}' ${stateDir}${prefix}chr${c}.state >> \
        ${stateBedDir}${sample}_state.bed
    done
    # sort file, for intersectBed
    sort -k1,1 -k2,2n ${stateBedDir}${sample}_state.bed > ${stateBedDir}${sample}_state_sorted.bed
done

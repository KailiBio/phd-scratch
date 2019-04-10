#!/bin/bash

# -- Kaili
# This script is for pooling all states bins together.

cd /data/zusers/fankaili/ideas/dhs_ctcf/
for i in {12..46}
do
    echo ${i}
    if [ -f ./pool-bins/dhs_ctcf_9-66_state_${i}_seperate_bins.bed ];then rm ./pool-bins/dhs_ctcf_9-66_state_${i}_seperate_bins.bed; fi
    while read sample
    do
        echo ${sample}
        awk -v i="$i" '{FS=OFS="\t"}{if($5==i){print $0}}' ./state_bed/${sample}_state.bed >> ./pool-bins/dhs_ctcf_9-66_state_${i}_seperate_bins.bed
    done < all_66_sample.txt
    sort -u ./pool-bins/dhs_ctcf_9-66_state_${i}_seperate_bins.bed | sort -k1,1 -k2,2n > ./pool-bins/dhs_ctcf_9-66_state_${i}_seperate_bins_sorted.bed
    bedtools merge -i ./pool-bins/dhs_ctcf_9-66_state_${i}_seperate_bins_sorted.bed > ./pool-bins/dhs_ctcf_9-66_state_${i}_pool-bins.bed
    awk -v i="$i" '{FS=OFS="\t"}{print $0,"state"i"_"NR}' ./pool-bins/dhs_ctcf_9-66_state_${i}_pool-bins.bed > ./pool-bins/dhs_ctcf_9-66_state_${i}_pool-bins_withID.bed
done

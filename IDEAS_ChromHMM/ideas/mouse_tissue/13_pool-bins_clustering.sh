#!/bin/bash

# -- Kaili
# This script is for clustering pool bins signal.
# use results from: /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_impute_result/

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# 1. pool all state bins in all samples
nohup bash ${scriptDir}pool_bins_each_state.sh > ./nohup/nohup.pool_bins_each_state.out 2>&1&

# 2. get signal of pool-bins
for state in {0..46}
do
    echo ${state}
    nohup bash ${scriptDir}get_pool-bins_signal.sh ${state} > ./nohup/nohup.get_pool-bins_signal_${state}.out 2>&1&
done

# 3. clustering
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3
do
    nohup bash ${scriptDir}make_cluster_for_pool-bins.sh ${mark} > ./nohup/nohup.make_cluster_for_pool-bins_${mark}.out 2>&1&
done


###############
# Apr 23
# test Arjan's 260,941 enhancer regions
awk '{FS=OFS="\t"}{print $1,$2,$3}' ./ss/Enh-embryonic_facial_prominence_12.5_H3K27ac | sort -k1,1 -k2,2n > ./ss/enhancer_list.bed
## 1) percentage of bins belong to enhancer_list in each state
if [ -f bins_percentage_in_enhancer_list.txt ];then rm bins_percentage_in_enhancer_list.txt; fi
for state in {0..46}
do
    echo ${state}
    #
    num=`wc -l ./pool-bins/dhs_ctcf_9-66_state_${state}_seperate_bins_sorted.bed | awk '{print $1}'`
    intersectBed -a ./pool-bins/dhs_ctcf_9-66_state_${state}_seperate_bins_sorted.bed -b ./ss/enhancer_list.bed -wa -wb | cut -f 4 | sort -u | wc -l | awk -v num="$num" -v state="$state" '{print state,$1,num,$1/num}' >>  bins_percentage_in_enhancer_list.txt
done
# Rscript make_enhancer_percentage_barplot.R
## 2) merge state 18,22,35,41,44
cat ./pool-bins/dhs_ctcf_9-66_state_18_pool-bins.bed ./pool-bins/dhs_ctcf_9-66_state_22_pool-bins.bed ./pool-bins/dhs_ctcf_9-66_state_35_pool-bins.bed ./pool-bins/dhs_ctcf_9-66_state_41_pool-bins.bed ./pool-bins/dhs_ctcf_9-66_state_44_pool-bins.bed | sort -k1,1 -k2,2n > tmp.bed
bedtools merge -i tmp.bed | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $0,"E"NR}' > ./pool-bins/dhs_ctcf_9-66_state_E_pool-bins_withID.bed
## 3) calculate signal for pool bins
nohup bash ${scriptDir}get_pool-bins_signal.sh "E" > ./nohup/nohup.get_pool-bins_signal_${state}.out 2>&1&
# 50781
## 4) do cluster
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3
do
    Rscript ${scriptDir}do_hclust_for_pool-bins.R "E" ${mark}
done
#
cd /data/zusers/fankaili/ideas/dhs_ctcf/pool-bins/hclust/
pdfjam --outfile stateE_hlcust.pdf stateE_H3K27ac_hlcust.pdf stateE_H3K4me3_hlcust.pdf stateE_H3K27me3_hlcust.pdf stateE_H3K36me3_hlcust.pdf stateE_H3K4me1_hlcust.pdf stateE_H3K4me2_hlcust.pdf stateE_H3K9ac_hlcust.pdf stateE_H3K9me3_hlcust.pdf

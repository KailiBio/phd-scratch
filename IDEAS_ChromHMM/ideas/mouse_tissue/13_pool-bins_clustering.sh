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

for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3
do
    pdfjam --outfile ${mark}_hlcust.pdf state0_${mark}_hlcust.pdf state1_${mark}_hlcust.pdf state2_${mark}_hlcust.pdf state3_${mark}_hlcust.pdf state4_${mark}_hlcust.pdf state5_${mark}_hlcust.pdf state6_${mark}_hlcust.pdf state7_${mark}_hlcust.pdf state8_${mark}_hlcust.pdf state9_${mark}_hlcust.pdf state10_${mark}_hlcust.pdf state11_${mark}_hlcust.pdf state12_${mark}_hlcust.pdf state13_${mark}_hlcust.pdf state14_${mark}_hlcust.pdf state15_${mark}_hlcust.pdf state16_${mark}_hlcust.pdf state17_${mark}_hlcust.pdf state18_${mark}_hlcust.pdf state19_${mark}_hlcust.pdf state20_${mark}_hlcust.pdf state21_${mark}_hlcust.pdf state22_${mark}_hlcust.pdf state23_${mark}_hlcust.pdf state24_${mark}_hlcust.pdf state25_${mark}_hlcust.pdf state26_${mark}_hlcust.pdf state27_${mark}_hlcust.pdf state28_${mark}_hlcust.pdf state29_${mark}_hlcust.pdf state30_${mark}_hlcust.pdf state31_${mark}_hlcust.pdf state32_${mark}_hlcust.pdf state33_${mark}_hlcust.pdf state34_${mark}_hlcust.pdf state35_${mark}_hlcust.pdf state36_${mark}_hlcust.pdf state37_${mark}_hlcust.pdf state38_${mark}_hlcust.pdf state39_${mark}_hlcust.pdf state40_${mark}_hlcust.pdf state41_${mark}_hlcust.pdf state42_${mark}_hlcust.pdf state43_${mark}_hlcust.pdf state44_${mark}_hlcust.pdf state45_${mark}_hlcust.pdf state46_${mark}_hlcust.pdf
done

#!/bin/bash

# -- Kaili
# This script is for getting pool bins cluster for DHS bins.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/"

cd ${workDir}

# 1. pool all state bins in all samples
nohup bash ${scriptDir}pool_bins_each_state.sh > ./nohup/nohup.pool_bins_each_state.out 2>&1&
#7439

# 2. get signal of pool-bins
for state in {26..43}
do
    echo ${state}
    nohup bash ${scriptDir}get_pool-bins_signal.sh ${state} > ./nohup/nohup.get_pool-bins_signal_${state}.out 2>&1&
done
# 0-10: z001 41810-41820
# 11-25: z010 58980-58994
# 26-43: z008 21295-21312


# 3. clustering
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3
do
    nohup bash ${scriptDir}make_cluster_for_pool-bins.sh ${mark} > ./nohup/nohup.make_cluster_for_pool-bins_${mark}.out 2>&1&
done

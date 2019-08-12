#!/bin/bash

# -- Kaili
# This script is for validating imputation-only result.


workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"

cd ${workDir}

# 1. 100X
prefix="9impute11-100"
mkdir state_bed_${prefix}
bash ${scriptDir}make_each_sample_state_bed_tmp.sh /data/zusers/fankaili/ideas/dhs_ctcf/9-impute-11_only_100_result/9-impute-11_only_100.tmp.1.state /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_${prefix}/ 2
#
bash ${scriptDir}make_PR_curve.sh ${prefix}

# 2. 10k
prefix="9impute11-10k"
mkdir state_bed_${prefix}
bash ${scriptDir}make_each_sample_state_bed_tmp.sh /data/zusers/fankaili/ideas/dhs_ctcf/9-impute-11_only_10k_result/9-impute-11_only_10k.tmp.1.state /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_${prefix}/ 2
#
bash ${scriptDir}make_PR_curve.sh ${prefix}

# 3. 1m
prefix="9impute11-1m"
mkdir state_bed_${prefix}
bash ${scriptDir}make_each_sample_state_bed_tmp.sh /data/zusers/fankaili/ideas/dhs_ctcf/9-impute-11_only_1m_result/9-impute-11_only_1m.tmp.1.state /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_${prefix}/ 2
#
bash ${scriptDir}make_PR_curve.sh ${prefix}


# 4. 100X + one-cluster
prefix="9impute11-100-1"
mkdir state_bed_${prefix}
bash ${scriptDir}make_each_sample_state_bed_tmp.sh /data/zusers/fankaili/ideas/dhs_ctcf/9-impute-11_only_100_result/9-impute-11_only_100.tmp.2.state /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_${prefix}/ 2
#
bash ${scriptDir}make_PR_curve.sh ${prefix}

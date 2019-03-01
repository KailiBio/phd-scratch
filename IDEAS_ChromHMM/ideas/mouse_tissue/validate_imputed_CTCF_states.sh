#!/bin/bash

# -- Kaili
# This script is for validating the inputed CTCF states using liver&lung e14.5 data.
# 1. assign state types
# 2. state conservation
# 3. make track hub
# 4. state CTCF signal


scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# 1. assign state types
Rscript ${scriptDir}cluster_state_IDEAS.R

# 2. state conservation
head -1 ./ctcf_9sample_impute_result/ctcf_9sample_impute.chr1.state > ctcf_9sample_impute.state
for i in {1..19} X Y
do
    echo ${i}
    awk '{if(NR>1){print $0}}' ./ctcf_9sample_impute_result/ctcf_9sample_impute.chr${i}.state >> ctcf_9sample_impute.state
done
#
# state_conservation_dhs_ctcf.sh
# z018 8288840


# 4. state CTCF signal
## 1) get state bed file in each sample
mkdir state_bed
#
stateDir="/data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_impute_result/"
prefix="ctcf_9sample_impute."
stateBedDir="/data/zusers/fankaili/ideas/dhs_ctcf/state_bed/"
nohup bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${stateDir} ${prefix} ${stateBedDir} 66 > /data/zusers/fankaili/ideas/dhs_ctcf/nohup/nohup.make_each_sample_state_bed.out 2>&1&
# 35115
## 2) get CTCF singal in each state for all samples
# mkdir state_CTCF_signal
# #
# for ((state=0; state<=46; state++))
# do
#     echo ${state}
#     bash ${scriptDir}get_CTCF_signal_dhs_ctcf.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state}
# done

## 3) get CTCF singal in each state for liver14.5, lung14.5, liver0, lung0, forebrain_0
# liver 14.5
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "liver_14.5"
done
# liver 0
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "liver_0"
done
# lung 14.5
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "lung_14.5"
done
# lung 0
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "lung_0"
done
# forebrain 0
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "forebrain_0"
done
# midbrain 0
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "midbrain_0"
done
# hindbrain 0
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "hindbrain_0"
done
# heart 0
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "heart_0"
done
# kidney 0
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "kidney_0"
done
# stomach 0
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "stomach_0"
done
# intestine 0
for ((state=0; state<=46; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_dhs_ctcf_given_sample.sh /data/zusers/fankaili/ideas/dhs_ctcf/ ${state} "intestine_0"
done
## get mean
if [ -f state_CTCF_signal_mean.txt ];then rm state_CTCF_signal_mean.txt;fi
#
for file in `ls /data/zusers/fankaili/ideas/dhs_ctcf/state_CTCF_signal/`
do
    echo $file
    #
    awk 'BEGIN{FS=OFS="\t";state=0;sample="";sum=0}{state=$1;sample=$2;sum+=$4}END{print state,sample,sum/NR}' ./state_CTCF_signal/${file} >> state_CTCF_signal_mean.txt
done

# Rscript make_barplot_mean_CTCF_signal_dhs-ctcf.R

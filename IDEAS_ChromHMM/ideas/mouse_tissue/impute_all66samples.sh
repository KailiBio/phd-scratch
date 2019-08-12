#!/bin/bash

# -- Kaili
# This script is for imputing CTCF states in all 66 samples.

workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"

cd ${workDir}

# 1. run 9impute66
cp 9impute11_average_motif_model.sh 9impute66_average_motif_model.sh
vim 9impute66_average_motif_model.sh
cp 9impute11_average_motif_model.parafile 9impute66_average_motif_model.parafile
vim 9impute66_average_motif_model.parafile
cp ctcf_9sample_impute.input 9impute66_average_motif_model.input
#
while read sample
do
    echo ${sample}" motif /data/zusers/fankaili/ideas/dhs_ctcf/CTCFmotif_signal_OCR.txt" >> 9impute66_average_motif_model.input
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/dhs_ctcf/9sample_averageCTCFsignal.txt" >> 9impute66_average_motif_model.input
done < all_66_sample.txt
#
nohup bash 9impute66_average_motif_model.sh > ./nohup/nohup.9impute66_average_motif_model.out 2>&1&
# z003 56240


# 2. only run IDEAS on 66 samples
cp 9impute66_average_motif_model.sh 66sample_average_motif_model.sh
vim 66sample_average_motif_model.sh
cp 9impute66_average_motif_model.parafile 66sample_average_motif_model.parafile
vim 66sample_average_motif_model.parafile
cp 9impute66_average_motif_model.input 66sample_average_motif_model.input
#
nohup bash 66sample_average_motif_model.sh > ./nohup/nohup.66sample_average_motif_model.out 2>&1&
# z010 30389


# 3. calculate PRAU for this two runs
## 1) 66sample_average_motif_model
mkdir state_bed_66sample_average_motif_model
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/dhs_ctcf/66sample_average_motif_model_result/ 66sample_average_motif_model. /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_66sample_average_motif_model/ 66
#
prefix="66sample_average_motif_model"
# bash ${scriptDir}make_CTCF_signal_based_curve.sh ${prefix}
#
bash ${scriptDir}make_PRcurve_peak_state-ranking3.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/66sample_average_motif_model_result/66sample_average_motif_model.para0
## 2) 66sample_average_motif_model
mkdir state_bed_9impute66_average_motif_model
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/dhs_ctcf/9impute66_average_motif_model_result/ 9impute66_average_motif_model. /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_9impute66_average_motif_model/ 66
#
prefix="9impute66_average_motif_model"
# bash ${scriptDir}make_CTCF_signal_based_curve.sh ${prefix}
#
bash ${scriptDir}make_PRcurve_peak_state-ranking3.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9impute66_average_motif_model_result/9impute66_average_motif_model.para0

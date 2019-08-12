#!/bin/bash

# -- Kaili
# This script is for using new ways to impuater CTCF states for new tissues.
## 8 tissues impute 3 brain tissues

# 1. 10 marks
# 2. 10 marks + CTCF motif
# 3. 10 marks + average CTCF signal
# 4. 10 makrs + CTCG motif + average CTCF signal
# 5. calculate PRAU, make barplot

workDir="/data/zusers/fankaili/ideas/dhs_ctcf2/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"

cd ${workDir}

# 0.
cp /data/zusers/fankaili/ideas/dhs_ctcf/all_ctcf_sample.txt 8samples.txt
cp /data/zusers/fankaili/ideas/dhs_ctcf/mm10_OCR-center_bins_v3_signal_based_space.bed ./
cp /data/zusers/fankaili/ideas/dhs_ctcf/mm10_OCR-center_bins_v3_signal_based_space.bed.inv ./
#
cp -r /data/zusers/fankaili/ideas/dhs_ctcf/bin ./
cp -r /data/zusers/fankaili/ideas/dhs_ctcf/data ./
## get 8 sample average signal
awk '{FS=OFS="\t"}{print $4,0}' /data/zusers/fankaili/ideas/dhs_ctcf/mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed > 8sample_sumCTCFsignal.txt
#
while read sample
do
    echo $sample
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $1,$2+a[$1]}}' ${signalDir}${sample}_CTCF_dhs.tab 8sample_sumCTCFsignal.txt > tmp.txt
    mv tmp.txt 8sample_sumCTCFsignal.txt
done < 8samples.txt
#
awk '{FS=OFS="\t"}{print $2/8}' 8sample_sumCTCFsignal.txt > 8sample_averageCTCFsignal.txt


# 1. 10 marks


# 2. 10 marks + CTCF motif
## 1) train model
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model.sh model_8sample_10marks_CTCFmotif.sh
vim model_8sample_10marks_CTCFmotif.sh
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model.parafile model_8sample_10marks_CTCFmotif.parafile
vim model_8sample_10marks_CTCFmotif.parafile
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model.input model_8sample_10marks_CTCFmotif.input
#
while read sample
do
    echo ${sample}" motif /data/zusers/fankaili/ideas/dhs_ctcf/CTCFmotif_signal_OCR.txt" >> model_8sample_10marks_CTCFmotif.input
done < 8samples.txt
#
nohup bash model_8sample_10marks_CTCFmotif.sh > ./nohup/nohup.model_8sample_10marks_CTCFmotif.out 2>&1&
# z003 25972
## 2) imputation
awk '{FS=OFS}{if(NR==1){print $0}else{printf 100*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 100*$i};printf "\n"}}' /data/zusers/fankaili/ideas/dhs_ctcf2/model_8sample_10marks_CTCFmotif_result/model_8sample_10marks_CTCFmotif.para0 > /data/zusers/fankaili/ideas/dhs_ctcf2/model_8sample_10marks_CTCFmotif_result/model_8sample_10marks_CTCFmotif_100fold.para0
#
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_impute_11sample.sh 8impute11_10marks_CTCFmotif.sh
vim 8impute11_10marks_CTCFmotif.sh
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_impute_11sample.parafile 8impute11_10marks_CTCFmotif.parafile
vim 8impute11_10marks_CTCFmotif.parafile
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_impute_11sample.input 8impute11_10marks_CTCFmotif.input
#
while read sample
do
    echo ${sample}" motif /data/zusers/fankaili/ideas/dhs_ctcf/CTCFmotif_signal_OCR.txt" >> 8impute11_10marks_CTCFmotif.input
done < /data/zusers/fankaili/ideas/dhs_ctcf/all_ctcf_sample.txt
#
nohup bash 8impute11_10marks_CTCFmotif.sh > ./nohup/nohup.8impute11_10marks_CTCFmotif.out 2>&1&
# z008 15886

# 3. 10 marks + average CTCF signal
## 1) train model
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model.sh model_8sample_10marks_CTCFaverage.sh
vim model_8sample_10marks_CTCFaverage.sh
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model.parafile model_8sample_10marks_CTCFaverage.parafile
vim model_8sample_10marks_CTCFaverage.parafile
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model.input model_8sample_10marks_CTCFaverage.input
while read sample
do
    echo ${sample}
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/dhs_ctcf2/8sample_averageCTCFsignal.txt" >> model_8sample_10marks_CTCFaverage.input
done < 8samples.txt
#
nohup bash model_8sample_10marks_CTCFaverage.sh > ./nohup/nohup.model_8sample_10marks_CTCFaverage.out 2>&1&
# z003 26316
## 2) imputation
awk '{FS=OFS}{if(NR==1){print $0}else{printf 100*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 100*$i};printf "\n"}}' /data/zusers/fankaili/ideas/dhs_ctcf2/model_8sample_10marks_CTCFaverage_result/model_8sample_10marks_CTCFaverage.para0 > /data/zusers/fankaili/ideas/dhs_ctcf2/model_8sample_10marks_CTCFaverage_result/model_8sample_10marks_CTCFaverage_100fold.para0
#
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_impute_11sample.sh 8impute11_10marks_CTCFaverage.sh
vim 8impute11_10marks_CTCFaverage.sh
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_impute_11sample.parafile 8impute11_10marks_CTCFaverage.parafile
vim 8impute11_10marks_CTCFaverage.parafile
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_impute_11sample.input 8impute11_10marks_CTCFaverage.input
#
while read sample
do
    echo ${sample}
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/dhs_ctcf2/8sample_averageCTCFsignal.txt" >> 8impute11_10marks_CTCFaverage.input
done < /data/zusers/fankaili/ideas/dhs_ctcf/all_ctcf_sample.txt
#
nohup bash 8impute11_10marks_CTCFaverage.sh > ./nohup/nohup.8impute11_10marks_CTCFaverage.out 2>&1&
# z008 16128


# 4. 10 makrs + CTCG motif + average CTCF signal
## 1) train model
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model.sh model_8sample_10marks_CTCFmotif_CTCFaverage.sh
vim model_8sample_10marks_CTCFmotif_CTCFaverage.sh
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model.parafile model_8sample_10marks_CTCFmotif_CTCFaverage.parafile
vim model_8sample_10marks_CTCFmotif_CTCFaverage.parafile
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_model.input model_8sample_10marks_CTCFmotif_CTCFaverage.input
while read sample
do
    echo ${sample}" motif /data/zusers/fankaili/ideas/dhs_ctcf/CTCFmotif_signal_OCR.txt" >> model_8sample_10marks_CTCFmotif_CTCFaverage.input
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/dhs_ctcf2/8sample_averageCTCFsignal.txt" >> model_8sample_10marks_CTCFmotif_CTCFaverage.input
done < 8samples.txt
#
nohup bash model_8sample_10marks_CTCFmotif_CTCFaverage.sh > ./nohup/nohup.model_8sample_10marks_CTCFmotif_CTCFaverage.out 2>&1&
# z008 21196
## 2) imputation
awk '{FS=OFS}{if(NR==1){print $0}else{printf 100*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 100*$i};printf "\n"}}' /data/zusers/fankaili/ideas/dhs_ctcf2/model_8sample_10marks_CTCFmotif_CTCFaverage_result/model_8sample_10marks_CTCFmotif_CTCFaverage.para0 > /data/zusers/fankaili/ideas/dhs_ctcf2/model_8sample_10marks_CTCFmotif_CTCFaverage_result/model_8sample_10marks_CTCFmotif_CTCFaverage_100fold.para0
#
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_impute_11sample.sh 8impute11_10marks_CTCFmotif_CTCFaverage.sh
vim 8impute11_10marks_CTCFmotif_CTCFaverage.sh
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_impute_11sample.parafile 8impute11_10marks_CTCFmotif_CTCFaverage.parafile
vim 8impute11_10marks_CTCFmotif_CTCFaverage.parafile
cp /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_8sample_impute_11sample.input 8impute11_10marks_CTCFmotif_CTCFaverage.input
#
while read sample
do
    echo ${sample}
    #
    echo ${sample}" motif /data/zusers/fankaili/ideas/dhs_ctcf/CTCFmotif_signal_OCR.txt" >> 8impute11_10marks_CTCFmotif_CTCFaverage.input
    #
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/dhs_ctcf2/8sample_averageCTCFsignal.txt" >> 8impute11_10marks_CTCFmotif_CTCFaverage.input
done < /data/zusers/fankaili/ideas/dhs_ctcf/all_ctcf_sample.txt
#
nohup bash 8impute11_10marks_CTCFmotif_CTCFaverage.sh > ./nohup/nohup.8impute11_10marks_CTCFmotif_CTCFaverage.out 2>&1&
# z003 31739

# 5. calculate PRAU, make barplot
## 1) 10 marks
prefix="ctcf_8sample_impute_11sample"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}

## 2) 10 marks + CTCF motif
prefix="8impute11_10marks_CTCFmotif"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}

## 3) 10 marks + average CTCF signal
prefix="8impute11_10marks_CTCFaverage"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}

## 4) 10 marks + CTCF motif + average CTCF signal
prefix="8impute11_10marks_CTCFmotif_CTCFaverage"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}

## 5) with real data
prefix="ctcf_8sample_to_11sample"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}

## 6) make barplot
# Rscript get_PRAU_barplot.R

#!/bin/bash

# -- Kaili
# This script is for defining the peak cut-off for RAMPAGE peaks.
## Conclusion: choose RMA>1, remove ENCSR026WQO dataset.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/"
promoterShapeResult="/data/zusers/zhangx/projects/rampage/0_rampage_peak/peak_uniq/"
peakDir="/data/zusers/zhangx/projects/RAMPAGE_peaks/hg38/"

cd ${workDir}

# 1. make RPM histogram for each sample
nohup Rscript ${scriptDir}make_RAMPAGE_RPM_histogram.R > ./nohup/nohup.make_RAMPAGE_RPM_histogram.out 2>&1&
# z001 57101

# 2. count number of peaks with different cut-off
if [ -f num_peak_with_diff_cutoff.txt ];then rm num_peak_with_diff_cutoff.txt;fi
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    sample=`awk '{print $4}' <<< ${line}`
    echo ${expID}
    #
    num_RPM1=`awk '{if($13>1){print $0}}' ${peakDir}${expID}_rampage_peaks.txt | wc -l`
    num_RPM2=`awk '{if($13>2){print $0}}' ${peakDir}${expID}_rampage_peaks.txt | wc -l`
    num_RPM10=`awk '{if($13>10){print $0}}' ${peakDir}${expID}_rampage_peaks.txt | wc -l`
    #
    echo -e ${line}"\t"${num_RPM1}"\t"${num_RPM2}"\t"${num_RPM10} >> num_peak_with_diff_cutoff.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt
# Rscript make_figs_for_choose_RMAPAGE_peak.R

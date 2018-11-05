#!/bin/bash

# -- Kaili
# This script is for doing PCA for 10 marks to call states.
# 1. pick chr1 as example, get signals & calculate z-score

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mark_pca/"
workDir="/data/zusers/fankaili/ideas/mark_pca/"

cd ${workDir}

# 1. pick chr1 as example, get signals & calculate z-score

## 1) get chr1 100bp-bins
grep "chr1 " /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed > mm10_chr1_space.txt
grep "chr1 " /data/zusers/fankaili/ideas/run_ideas_p_value/mm10_tab.bed > mm10_chr1_tab.bed

## 2) get signals & calculate z-score
mkdir chr1_signal
mkdir chr1_signal_p
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    mark=`awk '{print $2}' <<< ${line}`
    # calculate z-score
    if [ mark == "DNAme" ]
    then
        python ${scriptDir}zscore_normalization_noLog.py /data/zusers/fankaili/ideas/run_ideas_p_value/signal/${sample}_${mark}.tab 6 > ${workDir}chr1_signal/${sample}_${mark}.tab
    else
        python ${scriptDir}zscore_normalization_noLog.py /data/zusers/fankaili/ideas/run_ideas_p_value/signal/${sample}_${mark}.tab 5 > ${workDir}chr1_signal/${sample}_${mark}.tab
    fi
    # get chr1 z-score
    head -959499 ${workDir}chr1_signal/${sample}_${mark}.tab | awk '{print $2}' > ${workDir}chr1_signal/${sample}_${mark}_chr1.txt
    # get chr1 signal
    head -959499 /data/zusers/fankaili/ideas/run_ideas_p_value/signal/${sample}_${mark}.txt > ${workDir}chr1_signal_p/${sample}_${mark}_chr1_pvalue.txt
done < /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input

# 2.

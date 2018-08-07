#!/bin/bash

# -- Kaili
# This script is for making histogram of CTCF signal/z-score
# 1) make histogram to check bimodal distribution
# 2) classification using EM
# 3) get new matrix divide ubi-rDHS based on CTCF signal (1: higher group; 0: lower)
# 4) another matrix divide ubi-rDHS using zscore threshold 1.64
# INPUT: matrix of signal/z-score
# OUTPUT: histogram and corresponding threshold: /data/zusers/fankaili/ccre/tf/figs/

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scipts/"

# 1. make histogram
## 1) ubi-rDHS log10(signal+0.01)
Rscript ${scriptDir}make_histogram.R /data/zusers/fankaili/ccre/tf/matrix/hg19_ubi-rDHS_CTCF_signal_matrix.txt /data/zusers/fankaili/ccre/tf/figs/ hg19_ubi-rDHS_CTCF_signal_log10_histogram.pdf log10

## 2) ubi-rDHS zscore
Rscript ${scriptDir}make_histogram.R /data/zusers/fankaili/ccre/tf/matrix/hg19_ubi-rDHS_CTCF_zscore_matrix.txt /data/zusers/fankaili/ccre/tf/figs/ hg19_ubi-rDHS_CTCF_zscore_histogram.pdf zscore


# 2. classification using EM
## 1) ubi-rDHS log10(signal+0.01)
Rscript ${scriptDir}classify_bimodal_EM.R /data/zusers/fankaili/ccre/tf/matrix/hg19_ubi-rDHS_CTCF_signal_matrix.txt /data/zusers/fankaili/ccre/tf/matrix/ hg19_ubi-rDHS_CTCF_signal_log10_classification.txt \
/data/zusers/fankaili/ccre/tf/figs/hg19_ubi-rDHS_CTCF_signal_log10_classification.pdf log10

## 2) ubi-rDHS zscore
Rscript ${scriptDir}classify_bimodal_EM.R /data/zusers/fankaili/ccre/tf/matrix/hg19_ubi-rDHS_CTCF_zscore_matrix.txt /data/zusers/fankaili/ccre/tf/matrix/ hg19_ubi-rDHS_CTCF_zscore_classification.txt \
/data/zusers/fankaili/ccre/tf/figs/hg19_ubi-rDHS_CTCF_zscore_classification.pdf zscore

# 3. classification using zscore>1.64
awk '{FS=OFS="\t"}{if($1=="id"){print $0}else{printf $1;for(i=2;i<=NF;i++){if($i>1.64){printf "\t"1}else{printf "\t"0}};printf "\n"}}' /data/zusers/fankaili/ccre/tf/matrix/hg19_rDHS_CTCF_zscore_matrix.txt > \
/data/zusers/fankaili/ccre/tf/matrix/hg19_ubi-rDHS_CTCF_zscore_1.64_classification.txt

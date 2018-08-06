#!/bin/bash

# -- Kaili
# This script is for making histogram of CTCF signal/z-score
# 1) make histogram to check bimodal distribution
# 2) get threshold
# 3) remake histogram with threshold
# INPUT: matrix of signal/z-score
# OUTPUT: histogram and corresponding threshold: /data/zusers/fankaili/ccre/tf/figs/

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scipts/"

# 1. make histogram
## 1) ubi-rDHS log10(signal+0.01)
Rscript ${scriptDir}make_histogram.R /data/zusers/fankaili/ccre/tf/matrix/hg19_ubi-rDHS_CTCF_signal_matrix.txt /data/zusers/fankaili/ccre/tf/figs/ hg19_ubi-rDHS_CTCF_signal_log10_histogram.pdf log10

## 2) ubi-rDHS zscore
Rscript ${scriptDir}make_histogram.R /data/zusers/fankaili/ccre/tf/matrix/hg19_ubi-rDHS_CTCF_zscore_matrix.txt /data/zusers/fankaili/ccre/tf/figs/ hg19_ubi-rDHS_CTCF_zscore_histogram.pdf zscore


# 2. get threshold using EM

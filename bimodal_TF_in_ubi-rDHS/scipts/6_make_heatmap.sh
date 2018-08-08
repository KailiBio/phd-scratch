#!/bin/bash

# -- Kaili
# This script is for making heatmaps.
# INPUT
# OUTPUT


# 1. make sure all the matrix are using the ccRE ID
cd /data/zusers/fankaili/ccre/tf/matrix/
mv hg19_ubi-rDHS_CTCF_signal_log10_classification.txt hg19_ubi-rDHS_CTCF_signal_log10_classification_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_CTCF_zscore_1.64_classification.txt > \
hg19_ubi-rDHS_CTCF_zscore_1.64_classification_ccreid.txt
mv hg19_ubi-rDHS_CTCF_zscore_classification.txt hg19_ubi-rDHS_CTCF_zscore_classification_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification.txt > \
hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification.txt > \
hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification_ccreid.txt

# 2. get PLS and non-PLS name list
awk '{if($9=="255,0,0"){print $4}}' /data/zusers/fankaili/ccre/hg19-cREs.bed > /data/zusers/fankaili/ccre/hg19_PLS_list.txt
awk '{if($9!="255,0,0"){print $4}}' /data/zusers/fankaili/ccre/hg19-cREs.bed > /data/zusers/fankaili/ccre/hg19_non_PLS_list.txt

# 3. run Rscript to make heatmap
# run locally
Rscript make_heatmap_jaccard.R

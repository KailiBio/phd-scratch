#!/bin/bash

# -- Kaili
# This script is for making heatmaps.
# INPUT
# OUTPUT


# 1. make sure all the matrix are using the ccRE ID & sort by row
cd /data/zusers/fankaili/ccre/tf/matrix/
mv hg19_ubi-rDHS_CTCF_signal_log10_classification.txt hg19_ubi-rDHS_CTCF_signal_log10_classification_ccreid.txt
# zscore 1.64
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt hg19_rDHS_CTCF_zscore_1.64_classification.txt > hg19_ubi-rDHS_CTCF_zscore_1.64_classification.txt
head -1 hg19_rDHS_CTCF_zscore_1.64_classification.txt > hg19_ubi-rDHS_CTCF_zscore_1.64_classification_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_CTCF_zscore_1.64_classification.txt | sort -k1n >> hg19_ubi-rDHS_CTCF_zscore_1.64_classification_ccreid.txt
# zscore em
head -1 hg19_ubi-rDHS_CTCF_zscore_classification.txt > hg19_ubi-rDHS_CTCF_zscore_em_classification_ccreid.txt
awk 'NR>1' hg19_ubi-rDHS_CTCF_zscore_classification.txt | sort -k1n >> hg19_ubi-rDHS_CTCF_zscore_em_classification_ccreid.txt
# H3K27ac
head -1 hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification.txt > hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification.txt | \
sort -k1n >> hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification_ccreid.txt
# H3K4me3
head -1 hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification.txt > hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification.txt | \
sort -k1n >> hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification_ccreid.txt

# 2. get PLS and non-PLS name list
## ubi-rDHS PLS
awk '{if($9=="255,0,0"){print $4}}' /data/zusers/fankaili/ccre/hg19-cREs.bed > /data/zusers/fankaili/ccre/hg19_PLS_list.txt
cat /data/zusers/fankaili/ccre/hg19_PLS_list.txt /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt | sort | uniq -d > /data/zusers/fankaili/ccre/hg19_ubi-rDHS_PLS_list.txt
# ubi-rDHS ELS & CTCF-only
awk '{if($9!="255,0,0"){print $4}}' /data/zusers/fankaili/ccre/hg19-cREs.bed > /data/zusers/fankaili/ccre/hg19_non_PLS_list.txt
cat /data/zusers/fankaili/ccre/hg19_non_PLS_list.txt /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt | sort | uniq -d > /data/zusers/fankaili/ccre/hg19_ubi-rDHS_non_PLS_list.txt

# 3. run Rscript to make heatmap
# run locally
Rscript make_heatmap_jaccard.R

# 4. How many cellline are using
## CTCF
head -1 hg19_ubi-rDHS_CTCF_zscore_em_classification_ccreid.txt | awk '{for(i=1;i<=NF;i++){split($i,a,"_");print a[1]}}' | sort -u > hg19_ubi-rDHS_CTCF_cellline_list.txt
head -1 hg19_ubi-rDHS_CTCF_zscore_em_classification_ccreid.txt | awk '{for(i=1;i<=NF;i++){split($i,a,"_");print a[1]}}' | sort -u | wc -l
# 35
head -1 hg19_ubi-rDHS_CTCF_zscore_em_classification_ccreid.txt | awk '{print NF}'
# 52
## H3K4me3
head -1 hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification_ccreid.txt | awk '{print NF}'
#22
## H3K27ac
head -1 hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification_ccreid.txt | awk '{print NF}'
# 17

#!/bin/bash

# -- Kaili
# This script is for getting bimodal matrix of SMC3 and RAD21.
# INPUT:
# OUTPUT:

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scipts/"

#####################
## SMC3

# 1. calculate zscore
mkdir /data/zusers/fankaili/ccre/tf/zscore_smc3/
bash ${scriptDir}calculate_zscore.sh -l /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF_using.txt -d /data/zusers/fankaili/ccre/tf/zscore_smc3/ -t SMC3

# 2. make zscore matrix
bash ${scriptDir}make_tf_zscore_matrix.sh -l /data/zusers/moorej3/ENCODE-Registry/hg19/V4/hg19-rDHSs.bed -d /data/zusers/fankaili/ccre/tf/matrix/ \
 -m hg19_rDHS_SMC3_zscore_matrix.txt -t SMC3 -f /data/zusers/fankaili/ccre/tf/zscore_smc3/

## get ubi-rDHS signal matrix
cd /data/zusers/fankaili/ccre/tf/matrix/
awk '{FS=OFS="\t"}{if($1=="id"){print $0}}' hg19_rDHS_SMC3_zscore_matrix.txt > hg19_ubi-rDHS_SMC3_zscore_matrix_DHSID.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt hg19_rDHS_SMC3_zscore_matrix.txt >> hg19_ubi-rDHS_SMC3_zscore_matrix_DHSID.txt
awk '{FS=OFS="\t"}{if($1=="id"){print $0}}' hg19_ubi-rDHS_SMC3_zscore_matrix_DHSID.txt > hg19_ubi-rDHS_SMC3_zscore_matrix.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_SMC3_zscore_matrix_DHSID.txt >> hg19_ubi-rDHS_SMC3_zscore_matrix.txt

# 3. classification by EM
Rscript ${scriptDir}classify_bimodal_EM.R /data/zusers/fankaili/ccre/tf/matrix/ SMC3 zscore

cd /data/zusers/fankaili/ccre/tf/matrix/
head -1 hg19_ubi-rDHS_SMC3_zscore_classification.txt > hg19_ubi-rDHS_SMC3_zscore_em_classification_ccreid.txt
awk 'NR>1' hg19_ubi-rDHS_SMC3_zscore_classification.txt | sort -k1n >> hg19_ubi-rDHS_SMC3_zscore_em_classification_ccreid.txt




#####################
## RAD21

# 1. calculate zscore
mkdir /data/zusers/fankaili/ccre/tf/zscore_smc3/
bash ${scriptDir}calculate_zscore.sh -l /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF_using.txt -d /data/zusers/fankaili/ccre/tf/zscore_rad21/ -t RAD21

# 2. make zscore matrix
bash ${scriptDir}make_tf_zscore_matrix.sh -l /data/zusers/moorej3/ENCODE-Registry/hg19/V4/hg19-rDHSs.bed -d /data/zusers/fankaili/ccre/tf/matrix/ \
 -m hg19_rDHS_RAD21_zscore_matrix.txt -t RAD21 -f /data/zusers/fankaili/ccre/tf/zscore_rad21/

## get ubi-rDHS signal matrix
cd /data/zusers/fankaili/ccre/tf/matrix/
awk '{FS=OFS="\t"}{if($1=="id"){print $0}}' hg19_rDHS_RAD21_zscore_matrix.txt > hg19_ubi-rDHS_RAD21_zscore_matrix_DHSID.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt hg19_rDHS_RAD21_zscore_matrix.txt >> hg19_ubi-rDHS_RAD21_zscore_matrix_DHSID.txt
awk '{FS=OFS="\t"}{if($1=="id"){print $0}}' hg19_ubi-rDHS_RAD21_zscore_matrix_DHSID.txt > hg19_ubi-rDHS_RAD21_zscore_matrix.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_RAD21_zscore_matrix_DHSID.txt >> hg19_ubi-rDHS_RAD21_zscore_matrix.txt

# 3. classification by EM
Rscript ${scriptDir}classify_bimodal_EM.R /data/zusers/fankaili/ccre/tf/matrix/ RAD21 zscore

cd /data/zusers/fankaili/ccre/tf/matrix/
head -1 hg19_ubi-rDHS_RAD21_zscore_classification.txt > hg19_ubi-rDHS_RAD21_zscore_em_classification_ccreid.txt
awk 'NR>1' hg19_ubi-rDHS_RAD21_zscore_classification.txt | sort -k1n >> hg19_ubi-rDHS_RAD21_zscore_em_classification_ccreid.txt




#####################
#####################
## heatmap

# 1. merge all TF into a big matrix
cd /data/zusers/fankaili/ccre/tf/matrix/

awk '{FS=OFS="\t"}{if(NR==1){for(i=1;i<=NF;i++){printf "CTCF_"$i"\t"};printf "\n"}else{print $0}}' hg19_ubi-rDHS_CTCF_zscore_em_classification_ccreid.txt > hg19_ubi-rDHS_TF_zscore_em_classification_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){if(NR==1){for(i=1;i<=NF;i++){printf $i"\t"}}else{a[$1]=$0}}else{if(FNR==1){for(i=1;i<=NF;i++){printf "SMC3_"$i"\t"};printf "\n"}else{if(a[$1]){b=$1;$1="";print a[b],$0}}}}' hg19_ubi-rDHS_TF_zscore_em_classification_ccreid.txt hg19_ubi-rDHS_SMC3_zscore_em_classification_ccreid.txt > temp_smc3.txt
awk '{FS=OFS="\t"}{if(NR==FNR){if(NR==1){for(i=1;i<=NF;i++){printf $i"\t"}}else{a[$1]=$0}}else{if(FNR==1){for(i=1;i<=NF;i++){printf "RAD21_"$i"\t"};printf "\n"}else{if(a[$1]){b=$1;$1="";print a[b],$0}}}}' temp_smc3.txt hg19_ubi-rDHS_RAD21_zscore_em_classification_ccreid.txt > hg19_ubi-rDHS_TF_zscore_em_classification_ccreid.txt
rm temp_smc3.txt


# 2. make heatmap
# run locally
Rscript make_tf_heatmap.R

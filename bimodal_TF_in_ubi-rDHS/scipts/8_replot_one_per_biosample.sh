#!/bin/bash

# -- Kaili
# This script is for re-plotting all the heatmaps only used Jill's data.
# one exp per biosample

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scipts/"

# 1. get CTCF z-score matrix
## 1)
## from 5_get_H3K4me3_H3K27ac_zscore_matrix.sh, got /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_data_master_file.txt

## 2) get master file list
python ${scriptDir}get_zscore_file_list_from_master.py /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_data_master_file.txt CTCF /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_CTCF_zscore_filelist.txt

## 3) get zscore of ubi-rDHS
python ${scriptDir}get_zscore.py /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_CTCF_zscore_filelist.txt /data/zusers/fankaili/ccre/tf/zscore_ctcf2/ CTCF

## 4) merge zscore into matrix, then do classification
cd /data/zusers/fankaili/ccre/tf/matrix/
echo "id" > hg19_ubi-rDHS_CTCF_2_zscore_matrix.txt
cat /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt >> hg19_ubi-rDHS_CTCF_2_zscore_matrix.txt
#
for file in `ls /data/zusers/fankaili/ccre/tf/zscore_ctcf2/`
do
    filename0=${file%_CTCF_zscore.txt} ;
    filename=${filename0#hg19_ubi-rDHS_};
    echo ${filename};
    echo -e "id\t"${filename} > temp.txt ;
    awk '{FS=" ";OFS="\t"}{print $1,$2}' /data/zusers/fankaili/ccre/tf/zscore_ctcf2/${file} >> temp.txt ;
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print $0,b[$1]}}}' temp.txt hg19_ubi-rDHS_CTCF_2_zscore_matrix.txt > temp_matrix.txt ;
    mv temp_matrix.txt hg19_ubi-rDHS_CTCF_2_zscore_matrix.txt ;
    rm temp*.txt ;
done
# change ID
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(FNR==1){print $0}else if(a[$1]){b=$1;$1=a[b];print $0}}}' \
/data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_CTCF_2_zscore_matrix.txt > hg19_ubi-rDHS_CTCF_2_zscore_matrix_ccreid.txt

# 2. classification
## 1) z-score > 1.64
cd /data/zusers/fankaili/ccre/tf/matrix/
awk '{FS=OFS="\t"}{if($1=="id"){print $0}else{printf $1;for(i=2;i<=NF;i++){if($i>1.64){printf "\t"1}else{printf "\t"0}};printf "\n"}}' hg19_ubi-rDHS_CTCF_2_zscore_matrix_ccreid.txt > hg19_ubi-rDHS_CTCF_2_zscore_1.64_classification.txt
head -1 hg19_ubi-rDHS_CTCF_2_zscore_1.64_classification.txt > hg19_ubi-rDHS_CTCF_2_zscore_1.64_classification_ccreid.txt
awk 'NR>1' hg19_ubi-rDHS_CTCF_2_zscore_1.64_classification.txt | sort -k1n >> hg19_ubi-rDHS_CTCF_2_zscore_1.64_classification_ccreid.txt

## 2) EM
### run locally
Rscript classify_bimodal_EM_1.R
#
cd /data/zusers/fankaili/ccre/tf/matrix/
head -1 hg19_ubi-rDHS_CTCF_2_zscore_classification.txt > hg19_ubi-rDHS_CTCF_2_zscore_em_classification_ccreid.txt
awk 'NR>1' hg19_ubi-rDHS_CTCF_2_zscore_classification.txt | sort -k1n >> hg19_ubi-rDHS_CTCF_2_zscore_em_classification_ccreid.txt

# 3. make heatmap
### run locally
Rscript make_heatmap_jaccard_1.R



#####################################
# 4. classify H3K4me3&H3K27ac by EM
Rscript classify_bimodal_EM_3.R
Rscript classify_bimodal_EM_4.R

# 5. make new heatmaps for all using EM classification
## H3K27ac
head -1 hg19_ubi-rDHS_H3K27ac_zscore_classification.txt > hg19_ubi-rDHS_H3K27ac_zscore_em_classification_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_H3K27ac_zscore_classification.txt | \
sort -k1n >> hg19_ubi-rDHS_H3K27ac_zscore_em_classification_ccreid.txt
## H3K4me3
head -1 hg19_ubi-rDHS_H3K4me3_zscore_classification.txt > hg19_ubi-rDHS_H3K4me3_zscore_em_classification_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_H3K4me3_zscore_classification.txt | \
sort -k1n >> hg19_ubi-rDHS_H3K4me3_zscore_em_classification_ccreid.txt

#!/bin/bash

# -- Kaili
# This script is for getting zscore matrix of H3K4me3 & H3K27ac
# 1) cell line list for the cell line we used for CTCF
# 2) H3K4me3 zscore matrix
# 3) H3K27ac zscore matrix
# INPUT:
# OUTPUT:

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scipts/"

# 1. get cell line list
head -1 /data/zusers/fankaili/ccre/tf/matrix/hg19_ccRE_CTCF_signal_matrix.txt | awk '{FS=OFS="\t"}{for(i=2;i<=NF;i++){split($i,a,"_");print a[1]}}' | sort -u > /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF_using.txt

# 2. get H3K4me3 zscore matrix
## 1) get master file line
for cellline in `cat /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF_using.txt`
do
    grep $cellline /data/zusers/moorej3/Registry-of-ccREs/hg19/V4/Cell-Type-Specific/Master-Cell-List.txt ;
done > /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_data_master_file.txt
#
vim /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_data_master_file.txt # modifying

## 2) get master file list
python ${scriptDir}get_zscore_file_list_from_master.py /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_data_master_file.txt H3K4me3 /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_H3K4me3_zscore_filelist.txt

## 3) get zscore of ubi-rDHS
python ${scriptDir}get_zscore.py /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_H3K4me3_zscore_filelist.txt /data/zusers/fankaili/ccre/tf/zscore_h3k4me3/ H3K4me3

## 4) merge zscore into matrix, then do classification
cd /data/zusers/fankaili/ccre/tf/matrix/
echo "id" > hg19_ubi-rDHS_H3K4me3_zscore_matrix.txt
cat /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt >> hg19_ubi-rDHS_H3K4me3_zscore_matrix.txt
#
for file in `ls /data/zusers/fankaili/ccre/tf/zscore_h3k4me3/`
do
    filename0=${file%_H3K4me3_zscore.txt} ;
    filename=${filename0#hg19_ubi-rDHS_};
    echo ${filename};
    echo -e "id\t"${filename} > temp.txt ;
    awk '{FS=" ";OFS="\t"}{print $1,$2}' /data/zusers/fankaili/ccre/tf/zscore_h3k4me3/${file} >> temp.txt ;
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print $0,b[$1]}}}' temp.txt hg19_ubi-rDHS_H3K4me3_zscore_matrix.txt > temp_matrix.txt ;
    mv temp_matrix.txt hg19_ubi-rDHS_H3K4me3_zscore_matrix.txt ;
    rm temp*.txt ;
done
#
awk '{FS=OFS="\t"}{if($1=="id"){print $0}else{printf $1;for(i=2;i<=NF;i++){if($i>1.64){printf "\t"1}else{printf "\t"0}};printf "\n"}}' hg19_ubi-rDHS_H3K4me3_zscore_matrix.txt > hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification.txt

# 3. get H3K27ac zscore matrix

## 2) get master file list
python ${scriptDir}get_zscore_file_list_from_master.py /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_data_master_file.txt H3K27ac /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_H3K27ac_zscore_filelist.txt

## 3) get zscore of ubi-rDHS
python ${scriptDir}get_zscore.py /data/zusers/fankaili/ccre/tf/hg19_cellline_with_TF_H3K27ac_zscore_filelist.txt /data/zusers/fankaili/ccre/tf/zscore_h3k27ac/ H3K27ac

## 4) merge zscore into matrix, then do classification
cd /data/zusers/fankaili/ccre/tf/matrix/
echo "id" > hg19_ubi-rDHS_H3K27ac_zscore_matrix.txt
cat /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt >> hg19_ubi-rDHS_H3K27ac_zscore_matrix.txt
#
for file in `ls /data/zusers/fankaili/ccre/tf/zscore_h3k27ac/`
do
    filename0=${file%_H3K27ac_zscore.txt} ;
    filename=${filename0#hg19_ubi-rDHS_};
    echo ${filename};
    echo -e "id\t"${filename} > temp.txt ;
    awk '{FS=" ";OFS="\t"}{print $1,$2}' /data/zusers/fankaili/ccre/tf/zscore_h3k27ac/${file} >> temp.txt ;
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print $0,b[$1]}}}' temp.txt hg19_ubi-rDHS_H3K27ac_zscore_matrix.txt > temp_matrix.txt ;
    mv temp_matrix.txt hg19_ubi-rDHS_H3K27ac_zscore_matrix.txt ;
    rm temp*.txt ;
done
#
awk '{FS=OFS="\t"}{if($1=="id"){print $0}else{printf $1;for(i=2;i<=NF;i++){if($i>1.64){printf "\t"1}else{printf "\t"0}};printf "\n"}}' hg19_ubi-rDHS_H3K27ac_zscore_matrix.txt > hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification.txt

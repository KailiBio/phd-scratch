#!/bin/bash

# -- Kaili
# This script is for calculating tissue-specificity index for gene/TSS.


scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}
cp /data/zusers/zhangx/seq/mouse_ccRE/ts_all.py /data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scripts/
#--------------------------------------------------------------------------------
# 1. get tissue-specificity index for gene
# 113/183

## 1) get normal tissues expression matrix
# 46 tissues in total 113 samples
python ${scriptDir}get_GRCh38_tissue_RNA_exp.py
#
echo "gene_id" > hg38_tissue_gene_exp_matrix.txt
cut -f 1 ./all_gene_exp/ENCSR023ZXN.txt >> hg38_tissue_gene_exp_matrix.txt
while read line
do
    expID=`awk '{FS=OFS="\t"}{print $1}' <<< $line` ;
    echo ${expID} ;
    echo ${expID} > temp.txt;
    awk '{print $2}' ./all_gene_exp/${expID}.txt >> temp.txt ;
    paste hg38_tissue_gene_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_tissue_gene_exp_matrix.txt;
done < hg38_tissue_gene_exp_list.txt
rm temp.txt


## 2) calculate tissue specificity index
awk '{if(NR>1){print $0}}' hg38_tissue_gene_exp_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_gene_exp_TSscore.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_gene_exp_TSscore.txt
rm tmp.txt
# #
# awk '{if(NR>1){print $0}}' hg38_all_gene_exp_matrix.txt > tmp.txt
# python ${scriptDir}ts_all.py tmp.txt hg38_all_gene_exp_TSscore.txt
# sed -i 's/-0.100000/NA/g' hg38_all_gene_exp_TSscore.txt
# rm tmp.txt

## 3) histogram and get cutoff
# Rscript find_TSscore_cutoff.R


#--------------------------------------------------------------------------------
# 2. get tissue-specificity index for TSS
# 104/155

## 1) get normal tissues expression matrix
# 45 tissues in total 104 samples
python ${scriptDir}get_GRCh38_tissue_RAMPAGE_exp.py
#
echo "TSS_id" > hg38_tissue_TSS_exp_matrix.txt
cut -f 1 ./rampage/ENCFF794RVT.tab >> hg38_tissue_TSS_exp_matrix.txt
while read line
do
    plus=`awk '{FS=OFS="\t"}{print $2}' <<< $line` ;
    echo ${plus} ;
    echo ${plus} > temp.txt;
    awk '{print $5}' ./rampage/${plus}.tab >> temp.txt ;
    paste hg38_tissue_TSS_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_tissue_TSS_exp_matrix.txt;
    minus=`awk '{FS=OFS="\t"}{print $3}' <<< $line` ;
    echo ${minus} ;
    echo ${minus} > temp.txt;
    awk '{print $5}' ./rampage/${minus}.tab >> temp.txt ;
    paste hg38_tissue_TSS_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_tissue_TSS_exp_matrix.txt;
done < hg38_tissue_TSS_exp_list.txt
rm temp.txt

## 2) calculate tissue specificity index
awk '{if(NR>1){print $0}}' hg38_tissue_TSS_exp_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_TSS_exp_TSscore.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_TSS_exp_TSscore.txt
rm tmp.txt

## 3) histogram and get cutoff
# Rscript find_TSscore_cutoff.R

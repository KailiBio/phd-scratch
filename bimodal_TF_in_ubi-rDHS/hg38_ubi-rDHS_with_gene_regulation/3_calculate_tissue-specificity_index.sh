#!/bin/bash

# -- Kaili
# This script is for calculating tissue-specificity index for gene/TSS.
# 1. get tissue-specificity index for gene
# 2. get tissue-specificity index for TSS
# 3. comparing tissue-specificity index between gene and TSS


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

################################################
# Oct04
## get donor ID, only used the matched RNA-seq & RAMPAGE tissue sample to calculate tissue-specificity index
awk '{FS=OFS="\t"}{print $1,$2,$4"\n"$1,$3,$4}' hg38_RAMPAGE_list.txt > hg38_RAMPAGE_list_duplicate.txt
python ${scriptDir}get_donorID.py /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list_duplicate.txt \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list_duplicate_donor.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=1;b[$6]=1;c[$5]=$0;d[$6]=$0}else{if(a[$5] && ($5!="---")){print $0,c[$5]}else if(a[$6] && ($6!="---")){print $0,c[$6]}else if(b[$5] && ($5!="---")){print $0,d[$5]}else if(b[$6] && ($6!="---")){print $0,d[$6]}}}' \
hg38_closest_gene_exp_list_donor.txt hg38_RAMPAGE_list_duplicate_donor.txt > hg38_RNA_RAMPAGE_matched_list.txt
# get tissue RNAseq list with matched donor: 103
awk '{FS=OFS="\t"}{if(NR==FNR){a[$7]=1}else{if(a[$1]){print $0}}}' hg38_RNA_RAMPAGE_matched_list.txt hg38_tissue_gene_exp_list.txt \
> hg38_tissue_gene_exp_list2.txt
# get tissue RAMPAGE list with matched donor: 103
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' hg38_RNA_RAMPAGE_matched_list.txt hg38_tissue_TSS_exp_list.txt \
> hg38_tissue_TSS_exp_list2.txt

################################################

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
done < hg38_tissue_gene_exp_list2.txt
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

### Dec 09
# get TS-index for quantile normalized signal
awk '{if(NR>1){print $0}}' hg38_tissue_gene_exp_matrix_quantile.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_gene_exp_TSscore_quantile.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_gene_exp_TSscore_quantile.txt
rm tmp.txt

awk '{FS=OFS="\t"}{if(NR==1){print $0}else{printf $1;for(i=2;i<=NF;i++){if($i>100){$i=100};printf "\t"$i};printf "\n"}}' \
hg38_tissue_gene_exp_matrix_quantile.txt > hg38_tissue_gene_exp_matrix_quantile_max100.txt
awk '{if(NR>1){print $0}}' hg38_tissue_gene_exp_matrix_quantile_max100.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_gene_exp_TSscore_quantile_max100.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_gene_exp_TSscore_quantile_max100.txt
rm tmp.txt

awk '{FS=OFS="\t"}{if(NR==1){print $0}else{printf $1;for(i=2;i<=NF;i++){if($i>20){$i=20};printf "\t"$i};printf "\n"}}' \
hg38_tissue_gene_exp_matrix_quantile.txt > hg38_tissue_gene_exp_matrix_quantile_max20.txt
awk '{if(NR>1){print $0}}' hg38_tissue_gene_exp_matrix_quantile_max20.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_gene_exp_TSscore_quantile_max20.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_gene_exp_TSscore_quantile_max20.txt
rm tmp.txt

## only maximum
awk '{FS=OFS="\t"}{if(NR==1){print $0}else{printf $1;for(i=2;i<=NF;i++){if($i>100){$i=100};printf "\t"$i};printf "\n"}}' \
hg38_tissue_gene_exp_matrix.txt > hg38_tissue_gene_exp_matrix_max100.txt
awk '{if(NR>1){print $0}}' hg38_tissue_gene_exp_matrix_max100.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_gene_exp_TSscore_max100.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_gene_exp_TSscore_max100.txt
rm tmp.txt

awk '{FS=OFS="\t"}{if(NR==1){print $0}else{printf $1;for(i=2;i<=NF;i++){if($i>20){$i=20};printf "\t"$i};printf "\n"}}' \
hg38_tissue_gene_exp_matrix.txt > hg38_tissue_gene_exp_matrix_max20.txt
awk '{if(NR>1){print $0}}' hg38_tissue_gene_exp_matrix_max20.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_gene_exp_TSscore_max20.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_gene_exp_TSscore_max20.txt
rm tmp.txt

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
    awk '{print $4}' ./rampage/${plus}.tab >> temp.txt ;
    paste hg38_tissue_TSS_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_tissue_TSS_exp_matrix.txt;
    minus=`awk '{FS=OFS="\t"}{print $3}' <<< $line` ;
    echo ${minus} ;
    echo ${minus} > temp.txt;
    awk '{print $4}' ./rampage/${minus}.tab >> temp.txt ;
    paste hg38_tissue_TSS_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_tissue_TSS_exp_matrix.txt;
done < hg38_tissue_TSS_exp_list2.txt
rm temp.txt

## 2) calculate tissue specificity index
awk '{if(NR>1){print $0}}' hg38_tissue_TSS_exp_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_TSS_exp_TSscore0.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_TSS_exp_TSscore0.txt
rm tmp.txt
### remove duplicate TSSs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(FNR==1){print $0}else if(a[$1]){print $0}}}' TSS.Filtered.uniq.bed \
hg38_tissue_TSS_exp_TSscore0.txt > hg38_tissue_TSS_exp_TSscore.txt
# change to TSS uniqID
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $8,a[$4]}}' hg38_tissue_TSS_exp_TSscore0.txt TSS.Filtered.uniqID.bed | sort -u > \
hg38_tissue_TSS_exp_TSscore_uniqID.txt

## 3) histogram and get cutoff
# Rscript find_TSscore_cutoff.R


#--------------------------------------------------------------------------------
# 3. comparing tissue-specificity index between gene and TSS

## 1) correlated tissue-specificity index of gene and TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7}else{print $1,$2,a[$1]}}' TSS.Filtered.bed hg38_tissue_TSS_exp_TSscore.txt > temp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$3,a[$3]}}' hg38_tissue_gene_exp_TSscore.txt temp.txt > hg38_TSindex_TSS_gene.txt
sed -i '1c TSS\tts1\tgene\tts2' hg38_TSindex_TSS_gene.txt
rm temp.txt

## 2) plot
# Rscript compare_TSindex_gene_TSS.R

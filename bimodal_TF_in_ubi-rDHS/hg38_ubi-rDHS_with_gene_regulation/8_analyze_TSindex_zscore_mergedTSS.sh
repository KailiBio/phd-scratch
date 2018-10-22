#!/bin/bash

# -- Kaili
# This script is for analyzing tissue-specificity index (using RAMPAGE/RNAseq signal z-sxcore) for merged-TSS.
# 0. get tissue-specificity index from signal z-score


#######
## Since z-score have negValue,can't use z-score for tissue-specificity index calculation.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}


# 0. get tissue-specificity index from signal z-score
## calculate z-score
mkdir rampage_zscore
for file in `ls /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/`
do
    python ${scriptDir}zscore-normalization.py /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/${file} > \
    /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_zscore/${file} 4
done
#
mkdir rampage_mergedTSS_zscore
for file in `ls /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_mergedTSS/`
do
    python ${scriptDir}zscore-normalization.py /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_mergedTSS/${file} > \
    /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_mergedTSS_zscore/${file} 4
done
#
mkdir rna_zscore
for file in `ls /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/`
do
    python ${scriptDir}zscore-normalization.py /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/${file} > \
    /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rna_zscore/${file} 2
done

## calculate tissue-specificity index
## function for getting RAMPAGE signal
getRampageSignal(){
    sed -i 's/ \t/\t/g' ${rampage_signal_dir}$1.tab
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $8,a[$4]}}' ${rampage_signal_dir}$1.tab ${tss_plus} | sort -u > ${rampage_signal_tissue_dir}$3_rampage.txt
    sed -i 's/ \t/\t/g' ${rampage_signal_dir}$2.tab
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $8,a[$4]}}' ${rampage_signal_dir}$2.tab ${tss_minus} | sort -u >> ${rampage_signal_tissue_dir}$3_rampage.txt
}

## 1) TSS
rampage_signal_dir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_zscore/"
tss_plus="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered.uniqID_plus.bed"
tss_minus="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered.uniqID_minus.bed"
rampage_signal_tissue_dir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_zscore/"

echo "TSS_id" > hg38_tissue_TSS_exp_zscore_matrix.txt
cut -f 1 /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_zscore/adrenal_gland_37_year_rampage.txt >> hg38_tissue_TSS_exp_zscore_matrix.txt
#
while read line
do
    id=`awk '{print $1}' <<< ${line}`
    plus=`awk '{print $2}' <<< ${line}`
    minus=`awk '{print $3}' <<< ${line}`
    biosample=`awk '{print $4}' <<< ${line}`
    #
    getRampageSignal ${plus} ${minus} ${biosample}
    #
    echo ${id} > temp.txt
    awk '{print $2}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_zscore/${biosample}_rampage.txt >> temp.txt
    paste hg38_tissue_TSS_exp_zscore_matrix.txt temp.txt > temp2.txt
    mv temp2.txt hg38_tissue_TSS_exp_zscore_matrix.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_TSS_exp_list2.txt
rm temp.txt

### calculate tissue-specificity index
awk '{if(NR>1){print $0}}' hg38_tissue_TSS_exp_zscore_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_TSS_exp_zscore_TSscore.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_TSS_exp_zscore_TSscore.txt
rm tmp.txt

## 2) mergedTSS
rampage_signal_dir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_mergedTSS_zscore/"
tss_plus="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/GRCh38_merged-TSS_gene_length_plus.txt"
tss_minus="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/GRCh38_merged-TSS_gene_length_minus.txt"
rampage_signal_tissue_dir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS_zscore/"

echo "TSS_id" > hg38_tissue_mergedTSS_exp_zscore_matrix.txt
cut -f 1 /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS_zscore/A172_rampage.txt >> hg38_tissue_mergedTSS_exp_zscore_matrix.txt
#
while read line
do
    id=`awk '{print $1}' <<< ${line}`
    plus=`awk '{print $2}' <<< ${line}`
    minus=`awk '{print $3}' <<< ${line}`
    biosample=`awk '{print $4}' <<< ${line}`
    #
    getRampageSignal ${plus} ${minus} ${biosample}
    #
    echo ${id} > temp.txt
    awk '{print $2}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS_zscore/${biosample}_rampage.txt >> temp.txt
    paste hg38_tissue_mergedTSS_exp_zscore_matrix.txt temp.txt > temp2.txt
    mv temp2.txt hg38_tissue_mergedTSS_exp_zscore_matrix.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_TSS_exp_list2.txt
rm temp.txt

### calculate tissue-specificity index
awk '{if(NR>1){print $0}}' hg38_tissue_mergedTSS_exp_zscore_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_mergedTSS_exp_zscore_TSscore.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_mergedTSS_exp_zscore_TSscore.txt
rm tmp.txt

## 3) gene
echo "gene" > hg38_tissue_gene_exp_zscore_matrix.txt
cut -f 1 /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rna_zscore/ENCSR257NIR.txt >> hg38_tissue_gene_exp_zscore_matrix.txt
#
while read line
do
    id=`awk '{print $1}' <<< ${line}`
    #
    echo ${id} > temp.txt
    awk '{print $2}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rna_zscore/${biosample}_rampage.txt >> temp.txt
    paste hg38_tissue_gene_exp_zscore_matrix.txt temp.txt > temp2.txt
    mv temp2.txt hg38_tissue_gene_exp_zscore_matrix.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_tissue_gene_exp_list2.txt
rm temp.txt

### calculate tissue-specificity index
awk '{if(NR>1){print $0}}' hg38_tissue_gene_exp_zscore_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_gene_exp_zscore_TSscore.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_gene_exp_zscore_TSscore.txt
rm tmp.txt

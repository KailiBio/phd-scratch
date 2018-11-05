#!/bin/bash

# -- Kaili
# This scirpt is re-calculate TSindex using all the samples, then compare the results.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# 1. recalculate TS index
## 1) RNA-seq
awk '{if(NR>1){print $0}}' hg38_all_gene_exp_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_all_gene_exp_TSscore.txt
sed -i 's/-0.100000/NA/g' hg38_all_gene_exp_TSscore.txt
rm tmp.txt

## 2) RAMPAGE
echo "TSS_id" > hg38_all_TSS_exp_matrix.txt
cut -f 1 ./rampage/ENCFF794RVT.tab >> hg38_all_TSS_exp_matrix.txt
while read line
do
    plus=`awk '{FS=OFS="\t"}{print $2}' <<< $line` ;
    echo ${plus} ;
    echo ${plus} > temp.txt;
    awk '{print $4}' ./rampage/${plus}.tab >> temp.txt ;
    paste hg38_all_TSS_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_all_TSS_exp_matrix.txt;
    minus=`awk '{FS=OFS="\t"}{print $3}' <<< $line` ;
    echo ${minus} ;
    echo ${minus} > temp.txt;
    awk '{print $4}' ./rampage/${minus}.tab >> temp.txt ;
    paste hg38_all_TSS_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_all_TSS_exp_matrix.txt;
done < hg38_RAMPAGE_list.txt
rm temp.txt
##
awk '{if(NR>1){print $0}}' hg38_all_TSS_exp_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_all_TSS_exp_TSscore0.txt
sed -i 's/-0.100000/NA/g' hg38_all_TSS_exp_TSscore0.txt
rm tmp.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$8}else{if(FNR==1){print $0}else{print a[$1],$2}}}' TSS.Filtered.uniqID.bed \
hg38_all_TSS_exp_TSscore0.txt | sort -k1,1n | uniq > hg38_all_TSS_exp_TSscore.txt


# 2. make figures and comparison
# Rscript comparison_TSindex.R

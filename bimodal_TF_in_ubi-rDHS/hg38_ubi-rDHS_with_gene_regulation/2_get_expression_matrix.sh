#!/bin/bash

# -- Kaili
# This script is for getting expression data.
# 1. RNAseq: gene expression
# 2. RAMPAGE: TSS expression


### 0. get ref -- using GENCODE.v24
cp /data/zusers/moorej3/moorej.ghpcc.project/Reference/Human/hg38/GENCODE24/gencode.v24.annotation.gtf /home/fankaili/genome/
awk '{FS=" ";OFS="\t"}{if($3=="gene"){split($10,a,"\"");split($12,b,"\"");split($16,c,"\"");print a[2],b[2],c[2]}}' \
/home/fankaili/genome/gencode.v24.annotation.gtf > hg38_geneID_geneType_geneSymbol.txt
#
awk '{if($2=="protein_coding"){print $0}}' hg38_geneID_geneType_geneSymbol.txt > hg38_proteinCoding_geneID_geneType_geneSymbol.txt




scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
#--------------------------------------------------------------------------------
# 1. gene expression

mkdir /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/
cd /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/
python ${scriptDir}get_GRCh38_RNA_exp.py
#
echo "gene_id" > hg38_all_gene_exp_matrix.txt
cut -f 1 ./all_gene_exp/ENCSR023ZXN.txt >> hg38_all_gene_exp_matrix.txt
for file in `ls /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/`
do
    echo $file;
    echo ${file%.txt} > temp.txt;
    awk '{print $2}' ./all_gene_exp/${file} >> temp.txt ;
    paste hg38_all_gene_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_all_gene_exp_matrix.txt;
done
rm temp.txt


#--------------------------------------------------------------------------------
# 2. TSS expression

# sorted TSS file
## get RAMPAGE signal for TSS up/down stream 50bp
cp /data/zusers/moorej3/moorej.ghpcc.project/Reference/Human/hg38/GENCODE24/TSS.Filtered.bed ./
awk '{FS=OFS="\t"}{if(($2-500)>0){print $1,$2-500,$3+50,$4,1,$6}else{print $1,0,$3+50,$4,1,$6}}' TSS.Filtered.bed | sort -k1,1 -V -s > TSS.Filtered_sorted.bed

python ${scriptDir}get_hg38_RAMPAGE_signal.py
#
echo "id" > hg38_tss_rampage_signal_matrix.txt
cut -f 4 TSS.Filtered_sorted.bed >> hg38_tss_rampage_signal_matrix.txt
for file in `ls /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/`
do
    filename=${file%.tab};
    echo $filename;
    echo -e "id\t"${filename} > temp.txt;
    cut -f 1,4 /data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/${file} >> temp.txt;
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$1]){print $0,a[$1]}}}' temp.txt hg38_tss_rampage_signal_matrix.txt > temp2.txt;
    mv temp2.txt hg38_tss_rampage_signal_matrix.txt;
done
rm temp.txt

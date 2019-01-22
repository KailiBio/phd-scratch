#!/bin/bash

# -- Kaili
# This script is for making tissue-specificity index of RNA&RAMPAGE.
# fig2

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# 1. RNA-seq: ubi-rOCRs overlapped genes vs. remaining genes
## gene TSindex after quantile: hg38_tissue_gene_exp_TSscore_quantile.txt
## ubi-rOCRs overlapped gene list: GRCh38_ubi-rOCR_overlapped_gene_id.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"ubi-rOCRs_overlapped_genes"}else{if($1!="gene"){print $1,$2,"remaining_genes"}}}}' \
GRCh38_ubi-rOCR_overlapped_gene_id.txt hg38_tissue_gene_exp_TSscore_quantile.txt > \
GRCh38_gene_TSindex_quantile.txt

# 2. RAMPAGE: ubi-rOCRs overlapped TSSs vs. non_overlapped TSSs in overlapped genes vs. remaining TSSs
#### do quantile
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$8}else{printf a[$1];for(i=2;i<=NF;i++){printf "\t"$i};printf "\n"}}' \
TSS.Filtered.uniqID.bed hg38_tissue_TSS_exp_matrix.txt | sort -u > hg38_tissue_TSS_exp_matrix_uniqID.txt
## get after normalization expression matrix: hg38_tissue_TSS_exp_matrix_uniqID_quantile.txt
awk '{if(NR>1){print $0}}' hg38_tissue_TSS_exp_matrix_uniqID.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_TSS_exp_TSscore_uniqID_quantile.txt
sed -i 's/-0.100000/NA/g' hg38_tissue_TSS_exp_TSscore_uniqID_quantile.txt
rm tmp.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]){print $8}}}' \
GRCh38_ubi-rOCR_overlapped_gene_id.txt TSS.Filtered.uniqID.bed | sort -u > \
GRCh38_allTSS_ubi-rOCR_overlapped_gene.txt
##----------
## TSS TSindex after quantile: hg38_tissue_TSS_exp_TSscore_uniqID_quantile.txt
## ubi-rOCRs overlapped TSS list: GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt
## ubi-rOCR overlapped gene's TSS list: GRCh38_allTSS_ubi-rOCR_overlapped_gene.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$2,"other_TSSs_in_ubi-rOCRs_overlapped_genes"}else{if($1!="gene"){print $1,$2,"remaining_TSSs"}}}}' \
GRCh38_allTSS_ubi-rOCR_overlapped_gene.txt hg38_tissue_TSS_exp_TSscore_uniqID_quantile.txt > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$2,"ubi-rOCRs_overlapped_TSSs"}else{if($1!="gene"){print $1,$2,$3}}}}' \
GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt tmp.txt > GRCh38_TSS_TSindex_quantile.txt








#3. density line plot
#Rscript fig2.TSindex_density.R

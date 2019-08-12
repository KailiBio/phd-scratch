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
intersectBed -a TSS.Filtered.bed -b GRCh38-rOCRs.bed -wa -u | cut -f 7 | sort -u > \
GRCh38_gene_with_TSSs_overlapping_rOCRs_list.txt

awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR>1){if(a[$1]){print $1,$2,"other_genes_overlapping_rOCRs"}else{print $1,$2,"genes_not_overlapping_rOCRs"}}}}' \
GRCh38_gene_with_TSSs_overlapping_rOCRs_list.txt hg38_tissue_gene_exp_TSscore_quantile.txt > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"genes_overlapping_ubi-rOCRs"}else{print $0}}}' \
GRCh38_ubi-rOCR_overlapped_gene_id.txt tmp.txt > GRCh38_gene_TSindex_quantile.txt

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
intersectBed -a TSS.Filtered.uniq.bed -b GRCh38-rOCRs.bed -wa -u > GRCh38_TSS_rOCRs_overlapped.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$8]=1}else{if(a[$1]){print $1,$2,"other_TSSs_overlapping_rOCRs"}else{print $1,$2,"TSSs_not_overlapping_rOCRs"}}}' \
GRCh38_TSS_rOCRs_overlapped.txt hg38_tissue_TSS_exp_TSscore_uniqID_quantile.txt > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$2,"other_TSSs_of_genes_overlapping_ubi-rOCRs"}else{print $0}}}' \
GRCh38_allTSS_ubi-rOCR_overlapped_gene.txt tmp.txt > tmp2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$2,"TSSs_overlapping_ubi-rOCRs"}else{if($1!="gene"){print $1,$2,$3}}}}' \
GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt tmp2.txt > GRCh38_TSS_TSindex_quantile.txt

rm tmp*.txt







#3. density line plot
#Rscript fig2.TSindex_density.R

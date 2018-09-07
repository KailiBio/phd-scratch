#!/bin/bash

# -- Kaili
# This script is for biological analysis of ubi-rDHS.
# 1. get closest genes.
# 2. GO analysis
# 3. intersect with loops.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scripts/"
workDir="/data/zusers/fankaili/ccre/tf/closest_gene/"

cd ${workDir}

# 1. get closest gene list
## 1) get TSS file we used
cp /data/zusers/moorej3/moorej.ghpcc.project/Reference/Human/hg19/Gencode19/TSS.Filtered.bed ./

## 2) get geneID-geneSymbol file
awk '{FS=OFS="\t"}{if($3=="gene"){split($9,a,";");print a[1],a[3]}}' /home/fankaili/genome/Homo_sapiens.GRCh37.87.gtf | cut -d \" -f 2,4 | \
awk '{FS="\"";OFS="\t"}{print $1,$2}' > /home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt

## 3) for 10921 ubi-rDHS
### get closest gene IDs
bedtools closest -a /data/zusers/fankaili/ccre/hg19_ubi-rDHS.bed -b TSS.Filtered.bed -d -k 1 > hg19_ubi-rDHS_closest_gene.bed
cut -f 11 hg19_ubi-rDHS_closest_gene.bed | sort -u | awk '{split($1,a,".");print a[1]}' > hg19_ubi-rDHS_closest_gene_list.txt
### get gene symbol
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $2}}}' hg19_ubi-rDHS_closest_gene_list.txt \
/home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt > hg19_ubi-rDHS_closest_gene_symbol_list.txt
### make ccre-geneID-geneSymbol matrix
awk '{FS=OFS="\t"}{split($11,a,".");print $4,a[1]}' hg19_ubi-rDHS_closest_gene.bed | sort -u > hg19_ubi-rDHS_closest_geneID.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$2]}}' /home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt \
hg19_ubi-rDHS_closest_geneID.bed > hg19_ubi-rDHS_closest_geneID_symbol.bed

## 4) for 10921-ubi-rDHS (overlapped)
awk '{if($13>0){split($11,a,".");print $4,a[1]}}' hg19_ubi-rDHS_closest_gene.bed | sort -u > hg19_ubi-rDHS_overlapped_geneID.bed
awk '{FS=OFS}{if(NR==FNR){a[$1]=$2}else{print $0,a[$2]}}' /home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt \
hg19_ubi-rDHS_overlapped_geneID.bed > hg19_ubi-rDHS_overlapped_geneID_symbol.bed
awk '{print $2}' hg19_ubi-rDHS_overlapped_geneID_symbol.bed | sort -u > hg19_ubi-rDHS_overlapped_gene_list.txt

## 5) for 10,921 ubi-rDHS PLS, ELS/CTCF-only
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/hg19_ubi-rDHS_PLS_list.txt \
hg19_ubi-rDHS_closest_geneID_symbol.bed > hg19_ubi-rDHS_PLS_closest_geneID_symbol.bed
awk '{print $2}' hg19_ubi-rDHS_PLS_closest_geneID_symbol.bed | sort -u > hg19_ubi-rDHS_PLS_closest_gene_list.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/hg19_ubi-rDHS_non_PLS_list.txt \
hg19_ubi-rDHS_closest_geneID_symbol.bed > hg19_ubi-rDHS_nonPLS_closest_geneID_symbol.bed
awk '{print $2}' hg19_ubi-rDHS_nonPLS_closest_geneID_symbol.bed | sort -u > hg19_ubi-rDHS_nonPLS_closest_gene_list.txt

## 6) for 29,733 ubi-rDHS
### get closest gene IDs
bedtools closest -a /data/zusers/fankaili/ccre/hg19_ubi-rDHS_2.bed -b TSS.Filtered.bed -d -k 1 > hg19_ubi-rDHS_2_closest_gene.bed
cut -f 11 hg19_ubi-rDHS_2_closest_gene.bed | sort -u | awk '{split($1,a,".");print a[1]}' > hg19_ubi-rDHS_2_closest_gene_list.txt
### make ccre-geneID-geneSymbol matrix
awk '{FS=OFS="\t"}{split($11,a,".");print $4,a[1]}' hg19_ubi-rDHS_2_closest_gene.bed | sort -u > hg19_ubi-rDHS_2_closest_geneID.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$2]}}' /home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt \
hg19_ubi-rDHS_2_closest_geneID.bed > hg19_ubi-rDHS_2_closest_geneID_symbol.bed
### for overlappped
awk '{if($13>0){split($11,a,".");print $4,a[1]}}' hg19_ubi-rDHS_2_closest_gene.bed | sort -u > hg19_ubi-rDHS_2_overlapped_geneID.bed
awk '{FS=OFS}{if(NR==FNR){a[$1]=$2}else{print $0,a[$2]}}' /home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt \
hg19_ubi-rDHS_2_overlapped_geneID.bed > hg19_ubi-rDHS_2_overlapped_geneID_symbol.bed
awk '{print $2}' hg19_ubi-rDHS_2_overlapped_geneID_symbol.bed | sort -u > hg19_ubi-rDHS_2_overlapped_gene_list.txt



# 2. run PANTHER using hg19_ubi-rDHS_closest_gene_list.txt



# 3. make figures by PANTHER result
### locally
cd /Users/kaili/Dropbox\ \(UMass\ Medical\ School\)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_BP.txt > hg19_ubi-rDHS_BP_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_CC.txt > hg19_ubi-rDHS_CC_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_MF.txt > hg19_ubi-rDHS_MF_clean.txt
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_BP_clean.txt" \
"hg19_ubi-rDHS_BP_15" "10,921 ubi-rDHS (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_CC_clean.txt" \
"hg19_ubi-rDHS_CC_15" "10,921 ubi-rDHS (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_MF_clean.txt" \
"hg19_ubi-rDHS_MF_15" "10,921 ubi-rDHS (top 15)" 15

awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_overlapped_BP.txt > hg19_ubi-rDHS_overlapped_BP_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_overlapped_CC.txt > hg19_ubi-rDHS_overlapped_CC_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_overlapped_MF.txt > hg19_ubi-rDHS_overlapped_MF_clean.txt
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_overlapped_BP_clean.txt" \
"hg19_ubi-rDHS_overlapped_BP_15" "10,921 ubi-rDHS overlapped genes (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_overlapped_CC_clean.txt" \
"hg19_ubi-rDHS_overlapped_CC_15" "10,921 ubi-rDHS overlapped genes (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_overlapped_MF_clean.txt" \
"hg19_ubi-rDHS_overlapped_MF_15" "10,921 ubi-rDHS overlapped genes (top 15)" 15

awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_PLS_BP.txt > hg19_ubi-rDHS_PLS_BP_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_PLS_CC.txt > hg19_ubi-rDHS_PLS_CC_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_PLS_MF.txt > hg19_ubi-rDHS_PLS_MF_clean.txt
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_PLS_BP_clean.txt" \
"hg19_ubi-rDHS_PLS_BP_15" "9,009 ubi-rDHS PLS (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_PLS_CC_clean.txt" \
"hg19_ubi-rDHS_PLS_CC_15" "9,009 ubi-rDHS PLS (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_PLS_MF_clean.txt" \
"hg19_ubi-rDHS_PLS_MF_15" "9,009 ubi-rDHS PLS (top 15)" 15

awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_nonPLS_BP.txt > hg19_ubi-rDHS_nonPLS_BP_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_nonPLS_CC.txt > hg19_ubi-rDHS_nonPLS_CC_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_nonPLS_MF.txt > hg19_ubi-rDHS_nonPLS_MF_clean.txt
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_nonPLS_BP_clean.txt" \
"hg19_ubi-rDHS_nonPLS_BP_15" "1,912 ubi-rDHS ELS/CTCF-only (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_nonPLS_CC_clean.txt" \
"hg19_ubi-rDHS_nonPLS_CC_15" "1,912 ubi-rDHS ELS/CTCF-only (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_nonPLS_MF_clean.txt" \
"hg19_ubi-rDHS_nonPLS_MF_15" "1,912 ubi-rDHS ELS/CTCF-only (top 15)" 15

awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_2_BP.txt > hg19_ubi-rDHS_2_BP_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_2_CC.txt > hg19_ubi-rDHS_2_CC_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_2_MF.txt > hg19_ubi-rDHS_2_MF_clean.txt
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_2_BP_clean.txt" \
"hg19_ubi-rDHS_2_BP_15" "29,733 ubi-rDHS (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_2_CC_clean.txt" \
"hg19_ubi-rDHS_2_CC_15" "29,733 ubi-rDHS (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_2_MF_clean.txt" \
"hg19_ubi-rDHS_2_MF_15" "29,733 ubi-rDHS (top 15)" 15

awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_2_overlapped_BP.txt > hg19_ubi-rDHS_2_overlapped_BP_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_2_overlapped_CC.txt > hg19_ubi-rDHS_2_overlapped_CC_clean.txt
awk '{if(NR>10){print $0}}' hg19_ubi-rDHS_2_overlapped_MF.txt > hg19_ubi-rDHS_2_overlapped_MF_clean.txt
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_2_overlapped_BP_clean.txt" \
"hg19_ubi-rDHS_2_overlapped_BP_15" "29,733 ubi-rDHS overlapped genes (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_2_overlapped_CC_clean.txt" \
"hg19_ubi-rDHS_2_overlapped_CC_15" "29,733 ubi-rDHS overlapped genes (top 15)" 15
#
Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_2_overlapped_MF_clean.txt" \
"hg19_ubi-rDHS_2_overlapped_MF_15" "29,733 ubi-rDHS overlapped genes (top 15)" 15

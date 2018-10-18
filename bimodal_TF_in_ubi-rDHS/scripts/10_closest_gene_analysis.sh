#!/bin/bash

# -- Kaili
# This script is for biological analysis of ubi-rDHS.
# 1. get closest genes.
# 2. GO analysis
# 3. intersect with loops.
# 4. basic info
## 1) records in notebook
## 2) composition of all 11188 closest genes
## 3) how many TSS overlapped with ubi-rDHS?
## 4) for each gene overlapped with ubi-rDHS, how many TSSs are overlapped?
## 5) how many ubi-rDHS overlapped genes are bidirectional?
## 6) How many bidirectional pairs in genome?

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
cut -f 11 hg19_ubi-rDHS_closest_gene.bed | sort -u > hg19_ubi-rDHS_closest_gene_list0.txt

### get gene symbol
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $2}}}' hg19_ubi-rDHS_closest_gene_list.txt \
/home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt > hg19_ubi-rDHS_closest_gene_symbol_list.txt
### make ccre-geneID-geneSymbol matrix
awk '{FS=OFS="\t"}{split($11,a,".");print $4,a[1]}' hg19_ubi-rDHS_closest_gene.bed | sort -u > hg19_ubi-rDHS_closest_geneID.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$2]}}' /home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt \
hg19_ubi-rDHS_closest_geneID.bed > hg19_ubi-rDHS_closest_geneID_symbol.bed

## 4) for 10921-ubi-rDHS (overlapped)
awk '{if($13==0){split($11,a,".");print $4,a[1]}}' hg19_ubi-rDHS_closest_gene.bed | sort -u > hg19_ubi-rDHS_overlapped_geneID.bed
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

# awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $11}}}' /data/zusers/fankaili/ccre/hg19_ubi-rDHS_non_PLS_list.txt \
# hg19_ubi-rDHS_closest_gene.bed | sort -u | wc -l


## 6) for 29,733 ubi-rDHS
### get closest gene IDs
bedtools closest -a /data/zusers/fankaili/ccre/hg19_ubi-rDHS_2.bed -b TSS.Filtered.bed -d -k 1 > hg19_ubi-rDHS_2_closest_gene.bed
cut -f 11 hg19_ubi-rDHS_2_closest_gene.bed | sort -u | awk '{split($1,a,".");print a[1]}' > hg19_ubi-rDHS_2_closest_gene_list.txt
### make ccre-geneID-geneSymbol matrix
awk '{FS=OFS="\t"}{split($11,a,".");print $4,a[1]}' hg19_ubi-rDHS_2_closest_gene.bed | sort -u > hg19_ubi-rDHS_2_closest_geneID.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$2]}}' /home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt \
hg19_ubi-rDHS_2_closest_geneID.bed > hg19_ubi-rDHS_2_closest_geneID_symbol.bed
### for overlappped
awk '{if($13==0){split($11,a,".");print $4,a[1]}}' hg19_ubi-rDHS_2_closest_gene.bed | sort -u > hg19_ubi-rDHS_2_overlapped_geneID.bed
awk '{FS=OFS}{if(NR==FNR){a[$1]=$2}else{print $0,a[$2]}}' /home/fankaili/genome/Homo_sapiens.GRCh37.ensemblID_geneSymbol.txt \
hg19_ubi-rDHS_2_overlapped_geneID.bed > hg19_ubi-rDHS_2_overlapped_geneID_symbol.bed
awk '{print $2}' hg19_ubi-rDHS_2_overlapped_geneID_symbol.bed | sort -u > hg19_ubi-rDHS_2_overlapped_gene_list.txt

# bedtools intersect -u -a TSS.Filtered.bed -b /data/zusers/fankaili/ccre/hg19_ubi-rDHS_2.bed | awk '{print $7}' | sort -u | wc -l

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




# 4. basic info
## get basic info from /data/zusers/fankaili/ccre/tf/closest_gene/hg19_ubi-rDHS_closest_gene.bed
## 1) records in notebook

Rscript make_closes_gene_info_figs.R

## 2) composition of all 11188 closest genes
cd /data/zusers/fankaili/ccre/tf/closest_gene/
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if(a[$1]){print $1,a[$1]}}}' /home/fankaili/genome/hg19_geneID_geneType_geneSymbol.txt \
hg19_ubi-rDHS_closest_gene_list0.txt > hg19_ubi-rDHS_closest_gene_type.txt
awk '{FS=OFS="\t"}{sum[$2]++}END{for(i in sum){print i,sum[i]}}' hg19_ubi-rDHS_closest_gene_type.txt > hg19_ubi-rDHS_closest_gene_type_count.txt
#
Rscript make_closes_gene_info_figs.R


## 3) how many TSS overlapped with ubi-rDHS?
awk '{FS=OFS="\t"}{if($13==0){print $5,$6,$7,$10}}' hg19_ubi-rDHS_closest_gene.bed | sort -u > hg19_ubi-rDHS_overlapped_TSS.txt

## 4) for each gene overlapped with ubi-rDHS, how many TSSs are overlapped?
awk '{FS=OFS="\t"}{if($13==0){print $0}}' hg19_ubi-rDHS_closest_gene.bed | cut -f 5,6,7,10,11 | sort -u | cut -f 5 | sort | uniq -c > hg19_ubi-rDHS_overlaped_TSS_count.txt
awk '{FS=OFS="\t"}{if(NR==FNR){if($13==0){a[$11]=1}}else{if(a[$7]){print $1,$2,$3,$6,$7}}}' hg19_ubi-rDHS_closest_gene.bed TSS.Filtered.bed | sort -u | cut -f 5 | sort | uniq -c > hg19_ubi-rDHS_overlaped_TSS_total.txt

## 5) how many ubi-rDHS overlapped genes are bidirectional?
awk '{FS=OFS="\t"}{if($13==0 && $10=="+"){print $5,$6,$7,$11}}' hg19_ubi-rDHS_closest_gene.bed | sort -u | sort -k1,1 -k2,2n \
> hg19_ubi-rDHS_TSS_plus.bed
#
awk '{FS=OFS="\t"}{if($13==0 && $10=="-"){print $5,$6,$7,$11}}' hg19_ubi-rDHS_closest_gene.bed | sort -u | sort -k1,1 -k2,2n \
> hg19_ubi-rDHS_TSS_minus.bed
#
bedtools window -a hg19_ubi-rDHS_TSS_plus.bed -b hg19_ubi-rDHS_TSS_minus.bed -l 1000 -r 300 | awk '{FS=OFS="\t"}{print $4,$8}' | sort -u \
| sort -k1,1 -k2,2 > hg19_ubi-rDHS_bidireactional_gene.txt
### 2192 bidirectional pairs, 4189 genes

## 6) How many bidirectional pairs in genome?
awk '{FS=OFS="\t"}{if($6=="+"){print $1,$2,$3,$7}}' TSS.Filtered.bed | sort -u | sort -k1,1 -k2,2n > TSS.Filtered_plus.bed
awk '{FS=OFS="\t"}{if($6=="-"){print $1,$2,$3,$7}}' TSS.Filtered.bed | sort -u | sort -k1,1 -k2,2n > TSS.Filtered_minus.bed
bedtools window -a TSS.Filtered_plus.bed -b TSS.Filtered_minus.bed -l 1000 -r 300 | awk '{FS=OFS="\t"}{print $4,$8}' | sort -u \
| sort -k1,1 -k2,2 > hg19_bidireactional_gene.txt
### 4632 bidirectional pairs, 8624 genes
cat hg19_ubi-rDHS_bidireactional_gene.txt hg19_bidireactional_gene.txt | sort | uniq -u | sort -k1,1 -k2,2 \
> hg19_non_ubi-rDHS_bidirectional_gene.txt

#!/bin/bash

# -- Kaili
# This script is for getting ubi-rOCRs, and do basic analysis.

EDGE_hg38="/data/projects/psychencode/Registry/V1/GRCh38/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# 1. get ubi-rOCRs
cp ${EDGE_hg38}GRCh38-rOCRs.bed ./
cp ${EDGE_hg38}GRCh38-EDGEs.bed ./
awk '{FS=OFS="\t"}{if(NR==FNR && $2>=580){a[$1]=1}else{if(a[$4]){print $0}}}' \
/data/projects/psychencode/Registry/V1/GRCh38/Biosample-Counts/EDGE-DNase-Biosample-Counts.txt GRCh38-rOCRs.bed > GRCh38_ubi-rOCRs.bed

## if all ubi-rOCRs are EDGEs?
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$0}else{if(a[$4]){print b[$4]}}}' GRCh38-EDGEs.bed GRCh38_ubi-rOCRs.bed > GRCh38_ubi-rOCRs_EDGEid.bed



# 2. get nearest genes
mkdir closest_gene
bedtools closest -a GRCh38_ubi-rOCRs_EDGEid.bed -b TSS.Filtered.bed -d -k 1 > ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed

## how many genes
awk '{if($15==0){print $14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 8962
awk '{if($15<=2000){print $14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 9814
awk '{if($15<=5000){print $14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 9988
awk '{if($15<=10000){print $14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 10124
awk '{print $14}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 10348

## how many TSSs
awk '{if($15==0){print $8,$9,$10}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 31791
awk '{if($15<=2000){print $8,$9,$10}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 33078
awk '{if($15<=5000){print $8,$9,$10}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 33293
awk '{if($15<=10000){print $8,$9,$10}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 33471
awk '{print $8,$9,$10}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 33745

## how many ubi-rOCRs
awk '{if($15==0){print $4}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 8272
awk '{if($15<=2000){print $4}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 10959
awk '{if($15<=5000){print $4}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 11254
awk '{if($15<=10000){print $4}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 11494
awk '{print $4}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 11888

# Rscript make_basic_figures.R


# 3. how many TSSs overlapped with each ubi-rOCRs? genes?
awk '{FS=OFS="\t"}{if($15==0){print $4,$8,$9,$10}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | cut -f 1 | sort \
| uniq -c > ./closest_gene/GRCh38_ubi-rOCR_overlapped_TSS_count.txt
#
awk '{FS=OFS="\t"}{if($15==0){print $4,$14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | cut -f 1 | sort \
| uniq -c > ./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_count.txt

# Rscript make_basic_figures.R


# 4. GO analysis for ubi-rOCRs overlapped genes
awk '{if($15==0){print $14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u > GRCh38_ubi-rOCR_overlapped_gene_id.txt

## run PANTHER using file

awk '{if(NR>10){print $0}}' GRCh38_ubi-rOCR_overlapped_BP.txt > GRCh38_ubi-rOCR_overlapped_BP_clean.txt
awk '{if(NR>10){print $0}}' GRCh38_ubi-rOCR_overlapped_CC.txt > GRCh38_ubi-rOCR_overlapped_CC_clean.txt
awk '{if(NR>10){print $0}}' GRCh38_ubi-rOCR_overlapped_MF.txt > GRCh38_ubi-rOCR_overlapped_MF_clean.txt

Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCR_overlapped_BP_clean.txt" \
"GRCh38_ubi-rOCR_BP_15" "11,897 ubi-rDHS (top 15)" 15

Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCR_overlapped_CC_clean.txt" \
"GRCh38_ubi-rOCR_CC_15" "11,897 ubi-rDHS (top 15)" 15

Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCR_overlapped_MF_clean.txt" \
"GRCh38_ubi-rOCR_MF_15" "11,897 ubi-rDHS (top 15)" 15

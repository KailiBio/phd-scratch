#!/bin/bash

# -- Kaili
# This script is for analyzing bidirectional promoters.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

#Rscript make_bidirectional_figure.R

# 1. get bidireactional gene pairs
awk '{if($6=="+"){print $0}}' /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed > ./basic_annotation/hg38_v28_basic_TSS_plus.bed
awk '{if($6=="-"){print $0}}' /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed > ./basic_annotation/hg38_v28_basic_TSS_minus.bed
# calculate bidirectional
bedtools window -a ./basic_annotation/hg38_v28_basic_TSS_plus.bed -b ./basic_annotation/hg38_v28_basic_TSS_minus.bed -l 1000 -r 300 > ./bidirectional/bidirectional_pair.bed
# 15849 items
### basic calculation
cut -f 7,14 ./bidirectional/bidirectional_pair.bed | sort -u | wc -l
# 4039 gene pairs
cut -f 7 ./bidirectional/bidirectional_pair.bed | sort -u | wc -l
cut -f 14 ./bidirectional/bidirectional_pair.bed | sort -u | wc -l
# 3869 plus strand genes, 3801 minuts strand genes
# 7670 genes in total
# ~14% genes are bidirectional genes

# 2. bidireactional TSSs overlapping ubi-rOCRs
cut -f 1-7  ./bidirectional/bidirectional_pair.bed | sort -u > ./bidirectional/bidirectional_plus_TSS.bed
cut -f 8-14  ./bidirectional/bidirectional_pair.bed | sort -u > ./bidirectional/bidirectional_minus_TSS.bed
#
intersectBed -a ./bidirectional/bidirectional_plus_TSS.bed -b GRCh38_ubi-rOCRs_EDGEid.bed -u > ./bidirectional/bidirectional_plus_TSS_ubi-rOCRs.bed
wc -l ./bidirectional/bidirectional_plus_TSS_ubi-rOCRs.bed
intersectBed -a ./bidirectional/bidirectional_minus_TSS.bed -b GRCh38_ubi-rOCRs_EDGEid.bed -wa -u > ./bidirectional/bidirectional_minus_TSS_ubi-rOCRs.bed
wc -l ./bidirectional/bidirectional_minus_TSS_ubi-rOCRs.bed
###
intersectBed -a ./basic_annotation/hg38_v28_basic_TSS_plus.bed -b GRCh38_ubi-rOCRs_EDGEid.bed -u | cut -f 1-3,6 | sort -u | wc -l
intersectBed -a ./basic_annotation/hg38_v28_basic_TSS_minus.bed -b GRCh38_ubi-rOCRs_EDGEid.bed -u | cut -f 1-3,6 | sort -u | wc -l

# 3. How many pairs of bidirectional promoters have their TSSs overlap two ubi-rDHSs? One uib-rDHS? Only none-ubi-rDHSs?
intersectBed -a /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed -b GRCh38_ubi-rOCRs_EDGEid.bed -u > ./basic_annotation/hg38_v28_basic_TSS_overlapping_ubi-rOCRs.bed
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4] && a[$11]){print $0}}}' ./basic_annotation/hg38_v28_basic_TSS_overlapping_ubi-rOCRs.bed ./bidirectional/bidirectional_pair.bed | sort -u | wc -l
# 7204 both end overlap ubi-rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]==1 && a[$11]!=1){print $0}}}' ./basic_annotation/hg38_v28_basic_TSS_overlapping_ubi-rOCRs.bed ./bidirectional/bidirectional_pair.bed | sort -u | wc -l
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1 && a[$11]==1){print $0}}}' ./basic_annotation/hg38_v28_basic_TSS_overlapping_ubi-rOCRs.bed ./bidirectional/bidirectional_pair.bed | sort -u | wc -l
# 3728 just one end overlapping ubi-rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1 && a[$11]!=1){print $0}}}' ./basic_annotation/hg38_v28_basic_TSS_overlapping_ubi-rOCRs.bed ./bidirectional/bidirectional_pair.bed | sort -u | wc -l
# 4917 both end not overlapping ubi-rOCRs


# 4. for those both end overlapping ubi-rOCRs, do they overlapping the same ubi-rOCRs?
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4] && a[$11]){print $0}}}' ./basic_annotation/hg38_v28_basic_TSS_overlapping_ubi-rOCRs.bed ./bidirectional/bidirectional_pair.bed | sort -u > ./bidirectional/bidirectional_both_end_ubi-rOCRs.bed
awk '{FS=OFS="\t"}{if($2<$9){print $1,$2,$9,"bi_"NR}else{print $1,$9,$2,"bi_"NR}}' ./bidirectional/bidirectional_both_end_ubi-rOCRs.bed > ./bidirectional/bidirectional_both_end_ubi-rOCRs_region.bed
intersectBed -a ./bidirectional/bidirectional_both_end_ubi-rOCRs_region.bed -b GRCh38_ubi-rOCRs_EDGEid.bed -f 1 -u | wc -l
# 5894 pairs in the same ubi-rOCRs.

# 5. Do these tend to be protein-coding or lncRNAs?
## 1) both end
cut -f 7,14 ./bidirectional/bidirectional_both_end_ubi-rOCRs.bed > ./bidirectional/bidirectional_both_end_ubi-rOCRs_gene_list.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$6}else{print $1,a[$1],$2,a[$2]}}' /home/fankaili/genome/hg38_v28_basic_gene_filtered.txt ./bidirectional/bidirectional_both_end_ubi-rOCRs_gene_list.txt | sort -u > ./bidirectional/bidirectional_both_end_ubi-rOCRs_gene_type.txt
cut -f 1-2 ./bidirectional/bidirectional_both_end_ubi-rOCRs_gene_type.txt > tmp.txt
cut -f 3-4 ./bidirectional/bidirectional_both_end_ubi-rOCRs_gene_type.txt >> tmp.txt
sort -u tmp.txt | cut -f 2 | sort | uniq -c
## 2) all bidirectiona genes that overlapping ubi-rOCRs
cut -f 1-7 ./bidirectional/bidirectional_pair.bed > tmp.bed
cut -f 8-14 ./bidirectional/bidirectional_pair.bed >> tmp.bed
sort -u tmp.bed | sort -k1,1 -k2,2n > tmp2.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0}}}' ./basic_annotation/hg38_v28_basic_TSS_overlapping_ubi-rOCRs.bed tmp2.bed > tmp3.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$6}else{print $7,a[$7]}}' /home/fankaili/genome/hg38_v28_basic_gene_filtered.txt tmp3.txt | sort -u | cut -f 2 | sort | uniq -c

#!/bin/bash

# -- Kaili
# This script is for getting ubi-rOCRs, and do basic analysis.
# 1. get ubi-rOCRs
# 2. get nearest genes
# 3. how many TSSs overlapped with each ubi-rOCRs? genes?
# 4. GO analysis for ubi-rOCRs overlapped genes
# 5. get ubi-rORCs overlapped TSSs list and loci
# 6. for ubi-rOCRs overlapped genes, the percentage of overlapped TSSs.
# 7. how many ubi-rDHS overlapped genes are bidirectional?

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
### no


# 2. get nearest genes
mkdir closest_gene
bedtools closest -a GRCh38_ubi-rOCRs_EDGEid.bed -b TSS.Filtered.bed -d -k 1 > ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed
awk '{if($15==0){print $14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u > ./closest_gene/GRCh38_ubi-rOCR_closest_gene_list.bed
awk '{if($15==0){print $8,$9,$10,$13,$11}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -k1,4 -u | cut -d " " -f 5 > \
./closest_gene/GRCh38_ubi-rOCR_closest_TSS_list.bed


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
awk '{if($15==0){print $8,$9,$10,$13}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 31800
awk '{if($15<=2000){print $8,$9,$10,$13}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 33087
awk '{if($15<=5000){print $8,$9,$10,$13}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 33302
awk '{if($15<=10000){print $8,$9,$10,$13}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 33480
awk '{print $8,$9,$10,$13}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | wc -l
# 33754

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
awk '{FS=OFS="\t"}{if($15==0){print $4,$8,$9,$10,$13}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | cut -f 1 | sort \
| uniq -c > ./closest_gene/GRCh38_ubi-rOCR_overlapped_TSS_count.txt
#
# awk '{FS=OFS="\t"}{if($15==0){print $4,$14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u | cut -f 1 | sort \
# | uniq -c > ./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_count.txt
# add control
bedtools closest -a GRCh38-rOCRs.bed -b TSS.Filtered.bed -d -k 1 | awk '{FS=OFS="\t"}{if($12==0){print $4,$5,$6,$7,$10}}' | sort -u \
| cut -f 1 | sort | uniq -c > ./closest_gene/GRCh38_rOCR_overlapped_TSS_count.txt

bedtools closest -a GRCh38_ubi-rOCRs_EDGEid.bed -b TSS.Filtered.bed -d -k 1 > ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed
# Rscript make_basic_figures.R


# 4. GO analysis for ubi-rOCRs overlapped genes
awk '{if($15==0){print $14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u > GRCh38_ubi-rOCR_overlapped_gene_id.txt

## run PANTHER using file

awk '{if(NR>10){print $0}}' GRCh38_ubi-rOCR_overlapped_BP.txt > GRCh38_ubi-rOCR_overlapped_BP_clean.txt
awk '{if(NR>10){print $0}}' GRCh38_ubi-rOCR_overlapped_CC.txt > GRCh38_ubi-rOCR_overlapped_CC_clean.txt
awk '{if(NR>10){print $0}}' GRCh38_ubi-rOCR_overlapped_MF.txt > GRCh38_ubi-rOCR_overlapped_MF_clean.txt

Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCR_overlapped_BP_clean.txt" \
"GRCh38_ubi-rOCR_BP_15" "11,888 ubi-rDHS (top 15)" 15

Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCR_overlapped_CC_clean.txt" \
"GRCh38_ubi-rOCR_CC_15" "11,888 ubi-rDHS (top 15)" 15

Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/" \
"/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCR_overlapped_MF_clean.txt" \
"GRCh38_ubi-rOCR_MF_15" "11,888 ubi-rDHS (top 15)" 15


# 5. get ubi-rORCs overlapped TSSs list and loci
awk '{if($15==0){print $8,$9,$10}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u > GRCh38_ubi-rOCR_overlapped_TSS_loci.txt
awk '{if($15==0){print $8,$9,$10}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u > GRCh38_ubi-rOCR_overlapped_TSS_loci.txt


# 6. for ubi-rOCRs overlapped genes, the percentage of overlapped TSSs.

sort -k7,7 -k1,1 -k2,2n TSS.Filtered.bed > TSS.Filtered_sorted_gene.bed
## 1) merge TSSs within 50bp
python ${scriptDir}merge_TSS_50bp.py

## 2) get ubi-rOCR overlapped
sort -k1,1 -k2,2n hg38_merged_TSS_gene.bed > hg38_merged_TSS_gene_sorted.bed
bedtools closest -a GRCh38_ubi-rOCRs_EDGEid.bed -b hg38_merged_TSS_gene_sorted.bed -d -k 1 > ./closest_gene/GRCh38_ubi-rOCR_closest_gene_merged.bed
#
cut -f 11,14 ./closest_gene/GRCh38_ubi-rOCR_closest_gene_merged.bed | sort -u | cut -f 2 | sort | uniq -c > \
./closest_gene/hg38_gene_overlapped_TSS_count_merged.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$14]=1}else{if(a[$7]){print $4,$7}}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene_merged.bed \
hg38_merged_TSS_gene.bed | sort -u | cut -f 2 | sort | uniq -c > ./closest_gene/hg38_gene_all_TSS_count_merged.txt



# 7. how many ubi-rDHS overlapped genes are bidirectional?
## 1) GRCh38 overlapped
## get plus&minus gene with TSS.
awk '{FS=OFS="\t"}{if($15==0){print $11,$13,$14}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed | sort -u > ss.txt
awk '{FS=OFS="\t"}{if(NR==FNR){if($2=="+"){b[$1]=1}}else{if(b[$4]){print $1,$2,$3,$8,$5,$6,$7}}}' ss.txt TSS.Filtered.uniqID.bed | \
sort -u | sort -k1,1 -k2,2n > GRCh38_ubi-rOCR_overlapped_plus_TSSID_gene.bed
awk '{FS=OFS="\t"}{if(NR==FNR){if($2=="-"){b[$1]=1}}else{if(b[$4]){print $1,$2,$3,$8,$5,$6,$7}}}' ss.txt TSS.Filtered.uniqID.bed | \
sort -u | sort -k1,1 -k2,2n > GRCh38_ubi-rOCR_overlapped_minus_TSSID_gene.bed
rm ss.txt

## find bidirectional
bedtools window -a GRCh38_ubi-rOCR_overlapped_plus_TSSID_gene.bed -b GRCh38_ubi-rOCR_overlapped_minus_TSSID_gene.bed -l 1000 -r 300 | \
awk '{FS=OFS="\t"}{print $7,$14}' | sort -u | sort -k1,1 -k2,2 > GRCh38_ubi-rOCR_overlapped_bidirectional_gene.txt

## 2) How many bidirectional pairs in genome?
awk '{FS=OFS="\t"}{if($6=="+"){print $1,$2,$3,$8,$7}}' TSS.Filtered.uniqID.bed | sort -u | sort -k1,1 -k2,2n > TSS.Filtered.uniqID_plus.bed
awk '{FS=OFS="\t"}{if($6=="-"){print $1,$2,$3,$8,$7}}' TSS.Filtered.uniqID.bed | sort -u | sort -k1,1 -k2,2n > TSS.Filtered.uniqID_minus.bed
bedtools window -a TSS.Filtered.uniqID_plus.bed -b TSS.Filtered.uniqID_minus.bed -l 1000 -r 300 | awk '{FS=OFS="\t"}{print $5,$10}' | sort -u \
| sort -k1,1 -k2,2 > GRCh38_bidireactional_gene.txt

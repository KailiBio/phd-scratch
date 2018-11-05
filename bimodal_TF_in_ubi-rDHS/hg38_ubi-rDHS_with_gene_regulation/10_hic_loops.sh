#!/bin/bash

# -- Kaili
# This script is for analyzing 9,448 Hi-C loops.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/"

# 0. pre-processing
cd /data/zusers/fankaili/ccre/tf/loop/
## 1) get loops
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary%2Breplicate_HiCCUPS_looplist.txt.gz
--2018-10-22 09:55:07--  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary%2Breplicate_HiCCUPS_looplist.txt.gz
#
awk '{FS=OFS="\t"}{if(NR>1){print "chr"$1,$2,$3,"loop_"(NR-1)"_a","loop_"(NR-1)}}' GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt > GM12878_loop_a_hg19.txt
awk '{FS=OFS="\t"}{if(NR>1){print "chr"$4,$5,$6,"loop_"(NR-1)"_b","loop_"(NR-1)}}' GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt > GM12878_loop_b_hg19.txt

## 2) liftover to hg38
liftOver GM12878_loop_a_hg19.txt /home/fankaili/genome/hg19ToHg38.over.chain \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/GM12878_loop_a_hg38.txt \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/GM12878_loop_a_hg38_unMapped.txt
#
liftOver GM12878_loop_b_hg19.txt /home/fankaili/genome/hg19ToHg38.over.chain \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/GM12878_loop_b_hg38.txt \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/GM12878_loop_b_hg38_unMapped.txt
# merge
cd ${workDir}
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$0;b[$5]=1}else{if(b[$5]){print a[$5],$0}}}' GM12878_loop_a_hg38.txt GM12878_loop_b_hg38.txt > \
hg38_GM12878_loop.txt
### get 9441 hic loops
cut -f 1-5 hg38_GM12878_loop.txt > hg38_GM12878_loop_a.bed
cut -f 6-10 hg38_GM12878_loop.txt > hg38_GM12878_loop_b.bed

# 1. loops with ubi-rOCRs
intersectBed -a hg38_GM12878_loop_a.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > GRCh38_loop_a_overlapped_ubi-rOCRs.txt
intersectBed -a hg38_GM12878_loop_b.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > GRCh38_loop_b_overlapped_ubi-rOCRs.txt
#
cat GRCh38_loop_a_overlapped_ubi-rOCRs.txt GRCh38_loop_b_overlapped_ubi-rOCRs.txt > GRCh38_loop_overlapped_ubi-rOCRs.txt
cut -f 5 GRCh38_loop_overlapped_ubi-rOCRs.txt | sort -u | wc -l
cut -f 10,11 GRCh38_loop_overlapped_ubi-rOCRs.txt | sort -u | wc -l
cut -f 10,11 GRCh38_loop_overlapped_ubi-rOCRs.txt | sort -u | cut -f 2 | sort | uniq -c
#
cut -f 4,5 GRCh38_loop_overlapped_ubi-rOCRs.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l

# 2. loops with rOCRs
intersectBed -a hg38_GM12878_loop_a.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -wa -wb > GRCh38_loop_a_overlapped_rOCRs.txt
intersectBed -a hg38_GM12878_loop_b.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -wa -wb > GRCh38_loop_b_overlapped_rOCRs.txt

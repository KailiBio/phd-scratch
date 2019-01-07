#!/bin/bash

# -- Kaili
# This script is for basic analysis on ChIA-PET data.
# 0. pre-processing
# 1. tag overlapped with ubi-rOCRs
# 2. tag overlapped with GM-active rOCRs
# 3. tag overlapped with rOCRs


mkdir /data/zusers/fankaili/ccre/hg38_ubi-rDHS/chiapet

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/chiapet/"

cd ${workDir}

# 0. pre-processing
## 1) GM12878
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872886/suppl/GSM1872886_GM12878_CTCF_PET_clusters.txt.gz
gzip -d GSM1872886_GM12878_CTCF_PET_clusters.txt.gz
#
cut -f 1-3 GSM1872886_GM12878_CTCF_PET_clusters.txt | awk '{FS=OFS="\t"}{print $0,"loop_"NR"_a","loop_"NR}' \
> hg19_GM12878_loop_a.txt
cut -f 4-6 GSM1872886_GM12878_CTCF_PET_clusters.txt | awk '{FS=OFS="\t"}{print $0,"loop_"NR"_b","loop_"NR}' \
> hg19_GM12878_loop_b.txt
# liftOver
liftOver hg19_GM12878_loop_a.txt /home/fankaili/genome/hg19ToHg38.over.chain GRCh38_GM12878_loop_a.txt \
GRCh38_GM12878_loop_a_unMapped.txt
liftOver hg19_GM12878_loop_b.txt /home/fankaili/genome/hg19ToHg38.over.chain GRCh38_GM12878_loop_b.txt \
GRCh38_GM12878_loop_b_unMapped.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=1;b[$5]=$0}else{if(a[$5]){print b[$5],$0}}}' GRCh38_GM12878_loop_a.txt \
GRCh38_GM12878_loop_b.txt > GRCh38_GM12878_loop.txt

## 2) HeLa
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872888/suppl/GSM1872888_HeLa_CTCF_PET_clusters.txt.gz
gzip -d GSM1872888_HeLa_CTCF_PET_clusters.txt.gz
#
cut -f 1-3 GSM1872888_HeLa_CTCF_PET_clusters.txt | awk '{FS=OFS="\t"}{print $0,"loop_"NR"_a","loop_"NR}' \
> hg19_HeLa_loop_a.txt
cut -f 4-6 GSM1872888_HeLa_CTCF_PET_clusters.txt | awk '{FS=OFS="\t"}{print $0,"loop_"NR"_b","loop_"NR}' \
> hg19_HeLa_loop_b.txt
# liftOver
liftOver hg19_HeLa_loop_a.txt /home/fankaili/genome/hg19ToHg38.over.chain GRCh38_HeLa_loop_a.txt \
GRCh38_HeLa_loop_a_unMapped.txt
liftOver hg19_HeLa_loop_b.txt /home/fankaili/genome/hg19ToHg38.over.chain GRCh38_HeLa_loop_b.txt \
GRCh38_HeLa_loop_b_unMapped.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=1;b[$5]=$0}else{if(a[$5]){print b[$5],$0}}}' GRCh38_HeLa_loop_a.txt \
GRCh38_HeLa_loop_b.txt > GRCh38_HeLa_loop.txt



# 1. tag overlapped with ubi-rOCRs
cut -f 1-5 GRCh38_GM12878_loop.txt | sort -k1,1 -k2,2n > GRCh38_GM12878_loop_a_matched.bed
cut -f 6-10 GRCh38_GM12878_loop.txt | sort -k1,1 -k2,2n > GRCh38_GM12878_loop_b_matched.bed
intersectBed -a GRCh38_GM12878_loop_a_matched.bed -b \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > \
GRCh38_GM12878_PET_ubi-rOCR_a.txt
intersectBed -a GRCh38_GM12878_loop_b_matched.bed -b \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > \
GRCh38_GM12878_PET_ubi-rOCR_b.txt
#
cat GRCh38_GM12878_PET_ubi-rOCR_a.txt GRCh38_GM12878_PET_ubi-rOCR_b.txt > \
GRCh38_GM12878_PET_ubi-rOCR.txt
cut -f 5 GRCh38_GM12878_PET_ubi-rOCR.txt | sort -u | wc -l
# 16,498 loops overlapped
cut -f 4,5 GRCh38_GM12878_PET_ubi-rOCR.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l
# 1730 loops both ends
cut -f 9 GRCh38_GM12878_PET_ubi-rOCR.txt | sort -u | wc -l
# 5628 ubi-rOCRs
cut -f 9,11 GRCh38_GM12878_PET_ubi-rOCR.txt | sort -u | cut -f 2 | sort | uniq -c
# 11 CTCF-only, 753 ELS, 4864 PLS


# 2. tag overlapped with GM-active rOCRs
intersectBed -a GRCh38_GM12878_loop_a_matched.bed -b ../loop/GRCh38_GM12878_active_rOCRs.bed \
-wa -wb > GRCh38_GM12878_PET_active-rOCR_a.txt
intersectBed -a GRCh38_GM12878_loop_b_matched.bed -b ../loop/GRCh38_GM12878_active_rOCRs.bed \
-wa -wb > GRCh38_GM12878_PET_active-rOCR_b.txt
#
cat GRCh38_GM12878_PET_active-rOCR_a.txt GRCh38_GM12878_PET_active-rOCR_b.txt > \
GRCh38_GM12878_PET_active-rOCR.txt
cut -f 5 GRCh38_GM12878_PET_active-rOCR.txt | sort -u | wc -l
# 71,742 loops
cut -f 4-5 GRCh38_GM12878_PET_active-rOCR.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l
# 31,478 loops both ends overlapped
cut -f 9 GRCh38_GM12878_PET_active-rOCR.txt | sort -u | wc -l
# 40,520 active rOCRs


# 3. tag overlapped with rOCRs
intersectBed -a GRCh38_GM12878_loop_a_matched.bed -b ../GRCh38-rOCRs.bed \
-wa -wb > GRCh38_GM12878_PET_rOCR_a.txt
intersectBed -a GRCh38_GM12878_loop_b_matched.bed -b ../GRCh38-rOCRs.bed \
-wa -wb > GRCh38_GM12878_PET_rOCR_b.txt
#
cat GRCh38_GM12878_PET_rOCR_a.txt GRCh38_GM12878_PET_rOCR_b.txt > \
GRCh38_GM12878_PET_rOCR.txt
cut -f 5 GRCh38_GM12878_PET_rOCR.txt | sort -u | wc -l
# 92,407 loops
cut -f 4-5 GRCh38_GM12878_PET_rOCR.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l
# 86,147 loops both ends overlapped
cut -f 9 GRCh38_GM12878_PET_rOCR.txt | sort -u | wc -l
# 152,534 rOCRs

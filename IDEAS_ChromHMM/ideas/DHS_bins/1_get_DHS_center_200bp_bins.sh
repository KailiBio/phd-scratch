#!/bin/bash/

# -- Kaili
# This script is for getting DHS center 200-bp bins.
# 1. check mm10 DHS sizes and DHS-gap sizes
# 2. version 1: 100-300bp bins
# 3. version 2: all regions

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/"

cd ${workDir}


# 1. check mm10 DHS sizes and DHS-gap sizes
## get genomic region Yu used.
awk '{FS=" ";OFS="\t"}{print $1,$2,$3}' /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed | sort -k1,1 -k2,2n > mm10.bed
bedtools merge -i mm10.bed > mm10_chrom.bed
grep -v "chrM" mm10_chrom.bed > mm10_chrom_clean.bed
# get gap region
sort -k1,1 -k2,2n /data/projects/psychencode/Registry/V1/mm10/mm10-rOCRs.bed | cut -f 1-3 > mm10-rOCRs.bed
bedtools subtract -a mm10_chrom_clean.bed -b mm10-rOCRs.bed > mm10_rOCR_gap.bed

# Rscript mm10_DHS_length_basic.R


# 2. version 1: 100-300bp bins
python ${scriptDir}DHS-center_bins_v1.py /data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed /data/zusers/fankaili/ideas/dhs_bins/mm10_rOCR_gap.bed \
/data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v1.bed


# 3. version 2: all regions
python ${scriptDir}DHS-center_bins_v2.py /data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed /data/zusers/fankaili/ideas/dhs_bins/mm10_rOCR_gap.bed \
/data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v2.bed

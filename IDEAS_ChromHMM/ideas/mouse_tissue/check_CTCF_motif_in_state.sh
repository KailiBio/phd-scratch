#!/bin/bash

# -- Kaili
# This script is for checking CTCF motif in CTCF states.

cd /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/

# 1. get the CTCF motif region
scanMotifGenomeWide.pl mm10_ctcf_motif_file.motif /home/fankaili/genome/mm10.fa -bed > mm10_CTCF_motif_region.bed

# 2. get CTCF state region in chr1
awk '{FS=" ";OFS="\t"}{s=0;for(i=5;i<16;i++){if($i==10 ||$i==36 || $i==21 || $i==28 || $i==40 || $i==37 || $i==24 || $i==22){s=1}};if(s){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' \
./ctcf_samples_result/ctcf_samples.chr1.state > mm10_chr1_ctcf_state.txt

# 3. get bins overlapped with CTCF motif
intersectBed -a /data/zusers/fankaili/ideas/run_ideas_p_value/mm10_tab.bed -b mm10_CTCF_motif_region.bed -wa -wb \
> mm10_CTCF_motif_bins.bed

# 4. CTCF state regions
intersectBed -a mm10_chr1_ctcf_state.txt -b mm10_CTCF_motif_region.bed -wa -wb \
> mm10_CTCF_motif_chr1_ctcf_state.bed

intersectBed -b mm10_chr1_ctcf_state.txt -a mm10_CTCF_motif_region.bed -wa -wb \
> mm10_CTCF_motif_chr1_ctcf_state_2.bed

############
awk '{FS=" ";OFS="\t"}{s=0;for(i=5;i<16;i++){if($i==10){s=1}};if(s){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' \
./ctcf_samples_result/ctcf_samples.chr1.state > mm10_chr1_ctcf_state10.txt
intersectBed -a mm10_chr1_ctcf_state10.txt -b mm10_CTCF_motif_region.bed -wa -wb \
> mm10_CTCF_motif_chr1_ctcf_state10.bed
intersectBed -a mm10_chr1_ctcf_state10.txt -b mm10_CTCF_motif_region.bed -v \
> mm10_CTCF_motif_chr1_ctcf_state10_no.bed

awk '{FS=" ";OFS="\t"}{s=0;for(i=5;i<16;i++){if($i==36){s=1}};if(s){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' \
./ctcf_samples_result/ctcf_samples.chr1.state > mm10_chr1_ctcf_state36.txt
intersectBed -a mm10_chr1_ctcf_state36.txt -b mm10_CTCF_motif_region.bed -wa -wb \
> mm10_CTCF_motif_chr1_ctcf_state36.bed
intersectBed -a mm10_chr1_ctcf_state36.txt -b mm10_CTCF_motif_region.bed -v \
> mm10_CTCF_motif_chr1_ctcf_state36_no.bed


awk '{FS=" ";OFS="\t"}{s=0;for(i=5;i<16;i++){if($i==22){s=1}};if(s){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' \
./ctcf_samples_result/ctcf_samples.chr1.state > mm10_chr1_ctcf_state36.txt
intersectBed -a mm10_chr1_ctcf_state36.txt -b mm10_CTCF_motif_region.bed -wa -wb \
> mm10_CTCF_motif_chr1_ctcf_state36.bed
wc -l mm10_chr1_ctcf_state36.txt
wc -l mm10_CTCF_motif_chr1_ctcf_state36.bed

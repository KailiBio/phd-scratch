#!/bin/bash

# -- Kaili
# This script is for analyzing CTCF states along with CTCF peaks&motifs.
# 1. get CTCF singal in CTCF peaks.
# 2. pick CTCF peaks with signal over 50, intersect with CTCF motif
# 3. CTCF state in chr1 overlapped with CTCF high-signal peaks

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 0.
bash ${scriptDir}make_CTCF_track_hub.sh
bash ${scriptDir} check_CTCF_motif_in_state.#!/bin/sh

######
##eg: liver_14.5_day

# 1. get CTCF singal in CTCF peaks.
bigBedToBed /data/projects/encode/data/ENCSR397RHW/ENCFF368MGB.bigBed liver_14.5_day_CTCF_peak.bed
sort -k1,1 -k2,2n liver_14.5_day_CTCF_peak.bed | awk '{FS=OFS="\t"}{print $1,$2,$3,"peak_"NR,1,$6}' \
> liver_14.5_day_CTCF_peak_sorted.bed
bigWigAverageOverBed /data/projects/encode/data/ENCSR397RHW/ENCFF191NHQ.bigWig \
liver_14.5_day_CTCF_peak_sorted.bed liver_14.5_day_CTCF_peak_signal.txt
# Rscript figs_statistic_CTCF.R

# get a cut-off 50
# 2. pick CTCF peaks with signal over 50, intersect with CTCF motif
awk '{FS=OFS="\t"}{if(NR==FNR && $5>50){a[$1]=1}else{if(a[$4]){print $0}}}' liver_14.5_day_CTCF_peak_signal.txt \
liver_14.5_day_CTCF_peak_sorted.bed > liver_14.5_day_CTCF_peak_over50.bed
intersectBed -a liver_14.5_day_CTCF_peak_over50.bed -b mm10_CTCF_motif_region.bed -wa > \
liver_14.5_day_CTCF_peak_over50_withMotif.txt
intersectBed -a liver_14.5_day_CTCF_peak_sorted.bed -b mm10_CTCF_motif_region.bed -wa > \
liver_14.5_day_CTCF_peak_withMotif.txt


# 3. CTCF state in chr1 overlapped with CTCF high-signal peaks
awk '{if($1=="chr1"){print $0}}' liver_14.5_day_CTCF_peak_over50.bed > liver_14.5_day_CTCF_peak_over50_chr1.bed
awk '{if($1=="chr1"){print $0}}' liver_14.5_day_CTCF_peak_sorted.bed > liver_14.5_day_CTCF_peak_sorted_chr1.bed

intersectBed -a liver_14.5_day_CTCF_peak_over50_chr1.bed -b mm10_chr1_ctcf_state.txt -wa | sort -u | wc -l
intersectBed -a liver_14.5_day_CTCF_peak_sorted_chr1.bed -b mm10_chr1_ctcf_state.txt -wa | sort -u | wc -l
awk '{if($1=="chr1"){print $0}}' liver_14.5_day_CTCF_peak_over50.bed | wc -l

intersectBed -a mm10_chr1_ctcf_state.txt -b liver_14.5_day_CTCF_peak_sorted_chr1.bed -wa | sort -u | wc -l
intersectBed -a mm10_chr1_ctcf_state.txt -b liver_14.5_day_CTCF_peak_over50_chr1.bed -wa | sort -u | wc -l
awk '{FS=" ";OFS="\t"}{print $2,$3,$4,$5}' ./ctcf_samples_result/ctcf_samples.chr1.state | sort -k1,1 -k2,2n > mm10_chr1_state.txt
intersectBed -a mm10_chr1_state.txt -b liver_14.5_day_CTCF_peak_sorted_chr1.bed -wa | sort -u | wc -l
intersectBed -a mm10_chr1_state.txt -b liver_14.5_day_CTCF_peak_over50_chr1.bed -wa | sort -u | wc -l

# Rscript figs_statistic_CTCF.R

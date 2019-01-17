#!/bin/bash

# -- Kaili
# This script is for redo CTCF state analysis using CTCF state as unit.
# 1. Count num of states compare with CTCF peaks
# 2. compare three-group CTCF get_signal
# 3. CTCF signal of all the CTCF states (histogram)
# 4. 3-way venn diagram




scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 1. Count num of states compare with CTCF peaks
intersectBed -a mm10_ctcf_state_liver14.5_sorted.bed -b liver_14.5_day_CTCF_peak_sorted.bed -wa -f 0.5 -F 0.5 -e | sort -u | wc -l
intersectBed -a mm10_tab.bed -b liver_14.5_day_CTCF_peak_sorted.bed -wa -f 0.5 -F 0.5 -e | sort -u | wc -l
wc -l mm10_tab.bed
wc -l mm10_ctcf_state_liver14.5_sorted.bed

## peak length distribution
awk '{FS=OFS="\t"}{print $0,$3-$2}' liver_14.5_day_CTCF_peak_sorted.bed > liver_14.5_day_CTCF_peak_length.bed
# Rscript figs_statistic_CTCF.R


# 2. compare three-group CTCF get_signal
rep2="/data/projects/encode/data/ENCSR397RHW/ENCFF191NHQ.bigWig"
rep1_2="/data/projects/encode/data/ENCSR397RHW/ENCFF900YUX.bigWig"
#
intersectBed -a mm10_ctcf_state_liver14.5_sorted.bed -b liver_14.5_day_CTCF_peak_sorted.bed -wa -f 0.5 -F 0.5 -e | \
sort -u | sort -k1,1 -k2,2n > mm10_liver14.5_ctcf_state_ctcf_peak.bed
intersectBed -a mm10_ctcf_state_liver14.5_sorted.bed -b liver_14.5_day_CTCF_peak_sorted.bed -wa -f 0.5 -F 0.5 -e -v | \
sort -u | sort -k1,1 -k2,2n > mm10_liver14.5_ctcf_state_non_ctcf_peak.bed
intersectBed -a mm10_tab.bed -b mm10_ctcf_state_liver14.5_sorted.bed -wa -v | sort -u > mm10_non_ctcf_state_liver14.5_sorted.bed
intersectBed -a mm10_non_ctcf_state_liver14.5_sorted.bed -b liver_14.5_day_CTCF_peak_sorted.bed -wa -f 0.5 -F 0.5 -e | \
sort -u | sort -k1,1 -k2,2n > mm10_liver14.5_non_ctcf_state_ctcf_peak.bed
#
bigWigAverageOverBed ${rep1_2} mm10_liver14.5_ctcf_state_ctcf_peak.bed tmp.txt
awk '{FS=OFS="\t"}{print $1,$5,"CTCFstate_CTCFpeak"}' tmp.txt > liver_14.5_day_ctcf_peak_rep1.2_signal.txt
bigWigAverageOverBed ${rep1_2} mm10_liver14.5_ctcf_state_non_ctcf_peak.bed tmp.txt
awk '{FS=OFS="\t"}{print $1,$5,"CTCFstate_nonCTCFpeak"}' tmp.txt >> liver_14.5_day_ctcf_peak_rep1.2_signal.txt
bigWigAverageOverBed ${rep1_2} mm10_liver14.5_non_ctcf_state_ctcf_peak.bed tmp.txt
awk '{FS=OFS="\t"}{print $1,$5,"nonCTCFstate_CTCFpeak"}' tmp.txt >> liver_14.5_day_ctcf_peak_rep1.2_signal.txt
#
bigWigAverageOverBed ${rep2} mm10_liver14.5_ctcf_state_ctcf_peak.bed tmp.txt
awk '{FS=OFS="\t"}{print $1,$5,"CTCFstate_CTCFpeak"}' tmp.txt > liver_14.5_day_ctcf_peak_rep2_signal.txt
bigWigAverageOverBed ${rep2} mm10_liver14.5_ctcf_state_non_ctcf_peak.bed tmp.txt
awk '{FS=OFS="\t"}{print $1,$5,"CTCFstate_nonCTCFpeak"}' tmp.txt >> liver_14.5_day_ctcf_peak_rep2_signal.txt
bigWigAverageOverBed ${rep2} mm10_liver14.5_non_ctcf_state_ctcf_peak.bed tmp.txt
awk '{FS=OFS="\t"}{print $1,$5,"nonCTCFstate_CTCFpeak"}' tmp.txt >> liver_14.5_day_ctcf_peak_rep2_signal.txt

# Rscript figs_statistic_CTCF.R


# 3. CTCF signal of all the CTCF states (histogram)
bigWigAverageOverBed ${rep2} mm10_ctcf_state_liver14.5_sorted.bed tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$4]}}' tmp.txt mm10_ctcf_state_liver14.5_sorted.bed \
> mm10_ctcf_state_liver14.5_signal.txt
rm tmp.txt
# Rscript figs_statistic_CTCF.R


# 4. 3-way venn diagram
## liver 14.5
ctcf_state_file="mm10_ctcf_state_liver14.5_sorted.bed"
ctcf_peak_file="liver_14.5_day_CTCF_peak_sorted.bed"
ctcf_motif_file="mm10_CTCF_motif_region.bed"
non_ctcf_state_file="mm10_non_ctcf_state_liver14.5_sorted.bed"
#
intersectBed -a ${ctcf_state_file} -b ${ctcf_peak_file} -wa -f 0.5 -F 0.5 -e | sort -u > liver14.5_ctcf_state_peak.bed
intersectBed -a liver14.5_ctcf_state_peak.bed -b ${ctcf_motif_file} -wa | sort -u > liver14.5_ctcf_state_peak_motif.bed
wc -l liver14.5_ctcf_state_peak_motif.bed
wc -l liver14.5_ctcf_state_peak.bed
intersectBed -a ${ctcf_state_file} -b ${ctcf_motif_file} -wa | sort -u > liver14.5_ctcf_state_motif.bed
wc -l liver14.5_ctcf_state_motif.bed
intersectBed -a ${non_ctcf_state_file} -b ${ctcf_peak_file} -wa -f 0.5 -F 0.5 -e | sort -u > liver14.5_ctcf_non_state_peak.bed
intersectBed -a liver14.5_ctcf_non_state_peak.bed -b ${ctcf_motif_file} -wa | sort -u > liver14.5_ctcf_non_state_peak_motif.bed
wc -l liver14.5_ctcf_non_state_peak_motif.bed
wc -l liver14.5_ctcf_non_state_peak.bed
intersectBed -a ${non_ctcf_state_file} -b ${ctcf_motif_file} -wa | sort -u > liver14.5_ctcf_non_state_motif.bed
wc -l liver14.5_ctcf_non_state_motif.bed
#
wc -l ${ctcf_state_file}

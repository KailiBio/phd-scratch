#!/bin/bash

# -- Kaili
# This script is for analyzing CTCF states along with CTCF peaks&motifs.
# 1. get CTCF singal in CTCF peaks.
# 2. pick CTCF peaks with signal over 50, intersect with CTCF motif
# 3. CTCF state in chr1 overlapped with CTCF high-signal peaks
# 4. CTCF state in chr1 with CTCF motifs
# 5. running cut-off for CTCF signal in CTCF peaks

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 0.
bash ${scriptDir}make_CTCF_track_hub.sh
bash ${scriptDir} check_CTCF_motif_in_state.sh

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
intersectBed -a liver_14.5_day_CTCF_peak_over50.bed -b mm10_CTCF_motif_region.bed -wa | sort -u > \
liver_14.5_day_CTCF_peak_over50_withMotif.txt
intersectBed -a liver_14.5_day_CTCF_peak_sorted.bed -b mm10_CTCF_motif_region.bed -wa | sort -u > \
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



# 4. CTCF state in chr1 with CTCF motifs
## 1) contingency table
awk '{if($1=="chr1"){print $0}}' mm10_CTCF_motif_region.bed > mm10_CTCF_motif_region_chr1.bed
awk '{FS=" ";OFS="\t"}{if(NR>1){print $2,$3,$4,$1}}' \
/data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/run_IDEAS_8hm_atac_dname_pvalue.chr1.state \
| sort -u | sort -k1,1 -k2,2n > mm10_chr1_state.bed
intersectBed -a mm10_chr1_state.bed -b mm10_CTCF_motif_region_chr1.bed -wa | sort -u | wc -l
intersectBed -a mm10_chr1_ctcf_state.txt -b mm10_CTCF_motif_region_chr1.bed -wa | sort -u | wc -l

## 2) CTCF signal bettween CTCF states with CTCF motif & CTCF states without CTCF motif & non-CTCF states with CTCF motif
### CTCF states with CTCF motif
intersectBed -a mm10_chr1_ctcf_state.txt -b mm10_CTCF_motif_region_chr1.bed -wa | cut -f 1-4 | sort -u | \
sort -k1,1 -k2,2n > mm10_chr1_ctcf_state_ctcf_motif.bed
### CTCF states without CTCF motif
intersectBed -a mm10_chr1_ctcf_state.txt -b mm10_CTCF_motif_region_chr1.bed -wa -v | cut -f 1-4 | sort -u | \
sort -k1,1 -k2,2n > mm10_chr1_ctcf_state_non_ctcf_motif.bed
### non-CTCF states with CTCF motif
intersectBed -a mm10_chr1_state.bed -b mm10_chr1_ctcf_state.txt -wa -v > \
mm10_chr1_non_ctcf_state.bed
intersectBed -a mm10_chr1_non_ctcf_state.bed -b mm10_CTCF_motif_region_chr1.bed -wa | sort -u | \
sort -k1,1 -k2,2n> mm10_chr1_non_ctcf_state_ctcf_motif.bed
### get signal
bigWigAverageOverBed /data/projects/encode/data/ENCSR397RHW/ENCFF191NHQ.bigWig \
mm10_chr1_ctcf_state_ctcf_motif.bed liver_14.5_day_ctcf_state_ctcf_motif_signal.txt
bigWigAverageOverBed /data/projects/encode/data/ENCSR397RHW/ENCFF191NHQ.bigWig \
mm10_chr1_ctcf_state_non_ctcf_motif.bed liver_14.5_day_ctcf_state_non_ctcf_motif_signal.txt
bigWigAverageOverBed /data/projects/encode/data/ENCSR397RHW/ENCFF191NHQ.bigWig \
mm10_chr1_non_ctcf_state_ctcf_motif.bed liver_14.5_day_non_ctcf_state_ctcf_motif_signal.txt
#
awk '{FS=OFS="\t"}{print $1,$5,"CTCFstate_CTCFmotif"}' liver_14.5_day_ctcf_state_ctcf_motif_signal.txt \
> liver_14.5_day_chr1_ctcf_signal.txt
awk '{FS=OFS="\t"}{print $1,$5,"CTCFstate_nonCTCFmotif"}' liver_14.5_day_ctcf_state_non_ctcf_motif_signal.txt \
>> liver_14.5_day_chr1_ctcf_signal.txt
awk '{FS=OFS="\t"}{print $1,$5,"nonCTCFstate_CTCFmotif"}' liver_14.5_day_non_ctcf_state_ctcf_motif_signal.txt \
>> liver_14.5_day_chr1_ctcf_signal.txt
# Rscript figs_statistic_CTCF.R



# 5. running cut-off for CTCF signal in CTCF peaks
bash ${scriptDir}running_cutoff_CTCF_peak_signal.sh

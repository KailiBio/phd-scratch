#!/bin/bash/

# -- Kaili
# This script is for getting different CTCF singal cut-off.
# 0. get all CTCF state regions
# 1. count how many high CTCF peaks (running cut-off) with CTCF motif or CTCF states
# 2. count how many high CTCF states (running cut-off) are CTCF peaks
# 3. make histogram

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 0. get all CTCF state regions
if [ -f mm10_ctcf_state.bed ];then
    rm mm10_ctcf_state.bed
fi
#
for i in {1..19} X Y M
do
    awk '{FS=" ";OFS="\t"}{s=0;for(i=5;i<16;i++){if($i==28 ||$i==40 || $i==37 || $i==24 || $i==41 || $i==21 || $i==10 || $i==22 || $i==36){s=1}};if(s){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' \
    /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.chr${i}.state >> mm10_ctcf_state.bed
done
cut -f 1-4 mm10_ctcf_state.bed | sort -u | sort -k1,1 -k2,2n > mm10_ctcf_state_sorted.bed

## CTCF state in liver_14.5
if [ -f mm10_ctcf_state_liver14.5.bed ];then
    rm mm10_ctcf_state_liver14.5.bed
fi
#
for i in {1..19} X Y M
do
    awk '{FS=" ";OFS="\t"}{s=0;if($6==28 ||$6==40 || $6==37 || $6==24 || $6==41 || $6==21 || $6==10 || $6==22 || $6==36){s=1};if(s){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' \
    /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.chr${i}.state >> mm10_ctcf_state_liver14.5.bed
done
cut -f 1-4 mm10_ctcf_state_liver14.5.bed | sort -u | sort -k1,1 -k2,2n > mm10_ctcf_state_liver14.5_sorted.bed


# 1. count how many high CTCF peaks (running cut-off) with CTCF motif or CTCF states
if [ -f mm10_ctcf_peak_cutoff_with_ctcf_motif_num.txt ]; then
    rm mm10_ctcf_peak_cutoff_with_ctcf_motif_num.txt
fi
#
if [ -f mm10_ctcf_peak_cutoff_with_ctcf_state_num.txt ]; then
    rm mm10_ctcf_peak_cutoff_with_ctcf_state_num.txt
fi
#
for i in {0..247}
do
    echo ${i}
    awk -v cutoff="$i" '{FS=OFS="\t"}{if(NR==FNR && $5>cutoff){a[$1]=1}else{if(a[$4]){print $0}}}' liver_14.5_day_CTCF_peak_signal.txt \
    liver_14.5_day_CTCF_peak_sorted.bed > tmp_peak.bed
    num=`wc -l tmp_peak.bed | awk '{print $1}'`
    # high CTCF peaks contains motif
    a=`intersectBed -a tmp_peak.bed -b mm10_CTCF_motif_region.bed -wa | sort -u | wc -l`
    echo -e ${i}"\t"${a}"\t"${num} >> mm10_ctcf_peak_cutoff_with_ctcf_motif_num.txt
    # high CTCF peaks in CTCF states
    b=`intersectBed -a tmp_peak.bed -b mm10_ctcf_state_liver14.5_sorted.bed -wa | sort -u | wc -l`
    echo -e ${i}"\t"${b}"\t"${num} >>  mm10_ctcf_peak_cutoff_with_ctcf_state_num.txt
done
rm tmp_peak.bed

# 2. count how many high CTCF states (running cut-off) are CTCF peaks
if [ -f mm10_ctcf_state_cutoff_with_ctcf_peak_num.txt ]; then
    rm mm10_ctcf_state_cutoff_with_ctcf_peak_num.txt
fi
#
for i in {0..243}
do
    echo ${i}
    awk -v cutoff="$i" '{FS=OFS="\t"}{if($5>cutoff){print $0}}' mm10_ctcf_state_liver14.5_signal.txt > tmp_state.bed
    num=`wc -l tmp_state.bed | awk '{print $1}'`
    # high CTCF states are CTCF peaks
    a=`intersectBed -a tmp_state.bed -b liver_14.5_day_CTCF_peak_sorted.bed -wa -f 0.5 -F 0.5 -e | sort -u | wc -l`
    echo -e ${i}"\t"${a}"\t"${num} >> mm10_ctcf_state_cutoff_with_ctcf_peak_num.txt
done
rm tmp_state.bed


# 3. make histogram
# Rscript figs_statistic_CTCF.R

#!/bin/bash

# -- Kaili
# This script is for comparing ratio of CTCF peaks in CTCF states between 9to11 & 9impute11.

# 1. get CTCF peak in all tissues, get peak signal
# 2. calculate number of peaks in each states
# 3. make figures

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

mkdir peaks_validation_9To11
# 1. get CTCF peak in all tissues, get peak signal
while read line
do
    sample=`awk '{split($1,a,"_");print a[1]"_"a[2]}' <<< $line`
    expID=`awk '{print $4}' <<< $line`
    signal_file=`awk '{print $5}' <<< $line`
    peak_file=`awk '{print $6}' <<< $line`
    #
    bigBedToBed /data/projects/encode/data/${expID}/${peak_file}.bigBed ./peaks_validation_9To11/${sample}_ctcf_peak.bed
    awk '{FS=OFS="\t"}{print $1,$2,$3,"peak_"NR}' ./peaks_validation_9To11/${sample}_ctcf_peak.bed | sort -u | sort -k1,1 -k2,2n > ./peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed
    bigWigAverageOverBed /data/projects/encode/data/${expID}/${signal_file}.bigWig ./peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed ./peaks_validation_9To11/${sample}_ctcf_peak_sorted.tab
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$4]}}' ./peaks_validation_9To11/${sample}_ctcf_peak_sorted.tab ./peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed > ./peaks_validation_9To11/${sample}_ctcf_peak_withSignal.bed
    awk '{FS=OFS="\t"}{if(NR==FNR){a=int(($3-$2)/2);print $1,$2+a,$2+a,$4,$5,"+",$2,$3}}' ./peaks_validation_9To11/${sample}_ctcf_peak_withSignal.bed | sort -k1,1 -k2,2n > ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed
done < /data/zusers/fankaili/ideas/mm10_tissue_used_list_CTCF_peak_list.txt

# 2. calculate number of peaks in each states
## 1) 9impute11
mkdir state_bed_9impute11
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_impute_11sample_result/ ctcf_9sample_impute_11sample. /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_9impute11/ 11
#
while read sample
do
    echo ${sample}
    #
    intersectBed -a ./state_bed_9impute11/${sample}_state_sorted.bed -b ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
    # percentage of peak in each states
    num_peak=`wc -l ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
    cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peaks_validation_9To11/peak-ratio-in-state_9impute11_${sample}.txt
    # percentage of states have peaks
    cut -f 5 ./state_bed_9impute11/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peaks_validation_9To11/peak-ratio-in-state_9impute11_${sample}.txt > ./peaks_validation_9To11/bin-ratio-have-peak_9impute11_${sample}.txt
done < all_ctcf_sample.txt
## 2) 9to11
mkdir state_bed_9to11
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_to_11sample_result/ ctcf_9sample_to_11sample. /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_9to11/ 11
#
while read sample
do
    echo ${sample}
    #
    intersectBed -a ./state_bed_9to11/${sample}_state_sorted.bed -b ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
    # percentage of peak in each states
    num_peak=`wc -l ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
    cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peaks_validation_9To11/peak-ratio-in-state_9to11_${sample}.txt
    # percentage of states have peaks
    cut -f 5 ./state_bed_9to11/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peaks_validation_9To11/peak-ratio-in-state_9to11_${sample}.txt > ./peaks_validation_9To11/bin-ratio-have-peak_9to11_${sample}.txt
done < all_ctcf_sample.txt
#
rm tmp.bed tmp.txt

### heart 0 do not have any CTCF peaks in state28, mannully add one line there.

# 3. make figures
# Rscript ${scriptDir}make_barplot_percentage_state_peak.R

# 4. calculate percentage

# 5. percentage in 9impute66
for sample in liver_14.5 lung_14.5
do
    echo ${sample}
    #
    intersectBed -a ./state_bed/${sample}_state_sorted.bed -b ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
    # percentage of peak in each states
    num_peak=`wc -l ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
    cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peaks_validation_9To11/peak-ratio-in-state_9impute66_${sample}.txt
    # percentage of states have peaks
    cut -f 5 ./state_bed_9impute11/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peaks_validation_9To11/peak-ratio-in-state_9impute11_${sample}.txt > ./peaks_validation_9To11/bin-ratio-have-peak_9impute66_${sample}.txt
done

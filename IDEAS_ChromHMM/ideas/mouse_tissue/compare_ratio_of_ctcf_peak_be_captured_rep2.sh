#!/bin/bash

# -- Kaili
# This script is for comparing ratio of CTCF peaks in CTCF states between rep1 & rep1-impute-rep2.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

ctcf_peak_folder="/data/zusers/fankaili/ideas/dhs_ctcf/peaks_validation_9To11/"

cd ${workDir}

# 1. rep1
mkdir state_bed_rep1
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ ctcf_samples. /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_bed_rep1/ 11
#
mkdir peak_validation_rep1-2
while read sample
do
    echo ${sample}
    #
    intersectBed -a ./state_bed_rep1/${sample}_state_sorted.bed -b ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
    # percentage of peak in each states
    num_peak=`wc -l ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
    cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peak_validation_rep1-2/peak-ratio-in-state_rep1_${sample}.txt
    # percentage of states have peaks
    cut -f 5 ./state_bed_rep1/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peak_validation_rep1-2/peak-ratio-in-state_rep1_${sample}.txt > ./peak_validation_rep1-2/bin-ratio-have-peak_rep1_${sample}.txt
done < CTCF_sample_list.txt

# 2. rep1-impute-rep2
mkdir state_bed_rep1-impute-rep2
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_rep1-impute-rep2_result/ ctcf_samples_rep1-impute-rep2. /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_bed_rep1-impute-rep2/ 11
#
while read sample
do
    echo ${sample}
    #
    intersectBed -a ./state_bed_rep1-impute-rep2/${sample}_state_sorted.bed -b ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
    # percentage of peak in each states
    num_peak=`wc -l ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
    cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep2_${sample}.txt
    # percentage of states have peaks
    cut -f 5 ./state_bed_rep1-impute-rep2/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep2_${sample}.txt > ./peak_validation_rep1-2/bin-ratio-have-peak_rep1-impute-rep2_${sample}.txt
done < CTCF_sample_list.txt

# 3. rep1-impute-rep1
mkdir state_bed_rep1-impute-rep1
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_rep1-impute-rep1_result/ ctcf_samples_rep1-impute-rep1. /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_bed_rep1-impute-rep1/ 11
#
while read sample
do
    echo ${sample}
    #
    intersectBed -a ./state_bed_rep1-impute-rep1/${sample}_state_sorted.bed -b ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
    # percentage of peak in each states
    num_peak=`wc -l ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
    cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep1_${sample}.txt
    # percentage of states have peaks
    cut -f 5 ./state_bed_rep1-impute-rep1/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep1_${sample}.txt > ./peak_validation_rep1-2/bin-ratio-have-peak_rep1-impute-rep1_${sample}.txt
done < CTCF_sample_list.txt

# 4. rep1-impute-rep1 lung_14.5
mkdir state_bed_rep1-impute-rep1_lung_14.5
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_rep1-impute-rep1_oneSample_result/ ctcf_samples_rep1-impute-rep1_oneSample. /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_bed_rep1-impute-rep1_lung_14.5/ 1
#
sample="lung_14.5"
intersectBed -a ./state_bed_rep1-impute-rep1_lung_14.5/${sample}_state_sorted.bed -b ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
# percentage of peak in each states
num_peak=`wc -l ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep1_oneSample_${sample}.txt
# percentage of states have peaks
cut -f 5 ./state_bed_rep1-impute-rep1/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep1_oneSample_${sample}.txt > ./peak_validation_rep1-2/bin-ratio-have-peak_rep1-impute-rep1_oneSample_${sample}.txt

# 5. rep1-impute-rep2 lung_14.5
mkdir state_bed_rep1-impute-rep2_lung_14.5
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_rep1-impute-rep2_oneSample_result/ ctcf_samples_rep1-impute-rep2_oneSample. /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_bed_rep1-impute-rep2_lung_14.5/ 1
#
sample="lung_14.5"
intersectBed -a ./state_bed_rep1-impute-rep2_lung_14.5/${sample}_state_sorted.bed -b ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
# percentage of peak in each states
num_peak=`wc -l ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep2_oneSample_${sample}.txt
# percentage of states have peaks
cut -f 5 ./state_bed_rep1-impute-rep2/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep2_oneSample_${sample}.txt > ./peak_validation_rep1-2/bin-ratio-have-peak_rep1-impute-rep2_oneSample_${sample}.txt

# 6. make figures
# Rscript ${scriptDir}make_barplot_percentage_state_peak.R

# 7. rep1-impute-rep1 given9samples
mkdir state_bed_rep1-impute-rep1_given9samples
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_rep1-impute-rep1_given9samples_result/ ctcf_samples_rep1-impute-rep1_given9samples. /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_bed_rep1-impute-rep1_given9samples/ 11
#
while read sample
do
    echo ${sample}
    #
    intersectBed -a ./state_bed_rep1-impute-rep1_given9samples/${sample}_state_sorted.bed -b ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
    # percentage of peak in each states
    num_peak=`wc -l ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
    cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep1_given9samples_${sample}.txt
    # percentage of states have peaks
    cut -f 5 ./state_bed_rep1-impute-rep1_given9samples/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep1_given9samples_${sample}.txt > ./peak_validation_rep1-2/bin-ratio-have-peak_rep1-impute-rep1_given9samples_${sample}.txt
done < CTCF_sample_list.txt

# 8. rep1-impute-rep2 given9samples
mkdir state_bed_rep1-impute-rep2_given9samples
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_rep1-impute-rep2_given9samples_result/ ctcf_samples_rep1-impute-rep2_given9samples. /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_bed_rep1-impute-rep2_given9samples/ 11
#
while read sample
do
    echo ${sample}
    #
    intersectBed -a ./state_bed_rep1-impute-rep2_given9samples/${sample}_state_sorted.bed -b ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
    # percentage of peak in each states
    num_peak=`wc -l ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
    cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep2_given9samples_${sample}.txt
    # percentage of states have peaks
    cut -f 5 ./state_bed_rep1-impute-rep2_given9samples/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep2_given9samples_${sample}.txt > ./peak_validation_rep1-2/bin-ratio-have-peak_rep1-impute-rep2_given9samples_${sample}.txt
done < CTCF_sample_list.txt

# 9. rep1-impute-rep2 given9samples, 10k
mkdir state_bed_rep1-impute-rep2_given9samples_10k
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_rep1-impute-rep2_given9samples_10k_result/ ctcf_samples_rep1-impute-rep2_given9samples_10k. /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_bed_rep1-impute-rep2_given9samples_10k/ 11
#
while read sample
do
    echo ${sample}
    #
    intersectBed -a ./state_bed_rep1-impute-rep2_given9samples_10k/${sample}_state_sorted.bed -b ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed -wa -wb > tmp.bed
    # percentage of peak in each states
    num_peak=`wc -l ${ctcf_peak_folder}${sample}_ctcf_peak_center_withSignal.bed | awk '{print $1}'`
    cut -f 5,9 tmp.bed | sort -u | cut -f 1 | sort | uniq -c | awk -v num="$num_peak" '{OFS="\t"}{print $2,$1,$1/num}' | sort -k1,1n > ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep2_given9samples_10k_${sample}.txt
    # percentage of states have peaks
    cut -f 5 ./state_bed_rep1-impute-rep2_given9samples_10k/${sample}_state_sorted.bed | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k1,1n > tmp.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$2/a[$1]}}' tmp.txt ./peak_validation_rep1-2/peak-ratio-in-state_rep1-impute-rep2_given9samples_10k_${sample}.txt > ./peak_validation_rep1-2/bin-ratio-have-peak_rep1-impute-rep2_given9samples_10k_${sample}.txt
done < CTCF_sample_list.txt

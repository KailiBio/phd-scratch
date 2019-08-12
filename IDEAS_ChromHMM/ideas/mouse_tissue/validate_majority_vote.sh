#!/bin/bash

# -- Kaili
# This script is for valuing imputation using majority vote.
# 1. count voted times for each potential CTCF bins.
# 2. recall ratio
# 3. Venn: 9to11, 9impute11, majority-vote
# 4. accumulated bar of peak-recall for counted-bins.
# 5. check the recall of liver14.5/lung14.5 CTCF states

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# Rscript make_figures_for_majority_vote.R

# 1. count voted times for each potential CTCF bins.
# CTCF states: 12, 24, 26, 27, 30, 34, 37, 38, 39, 41, 45, 46
## 1) get state bed file
mkdir state_bed_9sample
#
stateDir="/data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/"
prefix="dhs_ctcf."
stateBedDir="/data/zusers/fankaili/ideas/dhs_ctcf/state_bed_9sample/"
nohup bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${stateDir} ${prefix} ${stateBedDir} 9 > /data/zusers/fankaili/ideas/dhs_ctcf/nohup/nohup.make_each_sample_state_bed.out 2>&1&
## 2) counting voted times
if [ -f ./state_bed_9sample/tmp.bed ]; then rm ./state_bed_9sample/tmp.bed; fi
for file in `ls ./state_bed_9sample/*_state_sorted.bed`
do
    sample=`awk '{split($1,a,"/");split(a[3],b,"_");print b[1]"_"b[2]}' <<< ${file}`
    echo ${sample}
    #
    awk -v sample="$sample" '{FS=OFS="\t"}{if($5==12 || $5==24 || $5==26 || $5==27 || $5==30 || $5==34 || $5==37 || $5==38 || $5==39 || $5==41 || $5==45 || $5==46){print $0,sample}}'  ${file} > ./state_bed_9sample/${sample}_CTCFstates.bed
    #
    cut -f 1-4 ./state_bed_9sample/${sample}_CTCFstates.bed >> ./state_bed_9sample/tmp.bed
done
sort ./state_bed_9sample/tmp.bed | uniq -c | awk '{OFS="\t"}{print $2,$3,$4,$5,$1}' | sort -k5,5nr -k1,1 -k2,2n > ./state_bed_9sample/CTCFstates_count.txt
cut -f 5 ./state_bed_9sample/CTCFstates_count.txt | sort | uniq -c
### get voted CTCF states
awk '{FS=OFS="\t"}{if($5>=5){print $0}}' ./state_bed_9sample/CTCFstates_count.txt > ./state_bed_9sample/CTCFstates_voted.txt

# 2. recall ratio
# bash compare_ratio_of_ctcf_peak_be_captured.sh

# 3. Venn: 9to11, 9impute11, majority-vote
## 1) liver_14.5
sample="liver_14.5"
awk '{FS=OFS="\t"}{if($5==12 || $5==24 || $5==26 || $5==27 || $5==30 || $5==34 || $5==37 || $5==38 || $5==39 || $5==41 || $5==45 || $5==46){print $0}}' ./state_bed_9to11/${sample}_state_sorted.bed > ./peaks_validation_9To11/${sample}_9to11_CTCFstates.txt
intersectBed -a ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -b ./peaks_validation_9To11/${sample}_9to11_CTCFstates.txt -wa -u | cut -f 4 > ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt
#
awk '{FS=OFS="\t"}{if($5==12 || $5==24 || $5==26 || $5==27 || $5==30 || $5==34 || $5==37 || $5==38 || $5==39 || $5==41 || $5==45 || $5==46){print $0}}' ./state_bed_9impute11/${sample}_state_sorted.bed > ./peaks_validation_9To11/${sample}_9impute11_CTCFstates.txt
intersectBed -a ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -b ./peaks_validation_9To11/${sample}_9impute11_CTCFstates.txt -wa -u | cut -f 4 > ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt
#
intersectBed -a ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -b ./state_bed_9sample/CTCFstates_voted.txt -wa -u | cut -f 4 > ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt
# get Venn
wc -l ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt
wc -l ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt
wc -l ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt
#
cat ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt | sort | uniq -d | wc -l
cat ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt | sort | uniq -d | wc -l
cat ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt | sort | uniq -d | wc -l
cat ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt | sort | uniq -c | awk '{if($1==3){print $0}}' | wc -l
wc -l ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed
cat ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt | sort -u | wc -l
## 2) lung_14.5
sample="lung_14.5"
awk '{FS=OFS="\t"}{if($5==12 || $5==24 || $5==26 || $5==27 || $5==30 || $5==34 || $5==37 || $5==38 || $5==39 || $5==41 || $5==45 || $5==46){print $0}}' ./state_bed_9to11/${sample}_state_sorted.bed > ./peaks_validation_9To11/${sample}_9to11_CTCFstates.txt
intersectBed -a ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -b ./peaks_validation_9To11/${sample}_9to11_CTCFstates.txt -wa -u | cut -f 4 > ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt
#
awk '{FS=OFS="\t"}{if($5==12 || $5==24 || $5==26 || $5==27 || $5==30 || $5==34 || $5==37 || $5==38 || $5==39 || $5==41 || $5==45 || $5==46){print $0}}' ./state_bed_9impute11/${sample}_state_sorted.bed > ./peaks_validation_9To11/${sample}_9impute11_CTCFstates.txt
intersectBed -a ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -b ./peaks_validation_9To11/${sample}_9impute11_CTCFstates.txt -wa -u | cut -f 4 > ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt
#
intersectBed -a ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -b ./state_bed_9sample/CTCFstates_voted.txt -wa -u | cut -f 4 > ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt
# get Venn
wc -l ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt
wc -l ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt
wc -l ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt
#
cat ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt | sort | uniq -d | wc -l
cat ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt | sort | uniq -d | wc -l
cat ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt | sort | uniq -d | wc -l
cat ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt | sort | uniq -c | awk '{if($1==3){print $0}}' | wc -l
wc -l ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed
cat ./peaks_validation_9To11/${sample}_9to11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_9impute11_recalled_peaks.txt ./peaks_validation_9To11/${sample}_voted_recalled_peaks.txt | sort -u | wc -l

# 4. accumulated bar of peak-recall for counted-bins.
## 1) liver_14.5
sample="liver_14.5"
for i in {1..9}
do
    echo $i
    #
    awk -v i="$i" '{if($5==i){print $0}}' ./state_bed_9sample/CTCFstates_count.txt > tmp.voted_${i}.txt
    intersectBed -a ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -b tmp.voted_${i}.txt -wa -u | cut -f 4 | wc -l
done
wc -l ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed

for i in {1..9}
do
    echo $i
    #
    awk -v i="$i" '{if($5==i){print $0}}' ./state_bed_9sample/CTCFstates_count.txt > tmp.voted_${i}.txt
    intersectBed -a tmp.voted_${i}.txt -b ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -wa -u | cut -f 4 | wc -l
done

for i in {1..9}
do
    echo $i
    #
    awk -v i="$i" '{if($5==i){print $0}}' ./state_bed_9sample/CTCFstates_count.txt > tmp.voted_${i}.txt
    intersectBed -a tmp.voted_${i}.txt -b ./peaks_validation_9To11/${sample}_ctcf_peak.bed -wa -u | cut -f 4 | wc -l
done


## 2) lung_14.5
sample="lung_14.5"
for i in {1..9}
do
    echo $i
    #
    awk -v i="$i" '{if($5<=i){print $0}}' ./state_bed_9sample/CTCFstates_count.txt > tmp.voted_${i}.txt
    intersectBed -a ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -b tmp.voted_${i}.txt -wa -u | cut -f 4 | wc -l
done
wc -l ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed

for i in {1..9}
do
    echo $i
    #
    awk -v i="$i" '{if($5==i){print $0}}' ./state_bed_9sample/CTCFstates_count.txt > tmp.voted_${i}.txt
    intersectBed -a tmp.voted_${i}.txt -b ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -wa -u | cut -f 4 | wc -l
done

for i in {1..9}
do
    echo $i
    #
    awk -v i="$i" '{if($5==i){print $0}}' ./state_bed_9sample/CTCFstates_count.txt > tmp.voted_${i}.txt
    intersectBed -a tmp.voted_${i}.txt -b ./peaks_validation_9To11/${sample}_ctcf_peak.bed -wa -u | cut -f 4 | wc -l
done

rm tmp.voted_*.txt

# 5. check the recall of liver14.5/lung14.5 CTCF states
### havn't finished. Maybe we don't need this, based on Venn.
if [ -f ./state_bed_9to11/tmp.bed ]; then rm ./state_bed_9to11/tmp.bed; fi
for file in `ls ./state_bed_9to11/*_state_sorted.bed`
do
    sample=`awk '{split($1,a,"/");split(a[3],b,"_");print b[1]"_"b[2]}' <<< ${file}`
    echo ${sample}
    #
    awk -v sample="$sample" '{FS=OFS="\t"}{if($5==12 || $5==24 || $5==26 || $5==27 || $5==30 || $5==34 || $5==37 || $5==38 || $5==39 || $5==41 || $5==45 || $5==46){print $0,sample}}'  ${file} > ./state_bed_9to11/${sample}_CTCFstates.bed
    #
    cut -f 1-4 ./state_bed_9to11/${sample}_CTCFstates.bed >> ./state_bed_9to11/tmp.bed
done
sort ./state_bed_9to11/tmp.bed | uniq -c | awk '{OFS="\t"}{print $2,$3,$4,$5,$1}' | sort -k5,5nr -k1,1 -k2,2n > ./state_bed_9to11/CTCFstates_count.txt
#
awk '{if($5==1){print $0}}' ./state_bed_9to11/CTCFstates_count.txt > ./state_bed_9to11/CTCFstates_specific.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0}}}' ./state_bed_9to11/liver_14.5_CTCFstates.bed ./state_bed_9to11/CTCFstates_specific.txt > ./state_bed_9to11/CTCFstates_liver_14.5_specific.txt



awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0}}}' ./state_bed_9to11/liver_0_CTCFstates.bed ./state_bed_9to11/CTCFstates_specific.txt | wc -l

#!/bin/bash

# -- Kaili
# This script is for analyzing tissue-specific CTCF sites.

# 1. basic info about CTCF peaks
# 2. find tissue-specific peaks (P0)
# 3. use tissue-specific peaks as pos, non-peak as neg, calculate AUPR.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# 1. basic info about CTCF peaks
# Rscript make_CTCF_peak_length_histogram.R

# 2. find tissue-specific peaks (P0)
## 1) peak exclusive
for sample1 in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 lung_0 stomach_0 liver_14.5 lung_14.5
do
    for sample2 in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 lung_0  stomach_0 liver_14.5 lung_14.5
    do
        if [ "$sample1" != "$sample2" ];
        then
            echo $sample1" VS "$sample2;
            intersectBed -a ./peaks_validation_9To11/${sample1}_ctcf_peak_sorted.bed -b ./peaks_validation_9To11/${sample2}_ctcf_peak_sorted.bed -v > ./tissue_specific_peak/${sample1}_minus_${sample2}_CTCFpeaks.txt
        fi
    done
done
## 2) count specificity
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 lung_0 stomach_0 liver_14.5 lung_14.5
do
    echo $sample
    #
    cat ./tissue_specific_peak/${sample}_minus_*_CTCFpeaks.txt | cut -f 4 | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k2nr > ./tissue_specific_peak/${sample}_specific_11CTCFpeaks_count.txt
    awk '{if($2==8){print $1}}' ./tissue_specific_peak/${sample}_specific_11CTCFpeaks_count.txt > ./tissue_specific_peak/${sample}_specific_11CTCFpeaks.txt
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $4,11-a[$4]}}' ./tissue_specific_peak/${sample}_specific_11CTCFpeaks_count.txt ./peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed > ./tissue_specific_peak/${sample}_specific_11CTCFpeaks_count_all.txt
done
##### Aug12
# tissue-specific Only
for sample in liver_14.5 lung_14.5
do
    echo $sample
    #
    cat ./tissue_specific_peak/${sample}_minus_*_CTCFpeaks.txt | cut -f 4 | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' | sort -k2nr > ./tissue_specific_peak/${sample}_specific_10CTCFpeaks_count.txt
    awk '{if($2==8){print $1}}' ./tissue_specific_peak/${sample}_specific_10CTCFpeaks_count.txt > ./tissue_specific_peak/${sample}_specific_10CTCFpeaks.txt
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $4,10-a[$4]}}' ./tissue_specific_peak/${sample}_specific_10CTCFpeaks_count.txt ./peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed > ./tissue_specific_peak/${sample}_specific_10CTCFpeaks_count_all.txt
done

## 3) make histogram for each tissue
# Rscript make_CTCFpeak_consistency_barplot.R


# 3. use tissue-specific peaks as pos, non-peak as neg, calculate AUPR.
bash ${scriptDir}calculate_AUPR_specificCTCF.sh



# 4. considering peak in the same tissue
# use tissue-specific peaks

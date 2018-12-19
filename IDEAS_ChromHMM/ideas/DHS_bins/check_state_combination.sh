#!/bin/bash

# -- Kaili
# This is script for checking the state combination for two closest OCR bins.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/"

cd ${workDir}

# 1. get close OCR bins (<100bp), and matched normal bins.
## 1) get close OCR bins (<100bp)
sort -k1,1 -k2,2n mm10_OCR-center_bins_v3_signal_based_OCR.bed > mm10_OCR-center_bins_v3_signal_based_OCR_sorted.bed
closestBed -a mm10_OCR-center_bins_v3_signal_based_OCR_sorted.bed -b mm10_OCR-center_bins_v3_signal_based_OCR_sorted.bed \
-d -io > mm10_OCR-center_bins_v3_signal_based_OCR_distance.txt
awk '{FS=OFS="\t"}{if($9<100){print $0}}' mm10_OCR-center_bins_v3_signal_based_OCR_distance.txt | \
awk '{FS=OFS="\t"}{if(NR%2==1){print $0}}'> mm10_OCR-center_bins_v3_signal_based_OCR_nearest.txt
# 147,632

## 2) get close OCR bins matched normal bins
awk '{FS=OFS="\t"}{center1=int($2+($3-$2)/2);center2=int($6+($7-$6)/2);print $1,center1,center1,$4"\n"$5,center2,center2,$8}' \
mm10_OCR-center_bins_v3_signal_based_OCR_nearest.txt | sort -u > mm10_OCR-center_bins_v3_signal_based_OCR_nearest_center.bed
intersectBed -a /data/zusers/fankaili/ideas/run_ideas_p_value/mm10_tab.bed -b \
mm10_OCR-center_bins_v3_signal_based_OCR_nearest_center.bed -wa -wb | sort -u | sort -k1,1 -k2,2n \
> mm10_nearest_OCR-center_matched_normal.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$8]=$0}else{print a[$4],a[$8]}}' \
mm10_nearest_OCR-center_matched_normal.bed mm10_OCR-center_bins_v3_signal_based_OCR_nearest.txt > \
mm10_OCR-center_bins_OCR_nearest_pairs.txt
awk '{FS=OFS="\t"}{if($4!="" && $8!="" && $12!="" && $16!=""){print $8,$16,$4,$12}}' \
mm10_OCR-center_bins_OCR_nearest_pairs.txt > OCR-center_bins_nearest_pairs.txt
# 147,628

# 2. get the states for bins
## 1) close OCR bins
# get states in intestine_16.5
dhs_state_path="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp_result/"
if [ -f  dhs-bins_state_intestine_16.5.txt ];then
    rm  dhs-bins_state_intestine_16.5.txt
fi
for i in {1..19} X Y
do
    awk '{FS=" ";OFS="\t"}{if(NR>1){print $1,$5}}' ${dhs_state_path}DHS_v3_100-400bp.chr${i}.state >> dhs-bins_state_intestine_16.5.txt
done
# get states
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,a[$1],$2,a[$2],$3,$4}}' dhs-bins_state_intestine_16.5.txt \
OCR-center_bins_nearest_pairs.txt > OCR-center_bins_nearest_pairs_tmp.txt

## 2) close OCR bins matched normal bins
# get states in intestine_16.5
normal_state_path="/data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/"
if [ -f  normal-bins_state_intestine_16.5.txt ];then
    rm  normal-bins_state_intestine_16.5.txt
fi
for i in {1..19} X Y
do
    awk '{FS=" ";OFS="\t"}{if(NR>1){print $1,$5}}' ${normal_state_path}run_IDEAS_8hm_atac_dname_pvalue.chr${i}.state >> \
    normal-bins_state_intestine_16.5.txt
done
# get states
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$3,$4,$5,a[$5],$6,a[$6]}}' normal-bins_state_intestine_16.5.txt \
OCR-center_bins_nearest_pairs_tmp.txt > OCR-center_bins_nearest_pairs_state.txt


# 3. make histogram
## 1) count state combination number
awk 'BEGIN{FS=OFS="\t";sum1=0;sum2=0}{if($2==$4){sum1+=1};if($6==$8){sum2+=1}}END{print sum1,sum2}' \
OCR-center_bins_nearest_pairs_state.txt
# 77921 68384

## 2) get combination
awk '{FS=OFS="\t"}{if($2<=$4){printf $2"-"$4}else{printf $4"-"$2};printf "\t";if($6<=$8){printf $6"-"$8}else{printf $8"-"$6};printf "\n"}' \
OCR-center_bins_nearest_pairs_state.txt > OCR-center_bins_nearest_pairs_state_combine.txt


# 4. figures
# Rscript plot_state_combination.R

#!/bin/bash

# -- Kaili
# This script is for making DHS-bins marks aggregation plot (corresponding regions).
# 0. get random OCR corresponding normal bins and gap regions.
# 1. cut bins into 50 windows, also up/downstream 500bp regions.
# 2. get signal for each small windows
# 3. get signal matrix, average signal matrix
# 4. make aggregation plot

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/"

cd ${workDir}

# # 0. get random OCR corresponding normal bins and gap regions.
# ## 1) get OCR matched normal bins
# awk '{FS=OFS="\t"}{center=int($2+($3-$2)/2);print $1,center,center,$4}' mm10_OCR-center_bins_v3_signal_based_OCR_random.bed > \
# mm10_OCR-center_bins_v3_signal_based_OCR_random_center.bed
# intersectBed -a /data/zusers/fankaili/ideas/run_ideas_p_value/mm10_tab.bed -b \
# mm10_OCR-center_bins_v3_signal_based_OCR_random_center.bed -wa | sort -u | sort -k1,1 -k2,2n \
# > mm10_OCR-center_matched_normal_random.bed
# # 100,489
#
# ## 2) get OCR flanking gap bins
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$1;b[$3]=$1;c[$2]=1;d[$3]=1}else{if(c[$3] && a[$3]==$1){print $0}else if(d[$2] && b[$2]==$1){print $0}}}' \
# mm10_OCR-center_bins_v3_signal_based_OCR_random.bed mm10_OCR-center_bins_v3_signal_based_gap.bed | \
# sort -u | sort -k1,1 -k2,2n > mm10_OCR-center_flanking_gap_random.bed
# # 166,781
#
#
# # 1. cut bins into 50 windows, also up/downstream 500bp regions.
# ## 1) get OCR matched normal bins
# bins="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random.bed"
# bins_out_file="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_bins.bed"
# upstream_out_file="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_up.bed"
# downstream_out_file="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_down.bed"
# nohup python ${scriptDir}cut_bins_into_50.py ${bins} ${bins_out_file} ${upstream_out_file} ${downstream_out_file} \
# > ${workDir}nohup.cut_bins_into_50_dhs_matched_normal.out 2>&1&
#
# ## 2) get OCR flanking gap bins
# bins="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_flanking_gap_random.bed"
# bins_out_file="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_flanking_gap_random_50windows_bins.bed"
# upstream_out_file="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_flanking_gap_random_50windows_up.bed"
# downstream_out_file="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_flanking_gap_random_50windows_down.bed"
# nohup python ${scriptDir}cut_bins_into_50.py ${bins} ${bins_out_file} ${upstream_out_file} ${downstream_out_file} \
# > ${workDir}nohup.cut_bins_into_50_dhs_flanking_gap.out 2>&1&
#
#
# # 2. get signal for each small windows
# get_signal(){
#     bedFile=$1
#     prefix=$2
#     #
#     codeFile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/code_for_getting_signal_p_value_"${prefix}".txt"
#     signalPath="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_"${prefix}"/"
#     #
#     python ${scriptDir}get_signal_value_code.py ${bedFile} ${codeFile} ${signalPath}
#     mkdir ${signalPath}
#     #
#     awk '{FS=OFS="\t"}{if(NR<1*133){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_1.sh
#     awk '{FS=OFS="\t"}{if(NR>(1*132) && NR<(2*132+1)){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_2.sh
#     awk '{FS=OFS="\t"}{if(NR>(2*132) && NR<(3*132+1)){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_3.sh
#     awk '{FS=OFS="\t"}{if(NR>(3*132) && NR<(4*132+1)){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_4.sh
#     awk '{FS=OFS="\t"}{if(NR>(4*132) && NR<(5*132+1)){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_5.sh
#     awk '{FS=OFS="\t"}{if(NR>(5*132) && NR<(6*132+1)){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_6.sh
#     awk '{FS=OFS="\t"}{if(NR>(6*132) && NR<(7*132+1)){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_7.sh
#     awk '{FS=OFS="\t"}{if(NR>(7*132) && NR<(8*132+1)){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_8.sh
#     awk '{FS=OFS="\t"}{if(NR>(8*132) && NR<(9*132+1)){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_9.sh
#     awk '{FS=OFS="\t"}{if(NR>(9*132) && NR<(10*132+1)){print $0}}' ${codeFile} > ./code_bins/code_for_getting_signal_p_value${prefix}_10.sh
#     #
#     for i in {1..10}
#     do
#         nohup bash ./code_bins/code_for_getting_signal_p_value${prefix}_${i}.sh > ./code_bins/nohup.code_for_getting_signal_p_value${prefix}_${i}.out 2>&1&
#     done
# }
#
# ## 1) get OCR matched normal bins
# bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_bins.bed"
# prefix="matched_normal_bins"
# get_signal ${bedfile} ${prefix}
# bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_up.bed"
# prefix="matched_normal_up"
# get_signal ${bedfile} ${prefix}
# bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_down.bed"
# prefix="matched_normal_down"
# get_signal ${bedfile} ${prefix}
#
# ## 2) get OCR flanking gap bins
# bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_flanking_gap_random_50windows_bins.bed"
# prefix="flanking_gap_bins"
# get_signal ${bedfile} ${prefix}
# bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_flanking_gap_random_50windows_up.bed"
# prefix="flanking_gap_up"
# get_signal ${bedfile} ${prefix}
# bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_flanking_gap_random_50windows_down.bed"
# prefix="flanking_down_down"
# get_signal ${bedfile} ${prefix}


# 3. get signal matrix, average signal matrix
## 1) get OCR matched normal bins
for file in `ls /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_matched_normal_up/*.txt`
do
    a=${file##/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_matched_normal_up/}
    filename=${a%.txt}
    echo ${filename}
    #
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_matched_normal_up/${filename}.txt > temp7.txt
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_matched_normal_bins/${filename}.txt > temp8.txt
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_matched_normal_down/${filename}.txt > temp9.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/100489}}' temp7.txt > \
    ./aggregation_signal/${filename}_matched_normal.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/100489}}' temp8.txt >> \
    ./aggregation_signal/${filename}_matched_normal.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/100489}}' temp9.txt >> \
    ./aggregation_signal/${filename}_matched_normal.txt
done

## 2) get OCR flanking gap bins
for file in `ls /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_flanking_gap_up/*.txt`
do
    a=${file##/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_flanking_gap_up/}
    filename=${a%.txt}
    echo ${filename}
    #
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_flanking_gap_up/${filename}.txt > temp1.txt
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_flanking_gap_bins/${filename}.txt > temp2.txt
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_flanking_gap_down/${filename}.txt > temp3.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/166781}}' temp1.txt > \
    ./aggregation_signal/${filename}_flanking_gap.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/166781}}' temp2.txt >> \
    ./aggregation_signal/${filename}_flanking_gap.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/166781}}' temp3.txt >> \
    ./aggregation_signal/${filename}_flanking_gap.txt
done

rm temp*.txt


# 4. make aggregation plot
for file in `ls /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_normal_up/*.txt`
do
    a=${file##/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/signal_normal_up/}
    filename=${a%.txt}
    echo ${filename}
    #
    Rscript ${scriptDir}make_matched_bins_aggregation_plot.R ${filename} "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_signal/"
done

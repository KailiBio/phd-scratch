#!/bin/bash

# -- Kaili
# This script is for making aggregation matrix for aggregation plot.

# INPUT:signal bigWig file.
#             bed file for each window. (file should be sorted and in standard bed4 format.)
#                       (name should be like region1_1, region1_2,...region1_50 then region2_1...)
#             num of windows.
#             isDNAme data? (1 or 0)
#             the path&name of output file. (This is a matrix with a list of num, can be use to plot)
# OUTPUT: a list of num that is the average signal of each windows.
# EXP: bash make_aggregation_matrix.sh /data/projects/encode/data/ENCSR100DLL/ENCFF574ZQY.bigWig
#           /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_OCR_random_50windows_merged_sorted.bed
#           150 0 /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/intestine_16.5_H3K4me2_OCR.txt

bwFile=$1
bedFile=$2
num=$3
isDNAme=$4
outfile=$5

############
# get signal for each bin
bigWigAverageOverBed ${bwFile} ${bedFile} tmp_make_aggregation_matrix_${bedFile}.tab
if [ isDNAme == 1 ]
then
    awk '{FS=OFS="\t"}{if(NR==FNR){if($3==0){a[$1]=-10}else{a[$1]=$6}}else{print $4,a[$4]}}' \
    tmp_make_aggregation_matrix_${bedFile}.tab ${bedFile} > tmp_make_aggregation_matrix_${bedFile}.txt
else
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $4,a[$4]}}' tmp_make_aggregation_matrix_${bedFile}.tab \
    ${bedFile} > tmp_make_aggregation_matrix_${bedFile}.txt
fi
# get average signal
awk '{FS=OFS="\t"}{split($1,a,"_");print $1,a[1],a[2],$2}' tmp_make_aggregation_matrix_${bedFile}.txt | sort -k2,2 -k3,3n \
> tmp_make_aggregation_matrix_${bedFile}_2.txt
awk '{FS=OFS="\t"}{if(NR==1){printf $4}else if(NR%150!=1){printf "\t"$4}else{printf "\n"$4}}END{printf "\n"}' \
tmp_make_aggregation_matrix_${bedFile}_2.txt > tmp_make_aggregation_matrix_${bedFile}_3.txt
awk 'BEGIN{FS=OFS="\t";for(i=1;i<=150;i++){b[i]=0}}{for(i=1;i<=150;i++){if($i>=0){a[i]+=$i;b[i]+=1}}}END{for(i=1;i<=150;i++){print a[i]/b[i]}}' \
tmp_make_aggregation_matrix_${bedFile}_3.txt > ${outfile}
